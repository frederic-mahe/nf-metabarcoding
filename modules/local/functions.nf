// Shared Groovy helper functions for nf-metabarcoding.
//
// Extracted verbatim from main.nf so they can be `include`d by the
// entry workflow and (as the per-part split proceeds) by the process
// modules that reference normalize_path() in their publishDir / storeDir
// directives. No behaviour change — pure code movement.

def normalize_path(value) {
    // [S60]: expand a leading shell-style `~` in a path string. Nextflow's
    // `file()` does not perform tilde expansion (that is a shell
    // construct), so a quoted '~/...' on the CLI or a '~/...' read from
    // a -params-file would otherwise be joined to launchDir and yield a
    // dangling work-dir symlink. Handles four shapes:
    //   - null         -> null   (preserve "unset")
    //   - List         -> recurse element-wise
    //   - String       -> expand if it starts with `~` (or `~/`, `~user/`);
    //                     comma-separated strings are split, expanded
    //                     per-segment, and re-joined so the downstream
    //                     `split(',')` callers keep working.
    //   - other        -> coerce to String and apply the rules above
    if ( value == null ) {
        return null
    }
    if ( value instanceof List ) {
        return value.collect { normalize_path(it) }
    }
    def s = value.toString()
    if ( s.contains(',') ) {
        return s.split(',')
            .collect { it.trim() }
            .findAll { it }
            .collect { normalize_path(it) }
            .join(',')
    }
    if ( s == '~' ) {
        return System.getProperty('user.home')
    }
    if ( s.startsWith('~/') ) {
        return System.getProperty('user.home') + s.substring(1)
    }
    if ( s.startsWith('~') ) {
        // `~user` or `~user/<rest>` — best-effort lookup.
        def slash = s.indexOf('/')
        def user
        def rest
        if ( slash < 0 ) {
            user = s.substring(1)
            rest = ''
        } else {
            user = s.substring(1, slash)
            rest = s.substring(slash)
        }
        def home = lookup_user_home(user)
        if ( home ) {
            return home + rest
        }
        // unknown user — pass through unchanged
    }
    return s
}


def lookup_user_home(String user) {
    // [S60] helper: resolve `~user` by trying `getent passwd <user>`
    // first, then falling back to a /etc/passwd scan. Returns null
    // when the user can't be found.
    // tokenize(':') drops empty fields, which shifts the home-dir
    // index when the GECOS column is blank (a common configuration).
    // split(':', -1) preserves every column, so fields[5] is always the
    // home directory per the passwd(5) layout.
    try {
        def proc = ['getent', 'passwd', user].execute()
        proc.waitFor()
        if ( proc.exitValue() == 0 ) {
            def fields = proc.text.trim().split(':', -1)
            if ( fields.size() >= 6 ) {
                return fields[5]
            }
        }
    } catch (Exception ignored) {
        // getent unavailable (e.g. non-Linux host) — fall through.
    }
    def passwd = new File('/etc/passwd')
    if ( passwd.canRead() ) {
        def match = passwd.readLines()
            .collect { it.split(':', -1) }
            .find { it.size() >= 6 && it[0] == user }
        if ( match ) {
            return match[5]
        }
    }
    return null
}


def effective_outdir(outdir, results_folder) {
    // [S71]: resolve the output root. Precedence: --outdir, then the
    // deprecated --results_folder alias, then the default 'results'.
    // Pure (no params access) so it is unit-testable; resolve_outdir()
    // is the params-reading wrapper. A leading `~` is expanded ([S60]).
    if ( outdir ) {
        return normalize_path(outdir)
    }
    if ( results_folder ) {
        return normalize_path(results_folder)
    }
    return 'results'
}


def resolve_outdir() {
    // [S71]: params-reading wrapper over effective_outdir(). Silent —
    // the --results_folder deprecation warning is emitted once at
    // startup by validate_params(), not here (this runs per publishDir
    // evaluation).
    return effective_outdir(params.outdir, params.results_folder)
}


def publish_dir(part) {
    // [S71]: publishDir target for a given output subdirectory
    // ('per_sample' / 'occurrence_table' / 'pipeline_info').
    return "${resolve_outdir()}/${part}"
}


def log_dir(part) {
    // [S71]/D16: publishDir target for step logs, grouped under a
    // dedicated logs/ tree organised by the producing pipeline stage
    // ('part_a/per_sample' / 'part_a' / 'part_b' / 'part_c' — D16
    // superseded D15's data-mirroring layout). Data files keep
    // publish_dir(); logs go here so the two never interleave.
    return "${resolve_outdir()}/logs/${part}"
}


def hash_relabel_flag() {
    // [S65]: the vsearch relabel flag matching params.hash_function,
    // threaded into filter_and_convert_to_fasta as an explicit input.
    // params.hash_function is validated at startup, so only the two
    // expected values reach here.
    return ( params.hash_function == 'md5' ) ? '--relabel_md5' : '--relabel_sha1'
}


def hash_id_length() {
    // [S65]: hex-string width of a params.hash_function digest — SHA1
    // is 40 chars, MD5 is 32 — threaded into the .qual dedup steps as
    // an explicit input so they feed `uniq --check-chars` the right
    // width without assuming a fixed hash length. (Helper functions are
    // not visible inside process shell blocks, hence the input route.)
    return ( params.hash_function == 'md5' ) ? 32 : 40
}


def samplesheet_profile() {
    // [S70]: infer the --input samplesheet profile ('fastq' → Part A,
    // 'fasta' → Part B) from its header columns so the entry-point
    // dispatch can route. Reads only the header line of the
    // user-supplied --input file (a declared input, not a folder scan
    // or an optional-sibling probe), so the parse-time read is
    // deliberate and minimal. Returns null when --input is unset.
    if ( !params.input ) {
        return null
    }
    def header
    header = file(normalize_path(params.input)).withReader { reader ->
        read_bounded_line(reader)
    }
    def cols
    cols = header ? header.split(',').collect { it.trim() } : []
    if ( 'fastq_1' in cols ) {
        return 'fastq'
    }
    if ( 'fasta' in cols ) {
        return 'fasta'
    }
    throw new IllegalArgumentException(
        "--input samplesheet header must contain a 'fastq_1' (fastq " +
        "profile) or 'fasta' (fasta profile) column; got ${cols}"
    )
}


def check_primer_format(String name, value) {
    // [S74]: validate a primer is an IUPAC nucleotide string — the codes
    // A C G T U R Y S W K M B D H V N plus I (inosine), upper- or
    // lower-case, at least 3 nt — the same alphabet reverse_complement.sh
    // understands. Returns the value when valid; throws
    // IllegalArgumentException naming the parameter on a malformed value,
    // so a typo aborts at startup rather than as a confusing cutadapt
    // error mid-run.
    if ( !(value ==~ /(?i)^[ACGTURYSWKMBDHVNI]{3,}$/) ) {
        throw new IllegalArgumentException(
            "--${name} must be an IUPAC nucleotide string (A C G T U plus " +
            "the ambiguity codes R Y S W K M B D H V N and I), at least 3 " +
            "characters, got '${value}'")
    }
    return value
}


def check_accession_format(value) {
    // [S98]: --accession accepts only bioproject and study accessions:
    //   bioproject — PRJ + (E|D|N) + one uppercase letter + digits
    //                (PRJEB / PRJNA / PRJDB)
    //   study      — (E|D|S) + RP + at least six digits (ERP / DRP / SRP)
    // Both patterns are anchored and case-sensitive, so a run / sample /
    // experiment accession — or a value that merely contains a valid
    // accession — is rejected. Returns the value when valid; throws
    // IllegalArgumentException naming the offending value so a bad
    // --accession aborts at startup before any network access, rather
    // than surfacing as an obscure fastq-dl error mid-run.
    def bioproject
    bioproject = value ==~ /^PRJ(E|D|N)[A-Z][0-9]+$/
    def study
    study = value ==~ /^(E|D|S)RP[0-9]{6,}$/
    if ( !bioproject && !study ) {
        throw new IllegalArgumentException(
            "--accession must be a bioproject (PRJEB/PRJNA/PRJDB) or study " +
            "(ERP/DRP/SRP) accession, got '${value}'")
    }
    return value
}


def parse_accessions(value) {
    // [S97]: normalise --accession into a list of accession strings. The
    // CLI form is a single accession or a comma-separated list; a
    // nextflow.config may instead supply a Groovy list. Trims each entry
    // and drops empties (so a trailing comma is harmless). Returns [] for
    // an unset value.
    if ( value == null ) {
        return []
    }
    def raw
    raw = ( value instanceof List )
        ? value.collect { it.toString() }
        : value.toString().split(',') as List
    return raw.collect { it.trim() }.findAll { it }
}


def read_bounded_line(reader, int max_chars = 65536) {
    // [S94]: read at most `max_chars` characters from `reader`, stopping
    // at the first newline. Bounds memory — and, for a gzip stream, how
    // much of it is decompressed — so a crafted input whose first line is
    // gigabytes long, or a decompression bomb whose first "line" inflates
    // without bound, cannot OOM the launcher during the startup header
    // sniffs ([S70] samplesheet profile / [S73] reference format). The
    // 64 KiB cap dwarfs any real FASTA or CSV header, so a legitimate
    // header is never truncated. Returns null on immediate EOF (empty
    // input); otherwise the line without its trailing newline (a trailing
    // carriage return is stripped so a CRLF-terminated header reads the
    // same as an LF one), matching BufferedReader.readLine() for the short
    // headers these sniffs read.
    // A single bounded read pulls at most `max_chars` characters out of
    // the (decompressed) stream — never the whole thing — so no explicit
    // loop is needed and a huge first line or a decompression bomb cannot
    // be materialised. We then cut at the first newline.
    def cb = java.nio.CharBuffer.allocate(max_chars)
    def n = reader.read(cb)
    if ( n < 0 ) {
        return null
    }
    cb.flip()
    def chunk = cb.toString()
    def nl = chunk.indexOf('\n')
    def line = ( nl >= 0 ) ? chunk.substring(0, nl) : chunk
    if ( line.endsWith('\r') ) {  // CRLF — strip the carriage return too
        line = line.substring(0, line.length() - 1)
    }
    return line
}


def reference_first_header(ref_file, boolean gzipped) {
    // [S73] helper: return the first line of a reference fasta (its FASTA
    // header), transparently decompressing gzip. Reads a single line —
    // a FASTA header is line 1 — using the same withReader idiom as
    // samplesheet_profile(). bzip2 is not handled here (the caller skips
    // `.bz2`).
    def header
    if ( gzipped ) {
        def stream
        stream = new java.util.zip.GZIPInputStream(ref_file.newInputStream())
        header = stream.withReader { reader -> read_bounded_line(reader) }
    } else {
        header = ref_file.withReader { reader -> read_bounded_line(reader) }
    }
    return header
}


def check_reference_format(value, String expected) {
    // [S73]: sniff the first FASTA header of a supplied reference and
    // verify it matches `expected` ('stampa' | 'sintax'). Returns the
    // value unchanged when valid (or unsniffable); throws
    // IllegalArgumentException naming the flag on a format mismatch.
    //   - stampa: `>id <space-separated lineage>` (header has whitespace)
    //   - sintax: `>id;tax=d:...,p:...;`           (header carries `tax=`)
    // A missing file is left for the [S47] presence assert / file()
    // staging; a `.bz2` reference is skipped (no pure-Groovy bzip2
    // decompressor) with a warning.
    def flag
    flag = ( expected == 'sintax' ) ? 'reference_dataset_sintax' : 'reference_dataset'
    def path_str
    path_str = normalize_path(value).toString()
    def ref_file
    ref_file = file(path_str)
    if ( !ref_file.exists() ) {
        return value
    }
    if ( path_str.endsWith('.bz2') ) {
        System.err.println(
            "WARNING: --${flag} is bzip2-compressed; skipping the startup " +
            "format check (vsearch reads it at runtime).")
        return value
    }
    def header
    header = reference_first_header(ref_file, path_str.endsWith('.gz'))
    if ( header == null || !header.startsWith('>') ) {
        throw new IllegalArgumentException(
            "--${flag} does not look like a FASTA file (first non-blank " +
            "line is not a '>' header): ${value}")
    }
    def body
    body = header.substring(1)
    if ( expected == 'sintax' ) {
        if ( !body.contains('tax=') ) {
            throw new IllegalArgumentException(
                "--${flag} must be sintax-formatted (header " +
                "'>id;tax=d:...,p:...;'); its first header carries no " +
                "'tax=' annotation: ${header}")
        }
    } else if ( !(body =~ /\s/) ) {
        throw new IllegalArgumentException(
            "--${flag} must be stampa-formatted (header " +
            "'>id <space-separated lineage>'); its first header has no " +
            "lineage after the id: ${header}")
    }
    return value
}


def validate_params() {
    // Order-independent parameter guards, run once at the top of the
    // entry workflow (after the --help short-circuit) so an invalid value
    // aborts immediately with a clear message instead of surfacing as an
    // obscure tool error mid-pipeline. Mode-specific requirements
    // (fastq_folder / primers / reference datasets) stay inline in the
    // entry workflow because they depend on the selected run mode.

    // [S02]/[S70]/[S48]: the five input-mode selectors are mutually
    // exclusive. Setting more than one is ambiguous (the dispatcher would
    // silently pick one and drop the rest), so abort up-front naming the
    // ones that were set. Subsumes the earlier --input-vs-folder check
    // ([S70]).
    def mode_selectors = [
        occurrence_table   : params.occurrence_table,
        representatives_fasta: params.representatives_fasta,
        input              : params.input,
        fasta_folder       : params.fasta_folder,
        fastq_folder       : params.fastq_folder,
        accession          : params.accession,
    ]
    def selectors_set
    selectors_set = mode_selectors
        .findAll { _name, value -> value }
        .collect { name, _value -> "--${name}" }
    assert selectors_set.size() <= 1 :
        "the input-mode selectors are mutually exclusive; set at most one " +
        "of --occurrence_table / --representatives_fasta / --input / " +
        "--fasta_folder / --fastq_folder / --accession, got: ${selectors_set}"

    // [S98]: validate every --accession entry up-front (a pure param
    // check, no network) so a bad accession aborts at startup naming the
    // offending value rather than surfacing as an obscure fastq-dl error
    // mid-download.
    parse_accessions(params.accession).each { check_accession_format(it) }

    // [S71]: --results_folder is a deprecated alias for --outdir. Warn
    // once at startup when it is the active output root.
    if ( !params.outdir && params.results_folder ) {
        System.err.println(
            "WARNING: --results_folder is deprecated; use --outdir. " +
            "Outputs now follow the " +
            "<outdir>/{per_sample,occurrence_table,pipeline_info}/ layout."
        )
    }

    // [S58]/[S61]/[S65]/[S63]/[S72]: the per-value type / enum / range
    // checks (--publish_mode, --taxonomy_method, --hash_function,
    // --join_padding_length, and the numeric ranges) are declared in
    // nextflow_schema.json and enforced by validateParameters(), which
    // the entry workflow calls *before* this function. The guards that
    // remain here are cross-parameter or file-content checks the schema
    // cannot express.

    // [S66]: majority assignment recomputes a per-OTU taxonomy from the
    // reference accessions listed in the `references` column. Only the
    // stampa method populates that column ([S50] leaves it at the NA
    // placeholder), so --majority_assignment is incompatible with
    // --taxonomy_method=sintax. Fail fast before any process is wired.
    assert !(params.majority_assignment && params.taxonomy_method == 'sintax') :
        "--majority_assignment requires --taxonomy_method=stampa " +
        "(sintax leaves the references column unpopulated)"

    // [S48]/[S66]: majority assignment recomputes a per-OTU taxonomy from
    // the assigned occurrence table. Fasta-input Part C
    // (--representatives_fasta) produces no occurrence table, so there is
    // nothing to compute a majority on. Fail fast rather than silently
    // ignore the flag.
    assert !(params.majority_assignment && params.representatives_fasta) :
        "--majority_assignment cannot be combined with --representatives_fasta " +
        "(fasta-input Part C produces no occurrence table to compute a " +
        "per-OTU majority on)"

    // [S102]/D20: --recluster_id is the master switch for the optional
    // post-mumu re-clustering pass; --recluster_iddef only tunes that
    // pass. Setting --recluster_iddef to a non-default value while
    // --recluster_id is unset would silently do nothing, so abort naming
    // the switch. The `2` mirrors the recluster_iddef default declared in
    // nextflow.config's params block (a user who explicitly passes the
    // default value is harmless and not flagged).
    if ( params.recluster_id == null && params.recluster_iddef != 2 ) {
        throw new IllegalArgumentException(
            "--recluster_iddef requires --recluster_id (the re-clustering " +
            "master switch): set --recluster_id to a real in (0, 1] to " +
            "enable the optional post-mumu re-clustering pass, or leave " +
            "--recluster_iddef at its default.")
    }

    // [S73]: sniff the header of each supplied reference so a swapped /
    // mis-formatted file aborts now rather than producing empty or wrong
    // assignments mid-pipeline. Format-only — presence is mode-specific
    // ([S47]). Kept last (after the numeric guards) so a more specific
    // value error still wins.
    if ( params.reference_dataset ) {
        check_reference_format(params.reference_dataset, 'stampa')
    }
    if ( params.reference_dataset_sintax ) {
        check_reference_format(params.reference_dataset_sintax, 'sintax')
    }
}


def resource_size_warnings(profile, dataset_size_gb, reference_size_gb, boolean reference_in_use) {
    // [S79]: the slurm profile scales the memory-bound steps off
    // --dataset_size_gb (the dataset-bound steps) and --reference_size_gb
    // (the taxonomic-assignment steps); when those hints are unset the
    // steps fall back to fixed defaults that can be far too small for a
    // large run and OOM mid-pipeline. Return a startup warning per unset
    // hint so a forgotten flag surfaces immediately instead of as a kill
    // deep in the run. Pure (no params / no `workflow` access) so it is
    // unit-testable; the entry workflow reads `workflow.profile` + params
    // and prints whatever this returns. The resource tiers only exist
    // under `-profile slurm`, so outside it this is silent.
    def warnings = []
    if ( !profile || !profile.toString().tokenize(',').contains('slurm') ) {
        return warnings
    }
    if ( !dataset_size_gb ) {
        warnings << ("--dataset_size_gb is unset: the memory-bound steps " +
            "(global_clustering, global_dereplication, chimera detection, " +
            "the all-vs-all search, mumu, occurrence-table assembly) fall " +
            "back to fixed memory that may be too small for a large dataset " +
            "and OOM mid-run. Set --dataset_size_gb to the approximate " +
            "input size in GB (or override the step's memory in a " +
            "-c site.config).")
    }
    if ( reference_in_use && !reference_size_gb ) {
        warnings << ("--reference_size_gb is unset: the taxonomic-assignment " +
            "steps fall back to fixed memory that may be too small for a " +
            "large reference database. Set --reference_size_gb to the " +
            "approximate reference size in GB.")
    }
    return warnings
}


def slurm_account_requirement_error(require_account, slurm_account) {
    // [S92]: a cluster profile may set params.require_slurm_account=true
    // when the site rejects jobs that are not charged to a project (the
    // abims profile does). Return an error message naming --slurm_account
    // when the requirement is in force but no account is set, so the run
    // aborts at startup instead of letting every sbatch bounce mid-run;
    // return null when the requirement is off or the account is supplied.
    // Pure (no params / no `workflow` access) so it is unit-testable; the
    // entry workflow reads params.require_slurm_account + params.slurm_account
    // and throws whatever this returns.
    if ( !require_account ) {
        return null
    }
    if ( slurm_account && slurm_account.toString().trim() ) {
        return null
    }
    return "--slurm_account is required on this cluster (the active profile " +
        "sets require_slurm_account): pass your project, e.g. " +
        "--slurm_account jedi_meta."
}
