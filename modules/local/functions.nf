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
        reader.readLine()
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


def usage() {
    // [S57]: usage block printed on --help. Group flags by the part
    // that consumes them; the three entry-point flags up top double
    // as mode selectors.
    return """\
nf-metabarcoding — swarm-based metabarcoding pipeline

Usage:
  nextflow run main.nf --fastq_folder PATH      [Part A → B [→ C]]
  nextflow run main.nf --fasta_folder PATH      [Part B standalone]
  nextflow run main.nf --occurrence_table PATH  [Part C standalone]

Entry-point selection (set exactly one):
  --input PATH                validated samplesheet CSV (primary input);
                              fastq profile (sample,fastq_1,fastq_2,run)
                              runs Part A, fasta profile
                              (sample,fasta,qual,stats) runs Part B.
                              Mutually exclusive with the folder-scan
                              inputs below
  --fastq_folder PATH         input fastq directory (one or more,
                              comma-separated); runs Part A and, when
                              --project_name is also set, Part B
  --fasta_folder PATH         per-sample .fas + .qual + .stats
                              directory; runs Part B standalone
  --occurrence_table PATH     Part B's <basename>_table.tsv; runs
                              Part C standalone

Part A — fastq → per-sample fasta:
  --fastq_pattern GLOB        R1/R2 pair pattern (canonical patterns
                              auto-detected; default: ${params.fastq_pattern})
  --fastq_encoding INT        Phred offset, 33 or 64 (default: ${params.fastq_encoding})
  --hash_function NAME        hash used to rename amplicons: 'sha1'
                              (default) or 'md5' (default: ${params.hash_function})
  --forward_primer SEQ        IUPAC forward primer (required unless --no_trimming)
  --reverse_primer SEQ        IUPAC reverse primer (required unless --no_trimming)
  --no_trimming               skip cutadapt primer trimming (default: ${params.no_trimming})
  --stripright INT            3' nt trimmed in the shadow pipeline
                              before R1/R2 are joined (default: ${params.stripright})
  --join_padding_length INT   length of the A run inserted between R1
                              and R2 by the shadow pipeline's
                              vsearch --fastq_join. The same length is
                              used for the matching Phred-40 ('I')
                              quality string. The A run is artificial
                              padding, not biological sequence — keep
                              that in mind when interpreting shadow
                              output (default: ${params.join_padding_length})

Part B — per-sample fasta → occurrence table:
  --project_name NAME         basename prefix for Part B output
                              (required when Part B runs)
  --results_folder PATH       where Part B publishes its artefacts
                              (required when Part B runs)
  --chimera_minsize INT       minimum abundance for vsearch
                              --uchime_denovo (default: ${params.chimera_minsize})
  --fastidious BOOL           run swarm + cluster_cleaver.py with
                              --fastidious; default true. Set false on
                              very large datasets where the fastidious
                              pass is too memory- or time-expensive
                              (default: ${params.fastidious})
  --percentage FLOAT          cleaving threshold: a sub-seed candidate
                              must appear as a per-sample cluster seed in
                              at least this fraction of samples; real in
                              (0, 1] (default: ${params.percentage})

Part C — taxonomic assignment:
  --reference_dataset PATH    stampa-formatted reference fasta
                              (.gz / .bz2 OK). Required when
                              --taxonomy_method=stampa (default).
  --reference_dataset_sintax PATH
                              sintax-formatted reference fasta
                              (.gz / .bz2 OK). Required when
                              --taxonomy_method=sintax; gates the
                              shadow Part C path — when unset, shadow
                              Part C is silently skipped.
  --taxonomy_method NAME      regular-path method: 'stampa' (default)
                              or 'sintax'. The shadow path always uses
                              sintax regardless of this flag.
  --iddef INT                 vsearch --iddef value (default: ${params.iddef})
  --stampa_chunk_size INT     sequences per stampa scatter chunk;
                              0 disables the split (default: ${params.stampa_chunk_size})
  --stampa_maxrejects INT     vsearch --maxrejects for the stampa
                              scatter step; 0 = no limit
                              (default: ${params.stampa_maxrejects})
  --stampa_id FLOAT           vsearch --id threshold for the stampa
                              scatter step (default: ${params.stampa_id})
  --sintax_cutoff FLOAT       vsearch --sintax_cutoff threshold; used by
                              the shadow path and by --taxonomy_method
                              sintax on the regular path
                              (default: ${params.sintax_cutoff})
  --majority_assignment       opt-in final step: recompute a
                              majority-rule taxonomy per OTU from the
                              reference accessions and publish
                              <basename>_table_assigned_majority.tsv.
                              Requires --taxonomy_method=stampa
                              (default: ${params.majority_assignment})

Runtime:
  --threads INT               threads per process (default: ${params.threads})
  --publish_mode MODE         publishDir mode: copy, copyNoFollow, link,
                              move, rellink, symlink. Use 'copy' when the
                              results folder is on a different filesystem
                              than the work directory (default: ${params.publish_mode})
  --help                      show this message and exit

See README.md for examples and SPECIFICATIONS.md for the behaviour
contract ([Sxx] IDs).
""".stripIndent()
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
        header = stream.withReader { reader -> reader.readLine() }
    } else {
        header = ref_file.withReader { reader -> reader.readLine() }
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


def numeric_param_spec() {
    // [S72]: accepted range for each numeric parameter, pinned to the
    // vsearch option it feeds. Pure data (no params access) so
    // check_numeric_param() stays unit-testable. `minExclusive` marks a
    // half-open lower bound (used by `percentage`, which must be > 0).
    return [
        fastq_encoding   : [type: 'enum', values: [33, 64]],
        threads          : [type: 'int',  min: 1,   max: 256],
        percentage       : [type: 'real', min: 0.0, max: 1.0, minExclusive: true],
        chimera_minsize  : [type: 'int',  min: 1],
        stripright       : [type: 'int',  min: 0],
        iddef            : [type: 'int',  min: 0,   max: 4],
        stampa_chunk_size: [type: 'int',  min: 0],
        stampa_maxrejects: [type: 'int',  min: 0],
        stampa_id        : [type: 'real', min: 0.0, max: 1.0],
        sintax_cutoff    : [type: 'real', min: 0.0, max: 1.0],
    ]
}


def check_numeric_param(String name, value) {
    // [S72]: validate one numeric parameter against numeric_param_spec().
    // Returns the value unchanged when valid; throws
    // IllegalArgumentException with a message naming the parameter when
    // the value is the wrong type or out of range, so the run aborts at
    // startup instead of failing as an obscure vsearch error
    // mid-pipeline.
    def spec = numeric_param_spec()[name]
    if ( spec == null ) {
        throw new IllegalArgumentException(
            "no numeric range is defined for parameter '${name}'")
    }
    if ( spec.type == 'enum' ) {
        if ( !(value instanceof Number) || !((value as int) in spec.values) ) {
            throw new IllegalArgumentException(
                "--${name} must be one of ${spec.values}, got '${value}'")
        }
        return value
    }
    if ( !(value instanceof Number) ) {
        throw new IllegalArgumentException(
            "--${name} must be a number, got '${value}'")
    }
    if ( spec.type == 'int' ) {
        def i = value as int
        if ( i != value ) {
            throw new IllegalArgumentException(
                "--${name} must be an integer, got '${value}'")
        }
        def range = (spec.max != null)
            ? "in ${spec.min}..${spec.max}"
            : ">= ${spec.min}"
        if ( i < spec.min || (spec.max != null && i > spec.max) ) {
            throw new IllegalArgumentException(
                "--${name} must be an integer ${range}, got '${value}'")
        }
    } else {
        def d = value as double
        def lower_ok = spec.minExclusive ? (d > spec.min) : (d >= spec.min)
        def open = spec.minExclusive ? '(' : '['
        if ( !lower_ok || d > spec.max ) {
            throw new IllegalArgumentException(
                "--${name} must be a real number in ${open}${spec.min}, ${spec.max}], " +
                "got '${value}'")
        }
    }
    return value
}


def validate_numeric_params() {
    // [S72]: params-reading wrapper over check_numeric_param() — validate
    // every numeric parameter against numeric_param_spec(). Called from
    // validate_params() so an out-of-range value aborts before any
    // process is wired.
    numeric_param_spec().each { name, _spec ->
        check_numeric_param(name, params[name])
    }
}


def validate_params() {
    // Order-independent parameter guards, run once at the top of the
    // entry workflow (after the --help short-circuit) so an invalid value
    // aborts immediately with a clear message instead of surfacing as an
    // obscure tool error mid-pipeline. Mode-specific requirements
    // (fastq_folder / primers / reference datasets) stay inline in the
    // entry workflow because they depend on the selected run mode.

    // [S70]: --input (samplesheet) is mutually exclusive with the
    // folder-scan inputs; setting both is ambiguous, so abort up-front.
    assert !(params.input && (params.fastq_folder || params.fasta_folder)) :
        "--input is mutually exclusive with --fastq_folder / --fasta_folder"

    // [S71]: --results_folder is a deprecated alias for --outdir. Warn
    // once at startup when it is the active output root.
    if ( !params.outdir && params.results_folder ) {
        System.err.println(
            "WARNING: --results_folder is deprecated; use --outdir. " +
            "Outputs now follow the " +
            "<outdir>/{per_sample,occurrence_table,pipeline_info}/ layout."
        )
    }

    // [S58]: validate --publish_mode before any process is wired so a
    // typo aborts the run immediately with a clear message instead of
    // failing on the first PublishDir attempt.
    def allowed_modes = ['copy', 'copyNoFollow', 'link', 'move', 'rellink', 'symlink']
    assert params.publish_mode in allowed_modes :
        "--publish_mode must be one of ${allowed_modes}, got '${params.publish_mode}'"

    // [S61]: validate --taxonomy_method up-front. Only the regular
    // Part C path consults this flag; shadow Part C always uses sintax
    // ([S50]).
    def allowed_taxonomy_methods = ['stampa', 'sintax']
    assert params.taxonomy_method in allowed_taxonomy_methods :
        "--taxonomy_method must be one of ${allowed_taxonomy_methods}, got '${params.taxonomy_method}'"

    // [S66]: majority assignment recomputes a per-OTU taxonomy from the
    // reference accessions listed in the `references` column. Only the
    // stampa method populates that column ([S50] leaves it at the NA
    // placeholder), so --majority_assignment is incompatible with
    // --taxonomy_method=sintax. Fail fast before any process is wired.
    assert !(params.majority_assignment && params.taxonomy_method == 'sintax') :
        "--majority_assignment requires --taxonomy_method=stampa " +
        "(sintax leaves the references column unpopulated)"

    // [S65]: validate --hash_function up-front so an unsupported hash
    // aborts before any process is wired rather than surfacing as an
    // obscure vsearch --relabel error mid-pipeline.
    def allowed_hash_functions = ['sha1', 'md5']
    assert params.hash_function in allowed_hash_functions :
        "--hash_function must be one of ${allowed_hash_functions}, got '${params.hash_function}'"

    // [S63]: validate --join_padding_length as a positive integer
    // before any process is scheduled. Non-positive values, non-integers,
    // and other invalid inputs would otherwise surface as a confusing
    // vsearch error mid-pipeline.
    def jpl = params.join_padding_length
    assert (jpl instanceof Number) && (jpl as int) == jpl && (jpl as int) >= 1 :
        "--join_padding_length must be a positive integer, got '${jpl}'"

    // [S72]: range-validate every numeric parameter against
    // numeric_param_spec(). Kept after the mode/string guards above
    // (which produce more specific messages) so they fire first.
    validate_numeric_params()

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
