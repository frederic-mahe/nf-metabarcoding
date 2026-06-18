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
}
