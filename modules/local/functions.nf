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
