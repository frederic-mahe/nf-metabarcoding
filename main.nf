#!/usr/bin/env nextflow

// [S57]: --help short-circuits the workflow before any required-param
// assertion fires, so users can discover the interface without first
// having to supply a mode-specific flag.
params.help = false

// [S58]: publishDir mode applied to every publishDir directive in the
// workflow. Default 'link' (fast hard link) preserves prior behaviour
// on same-filesystem layouts; switch to 'copy' when --fastq_folder or
// --results_folder is on a different filesystem than the work dir
// (hard links cannot cross devices).
params.publish_mode = 'link'

// optional (project defaults)
params.fastq_pattern = "/*_1_{1,2}.fastq.gz"
params.fastq_encoding = 33
params.threads = 4
params.no_trimming = false
// [S24] number of nucleotides stripped from the 3' end of each
// not-merged R1/R2 read before they are joined (shadow pipeline).
// Set to 0 to disable.
params.stripright = 30
// [S63] length (in nt) of the run of `A`s that vsearch --fastq_join
// inserts between R1 and R2 in the shadow Part A path ([S04]). The
// same length is used for the matching --join_padgapq Phred-40 ('I')
// quality string. Increase to make the artificial join site more
// visible to a human reader at the cost of slightly larger fastas.
params.join_padding_length = 8

// [S65] hash vsearch uses to rename amplicons in
// filter_and_convert_to_fasta (--relabel_sha1 / --relabel_md5).
// 'sha1' (default) or 'md5'; validated at workflow startup. The .qual
// dedup steps derive their uniq --check-chars width from this so no
// step assumes a fixed hash-string length.
params.hash_function = 'sha1'

// [S65]: single source of truth derived from hash_function and read
// directly in the process shell blocks (script-level functions are
// not visible there; params are). These are computed, not user knobs
// — hash_function itself is the validated input ([S65] assert below).
//   - hash_relabel_flag: the vsearch flag (--relabel_sha1/--relabel_md5)
//   - hash_id_length:     the digest hex width (SHA1 → 40, MD5 → 32),
//     fed to `uniq --check-chars` so the .qual dedup steps collapse
//     rows by the whole amplicon name without assuming a fixed width.
params.hash_relabel_flag = ( params.hash_function == 'md5' ) ? '--relabel_md5' : '--relabel_sha1'
params.hash_id_length    = ( params.hash_function == 'md5' ) ? 32 : 40

// forward_primer, reverse_primer, fastq_folder are required and have
// no default; the workflow asserts them at startup (see [S18]).
// [S20]: when params.no_trimming is true, forward_primer and
// reverse_primer must be empty and the trim_primers step is skipped.

// Part B parameters — required only when Part B runs.
// [S25] / [S26] / [S27]: setting --fasta_folder switches the
// workflow into Part B standalone mode; --project_name and
// --results_folder become required, and Part A is skipped.
params.fasta_folder    = null
params.project_name    = null
params.results_folder  = null
// [S34]: minimum size threshold for chimera detection
// (vsearch --fastx_filter --minsize). Records below this size are
// dropped before uchime_denovo runs.
params.chimera_minsize = 2
// [S22]/[S32]: propagation lever for the swarm --fastidious flag.
// When true (default), global_clustering's swarm and cleaving's
// cluster_cleaver.py both run fastidious and every artefact carries
// the `_1f` suffix; when false, both processes drop the flag and the
// suffix swaps to `_1`. Exposed for very large datasets where the
// fastidious pass is too expensive.
params.fastidious = true

// Part C parameters — required only when Part C runs.
// [S47]: path to a reference dataset (fasta, plain / .gz / .bz2).
// [S48]: occurrence table input (Part C standalone). When Parts B
// and C run end-to-end the Part C workflow consumes Part B's
// <basename>_table.tsv directly and this stays null.
// [S49] / [S50]: 'stampa' (default) drives the primary path,
// 'sintax' switches to the shadow path.
// [S49]: identity-definition flag passed to vsearch (iddef).
// [S49]: chunk size for the stampa scatter-gather. Each chunk of
// this many sequences runs vsearch + stampa_merge.py as a
// separate process. The default (1000) is slurm-tuned; the
// `local` profile in nextflow.config overrides this to 0
// (sentinel: bypass splitFasta and feed the whole fasta into
// one task).
// [S47]/[S64]: two reference flags, one per assignment method.
// `--reference_dataset` carries the stampa-formatted file (header
// `>id space-separated lineage`); consumed by the regular path when
// taxonomy_method='stampa' (default). `--reference_dataset_sintax`
// carries the sintax-formatted file (header `>id;tax=d:...,p:...;`);
// consumed by `part_C_shadow` (always sintax) and by the regular
// path when taxonomy_method='sintax'. The split is a stop-gap until
// a single multi-format reference loader lands — until then, the
// two formats differ in incompatible ways and a single file cannot
// satisfy both methods in production use.
params.reference_dataset        = null
params.reference_dataset_sintax = null
params.occurrence_table  = null
params.taxonomy_method   = 'stampa'
params.iddef             = 1
params.stampa_chunk_size = 1000
// [S49]: vsearch knobs exposed for the stampa scatter step.
// stampa_maxrejects=0 is vsearch's "no limit" sentinel — exhaustive
// search, slower but every reference seq is considered as a potential
// top hit. stampa_id is the --id threshold below which a hit is
// dropped (stampa_merge.py then emits "No_hit" for the amplicon).
params.stampa_maxrejects = 0
params.stampa_id         = 0.5
// [S50]: --sintax_cutoff threshold for vsearch's sintax classifier.
// Used by the shadow Part C path (always) and by the regular Part C
// path when params.taxonomy_method == 'sintax'.
params.sintax_cutoff     = 0.9


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


process merge_fastq_pairs {
    // --fastqout_notmerged_fwd/_rev capture reads that fail to merge;
    // they feed the shadow Part A pipeline ([S04]). Fwd and rev are
    // kept in sync by vsearch.
    publishDir path: { normalize_path(params.fastq_folder) }, mode: params.publish_mode, pattern: "*.log",
        enabled: params.fastq_folder != null

    input:
    tuple val(sampleId), path(fastq_pair)

    output:
    val sampleId
    path "merged_fastq"
    path "${sampleId}_merging.log"
    path "notmerged_fwd"
    path "notmerged_rev"

    shell:
    '''
    #!/bin/bash

    vsearch \
        --fastq_mergepairs !{fastq_pair[0]} \
        --reverse !{fastq_pair[1]} \
        --threads !{params.threads} \
        --fastq_ascii !{params.fastq_encoding} \
        --fastq_allowmergestagger \
        --quiet \
        --log !{sampleId}_merging.log \
        --fastqout merged_fastq \
        --fastqout_notmerged_fwd notmerged_fwd \
        --fastqout_notmerged_rev notmerged_rev
    '''
}


process strip_reads {
    // [S24] shadow pipeline only: trim params.stripright nucleotides
    // from the 3' end of each R1 and R2 not-merged read. Used to
    // discard the low-quality tails before --fastq_join so the join
    // padding is surrounded by higher-quality bases. vsearch
    // --fastq_stripright 0 is a valid no-op pass-through.

    input:
    val sampleId
    path notmerged_fwd
    path notmerged_rev

    output:
    val sampleId
    path "stripped_fwd"
    path "stripped_rev"

    shell:
    '''
    #!/bin/bash

    vsearch \
        --fastx_filter !{notmerged_fwd} \
        --fastq_ascii !{params.fastq_encoding} \
        --fastq_stripright !{params.stripright} \
        --quiet \
        --fastqout stripped_fwd

    vsearch \
        --fastx_filter !{notmerged_rev} \
        --fastq_ascii !{params.fastq_encoding} \
        --fastq_stripright !{params.stripright} \
        --quiet \
        --fastqout stripped_rev
    '''
}


process join_notmerged {
    // Shadow pipeline ([S04]) entry point: concatenate fwd/rev reads
    // that failed --fastq_mergepairs with a run of As (length
    // params.join_padding_length, default 8 — see [S63]) so they can
    // be processed by the rest of Part A as a single fastq. `A` is
    // used instead of vsearch's default `N` so the joined sequence
    // carries only A/C/G/T and swarm accepts it as-is later in the
    // shadow Part B path ([S56]) — no mask/restore round-trip needed.
    //
    // [S04]: the shadow path has no merging step — by definition the
    // reads in this branch could not be merged — so no `_merging.log`
    // is produced or published. vsearch is invoked without --log; the
    // three remaining shadow per-step logs (trimming / dereplicating
    // / clustering) reach `params.fastq_folder` through the regular
    // downstream processes.

    input:
    val sampleId
    path notmerged_fwd
    path notmerged_rev

    output:
    val sampleId
    path "joined_fastq"

    shell:
    '''
    #!/bin/bash

    vsearch \
        --fastq_join !{notmerged_fwd} \
        --reverse !{notmerged_rev} \
        --fastq_ascii !{params.fastq_encoding} \
        --join_padgap  !{'A' * (params.join_padding_length as int)} \
        --join_padgapq !{'I' * (params.join_padding_length as int)} \
        --quiet \
        --fastqout joined_fastq
    '''
}


process trim_primers {
    // search forward primer in both normal and revcomp: now all reads
    // are in the same orientation. Matching leftmost is the default.
    // Length and N-count filtering are delegated to
    // filter_and_convert_to_fasta (vsearch --fastq_minlen / --fastq_maxns).
    publishDir path: { normalize_path(params.fastq_folder) }, mode: params.publish_mode, pattern: "*.log",
        enabled: params.fastq_folder != null

    input:
    val sampleId
    path merged_fastq

    output:
    val sampleId
    path "trimmed_fastq"
    path "${sampleId}_trimming.log"

    shell:
    '''
    #!/bin/bash

    readonly ERROR_RATE=0.1

    reverse_primer_revcomp=$(reverse_complement.sh !{params.reverse_primer})

    MIN_F=$(( !{params.forward_primer.length()} * 2 / 3 ))  # match is >= 2/3 of primer length
    MIN_R=$(( !{params.reverse_primer.length()} * 2 / 3 ))
    {
        cutadapt \
            --cores=!{params.threads} \
            --error-rate "${ERROR_RATE}" \
            --revcomp \
            --rename="{id}" \
            --front "!{params.forward_primer};rightmost" \
            --overlap "${MIN_F}" \
            --discard-untrimmed \
            !{merged_fastq} | \
            cutadapt \
                --cores=!{params.threads} \
                --error-rate "${ERROR_RATE}" \
                --adapter "${reverse_primer_revcomp}" \
                --overlap "${MIN_R}" \
                --discard-untrimmed \
                - > trimmed_fastq
    } 2> !{sampleId}_trimming.log
    '''
}


process filter_and_convert_to_fasta {
    // use the params.hash_function digest (SHA1 or MD5, see [S65]) as
    // sequence names, compute expected error
    // values (ee), and apply the minimum-length / max-N filter
    // (--fastq_maxns 0: every N drops the read). The shadow path
    // ([S04]) pads with A/C/G/T (see [S63]), so the same max-N=0
    // threshold serves both the regular and the shadow path.
    input:
    val sampleId
    path trimmed_fastq

    output:
    val sampleId
    path "filtered_fasta"

    shell:
    '''
    #!/bin/bash

    readonly MIN_LENGTH=32

    vsearch \
        --fastq_filter !{trimmed_fastq} \
        --fastq_minlen "${MIN_LENGTH}" \
        --fastq_maxns 0 \
        !{params.hash_relabel_flag} \
        --fastq_ascii !{params.fastq_encoding} \
        --quiet \
        --eeout \
        --lengthout \
        --fasta_width 0 \
        --fastaout - > filtered_fasta
    '''
}


process extract_expected_error_values {
    // extract ee for future quality filtering (keep the lowest
    // observed expected error value for each unique sequence)
    publishDir path: { normalize_path(params.fastq_folder) }, mode: params.publish_mode

    input:
    val sampleId
    path filtered_fasta

    output:
    val sampleId
    path "${sampleId}.qual"

    shell:
    '''
    length_of_sequence_IDs=!{params.hash_id_length}
    extract_ee.awk !{filtered_fasta} | \
        sort --key=3,3n --key=1,1d --key=2,2n | \
        uniq --check-chars=${length_of_sequence_IDs} > !{sampleId}.qual
    '''
}


process dereplicate_fasta {
    // dereplicate and discard expected error values (ee)
    publishDir path: { normalize_path(params.fastq_folder) }, mode: params.publish_mode

    input:
    val sampleId
    path filtered_fasta

    output:
    val sampleId
    path "${sampleId}.fas"
    path "${sampleId}_dereplicating.log"

    shell:
    '''
    vsearch \
        --fastx_uniques !{filtered_fasta} \
        --sizeout \
        --quiet \
        --log !{sampleId}_dereplicating.log \
        --fasta_width 0 \
        --xee \
        --xlength \
        --fastaout - > !{sampleId}.fas
    '''
}


// refactoring:
// breakdown into two functions
// - clusterize_sample
// - list_per_sample_clusters

process list_local_clusters {
    // retain only clusters with more than 2 reads
    // (do not use the fastidious option here)
    publishDir path: { normalize_path(params.fastq_folder) }, mode: params.publish_mode

    input:
    val sampleId
    path dereplicated_fasta

    output:
    val sampleId
    path "${sampleId}.stats"
    path "${sampleId}_clustering.log"

    shell:
    '''
    #!/bin/bash

    swarm \
        --threads !{params.threads} \
        --differences 1 \
        --usearch-abundance \
        --log !{sampleId}_clustering.log \
        --output-file /dev/null \
        --statistics-file - \
        !{dereplicated_fasta} | \
        filter_swarm_stats.awk > !{sampleId}.stats
    '''
}


process discover_inputs {
    // Walk every folder listed in params.fastq_folder (a single path
    // or a comma-separated list) and emit a TSV:
    //     sample_id<TAB>r1<TAB>r2   (r2 empty for single-end samples)
    // bin/discover_fastq.py is the single source of truth for the
    // canonical paired-end pattern table ([S11]); params.fastq_pattern,
    // when set, is the user override checked before that table.

    output:
    path "samples.tsv"

    script:
    def folders = (params.fastq_folder instanceof List)
        ? params.fastq_folder.collect { normalize_path(it) }
        : params.fastq_folder
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it }
            .collect { normalize_path(it) }
    def folder_args = folders.collect { "'${it}'" }.join(' ')
    def extra_arg = (params.fastq_pattern
        && !params.fastq_pattern.toString().isEmpty())
            ? "--extra-pattern '${params.fastq_pattern}'"
            : ""
    """
    discover_fastq.py ${extra_arg} ${folder_args} > samples.tsv
    """
}


// ============================================================================
// Part B — global clustering / occurrence table
// ============================================================================
// Each process below receives the per-sample artefacts collected by
// the Part B fasta channel ([S27]) and a project-wide `basename`
// string of the shape `<project_name>_<N>_samples` (see [S25]). The
// processes can run in parallel — they do not depend on each other.

process discover_part_b_fasta {
    // [S27]/[S56]: walk params.fasta_folder, assert unique sample IDs,
    // and emit two TSVs:
    //   - fastas.tsv         — regular samples (excludes *_notmerged.fas)
    //   - shadow_fastas.tsv  — only *_notmerged.fas samples (shadow Part B)
    // Each row is `sample_id<TAB>fasta_path`. Either file may be empty
    // (e.g. shadow_fastas.tsv is empty when no _notmerged.fas exists).

    output:
    path "fastas.tsv"
    path "shadow_fastas.tsv"

    script:
    def folders = (params.fasta_folder instanceof List)
        ? params.fasta_folder.collect { normalize_path(it) }
        : params.fasta_folder
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it }
            .collect { normalize_path(it) }
    def folder_args = folders.collect { "'${it}'" }.join(' ')
    """
    discover_fasta.py          ${folder_args} > fastas.tsv
    discover_fasta.py --shadow ${folder_args} > shadow_fastas.tsv
    """
}


process build_expected_error_file {
    // [S28]: merge every per-sample <sampleId>.qual into one
    // project-wide <basename>.qual. The input files are already sorted
    // by length / hash / ee (see extract_expected_error_values), so
    // `sort --merge` is a straight k-way merge; uniq --check-chars
    // (width from params.hash_function, see [S65]) keeps the lowest-ee
    // row per amplicon name.
    //
    // [S59]: the .qual is an internal intermediate consumed by
    // build_occurrence_table — not published.

    input:
    path quals
    val basename

    output:
    path "${basename}.qual"

    shell:
    '''
    sort --key=3,3n --key=1,1d --key=2,2n --merge !{quals} | \
        uniq --check-chars=!{params.hash_id_length} > !{basename}.qual
    '''
}


process build_distribution_file {
    // [S29]: scan FASTA headers of every input .fas and emit the
    // sequence-to-sample mapping as tab-separated rows
    // <sha1>\t<sampleId>\t<size>. The sample ID is derived from the
    // fasta basename so the channel can be built without an explicit
    // sample-ID side car.
    //
    // [S59]: the .distr is an internal intermediate consumed by
    // build_occurrence_table — not published.

    input:
    path fastas
    val basename

    output:
    path "${basename}.distr"

    shell:
    '''
    : > !{basename}.distr
    for f in !{fastas} ; do
        sample="$(basename "${f}" .fas)"
        # `|| true`: empty samples ([S09]/[S27]) have a zero-record
        # .fas — grep returns 1, which would otherwise abort the
        # process. The empty sample legitimately contributes no rows.
        grep "^>" "${f}" | \
            sed 's/^>// ; s/;size=/\t/ ; s/;$//' | \
            awk -v s="${sample}" 'BEGIN {OFS = "\t"} {print $1, s, $2}' \
            >> !{basename}.distr || true
    done
    '''
}


process list_all_cluster_seeds_of_size_greater_than_2 {
    // [S30]: concatenate every per-sample <sampleId>.stats (already
    // filtered to clusters > 2 reads by Part A's list_local_clusters)
    // into a single project-wide file. Each row is prefixed with the
    // sample ID derived from the .stats filename.
    //
    // [S59]: this aggregated per-sample-OTUs stats file is an
    // internal intermediate consumed by cleaving — not published.

    input:
    path stats_files
    val basename

    output:
    path "${basename}_per_sample_OTUs.stats"

    shell:
    '''
    for f in !{stats_files} ; do
        sample="$(basename "${f}" .stats)"
        awk -v s="${sample}" 'BEGIN {OFS = "\t"} {print s, $0}' "${f}"
    done > !{basename}_per_sample_OTUs.stats
    '''
}


process global_dereplication {
    // [S31]: cat every input .fas and pass it through
    // vsearch --fastx_uniques. --sizein/--sizeout preserve the
    // per-sample size annotations so vsearch sums abundances across
    // samples. The vsearch log lands at <basename>_dereplication.log
    // ([S45]).
    //
    // [S59]: only the log reaches the results folder; the
    // dereplicated .fas is an internal intermediate.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode, pattern: "*.log",
        enabled: params.results_folder != null

    input:
    path fastas
    val basename

    output:
    path "${basename}.fas"
    path "${basename}_dereplication.log"

    shell:
    '''
    cat !{fastas} | \
        vsearch \
            --fastx_uniques - \
            --sizein \
            --sizeout \
            --fasta_width 0 \
            --quiet \
            --log !{basename}_dereplication.log \
            --fastaout !{basename}.fas
    '''
}


process global_clustering {
    // [S32]: swarm on the globally-dereplicated fasta. Output
    // filenames carry the swarm-parameters suffix `_1f`
    // (--fastidious, default) or `_1` (--no-fastidious), driven by
    // params.fastidious. See [S22]'s propagation lever. The swarm
    // log lands at <basename>_clustering.log ([S45]) — suffix-
    // independent.
    //
    // [S59]: only the log reaches the results folder; the
    // .swarms / .stats / .struct / _representatives.fas are
    // internal intermediates.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode, pattern: "*.log",
        enabled: params.results_folder != null

    input:
    path global_fasta
    val basename

    output:
    path "${basename}_${params.fastidious ? '1f' : '1'}.swarms"
    path "${basename}_${params.fastidious ? '1f' : '1'}.stats"
    path "${basename}_${params.fastidious ? '1f' : '1'}.struct"
    path "${basename}_${params.fastidious ? '1f' : '1'}_representatives.fas"
    path "${basename}_clustering.log"

    shell:
    def sfx = params.fastidious ? '1f' : '1'
    def fastidious_flag = params.fastidious ? '--fastidious' : ''
    """
    swarm \\
        --threads ${params.threads} \\
        --differences 1 \\
        ${fastidious_flag} \\
        --usearch-abundance \\
        --internal-structure ${basename}_${sfx}.struct \\
        --output-file ${basename}_${sfx}.swarms \\
        --statistics-file ${basename}_${sfx}.stats \\
        --seeds ${basename}_${sfx}_representatives.fas \\
        ${global_fasta} 2> ${basename}_clustering.log
    """
}


process fake_taxonomic_assignment {
    // [S33]: emit a placeholder TSV that the occurrence-table builder
    // can consume before Part C lands. Each row mirrors the stampa
    // output shape: <amplicon>\t<size>\t<identity>\t<taxonomy>\t<refs>.
    //
    // [S59]: placeholder taxonomy is an internal intermediate
    // consumed by build_occurrence_table — not published.

    input:
    path representatives
    val basename

    output:
    path "${basename}_1f_representatives.results"

    shell:
    '''
    grep "^>" !{representatives} | \
        sed -r 's/^>//
                s/;size=/\t/
                s/;?$/\t0.0\tNA\tNA/' > !{basename}_1f_representatives.results
    '''
}


process chimera_detection {
    // [S34]: filter representatives down to abundance >= chimera_minsize
    // (default 2), then run vsearch --uchime_denovo. The .uchime
    // hit table can be empty when no chimeras are found; the stderr
    // log captures the run.
    //
    // [S45] keeps this pre-cleave .log internal — the canonical
    // <basename>_chimera_detection.log is published by
    // chimera_detection_post_cleave and contains the concatenation of both
    // runs' stderr.
    //
    // [S59]: the .uchime hit table is an internal intermediate
    // consumed by build_occurrence_table — not published.

    input:
    path representatives
    val basename

    output:
    path "${basename}_1f_representatives.uchime"
    path "${basename}_1f_representatives.log"

    shell:
    '''
    vsearch \
        --fastx_filter !{representatives} \
        --minsize !{params.chimera_minsize} \
        --quiet \
        --fastaout - | \
    vsearch \
        --uchime_denovo - \
        --uchimeout !{basename}_1f_representatives.uchime \
        2> !{basename}_1f_representatives.log
    '''
}


process fake_taxonomic_assignment2 {
    // [S36]: same transform as [S33] but operating on the cleaver's
    // representatives (cleaving's third output, `*_1f_representatives.fas2`).
    // An empty fas2 (no clusters got cleaved) is a legitimate input
    // — the output is then an empty results2 file.
    //
    // [S59]: post-cleave placeholder taxonomy is an internal
    // intermediate consumed by build_occurrence_table — not
    // published.

    input:
    path cleaved_representatives
    val basename

    output:
    path "${basename}_1f_representatives.results2"

    shell:
    '''
    : > !{basename}_1f_representatives.results2
    grep "^>" !{cleaved_representatives} | \
        sed -r 's/^>//
                s/;size=/\t/
                s/;?$/\t0.0\tNA\tNA/' \
        >> !{basename}_1f_representatives.results2 || true
    '''
}


process chimera_detection_post_cleave {
    // [S37]: re-run uchime_denovo on cat(pre-cleave, cleaved). The
    // --minsize threshold drops to the smallest size found in the
    // cleaved fas2 so newly cleaved low-abundance clusters are still
    // searched for chimeras, but it never goes below
    // params.chimera_minsize. An empty cleaved input falls back to
    // params.chimera_minsize.
    //
    // [S45]: the canonical <basename>_chimera_detection.log is the
    // concatenation of both chimera-detection runs' stderr (pre-cleave
    // chimera_detection followed by this post-cleave run). Each
    // fragment is preceded by a `=== <process_name> ===` section
    // header so the two runs remain distinguishable.
    //
    // [S59]: only the log reaches the results folder; the .uchime2
    // hit table is an internal intermediate consumed by
    // build_occurrence_table.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode, pattern: "*.log",
        enabled: params.results_folder != null

    input:
    path representatives           // pre-cleave: <basename>_1f_representatives.fas
    path cleaved_representatives   // cleaver:    <basename>_1f_representatives.fas2
    path pre_cleave_log            // chimera_detection's stderr (<basename>_1f_representatives.log)
    val basename

    output:
    path "${basename}_1f_representatives.uchime2"
    path "${basename}_chimera_detection.log"

    shell:
    '''
    #!/bin/bash
    set -euo pipefail

    lowest="$(sed -rn '/^>/ s/.*;size=([0-9]+);?/\\1/p' !{cleaved_representatives} \
                | sort -n | head -n 1)"
    lowest="${lowest:-0}"
    if (( lowest < !{params.chimera_minsize} )) ; then
        lowest=!{params.chimera_minsize}
    fi

    cat !{representatives} !{cleaved_representatives} | \
        vsearch \
            --sortbysize - \
            --sizein \
            --minsize "${lowest}" \
            --sizeout \
            --quiet \
            --output - | \
        vsearch \
            --uchime_denovo - \
            --uchimeout !{basename}_1f_representatives.uchime2 \
            2> chimera_detection_post_cleave.log

    {
        echo "=== chimera_detection ==="
        cat !{pre_cleave_log}
        echo "=== chimera_detection_post_cleave ==="
        cat chimera_detection_post_cleave.log
    } > !{basename}_chimera_detection.log
    '''
}


process build_occurrence_table {
    // [S35]: merge swarm / uchime / quality / stampa / distribution
    // into the filtered occurrence table. Representatives, stats,
    // swarms, chimera annotations, and taxonomic assignments are
    // each concatenated across the pre-cleave outputs and the
    // post-cleave ([S22]/[S36]/[S37]) outputs so cleaved clusters
    // get a row.
    //
    // Concat order matters for the taxonomic assignments: stampa_parse
    // overwrites by amplicon ID, so the file listed *last* on stdin
    // wins for an overlapping seed. Bash uses ``${TAX}{2,}`` (results2
    // first, results last) so the pre-cleave assignment wins on
    // overlap — that ordering is preserved here.
    //
    // [S46]: this is an intermediate OTU table — it is **not**
    // published; only the final <basename>_table.tsv reaches the
    // results folder.

    input:
    path representatives,   stageAs: 'reps_global.fas'
    path representatives_2, stageAs: 'reps_cleaved.fas2'
    path stats,             stageAs: 'stats_global'
    path stats_2,           stageAs: 'stats_cleaved'
    path swarms,            stageAs: 'swarms_global'
    path swarms_2,          stageAs: 'swarms_cleaved'
    path uchime,            stageAs: 'uchime_global'
    path uchime_2,          stageAs: 'uchime_cleaved'
    path quality
    path assignments,       stageAs: 'assignments_global'
    path assignments_2,     stageAs: 'assignments_cleaved'
    path distribution
    val basename
    val sample_ids   // [S09]: comma-separated sample IDs; empty samples
                     // contribute no .distr rows but still need a zero
                     // column in the occurrence table.

    output:
    path "${basename}.OTU.filtered.cleaved.table"

    shell:
    '''
    build_filtered_contingency_table.py \
        --representatives <(cat !{representatives} !{representatives_2}) \
        --stats           <(cat !{stats} !{stats_2}) \
        --swarms          <(cat !{swarms} !{swarms_2}) \
        --chimera         <(cat !{uchime} !{uchime_2}) \
        --quality         !{quality} \
        --assignments     <(cat !{assignments_2} !{assignments}) \
        --distribution    !{distribution} \
        --samples         '!{sample_ids}' \
        > !{basename}.OTU.filtered.cleaved.table
    '''
}


process cleaving {
    // [S22]: re-cleave global swarm clusters along sub-seed
    // boundaries. The script does all the work — this process is the
    // nextflow wrapper around bin/cluster_cleaver.py. The output
    // names follow the legacy `<input>2` / `<basename>_<sfx>_representatives.fas2`
    // convention, where `<sfx>` is `1f` (default, --fastidious) or
    // `1` (--no-fastidious) per [S22]'s propagation lever.
    //
    // [S45]: bin/cluster_cleaver.py uses python's logging module to
    // emit INFO-level progress to stderr; the redirect captures that
    // as the canonical cleaving log.
    //
    // [S59]: only the log reaches the results folder; the .stats2 /
    // .swarms2 / _representatives.fas2 cleaver outputs are internal
    // intermediates consumed by build_occurrence_table.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode, pattern: "*.log",
        enabled: params.results_folder != null

    input:
    path global_stats        // <basename>_<sfx>.stats
    path struct              // <basename>_<sfx>.struct
    path swarms              // <basename>_<sfx>.swarms
    path global_fasta        // <basename>.fas
    path per_sample_stats    // <basename>_per_sample_OTUs.stats
    val basename

    output:
    path "${basename}_${params.fastidious ? '1f' : '1'}.stats2"
    path "${basename}_${params.fastidious ? '1f' : '1'}.swarms2"
    path "${basename}_${params.fastidious ? '1f' : '1'}_representatives.fas2"
    path "${basename}_cleaving.log"

    shell:
    def fastidious_flag = params.fastidious ? '--fastidious' : '--no-fastidious'
    """
    cluster_cleaver.py \\
        --global_stats ${global_stats} \\
        --per_sample_stats ${per_sample_stats} \\
        --struct ${struct} \\
        --swarms ${swarms} \\
        --fasta ${global_fasta} \\
        ${fastidious_flag} \\
        2> ${basename}_cleaving.log
    """
}


process search_for_terminal_gaps {
    // [S38]: self-cluster the OTU table at id=1.0 to find OTU pairs
    // that are identical modulo terminal gaps (sub-strings or
    // super-strings). The output is the `H` lines of vsearch's
    // `--uc` stream — `merge_substring_otus` consumes those to
    // collapse pupil OTUs onto their masters.
    //
    // grep "^H" returns exit 1 when there are no hits; `|| true`
    // keeps the process green (an empty .uc is a legitimate outcome).
    //
    // [S45]: vsearch's --log captures the search step's stderr; the
    // merge_substring_otus process cats it together with the merge
    // step's stderr to produce <basename>_superstring_clustering.log.

    input:
    path otu_table

    output:
    path "${otu_table.baseName}.uc"
    path "search.log"

    shell:
    '''
    awk 'NR > 1 {printf ">"$1"\\n"$10"\\n"}' !{otu_table} | \
        vsearch \
            --threads !{params.threads} \
            --cluster_smallmem - \
            --id 1.0 \
            --qmask none \
            --usersort \
            --log search.log \
            --uc - | \
        grep "^H" > !{otu_table.baseName}.uc || true
    '''
}


process merge_substring_otus {
    // [S39]: wraps merge_sub_superstring_OTUs_with_larger_OTUs.py
    // + the legacy bash `sort_occurrence_table` step
    // + the read-count invariant check. Pupil OTUs are merged into
    // their masters (samples summed, spread/total/cloud updated),
    // the resulting table is sorted by the OTU column, and we
    // assert that the global read count is conserved.
    //
    // [S45]: cat the upstream vsearch search.log with the merge
    // step's stderr to produce the combined
    // <basename>_superstring_clustering.log. The merged OTU table
    // itself is **not** published ([S46]).
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode,
        pattern: "*_superstring_clustering.log",
        enabled: params.results_folder != null

    input:
    path otu_table
    path matches
    path search_log
    val basename

    output:
    path "${otu_table.baseName}.nosubstringOTUs.table"
    path "${basename}_superstring_clustering.log"

    shell:
    '''
    #!/bin/bash
    set -euo pipefail

    out="!{otu_table.baseName}.nosubstringOTUs.table"
    log="!{basename}_superstring_clustering.log"
    tmp_table="$(mktemp)"
    merge_stderr="$(mktemp)"
    trap 'rm -f "${tmp_table}" "${merge_stderr}"' EXIT

    {
        echo "=== search_for_terminal_gaps (vsearch --cluster_smallmem) ==="
        cat !{search_log}
        echo
        echo "=== merge_substring_otus + invariant check ==="
    } > "${log}"

    {
        merge_substring_otus.py \
            -t !{otu_table} \
            -m !{matches} \
            -o "${tmp_table}"

        (head -n 1 "${tmp_table}"
         tail -n +2 "${tmp_table}" | sort -k1,1n) > "${out}"

        before="$(awk 'NR > 1 {t += $2} END {print t + 0}' !{otu_table})"
        after="$(awk 'NR > 1 {t += $2} END {print t + 0}' "${out}")"
        if (( before != after )) ; then
            echo "merge_substring_otus: read count changed (${before} -> ${after})" >&2
            exit 1
        fi
        echo "merge_substring_otus: read count conserved (${before})" >&2
    } 2> "${merge_stderr}"
    cat "${merge_stderr}" >> "${log}"
    '''
}


process extract_otu_fasta {
    // [S40]: emit a FASTA from an OTU table (every data row).
    // Header `<amplicon>;size=<total>;`; column 4 is the amplicon,
    // column 2 the total abundance, column 10 the sequence.
    //
    // [S59]: pre-mumu fasta is an internal intermediate consumed by
    // find_similar_sequences (mumu match list) — not published.
    // The post-mumu sibling `extract_mumu_fasta` produces the one
    // FASTA that reaches the results folder.

    input:
    path table

    output:
    path "${table.baseName}.fas"

    shell:
    '''
    awk 'NR > 1 {printf ">"$4";size="$2";\\n"$10"\\n"}' !{table} \
        > !{table.baseName}.fas
    '''
}


process extract_mumu_fasta {
    // [S40]: post-mumu sibling of `extract_otu_fasta` — same column
    // layout, but skips rows whose `total == 0`. After the
    // size=0 → 1 awk hotfix in `rebuild_post_mumu_table` ([S44]) no
    // row carries `$2 == 0` anymore, so the filter is a no-op
    // safety net retained for byte parity with the legacy bash.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode,
        enabled: params.results_folder != null

    input:
    path table

    output:
    path "${table.baseName}.fas"

    shell:
    '''
    awk 'NR > 1 && $2 != 0 {printf ">"$4";size="$2";\\n"$10"\\n"}' !{table} \
        > !{table.baseName}.fas
    '''
}


process trim_metadata_for_mumu {
    // [S41]: keep the amplicon column (col 4) and every sample
    // column (cols 14+). mumu's --otu_table consumes this shape.

    input:
    path table

    output:
    path "${table.baseName}_reduced.table"

    shell:
    '''
    cut -f 4,14- !{table} > !{table.baseName}_reduced.table
    '''
}


process find_similar_sequences {
    // [S42]: vsearch --usearch_global self-search; lulu-recommended
    // parameters. The legacy bash strips `;size=N;` from every
    // column with a sed pass — that's the format mumu accepts.

    input:
    path otu_fasta

    output:
    path "${otu_fasta.baseName}.match_list"

    shell:
    '''
    vsearch \
        --usearch_global !{otu_fasta} \
        --db !{otu_fasta} \
        --self \
        --threads !{params.threads} \
        --id 0.84 \
        --iddef 1 \
        --userfields query+target+id \
        --maxaccepts 0 \
        --query_cov 0.9 \
        --maxhits 10 \
        --quiet \
        --userout - | \
        sed -r 's/;size=[0-9]+;//g' > !{otu_fasta.baseName}.match_list
    '''
}


process run_mumu {
    // [S43]: mumu (>=1.1.1) post-clustering filter. Inputs are the
    // reduced OTU table (amplicon + sample cols) and the self-search
    // match list; outputs are the new OTU table and the analysis log.
    //
    // [S45]: the mumu --log output is the canonical post-clustering
    // curation log. The intermediate _raw_mumu.table is **not**
    // published ([S46]); the publishDir pattern keeps the log only.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode, pattern: "*.log",
        enabled: params.results_folder != null

    input:
    path reduced_table
    path match_list
    val basename

    output:
    path "${reduced_table.baseName.replaceFirst(/_reduced$/, '_raw_mumu')}.table"
    path "${basename}_post_clustering_curation.log"

    shell:
    def new_table = "${reduced_table.baseName.replaceFirst(/_reduced$/, '_raw_mumu')}.table"
    """
    mumu \\
        --threads ${params.threads} \\
        --otu_table ${reduced_table} \\
        --match_list ${match_list} \\
        --new_otu_table ${new_table} \\
        --log ${basename}_post_clustering_curation.log
    """
}


process rebuild_post_mumu_table {
    // [S44]: wraps rebuild_table_after_mumu.py + the legacy
    // size=0 → 1 awk hotfix (so downstream `vsearch --sizein`
    // consumers don't choke on a zero-abundance row).
    //
    // [S46]: emits the final occurrence table as
    // <basename>_table.tsv.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode,
        enabled: params.results_folder != null

    input:
    path mumu_table
    path old_table
    val basename

    output:
    path "${basename}_table.tsv"

    shell:
    """
    rebuild_table_after_mumu.py \\
        --mumu_table ${mumu_table} \\
        --old_table  ${old_table} | \\
        awk 'BEGIN {FS = OFS = "\\t"} {if (\$2 == 0) {\$2 = 1} ; print \$0}' \\
        > ${basename}_table.tsv
    """
}


workflow part_B {
    // Runs the five Part B pre-cleaving processes on already-collected
    // lists of per-sample fasta / qual / stats files. Shared by the
    // standalone (`workflow part_b`) and end-to-end paths so the
    // process wiring lives in one place.

    take:
    fasta_list   // value channel: List<Path>, one .fas per sample
    qual_list    // value channel: List<Path>, one .qual per sample (same order)
    stats_list   // value channel: List<Path>, one .stats per sample (same order)

    main:
    def basename = fasta_list.map { files ->
        "${params.project_name}_${files.size()}_samples"
    }

    // [S09]: empty samples write no rows to .distr, so the column list
    // derived inside build_filtered_contingency_table from the .distr
    // alone misses them. Hand the authoritative sample list down via
    // build_occurrence_table's --samples so zero columns appear.
    def sample_ids = fasta_list.map { files ->
        files.collect { it.baseName }.toSorted().join(",")
    }

    build_expected_error_file(qual_list, basename)
    build_distribution_file(fasta_list, basename)
    list_all_cluster_seeds_of_size_greater_than_2(stats_list, basename)
    global_dereplication(fasta_list, basename)
    global_clustering(global_dereplication.out[0], basename)

    // Placeholder taxonomy + chimera flag run in parallel — both
    // consume only the swarm representatives ([S33]/[S34]).
    fake_taxonomic_assignment(global_clustering.out[3], basename)
    chimera_detection(global_clustering.out[3], basename)

    // Re-cleave clusters using the per-sample stats ([S22]).
    cleaving(
        global_clustering.out[1],            // .stats
        global_clustering.out[2],            // .struct
        global_clustering.out[0],            // .swarms
        global_dereplication.out[0],         // global .fas
        list_all_cluster_seeds_of_size_greater_than_2.out[0],
        basename
    )

    // [S36]/[S37]: post-cleave taxonomy + chimera annotations.
    fake_taxonomic_assignment2(cleaving.out[2], basename)
    chimera_detection_post_cleave(
        global_clustering.out[3],            // pre-cleave reps
        cleaving.out[2],                     // cleaved reps (fas2)
        chimera_detection.out[1],            // pre-cleave stderr ([S45])
        basename
    )

    // [S35]: assemble the filtered occurrence table.
    build_occurrence_table(
        global_clustering.out[3],            // _1f_representatives.fas
        cleaving.out[2],                     // _1f_representatives.fas2
        global_clustering.out[1],            // _1f.stats
        cleaving.out[0],                     // _1f.stats2
        global_clustering.out[0],            // _1f.swarms
        cleaving.out[1],                     // _1f.swarms2
        chimera_detection.out[0],            // _1f_representatives.uchime
        chimera_detection_post_cleave.out[0],           // _1f_representatives.uchime2
        build_expected_error_file.out[0],    // .qual
        fake_taxonomic_assignment.out[0],    // _1f_representatives.results
        fake_taxonomic_assignment2.out[0],   // _1f_representatives.results2
        build_distribution_file.out[0],      // .distr
        basename,
        sample_ids,                          // [S09]
    )

    // [S38]/[S39]: collapse sub- and super-string OTUs. The combined
    // vsearch + python step emits <basename>_superstring_clustering.log
    // ([S45]).
    search_for_terminal_gaps(build_occurrence_table.out[0])
    merge_substring_otus(
        build_occurrence_table.out[0],
        search_for_terminal_gaps.out[0],   // .uc hits
        search_for_terminal_gaps.out[1],   // vsearch search.log
        basename,
    )

    // [S40]–[S44]: mumu (ex-lulu) post-clustering filter pass.
    extract_otu_fasta(merge_substring_otus.out[0])
    trim_metadata_for_mumu(merge_substring_otus.out[0])
    find_similar_sequences(extract_otu_fasta.out[0])
    run_mumu(trim_metadata_for_mumu.out[0], find_similar_sequences.out[0], basename)
    rebuild_post_mumu_table(run_mumu.out[0], merge_substring_otus.out[0], basename)
    extract_mumu_fasta(rebuild_post_mumu_table.out[0])

    emit:
    // [S46] final occurrence table — consumed by part_C when Part B
    // and Part C run end-to-end.
    table = rebuild_post_mumu_table.out[0]
}


workflow part_B_shadow {
    // [S56]: shadow Part B. Runs the **same processes** as part_B on
    // the Part A shadow outputs (<sampleId>_notmerged.{fas,qual,stats}
    // — see [S04]). The A-padding emitted by Part A is composed
    // entirely of A/C/G/T (see [S63]), so swarm accepts the sequences
    // as-is and no mask/restore wrapping is needed around
    // global_clustering. The basename carries a trailing _notmerged
    // token so every published artefact is distinguishable from the
    // regular Part B output.

    take:
    fasta_list   // value channel: List<Path>, one _notmerged.fas per sample
    qual_list    // value channel: List<Path>, one _notmerged.qual per sample
    stats_list   // value channel: List<Path>, one _notmerged.stats per sample

    main:
    def basename = fasta_list.map { files ->
        "${params.project_name}_${files.size()}_samples_notmerged"
    }

    // [S09]: same column-list contract as the regular Part B path —
    // empty shadow samples must still appear as zero columns.
    def sample_ids = fasta_list.map { files ->
        files.collect { it.baseName }.toSorted().join(",")
    }

    build_expected_error_file(qual_list, basename)
    build_distribution_file(fasta_list, basename)
    list_all_cluster_seeds_of_size_greater_than_2(stats_list, basename)
    global_dereplication(fasta_list, basename)

    global_clustering(global_dereplication.out[0], basename)
    def representatives = global_clustering.out[3]

    fake_taxonomic_assignment(representatives, basename)
    chimera_detection(representatives, basename)

    cleaving(
        global_clustering.out[1],            // .stats
        global_clustering.out[2],            // .struct
        global_clustering.out[0],            // .swarms
        global_dereplication.out[0],         // global .fas
        list_all_cluster_seeds_of_size_greater_than_2.out[0],
        basename
    )

    fake_taxonomic_assignment2(cleaving.out[2], basename)
    chimera_detection_post_cleave(
        representatives,
        cleaving.out[2],
        chimera_detection.out[1],
        basename
    )

    build_occurrence_table(
        representatives,
        cleaving.out[2],
        global_clustering.out[1],
        cleaving.out[0],
        global_clustering.out[0],
        cleaving.out[1],
        chimera_detection.out[0],
        chimera_detection_post_cleave.out[0],
        build_expected_error_file.out[0],
        fake_taxonomic_assignment.out[0],
        fake_taxonomic_assignment2.out[0],
        build_distribution_file.out[0],
        basename,
        sample_ids,                          // [S09]
    )

    search_for_terminal_gaps(build_occurrence_table.out[0])
    merge_substring_otus(
        build_occurrence_table.out[0],
        search_for_terminal_gaps.out[0],
        search_for_terminal_gaps.out[1],
        basename,
    )

    extract_otu_fasta(merge_substring_otus.out[0])
    trim_metadata_for_mumu(merge_substring_otus.out[0])
    find_similar_sequences(extract_otu_fasta.out[0])
    run_mumu(trim_metadata_for_mumu.out[0], find_similar_sequences.out[0], basename)
    rebuild_post_mumu_table(run_mumu.out[0], merge_substring_otus.out[0], basename)
    extract_mumu_fasta(rebuild_post_mumu_table.out[0])

    emit:
    // [S46] final shadow occurrence table — exposed so callers can
    // splice taxonomy onto the shadow path via a standalone Part C
    // invocation. End-to-end wiring runs Part C on the regular path
    // only; calling it twice from the same scope is forbidden in
    // DSL2.
    table = rebuild_post_mumu_table.out[0]
}


// ============================================================================
// Part C — taxonomic assignment (stampa re-implementation)
// ============================================================================
// Skeleton phase: the processes are wired but the test coverage is
// tagged `pending` because D04 still needs to land before the
// stand-alone CLI / output policy is finalised.

process extract_fasta_sequences_from_occurrence_table {
    // [S48]: extract a representatives FASTA from the occurrence
    // table. Column 4 = amplicon ID, column 2 = abundance, column
    // 10 = sequence (same layout as [S40]'s extract_otu_fasta).

    input:
    path occurrence_table

    output:
    path "${occurrence_table.baseName}_representatives.fas"

    shell:
    '''
    awk 'NR > 1 {printf ">"$4";size="$2";\\n"$10"\\n"}' !{occurrence_table} \
        > !{occurrence_table.baseName}_representatives.fas
    '''
}


process assign_taxonomy_stampa {
    // [S49]: per-chunk step of the stampa scatter-gather. Run
    // vsearch --usearch_global against the reference dataset on
    // one fasta chunk, then collapse the top hits per amplicon
    // with bin/stampa_merge.py. The workflow uses splitFasta
    // upstream to produce chunks and collectFile downstream to
    // gather + sort the slices into the final
    // <basename>_taxonomy_stampa.tsv. --notrunclabels keeps the
    // ";size=N;" annotation in vsearch's --userout, which
    // stampa_merge.py parses directly.

    input:
    path chunk
    path reference_dataset

    output:
    path "stampa_chunk.tsv"

    shell:
    '''
    vsearch \
        --usearch_global !{chunk} \
        --threads !{params.threads} \
        --db !{reference_dataset} \
        --dbmask none \
        --qmask none \
        --rowlen 0 \
        --notrunclabels \
        --userfields query+id!{params.iddef}+target \
        --maxaccepts 0 \
        --maxrejects !{params.stampa_maxrejects} \
        --top_hits_only \
        --output_no_hits \
        --id !{params.stampa_id} \
        --iddef !{params.iddef} \
        --userout hits \
        2> vsearch.log

    stampa_merge.py hits > stampa_chunk.tsv
    '''
}


process assign_taxonomy_sintax {
    // [S50]: shadow Part C runs `vsearch --sintax` against the same
    // --reference_dataset ([S47]) and reshapes the tabbed output to
    // the canonical assignments TSV
    // (amplicon\tabundance\tidentity\ttaxonomy\treferences).
    //
    // Column mapping from vsearch --sintax --tabbedout:
    //   col 1 → amplicon (after stripping ;size=N;)
    //   col 2 → taxonomy (bootstrap-annotated lineage, e.g.
    //           `d:Bacteria(1.00),p:Proteobacteria(0.85)`)
    //   col 3 (strand) and col 4 (cutoff-filtered lineage) are ignored.
    // `identity` and `references` stay at placeholder values
    // ([S33]/[S46]) — sintax does not produce a percent identity, and
    // there is no single "reference accession" to record.
    //
    // The process is also used by the regular path when
    // params.taxonomy_method == 'sintax' ([S61]); the published
    // filename embeds `basename`, which carries the `_notmerged` token
    // on the shadow path.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode,
        enabled: params.results_folder != null

    input:
    path representatives
    path reference_dataset
    val basename

    output:
    path "${basename}_taxonomy_sintax.tsv"
    path "${basename}_taxonomy.log"

    shell:
    '''
    vsearch \
        --sintax !{representatives} \
        --threads !{params.threads} \
        --db !{reference_dataset} \
        --dbmask none \
        --sintax_cutoff !{params.sintax_cutoff} \
        --tabbedout raw_sintax.tsv \
        --log !{basename}_taxonomy.log

    awk 'BEGIN { FS = OFS = "\t" }
         {
             split($1, parts, ";size=")
             amplicon = parts[1]
             gsub(/;.*$/, "", parts[2])
             abundance = parts[2] + 0
             taxonomy = ($2 == "") ? "NA" : $2
             print amplicon, abundance, "0.0", taxonomy, "NA"
         }' raw_sintax.tsv > !{basename}_taxonomy_sintax.tsv
    '''
}


process update_occurrence_table {
    // [S51]: splice the taxonomy assignments back onto Part B's
    // <basename>_table.tsv. The output is published as a sibling
    // file `<basename>_table_assigned.tsv` so Part B's unannotated
    // table is preserved alongside the annotated one (D04 sub-q2).
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode,
        enabled: params.results_folder != null

    input:
    path occurrence_table
    path assignments
    val basename

    output:
    path "${basename}_table_assigned.tsv"

    shell:
    '''
    update_occurrence_table.py \
        --occurrence_table !{occurrence_table} \
        --assignments !{assignments} \
        > !{basename}_table_assigned.tsv
    '''
}


workflow part_C {
    // Re-usable Part C wiring shared by the standalone workflow and
    // the Part B → Part C end-to-end path.

    take:
    occurrence_table   // value channel: Path to <basename>_table.tsv
    basename           // value channel: <project>_<N>_samples

    main:
    // [S47]/[S64]: pick the reference matching the assignment method.
    // The two formats are not interchangeable; stampa expects a
    // space-separated lineage in the header and sintax expects a
    // `;tax=...;` annotation. The relevant assert in the entry point
    // already guarantees the matching flag is set before we land here.
    def reference = ( params.taxonomy_method == 'sintax' )
        ? file(normalize_path(params.reference_dataset_sintax))
        : file(normalize_path(params.reference_dataset))

    // [S48]: extract a representatives FASTA from the occurrence
    // table. (The fasta-input branch promised by [S48] is blocked
    // on D04 and is not wired in this skeleton.)
    extract_fasta_sequences_from_occurrence_table(occurrence_table)

    if ( params.taxonomy_method == 'sintax' ) {
        // [S50]: shadow path.
        assign_taxonomy_sintax(
            extract_fasta_sequences_from_occurrence_table.out[0],
            reference,
            basename,
        )
        update_occurrence_table(
            occurrence_table,
            assign_taxonomy_sintax.out[0],
            basename,
        )
    } else {
        // [S49]: stampa scatter-gather. Split the representatives
        // fasta into chunks (or pass it through unchanged when
        // params.stampa_chunk_size == 0), fan each chunk through
        // assign_taxonomy_stampa in parallel, then concatenate +
        // sort the slices via collectFile.
        //
        // Each chunk's TSV is paired with the basename value
        // channel via `combine` so the collectFile closure can
        // build "${basename}_taxonomy_stampa.tsv" per-item
        // (basename can't be interpolated into collectFile's
        // `name:` String at workflow-build time).
        //
        // Note on the sort closure: collectFile sorts in JVM heap,
        // which is fine for nf-metabarcoding's representative
        // counts (typically < a few million OTUs). If a real run
        // ever hits the OOM line, swap this for an external
        // `sort_taxonomy` process running `LC_ALL=C sort -k2,2nr
        // -k1,1d` on an unsorted collectFile output. See Plan B
        // (2026-05-19) for the fallback writeup.
        def reps_ch = extract_fasta_sequences_from_occurrence_table.out[0]
        def chunks = (params.stampa_chunk_size > 0)
            ? reps_ch.splitFasta(by: params.stampa_chunk_size, file: true)
            : reps_ch

        def merged = assign_taxonomy_stampa(chunks, reference)
            .combine(basename)
            .collectFile(
                storeDir: normalize_path(params.results_folder),
                sort: { line ->
                    // -k2,2nr -k1,1d == abundance desc, amplicon asc.
                    // collectFile's sort closure takes one argument
                    // (a line) and must return a Comparable key.
                    // ArrayLists aren't Comparable on the JVM side,
                    // so we build a single String key: width-13 inverted
                    // abundance (so descending becomes ascending under
                    // lexicographic ordering), then a tab, then the
                    // amplicon name (already ascending).
                    def f = line.tokenize('\t')
                    String.format(
                        "%013d\t%s",
                        9999999999999L - (f[1] as Long),
                        f[0],
                    )
                },
            ) { chunk_tsv, bn ->
                ["${bn}_taxonomy_stampa.tsv", chunk_tsv.text]
            }

        // [S51] empty-input contract: when the occurrence table has no
        // rows, the extracted representatives fasta is empty,
        // splitFasta emits zero chunks, and `merged` never fires —
        // which would leave update_occurrence_table dangling. Fall
        // back to an empty assignments file so update_occurrence_table
        // still runs and publishes a header-only
        // <basename>_table_assigned.tsv.
        def empty_assignments = file("${workflow.workDir}/empty_assignments.tsv")
        empty_assignments.text = ""

        update_occurrence_table(
            occurrence_table,
            merged.ifEmpty(empty_assignments),
            basename,
        )
    }

    emit:
    // [S51] annotated occurrence table — sibling of Part B's
    // <basename>_table.tsv, named <basename>_table_assigned.tsv.
    table = update_occurrence_table.out[0]
}


workflow part_C_shadow {
    // [S50]: shadow Part C. Runs `vsearch --sintax` (always, no
    // toggle) on the shadow occurrence table produced by shadow
    // Part B ([S56]) or supplied via [S62]'s standalone-mode probe.
    // The basename carries the `_notmerged` token so every
    // published artefact is distinguishable from the regular Part C
    // output. Shares `extract_fasta_sequences_from_occurrence_table`
    // and `update_occurrence_table` with `part_C` — DSL2 allows
    // these processes to be invoked from both sub-workflows because
    // each invocation lives in its own workflow scope.
    //
    // [S50] short-circuit: when the upstream produced no `_notmerged`
    // samples, `part_B_shadow` still emits a header-only table
    // (`_0_samples_notmerged_table.tsv`). The filter on `populated`
    // drops that case so no shadow artefacts are published when there
    // is nothing to assign. Standalone mode reaches this workflow only
    // when the sibling file exists ([S62]), but the filter is also a
    // safety net there (empty sibling → no published shadow artefact).

    take:
    occurrence_table   // value channel: Path to <basename>_notmerged_table.tsv
    basename           // value channel: <project>_<N>_samples_notmerged

    main:
    // [S64]: shadow Part C consumes the sintax-formatted reference.
    // Entry points gate this workflow on params.reference_dataset_sintax
    // being non-null, so the file() call below is safe by the time the
    // workflow is invoked.
    def reference = file(normalize_path(params.reference_dataset_sintax))

    def populated = occurrence_table.filter { tbl ->
        // any non-empty line after the header → there is at least one
        // amplicon to assign.
        tbl.toFile().readLines().drop(1).any { it.trim() }
    }

    extract_fasta_sequences_from_occurrence_table(populated)

    assign_taxonomy_sintax(
        extract_fasta_sequences_from_occurrence_table.out[0],
        reference,
        basename,
    )

    update_occurrence_table(
        populated,
        assign_taxonomy_sintax.out[0],
        basename,
    )

    emit:
    table = update_occurrence_table.out[0]
}


workflow part_c {
    // Standalone Part C ([S47]/[S48]): the user provides
    // --occurrence_table (Part B's <basename>_table.tsv) and a
    // reference whose format matches --taxonomy_method ([S47]/[S64]).
    // The fasta-input branch is blocked on D04 and is not exposed by
    // this skeleton.
    if ( params.taxonomy_method == 'sintax' ) {
        assert params.reference_dataset_sintax :
            "--reference_dataset_sintax must be set when --taxonomy_method=sintax"
    } else {
        assert params.reference_dataset :
            "--reference_dataset must be set when --taxonomy_method=stampa"
    }
    assert params.occurrence_table :
        "--occurrence_table must be set when running Part C standalone " +
        "(fasta input is blocked on DECISIONS.md D04)"
    assert params.results_folder :
        "--results_folder must be set when running Part C"

    def results_dir = new File(normalize_path(params.results_folder).toString())
    results_dir.mkdirs()

    def table_path = file(normalize_path(params.occurrence_table))
    def derived_basename = table_path.baseName.replaceFirst(/_table$/, '')
    def basename_ch = Channel.value(derived_basename)
    def table_ch    = Channel.fromPath(normalize_path(params.occurrence_table))

    part_C(table_ch, basename_ch)

    // [S62]/[S64]: probe for a shadow sibling next to the input
    // table. The file lookup is done at workflow-build time (not
    // inside a channel operator) so the toggle is "the file exists on
    // disk now". When the sibling is present **and**
    // --reference_dataset_sintax is set, route the sibling through
    // part_C_shadow alongside the regular part_C invocation — DSL2
    // allows two distinct sub-workflows in the same scope. When
    // either condition is unmet, the shadow workflow is simply not
    // invoked.
    def shadow_table_path = file(
        "${table_path.parent}/${derived_basename}_notmerged_table.tsv"
    )
    if ( shadow_table_path.exists() && params.reference_dataset_sintax ) {
        def shadow_basename_ch = Channel.value("${derived_basename}_notmerged")
        def shadow_table_ch    = Channel.fromPath(shadow_table_path.toString())
        part_C_shadow(shadow_table_ch, shadow_basename_ch)
    }
}


workflow part_b {
    // Standalone Part B ([S25]/[S26]/[S27]/[S56]): discover the
    // per-sample .fas files under params.fasta_folder and route them
    // to two parallel workflows:
    //   - part_B        for samples without _notmerged suffix
    //   - part_B_shadow for samples with the _notmerged suffix
    // .qual and .stats are derived as sibling files of each .fas, so
    // the same per-sample file-naming convention applies to both paths.
    //
    // When --reference_dataset is also set, Part C is wired onto the
    // regular path so the table_assigned.tsv lands alongside the Part B
    // outputs. Shadow taxonomy ([S50]) is wired via the sibling
    // workflow `part_C_shadow` on `part_B_shadow.out.table` — DSL2
    // allows two distinct sub-workflows in the same scope.
    assert params.project_name :
        "--project_name must be set when --fasta_folder is set"
    assert params.results_folder :
        "--results_folder must be set when --fasta_folder is set"

    // [S26]: create the results folder (and any missing parents)
    // before anything publishes into it.
    def results_dir = new File(normalize_path(params.results_folder).toString())
    results_dir.mkdirs()

    discover_part_b_fasta()

    def regular_samples_ch = discover_part_b_fasta.out[0]
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], file(row[1])) }
    def shadow_samples_ch = discover_part_b_fasta.out[1]
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], file(row[1])) }

    def fasta_list = regular_samples_ch.map { _id, f -> f }.collect()
    def qual_list  = regular_samples_ch
        .map { id, f -> file("${f.parent}/${id}.qual") }
        .collect()
    def stats_list = regular_samples_ch
        .map { id, f -> file("${f.parent}/${id}.stats") }
        .collect()
    part_B(fasta_list, qual_list, stats_list)

    def s_fasta_list = shadow_samples_ch.map { _id, f -> f }.collect()
    def s_qual_list  = shadow_samples_ch
        .map { id, f -> file("${f.parent}/${id}.qual") }
        .collect()
    def s_stats_list = shadow_samples_ch
        .map { id, f -> file("${f.parent}/${id}.stats") }
        .collect()
    part_B_shadow(s_fasta_list, s_qual_list, s_stats_list)

    // [S47]/[S64]: invoke Part C when the user supplied the reference
    // matching the selected method. Shadow Part C ([S50]) additionally
    // requires --reference_dataset_sintax — see [S64].
    def regular_ref_present = ( params.taxonomy_method == 'sintax' )
        ? params.reference_dataset_sintax
        : params.reference_dataset
    if ( regular_ref_present ) {
        def basename_ch = fasta_list.map { files ->
            "${params.project_name}_${files.size()}_samples"
        }
        part_C(part_B.out.table, basename_ch)
    }
    if ( params.reference_dataset_sintax ) {
        // [S50]/[S64]: shadow taxonomy on part_B_shadow's occurrence
        // table. The header-only gate inside part_C_shadow suppresses
        // output when the shadow side had no samples upstream.
        def shadow_basename_ch = s_fasta_list.map { files ->
            "${params.project_name}_${files.size()}_samples_notmerged"
        }
        part_C_shadow(part_B_shadow.out.table, shadow_basename_ch)
    }
}


workflow part_A {
    // Part A end-to-end: fastq → per-sample dereplicated fasta + qual
    // + stats. Wraps the regular and shadow ([S04]) pipelines so every
    // process tagged with this part shows up as `part_A:<process>` in
    // the Nextflow tty output.
    //
    // Emits three (sampleId, path) channels — one per per-sample
    // artefact — so the main workflow can split them into regular vs
    // shadow streams when feeding part_B / part_B_shadow.

    main:
    // [S10]/[S11]/[S12]/[S21]: pattern-driven discovery. Pairs go
    // through merge_fastq_pairs; single-end files skip the merging
    // step.
    discover_inputs()

    def branched = discover_inputs.out
        .splitCsv(sep: '\t')
        .map { row ->
            def sample = row[0]
            def r1 = file(row[1])
            def r2 = (row.size() > 2 && row[2]) ? file(row[2]) : null
            tuple(sample, r1, r2)
        }
        .branch { _sample, _r1, r2 ->
            paired:   r2 != null
            unpaired: r2 == null
        }

    def paired_ch = branched.paired.map { sample, r1, r2 ->
        tuple(sample, [r1, r2])
    }
    def unpaired_ch = branched.unpaired.map { sample, r1, _r2 ->
        tuple(sample, r1)
    }

    // discover pairs and merge
    merge_fastq_pairs(paired_ch)

    // [S04] shadow pipeline: optionally strip the low-quality 3' tails
    // ([S24]) and then join unmerged R1/R2 with A-padding (length
    // params.join_padding_length, default 8 — see [S63]). The shadow
    // sampleId is `<sampleId>_notmerged`, so all subsequent processes
    // naturally publish artefacts at `<sampleId>_notmerged.*` without
    // touching the regular pipeline. Because the padding uses A/C/G/T
    // only, the shadow and regular paths share the same downstream
    // filters (--fastq_maxns 0) and the same swarm invocation —
    // no per-path branching is needed after this point.
    strip_reads(
        merge_fastq_pairs.out[0],
        merge_fastq_pairs.out[3],
        merge_fastq_pairs.out[4]
    )
    def shadow_id = strip_reads.out[0].map { it + "_notmerged" }
    join_notmerged(shadow_id, strip_reads.out[1], strip_reads.out[2])

    // Build a unified (id, fastq) stream for the rest of Part A.
    def regular_ch = merge_fastq_pairs.out[0]
        .merge(merge_fastq_pairs.out[1])
        .mix(unpaired_ch)
    def shadow_ch = join_notmerged.out[0]
        .merge(join_notmerged.out[1])

    def to_process = regular_ch.mix(shadow_ch).multiMap { id, f ->
        id:   id
        file: f
    }

    // trim primers (skipped when --no_trimming is set)
    def sampleId_ch
    def fastq_ch
    if ( params.no_trimming ) {
        sampleId_ch = to_process.id
        fastq_ch    = to_process.file
    } else {
        trim_primers(to_process.id, to_process.file)
        sampleId_ch = trim_primers.out[0]
        fastq_ch    = trim_primers.out[1]
    }

    // convert to fasta with SHA1 + ee, apply min-length / max-N=0
    // filters. The A-padded shadow reads carry no Ns, so they pass
    // the same max-N=0 threshold as the regular path.
    filter_and_convert_to_fasta(sampleId_ch, fastq_ch)

    // set aside EE values
    extract_expected_error_values(
        filter_and_convert_to_fasta.out[0], filter_and_convert_to_fasta.out[1]
    )

    // dereplicate
    dereplicate_fasta(
        filter_and_convert_to_fasta.out[0], filter_and_convert_to_fasta.out[1]
    )

    list_local_clusters(dereplicate_fasta.out[0], dereplicate_fasta.out[1])

    emit:
    // (sampleId, path) tuples — shadow samples carry a _notmerged
    // suffix so the main workflow can filter them into the shadow
    // stream before calling part_B / part_B_shadow.
    fasta = dereplicate_fasta.out[0].merge(dereplicate_fasta.out[1])
    qual  = extract_expected_error_values.out[0]
        .merge(extract_expected_error_values.out[1])
    stats = list_local_clusters.out[0].merge(list_local_clusters.out[1])
}


workflow {
    // [S57]: --help short-circuits before any required-param assert.
    // print() (not log.info) so the usage block lands on stdout —
    // the conventional channel for help output and what nf-test
    // captures in workflow.stdout.
    if ( params.help ) {
        print usage()
        return
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

    // [S60]: every path-typed param is read through `normalize_path()`
    // at its use site (file(), publishDir, etc.) — see the helper at the
    // top of this file. Nextflow 25's `params` map is read-only from
    // inside a workflow body, so a one-shot mutation here is silently
    // dropped; the wrap-at-use-site pattern is the workaround.

    // [S47]/[S48]: --occurrence_table switches the pipeline into
    // Part C standalone mode (Parts A and B do not run). The
    // mode also requires --reference_dataset and --results_folder.
    if ( params.occurrence_table ) {
        part_c()
        return
    }

    // [S27]: --fasta_folder switches the pipeline into Part B
    // standalone mode. Part A's --fastq_folder requirement is lifted
    // in that mode (Part A does not run). When --reference_dataset is
    // also set, `workflow part_b` chains Part C onto Part B's regular
    // output.
    if ( params.fasta_folder ) {
        part_b()
        return
    }

    // required parameters (no default — supply via CLI or project config)
    assert params.fastq_folder : "--fastq_folder must be set (no default)"

    // [S18]/[S20]: primers and --no_trimming are mutually exclusive
    if ( params.no_trimming ) {
        assert !params.forward_primer : "--forward_primer must be empty when --no_trimming is set"
        assert !params.reverse_primer : "--reverse_primer must be empty when --no_trimming is set"
    } else {
        assert params.forward_primer : "--forward_primer must be set (no default)"
        assert params.reverse_primer : "--reverse_primer must be set (no default)"
    }

    part_A()

    // ----- Part A → Part B handoff ([S27]/[S56]) -----
    // When --project_name is set, Part B runs on the Part A outputs
    // (without going through disk discovery). Regular and shadow
    // samples ([S04]) feed two parallel workflows:
    //   - part_B        on samples whose ID does NOT end in _notmerged
    //   - part_B_shadow on samples whose ID DOES end in _notmerged
    // The shadow workflow's published artefacts carry a _notmerged
    // token in their basename so they sit alongside the regular
    // ones in --results_folder without colliding. When
    // --reference_dataset is also set, Part C runs on the regular
    // Part B table and `part_C_shadow` ([S50]) runs on the shadow
    // table — DSL2 allows two distinct sub-workflows in the same
    // scope.
    if ( params.project_name ) {
        assert params.results_folder :
            "--results_folder must be set when --project_name is set"
        def results_dir = new File(normalize_path(params.results_folder).toString())
        results_dir.mkdirs()

        def regular_fasta = part_A.out.fasta.filter { id, _f -> !id.endsWith("_notmerged") }
        def regular_qual  = part_A.out.qual .filter { id, _f -> !id.endsWith("_notmerged") }
        def regular_stats = part_A.out.stats.filter { id, _f -> !id.endsWith("_notmerged") }

        def shadow_fasta = part_A.out.fasta.filter { id, _f -> id.endsWith("_notmerged") }
        def shadow_qual  = part_A.out.qual .filter { id, _f -> id.endsWith("_notmerged") }
        def shadow_stats = part_A.out.stats.filter { id, _f -> id.endsWith("_notmerged") }

        // join the three streams on sample ID so the lists stay
        // aligned even when processes complete out of order.
        def joined_b = regular_fasta
            .join(regular_qual)
            .join(regular_stats)
        def b_fasta = joined_b.map { _id, f, _q, _s -> f }.collect()
        def b_qual  = joined_b.map { _id, _f, q, _s -> q }.collect()
        def b_stats = joined_b.map { _id, _f, _q, s -> s }.collect()
        part_B(b_fasta, b_qual, b_stats)

        def joined_shadow = shadow_fasta
            .join(shadow_qual)
            .join(shadow_stats)
        def s_fasta = joined_shadow.map { _id, f, _q, _s -> f }.collect()
        def s_qual  = joined_shadow.map { _id, _f, q, _s -> q }.collect()
        def s_stats = joined_shadow.map { _id, _f, _q, s -> s }.collect()
        part_B_shadow(s_fasta, s_qual, s_stats)

        // [S47]/[S64]: same routing as `workflow part_b`. Run regular
        // Part C when the reference matching --taxonomy_method is set;
        // run shadow Part C ([S50]) when --reference_dataset_sintax is
        // set. The two gates are independent: a stampa-only run still
        // produces regular taxonomy; supplying just the sintax flag
        // skips the regular part_C and produces shadow taxonomy only.
        def regular_ref_present = ( params.taxonomy_method == 'sintax' )
            ? params.reference_dataset_sintax
            : params.reference_dataset
        if ( regular_ref_present ) {
            def basename_ch = b_fasta.map { files ->
                "${params.project_name}_${files.size()}_samples"
            }
            part_C(part_B.out.table, basename_ch)
        }
        if ( params.reference_dataset_sintax ) {
            // [S50]/[S64]: shadow taxonomy on part_B_shadow's table.
            // The header-only gate inside part_C_shadow suppresses
            // output when the shadow side had zero _notmerged samples.
            def shadow_basename_ch = s_fasta.map { files ->
                "${params.project_name}_${files.size()}_samples_notmerged"
            }
            part_C_shadow(part_B_shadow.out.table, shadow_basename_ch)
        }
    }
}
