#!/usr/bin/env nextflow

// optional (project defaults)
params.fastq_pattern = "/*_1_{1,2}.fastq.gz"
params.fastq_encoding = 33
params.threads = 4
params.no_trimming = false
// [S24] number of nucleotides stripped from the 3' end of each
// not-merged R1/R2 read before they are joined (shadow pipeline).
// Set to 0 to disable.
params.stripright = 30

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


process merge_fastq_pairs {
    // --fastqout_notmerged_fwd/_rev capture reads that fail to merge;
    // they feed the shadow Part A pipeline ([S04]). Fwd and rev are
    // kept in sync by vsearch.
    publishDir path: { params.fastq_folder }, mode: 'link', pattern: "*.log",
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
    // that failed --fastq_mergepairs with the default 8-N padding so
    // they can be processed by the rest of Part A as a single fastq.
    // The sampleId already carries the `_notmerged` suffix, so the
    // published log lands at <sampleId>_notmerged_merging.log.
    publishDir path: { params.fastq_folder }, mode: 'link', pattern: "*.log",
        enabled: params.fastq_folder != null

    input:
    val sampleId
    path notmerged_fwd
    path notmerged_rev

    output:
    val sampleId
    path "joined_fastq"
    path "${sampleId}_merging.log"

    shell:
    '''
    #!/bin/bash

    vsearch \
        --fastq_join !{notmerged_fwd} \
        --reverse !{notmerged_rev} \
        --fastq_ascii !{params.fastq_encoding} \
        --quiet \
        --log !{sampleId}_merging.log \
        --fastqout joined_fastq
    '''
}


process mask_ns_for_swarm {
    // Shadow pipeline ([S04]) only: swarm rejects Ns, so every N in
    // sequence lines is rewritten to A just before clustering. Header
    // lines are left untouched so the SHA1 IDs computed by
    // filter_and_convert_to_fasta stay consistent across .fas, .qual,
    // and .stats. The masked fasta is intentionally NOT published.

    input:
    val sampleId
    path fasta

    output:
    val sampleId
    path "masked_fasta"

    shell:
    '''
    sed '/^>/!s/N/A/g' !{fasta} > masked_fasta
    '''
}


process trim_primers {
    // search forward primer in both normal and revcomp: now all reads
    // are in the same orientation. Matching leftmost is the default.
    // Length and N-count filtering are delegated to
    // filter_and_convert_to_fasta (vsearch --fastq_minlen / --fastq_maxns).
    publishDir path: { params.fastq_folder }, mode: 'link', pattern: "*.log",
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
    // use SHA1 values as sequence names, compute expected error
    // values (ee), and apply the minimum-length / max-N filters.
    // max_n is a caller-supplied input so the same process can serve
    // merged reads (max_n=0) and the [S04] unmerged-pair path (max_n
    // = size of the N-join insert).
    input:
    val sampleId
    path trimmed_fastq
    val max_n

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
        --fastq_maxns !{max_n} \
        --relabel_sha1 \
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
    publishDir params.fastq_folder, mode: 'link'

    input:
    val sampleId
    path filtered_fasta

    output:
    val sampleId
    path "${sampleId}.qual"

    shell:
    '''
    length_of_sequence_IDs=40
    extract_ee.awk !{filtered_fasta} | \
        sort --key=3,3n --key=1,1d --key=2,2n | \
        uniq --check-chars=${length_of_sequence_IDs} > !{sampleId}.qual
    '''
}


process dereplicate_fasta {
    // dereplicate and discard expected error values (ee)
    publishDir params.fastq_folder, mode: 'link'

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
    publishDir params.fastq_folder, mode: 'link'

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
        ? params.fastq_folder
        : params.fastq_folder
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it }
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
    // [S27]: walk params.fasta_folder, drop *_notmerged.fas, assert
    // unique sample IDs. Emits TSV `sample_id<TAB>fasta_path`.

    output:
    path "fastas.tsv"

    script:
    def folders = (params.fasta_folder instanceof List)
        ? params.fasta_folder
        : params.fasta_folder
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it }
    def folder_args = folders.collect { "'${it}'" }.join(' ')
    """
    discover_fasta.py ${folder_args} > fastas.tsv
    """
}


process build_expected_error_file {
    // [S28]: merge every per-sample <sampleId>.qual into one
    // project-wide <basename>.qual. The input files are already sorted
    // by length / SHA1 / ee (see extract_expected_error_values), so
    // `sort --merge` is a straight k-way merge; uniq --check-chars=40
    // keeps the lowest-ee row per SHA1.
    publishDir params.results_folder, mode: 'link',
        enabled: params.results_folder != null

    input:
    path quals
    val basename

    output:
    path "${basename}.qual"

    shell:
    '''
    sort --key=3,3n --key=1,1d --key=2,2n --merge !{quals} | \
        uniq --check-chars=40 > !{basename}.qual
    '''
}


process build_distribution_file {
    // [S29]: scan FASTA headers of every input .fas and emit the
    // sequence-to-sample mapping as tab-separated rows
    // <sha1>\t<sampleId>\t<size>. The sample ID is derived from the
    // fasta basename so the channel can be built without an explicit
    // sample-ID side car.
    publishDir params.results_folder, mode: 'link',
        enabled: params.results_folder != null

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
    publishDir params.results_folder, mode: 'link',
        enabled: params.results_folder != null

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
    // samples.
    publishDir params.results_folder, mode: 'link',
        enabled: params.results_folder != null

    input:
    path fastas
    val basename

    output:
    path "${basename}.fas"
    path "${basename}.log"

    shell:
    '''
    cat !{fastas} | \
        vsearch \
            --fastx_uniques - \
            --sizein \
            --sizeout \
            --fasta_width 0 \
            --quiet \
            --log !{basename}.log \
            --fastaout !{basename}.fas
    '''
}


process global_clustering {
    // [S32]: swarm on the globally-dereplicated fasta. Output
    // filenames carry the `_1f` suffix (resolution=1, --fastidious)
    // to mirror the bash-reference naming scheme.
    publishDir params.results_folder, mode: 'link',
        enabled: params.results_folder != null

    input:
    path global_fasta
    val basename

    output:
    path "${basename}_1f.swarms"
    path "${basename}_1f.stats"
    path "${basename}_1f.struct"
    path "${basename}_1f_representatives.fas"
    path "${basename}_1f.log"

    shell:
    '''
    swarm \
        --threads !{params.threads} \
        --differences 1 \
        --fastidious \
        --usearch-abundance \
        --internal-structure !{basename}_1f.struct \
        --output-file !{basename}_1f.swarms \
        --statistics-file !{basename}_1f.stats \
        --seeds !{basename}_1f_representatives.fas \
        !{global_fasta} 2> !{basename}_1f.log
    '''
}


process fake_taxonomic_assignment {
    // [S33]: emit a placeholder TSV that the occurrence-table builder
    // can consume before Part C lands. Each row mirrors the stampa
    // output shape: <amplicon>\t<size>\t<identity>\t<taxonomy>\t<refs>.
    publishDir params.results_folder, mode: 'link',
        enabled: params.results_folder != null

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
    publishDir params.results_folder, mode: 'link',
        enabled: params.results_folder != null

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
    publishDir params.results_folder, mode: 'link',
        enabled: params.results_folder != null

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


process chimera_detection2 {
    // [S37]: re-run uchime_denovo on cat(pre-cleave, cleaved). The
    // --minsize threshold drops to the smallest size found in the
    // cleaved fas2 so newly cleaved low-abundance clusters are still
    // searched for chimeras, but it never goes below
    // params.chimera_minsize. An empty cleaved input falls back to
    // params.chimera_minsize.
    publishDir params.results_folder, mode: 'link',
        enabled: params.results_folder != null

    input:
    path representatives           // pre-cleave: <basename>_1f_representatives.fas
    path cleaved_representatives   // cleaver:    <basename>_1f_representatives.fas2
    val basename

    output:
    path "${basename}_1f_representatives.uchime2"
    path "${basename}_1f_representatives.log2"

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
            2> !{basename}_1f_representatives.log2
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
    publishDir params.results_folder, mode: 'link',
        enabled: params.results_folder != null

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
        > !{basename}.OTU.filtered.cleaved.table
    '''
}


process cleaving {
    // [S22]: re-cleave global swarm clusters along sub-seed
    // boundaries. The script does all the work — this process is the
    // nextflow wrapper around bin/cluster_cleaver.py. The output
    // names follow the legacy `<input>2` / `<basename>_1f_representatives.fas2`
    // convention so downstream concatenations (used by the
    // occurrence-table builder) keep working unchanged.
    publishDir params.results_folder, mode: 'link',
        enabled: params.results_folder != null

    input:
    path global_stats        // <basename>_1f.stats
    path struct              // <basename>_1f.struct
    path swarms              // <basename>_1f.swarms
    path global_fasta        // <basename>.fas
    path per_sample_stats    // <basename>_per_sample_OTUs.stats
    val basename

    output:
    path "${basename}_1f.stats2"
    path "${basename}_1f.swarms2"
    path "${basename}_1f_representatives.fas2"

    shell:
    '''
    cluster_cleaver.py \
        --global_stats !{global_stats} \
        --per_sample_stats !{per_sample_stats} \
        --struct !{struct} \
        --swarms !{swarms} \
        --fasta !{global_fasta} \
        --fastidious
    '''
}


workflow part_b_processes {
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
    chimera_detection2(
        global_clustering.out[3],            // pre-cleave reps
        cleaving.out[2],                     // cleaved reps (fas2)
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
        chimera_detection2.out[0],           // _1f_representatives.uchime2
        build_expected_error_file.out[0],    // .qual
        fake_taxonomic_assignment.out[0],    // _1f_representatives.results
        fake_taxonomic_assignment2.out[0],   // _1f_representatives.results2
        build_distribution_file.out[0],      // .distr
        basename,
    )
}


workflow part_b {
    // Standalone Part B ([S25]/[S26]/[S27]): discover per-sample .fas
    // files (skipping shadow-pipeline _notmerged artefacts), assert
    // sample-ID uniqueness, and run the pre-cleaving processes. .qual
    // and .stats are derived as sibling files of each .fas.
    assert params.project_name :
        "--project_name must be set when --fasta_folder is set"
    assert params.results_folder :
        "--results_folder must be set when --fasta_folder is set"

    // [S26]: create the results folder (and any missing parents)
    // before anything publishes into it.
    def results_dir = new File(params.results_folder.toString())
    results_dir.mkdirs()

    discover_part_b_fasta()

    def samples_ch = discover_part_b_fasta.out
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], file(row[1])) }

    def fasta_list = samples_ch.map { _id, f -> f }.collect()
    def qual_list  = samples_ch
        .map { id, f -> file("${f.parent}/${id}.qual") }
        .collect()
    def stats_list = samples_ch
        .map { id, f -> file("${f.parent}/${id}.stats") }
        .collect()

    part_b_processes(fasta_list, qual_list, stats_list)
}


workflow {
    // [S27]: --fasta_folder switches the pipeline into Part B
    // standalone mode. Part A's --fastq_folder requirement is lifted
    // in that mode (Part A does not run).
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
    // ([S24]) and then join unmerged R1/R2 with 8-N padding. The
    // shadow sampleId is `<sampleId>_notmerged`, so all subsequent
    // processes naturally publish artefacts at `<sampleId>_notmerged.*`
    // without touching the regular pipeline.
    strip_reads(
        merge_fastq_pairs.out[0],
        merge_fastq_pairs.out[3],
        merge_fastq_pairs.out[4]
    )
    def shadow_id = strip_reads.out[0].map { it + "_notmerged" }
    join_notmerged(shadow_id, strip_reads.out[1], strip_reads.out[2])

    // Build a unified (id, fastq, max_n) stream for the rest of Part A.
    //   regular path uses max_n=0; shadow path uses max_n=8 (the
    //   vsearch --fastq_join default padding length).
    def regular_ch = merge_fastq_pairs.out[0]
        .merge(merge_fastq_pairs.out[1])
        .mix(unpaired_ch)
        .map { id, f -> tuple(id, f, 0) }
    def shadow_ch = join_notmerged.out[0]
        .merge(join_notmerged.out[1])
        .map { id, f -> tuple(id, f, 8) }

    def to_process = regular_ch.mix(shadow_ch).multiMap { id, f, max_n ->
        id:    id
        file:  f
        max_n: max_n
    }

    // trim primers (skipped when --no_trimming is set)
    def sampleId_ch
    def fastq_ch
    def max_n_ch
    if ( params.no_trimming ) {
        sampleId_ch = to_process.id
        fastq_ch    = to_process.file
        max_n_ch    = to_process.max_n
    } else {
        trim_primers(to_process.id, to_process.file)
        // trim_primers may reorder items (parallel execution), so
        // re-attach max_n by joining on sampleId.
        def id_max_n = to_process.id.merge(to_process.max_n)
        def joined = trim_primers.out[0]
            .merge(trim_primers.out[1])
            .join(id_max_n)
            .multiMap { id, f, m ->
                id:    id
                file:  f
                max_n: m
            }
        sampleId_ch = joined.id
        fastq_ch    = joined.file
        max_n_ch    = joined.max_n
    }

    // convert to fasta with SHA1 + ee, apply min-length / max-N
    // filters (max_n=0 for the regular path, max_n=8 for the shadow
    // path so the join-padding Ns survive)
    filter_and_convert_to_fasta(sampleId_ch, fastq_ch, max_n_ch)

    // set aside EE values
    extract_expected_error_values(
        filter_and_convert_to_fasta.out[0], filter_and_convert_to_fasta.out[1]
    )

    // dereplicate
    dereplicate_fasta(
        filter_and_convert_to_fasta.out[0], filter_and_convert_to_fasta.out[1]
    )

    // [S04] swarm rejects Ns, so the shadow path's dereplicated fasta
    // gets a transient N->A pass before clustering. Regular-path items
    // (no _notmerged suffix) bypass the mask step.
    def derep_branched = dereplicate_fasta.out[0]
        .merge(dereplicate_fasta.out[1])
        .branch { id, _f ->
            shadow:  id.endsWith("_notmerged")
            regular: true
        }
    mask_ns_for_swarm(
        derep_branched.shadow.map { id, _f -> id },
        derep_branched.shadow.map { _id, f -> f }
    )
    def ready_for_swarm = derep_branched.regular
        .mix(mask_ns_for_swarm.out[0].merge(mask_ns_for_swarm.out[1]))
        .multiMap { id, f ->
            id:   id
            file: f
        }

    list_local_clusters(ready_for_swarm.id, ready_for_swarm.file)

    // ----- Part A → Part B handoff ([S27]) -----
    // When --project_name is set, Part B runs on the Part A outputs
    // (without going through disk discovery). The shadow path's
    // outputs are dropped from the Part B channel: [S04] specifies
    // they belong to a separate downstream path, not the regular
    // Part B fasta channel. They remain published to params.fastq_folder.
    if ( params.project_name ) {
        assert params.results_folder :
            "--results_folder must be set when --project_name is set"
        def results_dir = new File(params.results_folder.toString())
        results_dir.mkdirs()

        def regular_fasta = dereplicate_fasta.out[0]
            .merge(dereplicate_fasta.out[1])
            .filter { id, _f -> !id.endsWith("_notmerged") }
        def regular_qual = extract_expected_error_values.out[0]
            .merge(extract_expected_error_values.out[1])
            .filter { id, _f -> !id.endsWith("_notmerged") }
        def regular_stats = list_local_clusters.out[0]
            .merge(list_local_clusters.out[1])
            .filter { id, _f -> !id.endsWith("_notmerged") }

        // join the three streams on sample ID so the lists stay
        // aligned even when processes complete out of order.
        def joined_b = regular_fasta
            .join(regular_qual)
            .join(regular_stats)

        def b_fasta = joined_b.map { _id, f, _q, _s -> f }.collect()
        def b_qual  = joined_b.map { _id, _f, q, _s -> q }.collect()
        def b_stats = joined_b.map { _id, _f, _q, s -> s }.collect()

        part_b_processes(b_fasta, b_qual, b_stats)
    }
}
