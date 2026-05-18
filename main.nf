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


workflow {
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
}
