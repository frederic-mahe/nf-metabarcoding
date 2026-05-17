#!/usr/bin/env nextflow

// optional (project defaults)
params.fastq_pattern = "/*_1_{1,2}.fastq.gz"
params.fastq_encoding = 33
params.threads = 4
params.no_trimming = false

// forward_primer, reverse_primer, fastq_folder are required and have
// no default; the workflow asserts them at startup (see [S18]).
// [S20]: when params.no_trimming is true, forward_primer and
// reverse_primer must be empty and the trim_primers step is skipped.


process merge_fastq_pairs {
    publishDir path: { params.fastq_folder }, mode: 'link', pattern: "*.log",
        enabled: params.fastq_folder != null

    input:
    tuple val(sampleId), path(fastq_pair)

    output:
    val sampleId
    path "merged_fastq"
    path "${sampleId}_merging.log"

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
        --fastqout merged_fastq
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

    // [S21]: collect every fastq file in fastq_folder. Files that
    // match --fastq_pattern form pairs (merged by merge_fastq_pairs);
    // anything else is processed as an unpaired single-end sample
    // that skips the merging step.
    def paired_ch = channel.fromFilePairs(params.fastq_folder + params.fastq_pattern)
    def paired_paths = paired_ch
        .flatMap { id, pair -> pair }
        .collect()
        .ifEmpty([])
        .map { it as Set }

    def unpaired_ch = channel
        .fromPath(params.fastq_folder + "/*.{fastq,fq}{,.gz,.bz2}")
        .combine(paired_paths)
        .filter { p, paired -> !paired.contains(p) }
        .map { p, paired ->
            def name = p.getFileName().toString()
            def sampleId = name.replaceFirst(/\.(fastq|fq)(\.(gz|bz2))?$/, '')
            tuple(sampleId, p)
        }

    // discover pairs and merge
    merge_fastq_pairs(paired_ch)

    // re-pair merge outputs into tuples, mix with unpaired files,
    // then split back into two synchronised channels for downstream
    def to_process = merge_fastq_pairs.out[0]
        .merge(merge_fastq_pairs.out[1])
        .mix(unpaired_ch)
        .multiMap { id, f ->
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

    // convert to fasta with SHA1 + ee, apply min-length / max-N
    // filters (max_n=0 for merged reads; the [S04] unmerged-pair
    // path will pass the N-join insert size when implemented)
    filter_and_convert_to_fasta(sampleId_ch, fastq_ch, 0)

    // set aside EE values
    extract_expected_error_values(
        filter_and_convert_to_fasta.out[0], filter_and_convert_to_fasta.out[1]
    )

    // dereplicate and clusterize
    dereplicate_fasta(
        filter_and_convert_to_fasta.out[0], filter_and_convert_to_fasta.out[1]
    )
    list_local_clusters(dereplicate_fasta.out[0], dereplicate_fasta.out[1])
}
