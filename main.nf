#!/usr/bin/env nextflow

// optional (project defaults)
params.fastq_pattern = "/*_1_{1,2}.fastq.gz"
params.fastq_encoding = 33
params.threads = 4

// forward_primer, reverse_primer, fastq_folder are required and have
// no default; the workflow asserts them at startup (see [S18]).


process merge_fastq_pairs {
    input:
    tuple val(sampleId), path(fastq_pair)

    output:
    val sampleId
    path "merged_fastq"

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
        --fastqout merged_fastq
    '''
}


process trim_primers {
    // search forward primer in both normal and revcomp: now all reads
    // are in the same orientation. Matching leftmost is the default.
    // max_n is a caller-supplied input so the same process can serve
    // merged reads (max_n=0) and the [S04] unmerged-pair path (max_n
    // = size of the N-join insert).
    input:
    val sampleId
    path merged_fastq
    val max_n

    output:
    val sampleId
    path "trimmed_fastq"

    shell:
    '''
    #!/bin/bash

    readonly MIN_LENGTH=32
    readonly ERROR_RATE=0.1

    reverse_primer_revcomp=$(reverse_complement.sh !{params.reverse_primer})

    MIN_F=$(( !{params.forward_primer.length()} * 2 / 3 ))  # match is >= 2/3 of primer length
    MIN_R=$(( !{params.reverse_primer.length()} * 2 / 3 ))
    cutadapt \
        --cores=!{params.threads} \
        --minimum-length "${MIN_LENGTH}" \
        --error-rate "${ERROR_RATE}" \
        --revcomp \
        --rename="{id}" \
        --front "!{params.forward_primer};rightmost" \
        --overlap "${MIN_F}" \
        --discard-untrimmed \
        !{merged_fastq} | \
        cutadapt \
            --cores=!{params.threads} \
            --minimum-length "${MIN_LENGTH}" \
            --error-rate "${ERROR_RATE}" \
            --adapter "${reverse_primer_revcomp}" \
            --overlap "${MIN_R}" \
            --discard-untrimmed \
            --max-n "!{max_n}" \
            - > trimmed_fastq
    '''
}


process convert_fastq_to_fasta {
    // use SHA1 values as sequence names,
    // compute expected error values (ee)
    input:
    val sampleId
    path trimmed_fastq

    output:
    val sampleId
    path "filtered_fasta"

    shell:
    '''
    #!/bin/bash

    vsearch \
        --fastq_filter !{trimmed_fastq} \
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
    publishDir params.fastq_folder

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
    publishDir params.fastq_folder

    input:
    val sampleId
    path filtered_fasta

    output:
    val sampleId
    path "${sampleId}.fas"

    shell:
    '''
    vsearch \
        --derep_fulllength !{filtered_fasta} \
        --sizeout \
        --quiet \
        --fasta_width 0 \
        --xee \
        --xlength \
        --output - > !{sampleId}.fas
    '''
}


// refactoring:
// breakdown into two functions
// - clusterize_sample
// - list_per_sample_clusters

process list_local_clusters {
    // retain only clusters with more than 2 reads
    // (do not use the fastidious option here)
    publishDir params.fastq_folder

    input:
    val sampleId
    path dereplicated_fasta

    output:
    val sampleId
    path "${sampleId}.stats"

    shell:
    '''
    #!/bin/bash

    swarm \
        --threads !{params.threads} \
        --differences 1 \
        --usearch-abundance \
        --log /dev/null \
        --output-file /dev/null \
        --statistics-file - \
        !{dereplicated_fasta} | \
        filter_swarm_stats.awk > !{sampleId}.stats
    '''
}


workflow {
    // required parameters (no default — supply via CLI or project config)
    assert params.forward_primer : "--forward_primer must be set (no default)"
    assert params.reverse_primer : "--reverse_primer must be set (no default)"
    assert params.fastq_folder   : "--fastq_folder must be set (no default)"

    // discover pairs and merge
    merge_fastq_pairs(channel.fromFilePairs(params.fastq_folder + params.fastq_pattern))

    // trim primers (max_n=0 for merged reads; the [S04] unmerged-pair
    // path will pass the N-join insert size when implemented)
    trim_primers(merge_fastq_pairs.out[0], merge_fastq_pairs.out[1], 0)

    // convert to fasta with SHA1 + ee, then fan out
    ch_filtered_fasta = trim_primers.out | convert_fastq_to_fasta

    // set aside EE values
    ch_filtered_fasta |
        extract_expected_error_values

    // dereplicate and clusterize
    ch_filtered_fasta |
        dereplicate_fasta |
        list_local_clusters
}
