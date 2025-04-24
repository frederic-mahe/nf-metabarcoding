#!/usr/bin/env nextflow

params.forward_primer = "CCAGCASCYGCGGTAATTCC"
params.reverse_primer = "ACTTTCGTTCTTGATYRA"  // should be TYRATCAAGAACGAAAGT
params.fastq_folder = "data"
params.fastq_pattern = "/*_1_{1,2}.fastq.gz"
params.fastq_encoding = 33
params.threads = 4


process generate_test_data_urls {
    output:
    path "test_data_urls.list"

    shell:
    '''
    #!/bin/bash

    URL="https://github.com/frederic-mahe/BIO9905MERG1_vsearch_swarm_pipeline/raw/main/data"

    (echo "${URL}/MD5SUM"
     for SAMPLE in {B,L}{010..100..10} ; do
         for READ in 1 2 ; do
             echo "${URL}/${SAMPLE}_1_${READ}.fastq.gz"
         done
     done
    ) > test_data_urls.list
    '''
}


process download_list_of_urls {
    input:
    path "urls"

    publishDir params.fastq_folder

    output:
    path "*.fastq.gz"

    shell:
    '''
    #!/bin/bash

    wget --continue --quiet --input-file="!{urls}"
    '''
}


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


// refactoring:
// extract function reverse_complement_sequence

process trim_primers {
    // search forward primer in both normal and revcomp: now all reads
    // are in the same orientation. Matching leftmost is the default.
    input:
    val sampleId
    path merged_fastq

    output:
    val sampleId
    path "trimmed_fastq"

    shell:
    '''
    #!/bin/bash

    nucleotides="acgturykmbdhvswACGTURYKMBDHVSW"
    complements="tgcaayrmkvhdbswTGCAAYRMKVHDBSW"
    reverse_primer_revcomp=$(tr "${nucleotides}" "${complements}" <<< !{params.reverse_primer} | rev)

    MIN_F=$(( !{params.forward_primer.length()} * 2 / 3 ))  # match is >= 2/3 of primer length
    MIN_R=$(( !{params.reverse_primer.length()} * 2 / 3 ))
    cutadapt \
        --revcomp \
        --front "!{params.forward_primer};rightmost" \
        --overlap "${MIN_F}" !{merged_fastq} | \
        cutadapt \
            --adapter "${reverse_primer_revcomp}" \
            --overlap "${MIN_R}" \
            --max-n 0 - > trimmed_fastq
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
        --fasta_width 0 \
        --fastaout - > filtered_fasta
    '''
}


process extract_expected_error_values {
    // extract ee for future quality filtering (keep the lowest
    // observed expected error value for each unique sequence)
    input:
    val sampleId
    path filtered_fasta

    publishDir params.fastq_folder

    output:
    val sampleId
    path "${sampleId}.qual"

    shell:
    '''
    length_of_sequence_IDs=40
    paste - - < !{filtered_fasta} | \
        awk 'BEGIN {FS = "[>;=\t]"} {print $2, $4, length($NF)}' | \
        sort --key=3,3n --key=1,1d --key=2,2n | \
        uniq --check-chars=${length_of_sequence_IDs} > !{sampleId}.qual
    '''
}


process dereplicate_fasta {
    // dereplicate and discard expected error values (ee)
    input:
    val sampleId
    path filtered_fasta

    publishDir params.fastq_folder

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
    input:
    val sampleId
    path dereplicated_fasta

    publishDir params.fastq_folder

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
        awk 'BEGIN {FS = OFS = "\t"} $2 > 2' > !{sampleId}.stats
    '''
}


workflow {
    // collect test data
    generate_test_data_urls |
        download_list_of_urls

    // merge, trim, convert
    ch_filtered_fasta = channel.fromFilePairs(params.fastq_folder + params.fastq_pattern) |
        merge_fastq_pairs |
        trim_primers |
        convert_fastq_to_fasta

    // set aside EE values
    extract_expected_error_values(ch_filtered_fasta)

    // dereplicate and clusterize
    dereplicate_fasta(ch_filtered_fasta) |
        list_local_clusters
}
