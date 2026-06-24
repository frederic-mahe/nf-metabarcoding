process filter_and_convert_to_fasta {
    // use the params.hash_function digest (SHA1 or MD5, see [S65]) as
    // sequence names, compute expected error
    // values (ee), and apply the minimum-length / max-N filter
    // (--fastq_maxns 0: every N drops the read). The shadow path
    // ([S04]) pads with A/C/G/T (see [S63]), so the same max-N=0
    // threshold serves both the regular and the shadow path.
    input:
    tuple val(sampleId), path(trimmed_fastq)
    val relabel_flag

    output:
    tuple val(sampleId), path("filtered_fasta"), emit: fasta

    shell:
    '''
    #!/bin/bash

    readonly MIN_LENGTH=!{params.fastq_minlen}

    vsearch \
        --fastq_filter !{trimmed_fastq} \
        --fastq_minlen "${MIN_LENGTH}" \
        --fastq_maxns 0 \
        !{relabel_flag} \
        --fastq_ascii !{params.fastq_encoding} \
        --quiet \
        --eeout \
        --lengthout \
        --fasta_width 0 \
        --fastaout - > filtered_fasta
    '''

    stub:
    """
    touch filtered_fasta
    """
}
