process strip_reads {
    // [S24] shadow pipeline only: trim params.stripright nucleotides
    // from the 3' end of each R1 and R2 not-merged read. Used to
    // discard the low-quality tails before --fastq_join so the join
    // padding is surrounded by higher-quality bases. vsearch
    // --fastq_stripright 0 is a valid no-op pass-through.

    input:
    tuple val(sampleId), path(notmerged_fwd), path(notmerged_rev)

    output:
    tuple val(sampleId), path("stripped_fwd"), path("stripped_rev"), emit: stripped

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

    stub:
    """
    touch stripped_fwd stripped_rev
    """
}
