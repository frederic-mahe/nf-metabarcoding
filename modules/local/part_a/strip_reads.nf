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

    readonly OFFSET=!{params.fastq_encoding}
    # [S106] since vsearch 2.32.0, --fastq_qmax defaults to the highest
    # quality the offset can represent (93 at offset 33, 62 at offset
    # 64), so the full range is accepted without a flag.

    vsearch \
        --fastx_filter !{notmerged_fwd} \
        --fastq_ascii "${OFFSET}" \
        --fastq_stripright !{params.stripright} \
        --quiet \
        --fastqout stripped_fwd

    vsearch \
        --fastx_filter !{notmerged_rev} \
        --fastq_ascii "${OFFSET}" \
        --fastq_stripright !{params.stripright} \
        --quiet \
        --fastqout stripped_rev
    '''

    stub:
    """
    touch stripped_fwd stripped_rev
    """
}
