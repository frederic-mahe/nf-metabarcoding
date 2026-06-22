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
    tuple val(sampleId), path(notmerged_fwd), path(notmerged_rev)

    output:
    tuple val(sampleId), path("joined_fastq"), emit: joined

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

    stub:
    """
    touch joined_fastq
    """
}
