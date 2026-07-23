include { log_dir } from '../functions.nf'


process merge_fastq_pairs {
    // --fastqout_notmerged_fwd/_rev capture reads that fail to merge;
    // they feed the shadow Part A pipeline ([S04]). Fwd and rev are
    // kept in sync by vsearch.
    publishDir path: { log_dir('part_a/per_sample') }, mode: params.publish_mode, pattern: "*.log"

    input:
    tuple val(sampleId), path(fastq_pair)

    output:
    tuple val(sampleId), path("merged_fastq"),                        emit: merged
    tuple val(sampleId), path("notmerged_fwd"), path("notmerged_rev"), emit: notmerged
    path "${sampleId}_merging.log",                                   emit: log

    shell:
    '''
    #!/bin/bash

    readonly OFFSET=!{params.fastq_encoding}
    # [S106] vsearch requires offset + qmax <= 126 (126 = last printable
    # ASCII); 126 - offset is the highest representable quality, so this
    # accepts the full range for either encoding (93 at offset 33, 62 at
    # offset 64). --fastq_qmaxout uses the same ceiling so the merged
    # region is written at up to that value rather than clamped to 41.
    readonly QMAX=$((126 - OFFSET))

    vsearch \
        --fastq_mergepairs !{fastq_pair[0]} \
        --reverse !{fastq_pair[1]} \
        --threads !{task.cpus} \
        --fastq_ascii "${OFFSET}" \
        --fastq_qmax "${QMAX}" \
        --fastq_qmaxout "${QMAX}" \
        --fastq_allowmergestagger \
        --quiet \
        --log !{sampleId}_merging.log \
        --fastqout merged_fastq \
        --fastqout_notmerged_fwd notmerged_fwd \
        --fastqout_notmerged_rev notmerged_rev
    '''

    stub:
    """
    touch merged_fastq notmerged_fwd notmerged_rev ${sampleId}_merging.log
    """
}
