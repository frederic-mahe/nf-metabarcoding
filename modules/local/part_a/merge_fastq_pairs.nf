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
    # [S106] since vsearch 2.32.0, --fastq_qmax defaults to the highest
    # quality the offset can represent (93 at offset 33, 62 at offset
    # 64; vsearch requires offset + qmax <= 126, the last printable
    # ASCII), so input needs no flag. --fastq_mergepairs still clamps
    # its output at 41 by default, so --fastq_qmaxout is raised to that
    # same ceiling and the merged region is not clamped to 41.
    readonly QMAXOUT=$((126 - OFFSET))

    vsearch \
        --fastq_mergepairs !{fastq_pair[0]} \
        --reverse !{fastq_pair[1]} \
        --threads !{task.cpus} \
        --fastq_ascii "${OFFSET}" \
        --fastq_qmaxout "${QMAXOUT}" \
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
