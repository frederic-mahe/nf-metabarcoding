include { publish_dir; log_dir } from '../functions.nf'


// refactoring:
// breakdown into two functions
// - clusterize_sample
// - list_per_sample_clusters

process list_local_clusters {
    // retain only clusters with more than 2 reads
    // (do not use the fastidious option here)
    // [D16]: data (.stats) to per_sample/, log to logs/part_a/per_sample/.
    publishDir path: { publish_dir('per_sample') }, mode: params.publish_mode, pattern: "*.stats"
    publishDir path: { log_dir('part_a/per_sample') }, mode: params.publish_mode, pattern: "*.log"

    input:
    tuple val(sampleId), path(dereplicated_fasta)

    output:
    tuple val(sampleId), path("${sampleId}.stats"), emit: stats
    path "${sampleId}_clustering.log",              emit: log

    shell:
    '''
    #!/bin/bash
    set -euo pipefail

    swarm \
        --threads !{task.cpus} \
        --differences 1 \
        --usearch-abundance \
        --log !{sampleId}_clustering.log \
        --output-file /dev/null \
        --statistics-file - \
        !{dereplicated_fasta} | \
        filter_swarm_stats.awk > !{sampleId}.stats
    '''

    stub:
    """
    touch ${sampleId}.stats ${sampleId}_clustering.log
    """
}
