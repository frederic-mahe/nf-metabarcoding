include { publish_dir } from '../functions.nf'


// refactoring:
// breakdown into two functions
// - clusterize_sample
// - list_per_sample_clusters

process list_local_clusters {
    // retain only clusters with more than 2 reads
    // (do not use the fastidious option here)
    publishDir path: { publish_dir('per_sample') }, mode: params.publish_mode

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
