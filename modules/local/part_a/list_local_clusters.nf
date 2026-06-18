include { normalize_path } from '../functions.nf'


// refactoring:
// breakdown into two functions
// - clusterize_sample
// - list_per_sample_clusters

process list_local_clusters {
    // retain only clusters with more than 2 reads
    // (do not use the fastidious option here)
    publishDir path: { normalize_path(params.fastq_folder) }, mode: params.publish_mode,
        enabled: params.fastq_folder != null

    input:
    tuple val(sampleId), path(dereplicated_fasta)

    output:
    tuple val(sampleId), path("${sampleId}.stats"), emit: stats
    path "${sampleId}_clustering.log",              emit: log

    shell:
    '''
    #!/bin/bash

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
}
