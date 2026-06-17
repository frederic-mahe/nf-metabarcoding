include { normalize_path } from '../functions.nf'


process global_clustering {
    // [S32]: swarm on the globally-dereplicated fasta. Output
    // filenames carry the swarm-parameters suffix `_1f`
    // (--fastidious, default) or `_1` (--no-fastidious), driven by
    // params.fastidious. See [S22]'s propagation lever. The swarm
    // log lands at <basename>_clustering.log ([S45]) — suffix-
    // independent.
    //
    // [S59]: only the log reaches the results folder; the
    // .swarms / .stats / .struct / _representatives.fas are
    // internal intermediates.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode, pattern: "*.log",
        enabled: params.results_folder != null

    input:
    path global_fasta
    val basename

    output:
    path "${basename}_${params.fastidious ? '1f' : '1'}.swarms"
    path "${basename}_${params.fastidious ? '1f' : '1'}.stats"
    path "${basename}_${params.fastidious ? '1f' : '1'}.struct"
    path "${basename}_${params.fastidious ? '1f' : '1'}_representatives.fas"
    path "${basename}_clustering.log"

    shell:
    def sfx = params.fastidious ? '1f' : '1'
    def fastidious_flag = params.fastidious ? '--fastidious' : ''
    """
    swarm \\
        --threads ${task.cpus} \\
        --differences 1 \\
        ${fastidious_flag} \\
        --usearch-abundance \\
        --internal-structure ${basename}_${sfx}.struct \\
        --output-file ${basename}_${sfx}.swarms \\
        --statistics-file ${basename}_${sfx}.stats \\
        --seeds ${basename}_${sfx}_representatives.fas \\
        ${global_fasta} 2> ${basename}_clustering.log
    """
}
