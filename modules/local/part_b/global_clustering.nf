include { log_dir } from '../functions.nf'


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
    // internal intermediates. [D16]: logs go to logs/part_b/.
    publishDir path: { log_dir('part_b') }, mode: params.publish_mode, pattern: "*.log"

    input:
    path global_fasta
    val basename

    output:
    path "${basename}_${params.fastidious ? '1f' : '1'}.swarms",              emit: swarms
    path "${basename}_${params.fastidious ? '1f' : '1'}.stats",               emit: stats
    path "${basename}_${params.fastidious ? '1f' : '1'}.struct",              emit: struct
    path "${basename}_${params.fastidious ? '1f' : '1'}_representatives.fas", emit: representatives
    path "${basename}_clustering.log",                                        emit: log

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

    stub:
    """
    touch ${basename}_${params.fastidious ? '1f' : '1'}.swarms ${basename}_${params.fastidious ? '1f' : '1'}.stats ${basename}_${params.fastidious ? '1f' : '1'}.struct ${basename}_${params.fastidious ? '1f' : '1'}_representatives.fas ${basename}_clustering.log
    """
}
