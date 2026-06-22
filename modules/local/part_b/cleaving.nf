include { log_dir } from '../functions.nf'


process cleaving {
    // [S22]: re-cleave global swarm clusters along sub-seed
    // boundaries. The script does all the work — this process is the
    // nextflow wrapper around bin/cluster_cleaver.py. The output
    // names follow the legacy `<input>2` / `<basename>_<sfx>_representatives.fas2`
    // convention, where `<sfx>` is `1f` (default, --fastidious) or
    // `1` (--no-fastidious) per [S22]'s propagation lever.
    //
    // [S45]: bin/cluster_cleaver.py uses python's logging module to
    // emit INFO-level progress to stderr; the redirect captures that
    // as the canonical cleaving log.
    //
    // [S59]: only the log reaches the results folder; the .stats2 /
    // .swarms2 / _representatives.fas2 cleaver outputs are internal
    // intermediates consumed by build_occurrence_table. [D16]: logs go
    // to logs/part_b/.
    publishDir path: { log_dir('part_b') }, mode: params.publish_mode, pattern: "*.log"

    input:
    path global_stats        // <basename>_<sfx>.stats
    path struct              // <basename>_<sfx>.struct
    path swarms              // <basename>_<sfx>.swarms
    path global_fasta        // <basename>.fas
    path per_sample_stats    // <basename>_per_sample_OTUs.stats
    val basename

    output:
    path "${basename}_${params.fastidious ? '1f' : '1'}.stats2",               emit: stats
    path "${basename}_${params.fastidious ? '1f' : '1'}.swarms2",              emit: swarms
    path "${basename}_${params.fastidious ? '1f' : '1'}_representatives.fas2", emit: representatives
    path "${basename}_cleaving.log",                                          emit: log

    shell:
    def fastidious_flag = params.fastidious ? '--fastidious' : '--no-fastidious'
    """
    cluster_cleaver.py \\
        --global_stats ${global_stats} \\
        --per_sample_stats ${per_sample_stats} \\
        --struct ${struct} \\
        --swarms ${swarms} \\
        --fasta ${global_fasta} \\
        --percentage ${params.percentage} \\
        ${fastidious_flag} \\
        2> ${basename}_cleaving.log
    """

    stub:
    """
    touch ${basename}_${params.fastidious ? '1f' : '1'}.stats2 ${basename}_${params.fastidious ? '1f' : '1'}.swarms2 ${basename}_${params.fastidious ? '1f' : '1'}_representatives.fas2 ${basename}_cleaving.log
    """
}
