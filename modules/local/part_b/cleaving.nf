include { normalize_path } from '../functions.nf'


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
    // intermediates consumed by build_occurrence_table.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode, pattern: "*.log",
        enabled: params.results_folder != null

    input:
    path global_stats        // <basename>_<sfx>.stats
    path struct              // <basename>_<sfx>.struct
    path swarms              // <basename>_<sfx>.swarms
    path global_fasta        // <basename>.fas
    path per_sample_stats    // <basename>_per_sample_OTUs.stats
    val basename

    output:
    path "${basename}_${params.fastidious ? '1f' : '1'}.stats2"
    path "${basename}_${params.fastidious ? '1f' : '1'}.swarms2"
    path "${basename}_${params.fastidious ? '1f' : '1'}_representatives.fas2"
    path "${basename}_cleaving.log"

    shell:
    def fastidious_flag = params.fastidious ? '--fastidious' : '--no-fastidious'
    """
    cluster_cleaver.py \\
        --global_stats ${global_stats} \\
        --per_sample_stats ${per_sample_stats} \\
        --struct ${struct} \\
        --swarms ${swarms} \\
        --fasta ${global_fasta} \\
        ${fastidious_flag} \\
        2> ${basename}_cleaving.log
    """
}
