process chimera_detection {
    // [S34]: filter representatives down to abundance >= chimera_minsize
    // (default 2), then run vsearch --uchime_denovo. The .uchime
    // hit table can be empty when no chimeras are found; the stderr
    // log captures the run.
    //
    // [S45] keeps this pre-cleave .log internal — the canonical
    // <basename>_chimera_detection.log is published by
    // chimera_detection_post_cleave and contains the concatenation of both
    // runs' stderr.
    //
    // [S59]: the .uchime hit table is an internal intermediate
    // consumed by build_occurrence_table — not published.

    input:
    path representatives
    val basename

    output:
    path "${basename}_1f_representatives.uchime", emit: uchime
    path "${basename}_1f_representatives.log",    emit: log

    shell:
    '''
    set -euo pipefail

    vsearch \
        --fastx_filter !{representatives} \
        --minsize !{params.chimera_minsize} \
        --quiet \
        --fastaout - | \
    vsearch \
        --uchime_denovo - \
        --uchimeout !{basename}_1f_representatives.uchime \
        2> !{basename}_1f_representatives.log
    '''

    stub:
    """
    touch ${basename}_1f_representatives.uchime ${basename}_1f_representatives.log
    """
}
