include { log_dir } from '../functions.nf'


process chimera_detection_post_cleave {
    // [S37]: re-run uchime_denovo on cat(pre-cleave, cleaved). The
    // --minsize threshold drops to the smallest size found in the
    // cleaved fas2 so newly cleaved low-abundance clusters are still
    // searched for chimeras, but it never goes below
    // params.chimera_minsize. An empty cleaved input falls back to
    // params.chimera_minsize.
    //
    // [S45]: the canonical <basename>_chimera_detection.log is the
    // concatenation of both chimera-detection runs' stderr (pre-cleave
    // chimera_detection followed by this post-cleave run). Each
    // fragment is preceded by a `=== <process_name> ===` section
    // header so the two runs remain distinguishable.
    //
    // [S59]: only the log reaches the results folder; the .uchime2
    // hit table is an internal intermediate consumed by
    // build_occurrence_table. [D15]: logs go to logs/occurrence_table/.
    publishDir path: { log_dir('occurrence_table') }, mode: params.publish_mode, pattern: "*.log"

    input:
    path representatives           // pre-cleave: <basename>_1f_representatives.fas
    path cleaved_representatives   // cleaver:    <basename>_1f_representatives.fas2
    path pre_cleave_log            // chimera_detection's stderr (<basename>_1f_representatives.log)
    val basename

    output:
    path "${basename}_1f_representatives.uchime2", emit: uchime
    path "${basename}_chimera_detection.log",      emit: log

    shell:
    '''
    #!/bin/bash
    set -euo pipefail

    # [S81]: `sed -n 1p` (not `head -n 1`) takes the smallest size. `head`
    # closes the pipe after one line, which SIGPIPEs `sort` (exit 141)
    # once the sorted output overflows the OS pipe buffer on a large
    # cleaved file — under `set -euo pipefail` that aborts the task.
    # `sed -n 1p` reads its input to completion, so the producer never
    # gets SIGPIPE.
    lowest="$(sed -rn '/^>/ s/.*;size=([0-9]+);?/\\1/p' !{cleaved_representatives} \
                | sort -n | sed -n '1p')"
    lowest="${lowest:-0}"
    if (( lowest < !{params.chimera_minsize} )) ; then
        lowest=!{params.chimera_minsize}
    fi

    cat !{representatives} !{cleaved_representatives} | \
        vsearch \
            --sortbysize - \
            --sizein \
            --minsize "${lowest}" \
            --sizeout \
            --quiet \
            --output - | \
        vsearch \
            --uchime_denovo - \
            --uchimeout !{basename}_1f_representatives.uchime2 \
            2> chimera_detection_post_cleave.log

    {
        echo "=== chimera_detection ==="
        cat !{pre_cleave_log}
        echo "=== chimera_detection_post_cleave ==="
        cat chimera_detection_post_cleave.log
    } > !{basename}_chimera_detection.log
    '''

    stub:
    """
    touch ${basename}_1f_representatives.uchime2 ${basename}_chimera_detection.log
    """
}
