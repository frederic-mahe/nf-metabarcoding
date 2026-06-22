include { log_dir } from '../functions.nf'


process trim_primers {
    // search forward primer in both normal and revcomp: now all reads
    // are in the same orientation. Matching leftmost is the default.
    // Length and N-count filtering are delegated to
    // filter_and_convert_to_fasta (vsearch --fastq_minlen / --fastq_maxns).
    publishDir path: { log_dir('part_a/per_sample') }, mode: params.publish_mode, pattern: "*.log"

    input:
    tuple val(sampleId), path(merged_fastq)

    output:
    tuple val(sampleId), path("trimmed_fastq"),     emit: trimmed
    path "${sampleId}_trimming_forward.log",        emit: log_forward
    path "${sampleId}_trimming_reverse.log",        emit: log_reverse

    shell:
    '''
    #!/bin/bash
    set -euo pipefail

    readonly ERROR_RATE=0.1

    reverse_primer_revcomp=$(reverse_complement.sh "!{params.reverse_primer}")

    MIN_F=$(( !{params.forward_primer.length()} * 2 / 3 ))  # match is >= 2/3 of primer length
    MIN_R=$(( !{params.reverse_primer.length()} * 2 / 3 ))
    # The two passes run concurrently through the pipe, so naming
    # task.cpus on each would request twice the reservation. Split the
    # budget instead: the forward pass searches both orientations
    # (--revcomp), roughly double the reverse pass' matching work, so it
    # gets ~2/3 of the cores. Both are floored at 1 (cutadapt treats
    # --cores=0 as "use every core on the machine").
    readonly FWD_CORES=$(( !{task.cpus} > 2 ? !{task.cpus} * 2 / 3 : 1 ))
    readonly REV_CORES=$(( !{task.cpus} - FWD_CORES > 0 ? !{task.cpus} - FWD_CORES : 1 ))
    # Two cutadapt passes (forward primer, then reverse primer) piped
    # together. Each pass writes its own report so the per-primer
    # trimming statistics stay separable ([S19]).
    cutadapt \
        --cores="${FWD_CORES}" \
        --error-rate "${ERROR_RATE}" \
        --revcomp \
        --rename="{id}" \
        --front "!{params.forward_primer};rightmost" \
        --overlap "${MIN_F}" \
        --discard-untrimmed \
        !{merged_fastq} 2> !{sampleId}_trimming_forward.log | \
        cutadapt \
            --cores="${REV_CORES}" \
            --error-rate "${ERROR_RATE}" \
            --adapter "${reverse_primer_revcomp}" \
            --overlap "${MIN_R}" \
            --discard-untrimmed \
            - > trimmed_fastq 2> !{sampleId}_trimming_reverse.log
    '''

    stub:
    """
    touch trimmed_fastq ${sampleId}_trimming_forward.log ${sampleId}_trimming_reverse.log
    """
}
