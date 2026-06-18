include { normalize_path } from '../functions.nf'


process trim_primers {
    // search forward primer in both normal and revcomp: now all reads
    // are in the same orientation. Matching leftmost is the default.
    // Length and N-count filtering are delegated to
    // filter_and_convert_to_fasta (vsearch --fastq_minlen / --fastq_maxns).
    publishDir path: { normalize_path(params.fastq_folder) }, mode: params.publish_mode, pattern: "*.log",
        enabled: params.fastq_folder != null

    input:
    tuple val(sampleId), path(merged_fastq)

    output:
    tuple val(sampleId), path("trimmed_fastq"), emit: trimmed
    path "${sampleId}_trimming.log",            emit: log

    shell:
    '''
    #!/bin/bash

    readonly ERROR_RATE=0.1

    reverse_primer_revcomp=$(reverse_complement.sh !{params.reverse_primer})

    MIN_F=$(( !{params.forward_primer.length()} * 2 / 3 ))  # match is >= 2/3 of primer length
    MIN_R=$(( !{params.reverse_primer.length()} * 2 / 3 ))
    {
        cutadapt \
            --cores=!{task.cpus} \
            --error-rate "${ERROR_RATE}" \
            --revcomp \
            --rename="{id}" \
            --front "!{params.forward_primer};rightmost" \
            --overlap "${MIN_F}" \
            --discard-untrimmed \
            !{merged_fastq} | \
            cutadapt \
                --cores=!{task.cpus} \
                --error-rate "${ERROR_RATE}" \
                --adapter "${reverse_primer_revcomp}" \
                --overlap "${MIN_R}" \
                --discard-untrimmed \
                - > trimmed_fastq
    } 2> !{sampleId}_trimming.log
    '''
}
