include { normalize_path } from '../functions.nf'


process dereplicate_fasta {
    // dereplicate and discard expected error values (ee)
    publishDir path: { normalize_path(params.fastq_folder) }, mode: params.publish_mode

    input:
    val sampleId
    path filtered_fasta

    output:
    val sampleId
    path "${sampleId}.fas"
    path "${sampleId}_dereplicating.log"

    shell:
    '''
    vsearch \
        --fastx_uniques !{filtered_fasta} \
        --sizeout \
        --quiet \
        --log !{sampleId}_dereplicating.log \
        --fasta_width 0 \
        --xee \
        --xlength \
        --fastaout - > !{sampleId}.fas
    '''
}
