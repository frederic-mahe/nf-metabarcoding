include { publish_dir } from '../functions.nf'


process extract_expected_error_values {
    // extract ee for future quality filtering (keep the lowest
    // observed expected error value for each unique sequence)
    publishDir path: { publish_dir('per_sample') }, mode: params.publish_mode

    input:
    tuple val(sampleId), path(filtered_fasta)
    val id_length

    output:
    tuple val(sampleId), path("${sampleId}.qual"), emit: qual

    shell:
    '''
    set -euo pipefail

    length_of_sequence_IDs=!{id_length}
    extract_ee.awk !{filtered_fasta} | \
        sort --key=3,3n --key=1,1d --key=2,2n | \
        uniq --check-chars=${length_of_sequence_IDs} > !{sampleId}.qual
    '''

    stub:
    """
    touch ${sampleId}.qual
    """
}
