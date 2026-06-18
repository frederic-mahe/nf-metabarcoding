process build_expected_error_file {
    // [S28]: merge every per-sample <sampleId>.qual into one
    // project-wide <basename>.qual. The input files are already sorted
    // by length / hash / ee (see extract_expected_error_values), so
    // `sort --merge` is a straight k-way merge; uniq --check-chars
    // (width from params.hash_function, see [S65]) keeps the lowest-ee
    // row per amplicon name.
    //
    // [S59]: the .qual is an internal intermediate consumed by
    // build_occurrence_table — not published.

    input:
    path quals
    val basename
    val id_length

    output:
    path "${basename}.qual", emit: qual

    shell:
    '''
    sort --key=3,3n --key=1,1d --key=2,2n --merge !{quals} | \
        uniq --check-chars=!{id_length} > !{basename}.qual
    '''
}
