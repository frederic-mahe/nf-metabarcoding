process check_reference {
    // [S114]-[S117] (D22): read every header of the reference once and
    // check it against the grammar of the method that reads it
    // (stampa | sintax). The task always succeeds, even when it finds
    // errors: the subworkflow publishes the report, then stops the run
    // before any assignment (reference_check_gate in functions.nf).
    // No publishDir here, for that reason.

    input:
    path reference_dataset
    val grammar      // 'stampa' | 'sintax'
    val parameter    // parameter name quoted in the messages
    val majority     // true: reject commas in accessions ([S111])
    val basename

    output:
    path "${basename}_reference_check.tsv", emit: report
    path "reference_check_summary.txt",     emit: summary

    script:
    def majority_flag = majority ? '--majority' : ''
    """
    check_reference.py \\
        --grammar ${grammar} \\
        --parameter ${parameter} \\
        ${majority_flag} \\
        --report ${basename}_reference_check.tsv \\
        --summary reference_check_summary.txt \\
        ${reference_dataset}
    """

    stub:
    """
    printf 'line\\taccession\\tissue\\n' > ${basename}_reference_check.tsv
    echo clean > reference_check_summary.txt
    """
}
