process build_distribution_file {
    // [S29]: scan FASTA headers of every input .fas and emit the
    // sequence-to-sample mapping as tab-separated rows
    // <sha1>\t<sampleId>\t<size>. The sample ID is derived from the
    // fasta basename so the channel can be built without an explicit
    // sample-ID side car.
    //
    // [S59]: the .distr is an internal intermediate consumed by
    // build_occurrence_table — not published.

    input:
    path fastas
    val basename

    output:
    path "${basename}.distr", emit: distr

    shell:
    '''
    : > !{basename}.distr
    for f in !{fastas} ; do
        sample="$(basename "${f}" .fas)"
        # `|| true`: empty samples ([S09]/[S27]) have a zero-record
        # .fas — grep returns 1, which would otherwise abort the
        # process. The empty sample legitimately contributes no rows.
        grep "^>" "${f}" | \
            sed 's/^>// ; s/;size=/\t/ ; s/;$//' | \
            awk -v s="${sample}" 'BEGIN {OFS = "\t"} {print $1, s, $2}' \
            >> !{basename}.distr || true
    done
    '''
}
