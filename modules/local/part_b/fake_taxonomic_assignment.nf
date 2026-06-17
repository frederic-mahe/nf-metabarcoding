process fake_taxonomic_assignment {
    // [S33]: emit a placeholder TSV that the occurrence-table builder
    // can consume before Part C lands. Each row mirrors the stampa
    // output shape: <amplicon>\t<size>\t<identity>\t<taxonomy>\t<refs>.
    //
    // [S59]: placeholder taxonomy is an internal intermediate
    // consumed by build_occurrence_table — not published.

    input:
    path representatives
    val basename

    output:
    path "${basename}_1f_representatives.results"

    shell:
    '''
    grep "^>" !{representatives} | \
        sed -r 's/^>//
                s/;size=/\t/
                s/;?$/\t0.0\tNA\tNA/' > !{basename}_1f_representatives.results
    '''
}
