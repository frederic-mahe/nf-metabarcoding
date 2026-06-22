process fake_taxonomic_assignment2 {
    // [S36]: same transform as [S33] but operating on the cleaver's
    // representatives (cleaving's third output, `*_1f_representatives.fas2`).
    // An empty fas2 (no clusters got cleaved) is a legitimate input
    // — the output is then an empty results2 file.
    //
    // [S59]: post-cleave placeholder taxonomy is an internal
    // intermediate consumed by build_occurrence_table — not
    // published.

    input:
    path cleaved_representatives
    val basename

    output:
    path "${basename}_1f_representatives.results2", emit: results

    shell:
    '''
    set -euo pipefail

    : > !{basename}_1f_representatives.results2
    grep "^>" !{cleaved_representatives} | \
        sed -r 's/^>//
                s/;size=/\t/
                s/;?$/\t0.0\tNA\tNA/' \
        >> !{basename}_1f_representatives.results2 || true
    '''

    stub:
    """
    touch ${basename}_1f_representatives.results2
    """
}
