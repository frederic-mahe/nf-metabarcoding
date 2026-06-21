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
    path "${basename}_1f_representatives.results", emit: results

    shell:
    '''
    set -euo pipefail

    # [S81]: `grep "^>"` exits 1 when the representatives fasta carries
    # no headers (an empty input), which would abort the task under
    # `set -euo pipefail`. Pre-create the output and trail the pipe with
    # `|| true` so an empty input yields an empty .results and succeeds
    # (mirrors fake_taxonomic_assignment2).
    : > !{basename}_1f_representatives.results
    grep "^>" !{representatives} | \
        sed -r 's/^>//
                s/;size=/\t/
                s/;?$/\t0.0\tNA\tNA/' \
        >> !{basename}_1f_representatives.results || true
    '''
}
