process find_similar_sequences {
    // [S42]: vsearch --usearch_global self-search; lulu-recommended
    // parameters. The legacy bash strips `;size=N;` from every
    // column with a sed pass — that's the format mumu accepts.

    input:
    path otu_fasta

    output:
    path "${otu_fasta.baseName}.match_list", emit: matches

    shell:
    '''
    set -euo pipefail

    vsearch \
        --usearch_global !{otu_fasta} \
        --db !{otu_fasta} \
        --self \
        --threads !{task.cpus} \
        --id !{params.similar_id} \
        --iddef 1 \
        --userfields query+target+id \
        --maxaccepts 0 \
        --query_cov !{params.similar_query_cov} \
        --maxhits !{params.similar_maxhits} \
        --quiet \
        --userout - | \
        sed -r 's/;size=[0-9]+;//g' > !{otu_fasta.baseName}.match_list
    '''

    stub:
    """
    touch ${otu_fasta.baseName}.match_list
    """
}
