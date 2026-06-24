process assign_taxonomy_stampa {
    // [S49]: per-chunk step of the stampa scatter-gather. Run
    // vsearch --usearch_global against the reference dataset on
    // one fasta chunk, then collapse the top hits per amplicon
    // with bin/stampa_merge.py. The workflow uses splitFasta
    // upstream to produce chunks and collectFile downstream to
    // gather + sort the slices into the final
    // <basename>_taxonomy_stampa.tsv. --notrunclabels keeps the
    // ";size=N;" annotation in vsearch's --userout, which
    // stampa_merge.py parses directly.

    input:
    path chunk
    path reference_dataset

    output:
    path "stampa_chunk.tsv", emit: taxonomy
    // [S49]/[S59]/D16: the per-chunk vsearch log is gathered by the
    // subworkflow into the published logs/part_c/<basename>_taxonomy.log
    // (the stampa counterpart of the sintax path's --log). Not published
    // here: a fixed-name per-chunk file would collide across the scatter.
    path "vsearch.log",      emit: log

    shell:
    '''
    vsearch \
        --usearch_global !{chunk} \
        --threads !{task.cpus} \
        --db !{reference_dataset} \
        --dbmask none \
        --qmask none \
        --rowlen 0 \
        --notrunclabels \
        --userfields query+id!{params.iddef}+target \
        --maxaccepts 0 \
        --maxrejects !{params.stampa_maxrejects} \
        --top_hits_only \
        --output_no_hits \
        --id !{params.stampa_id} \
        --iddef !{params.iddef} \
        --userout hits \
        2> vsearch.log

    stampa_merge.py hits > stampa_chunk.tsv
    '''

    stub:
    """
    touch stampa_chunk.tsv vsearch.log
    """
}
