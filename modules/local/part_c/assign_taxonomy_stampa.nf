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
    //
    // [S108]: --minwordmatches 0 turns off the k-mer pre-filter, so
    // that, with --maxaccepts 0 and --maxrejects 0, every reference is
    // aligned (vsearch's documented exhaustive search). Without it, a
    // reference sharing fewer than 12 distinct 8-mers with the query
    // is never considered, which drops the best hit of short or
    // divergent amplicons.
    //
    // TODO: refactor with vsearch --search_global once vsearch v2.33.0
    // is available on bioconda. That command aligns every query
    // against every reference with no pre-filter and no early stop,
    // and ignores --minwordmatches, --maxaccepts and --maxrejects, so
    // it replaces --usearch_global plus these three options (and makes
    // params.stampa_maxrejects moot).

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
        --minwordmatches 0 \
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
