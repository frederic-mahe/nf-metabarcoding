include { publish_dir; log_dir } from '../functions.nf'


process assign_taxonomy_sintax {
    // [S50]: shadow Part C runs `vsearch --sintax` against the same
    // --reference_dataset ([S47]) and reshapes the tabbed output to
    // the canonical assignments TSV
    // (amplicon\tabundance\tidentity\ttaxonomy\treferences).
    //
    // Column mapping from vsearch --sintax --tabbedout:
    //   col 1 → amplicon (after stripping ;size=N;)
    //   col 2 → taxonomy (bootstrap-annotated lineage, e.g.
    //           `d:Bacteria(1.00),p:Proteobacteria(0.85)`)
    //   col 3 (strand) and col 4 (cutoff-filtered lineage) are ignored.
    // `identity` and `references` stay at placeholder values
    // ([S33]/[S46]) — sintax does not produce a percent identity, and
    // there is no single "reference accession" to record.
    //
    // The process emits two TSVs:
    //   * <basename>_assignments_sintax.tsv — the canonical 5-column
    //     join format above; an intermediate consumed by
    //     update_occurrence_table that stays in the work dir
    //     (mirroring the stampa path's unpublished chunks).
    //   * <basename>_taxonomy_sintax.tsv — vsearch's verbatim 4-column
    //     --tabbedout (query\ttaxonomy\tstrand\tcutoff_taxonomy) with a
    //     header row, published to occurrence_table/ as the sintax
    //     counterpart of the stampa path's <basename>_taxonomy_stampa.tsv
    //     ([S59]/[S61]).
    //
    // The process is also used by the shadow path, whose `basename`
    // carries the `_notmerged` token. Per the developer's policy the
    // standalone sintax table is published only when the user explicitly
    // selects --taxonomy_method=sintax (regular path) — never for the
    // shadow run — so the occurrence_table publishDir skips any
    // `_notmerged` basename via its saveAs closure. The log is always
    // published ([D16]: under logs/part_c/).
    publishDir path: { log_dir('part_c') }, mode: params.publish_mode, pattern: "*.log"
    publishDir(
        path: { publish_dir('occurrence_table') },
        mode: params.publish_mode,
        pattern: "*_taxonomy_sintax.tsv",
        saveAs: { filename -> basename.toString().contains('_notmerged') ? null : filename },
    )

    input:
    path representatives
    path reference_dataset
    val basename

    output:
    path "${basename}_assignments_sintax.tsv", emit: taxonomy
    path "${basename}_taxonomy_sintax.tsv",    emit: published
    path "${basename}_taxonomy.log",           emit: log

    shell:
    '''
    vsearch \
        --sintax !{representatives} \
        --threads !{task.cpus} \
        --randseed !{params.sintax_randseed} \
        --db !{reference_dataset} \
        --dbmask none \
        --sintax_cutoff !{params.sintax_cutoff} \
        --tabbedout raw_sintax.tsv \
        --log !{basename}_taxonomy.log

    awk 'BEGIN { FS = OFS = "\t" }
         {
             split($1, parts, ";size=")
             amplicon = parts[1]
             gsub(/;.*$/, "", parts[2])
             abundance = parts[2] + 0
             taxonomy = ($2 == "") ? "NA" : $2
             print amplicon, abundance, "0.0", taxonomy, "NA"
         }' raw_sintax.tsv > !{basename}_assignments_sintax.tsv

    {
        printf 'query\\ttaxonomy\\tstrand\\tcutoff_taxonomy\\n'
        cat raw_sintax.tsv
    } > !{basename}_taxonomy_sintax.tsv
    '''

    stub:
    """
    touch ${basename}_assignments_sintax.tsv ${basename}_taxonomy_sintax.tsv ${basename}_taxonomy.log
    """
}
