include { normalize_path } from '../functions.nf'


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
    // The process is also used by the regular path when
    // params.taxonomy_method == 'sintax' ([S61]); the published
    // filename embeds `basename`, which carries the `_notmerged` token
    // on the shadow path.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode,
        enabled: params.results_folder != null

    input:
    path representatives
    path reference_dataset
    val basename

    output:
    path "${basename}_taxonomy_sintax.tsv"
    path "${basename}_taxonomy.log"

    shell:
    '''
    vsearch \
        --sintax !{representatives} \
        --threads !{task.cpus} \
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
         }' raw_sintax.tsv > !{basename}_taxonomy_sintax.tsv
    '''
}
