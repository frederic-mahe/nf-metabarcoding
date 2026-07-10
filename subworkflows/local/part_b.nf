
// ============================================================================
// Part B — global clustering / occurrence table
// ============================================================================
// Each process below receives the per-sample artefacts collected by
// the Part B fasta channel ([S27]) and a project-wide `basename`
// string of the shape `<project_name>_<N>_samples` (see [S25]). The
// processes can run in parallel — they do not depend on each other.

include { build_expected_error_file }                      from '../../modules/local/part_b/build_expected_error_file.nf'
include { build_distribution_file }                        from '../../modules/local/part_b/build_distribution_file.nf'
include { list_all_cluster_seeds_of_size_greater_than_2 }  from '../../modules/local/part_b/list_all_cluster_seeds_of_size_greater_than_2.nf'
include { global_dereplication }                           from '../../modules/local/part_b/global_dereplication.nf'
include { global_clustering }                              from '../../modules/local/part_b/global_clustering.nf'
include { fake_taxonomic_assignment }                      from '../../modules/local/part_b/fake_taxonomic_assignment.nf'
include { chimera_detection }                              from '../../modules/local/part_b/chimera_detection.nf'
include { fake_taxonomic_assignment2 }                     from '../../modules/local/part_b/fake_taxonomic_assignment2.nf'
include { chimera_detection_post_cleave }                  from '../../modules/local/part_b/chimera_detection_post_cleave.nf'
include { build_occurrence_table }                         from '../../modules/local/part_b/build_occurrence_table.nf'
include { cleaving }                                       from '../../modules/local/part_b/cleaving.nf'
include { search_for_terminal_gaps }                       from '../../modules/local/part_b/search_for_terminal_gaps.nf'
include { merge_substring_otus }                           from '../../modules/local/part_b/merge_substring_otus.nf'
include { extract_otu_fasta }                              from '../../modules/local/part_b/extract_otu_fasta.nf'
include { extract_mumu_fasta }                             from '../../modules/local/part_b/extract_mumu_fasta.nf'
include { trim_metadata_for_mumu }                         from '../../modules/local/part_b/trim_metadata_for_mumu.nf'
include { find_similar_sequences }                         from '../../modules/local/part_b/find_similar_sequences.nf'
include { run_mumu }                                       from '../../modules/local/part_b/run_mumu.nf'
include { rebuild_post_mumu_table }                        from '../../modules/local/part_b/rebuild_post_mumu_table.nf'
include { recluster_search }                               from '../../modules/local/part_b/recluster_search.nf'
include { recluster_merge }                                from '../../modules/local/part_b/recluster_merge.nf'
include { extract_recluster_fasta }                        from '../../modules/local/part_b/extract_recluster_fasta.nf'
include { hash_id_length }                                 from '../../modules/local/functions.nf'


workflow part_B {
    // Runs the five Part B pre-cleaving processes on already-collected
    // lists of per-sample fasta / qual / stats files. Shared by the
    // standalone (`workflow part_b`) and end-to-end paths so the
    // process wiring lives in one place.

    take:
    fasta_list   // value channel: List<Path>, one .fas per sample
    qual_list    // value channel: List<Path>, one .qual per sample (same order)
    stats_list   // value channel: List<Path>, one .stats per sample (same order)

    main:
    def basename = fasta_list.map { files ->
        "${params.project_name}_${files.size()}_samples"
    }

    // [S09]: empty samples write no rows to .distr, so the column list
    // derived inside build_filtered_contingency_table from the .distr
    // alone misses them. Hand the authoritative sample list down via
    // build_occurrence_table's --samples so zero columns appear.
    def sample_ids = fasta_list.map { files ->
        files.collect { it.baseName }.toSorted().join(",")
    }

    build_expected_error_file(qual_list, basename, hash_id_length())
    build_distribution_file(fasta_list, basename)
    list_all_cluster_seeds_of_size_greater_than_2(stats_list, basename)
    global_dereplication(fasta_list, basename)
    global_clustering(global_dereplication.out.fasta, basename)

    // Placeholder taxonomy + chimera flag run in parallel — both
    // consume only the swarm representatives ([S33]/[S34]).
    fake_taxonomic_assignment(global_clustering.out.representatives, basename)
    chimera_detection(global_clustering.out.representatives, basename)

    // Re-cleave clusters using the per-sample stats ([S22]).
    cleaving(
        global_clustering.out.stats,
        global_clustering.out.struct,
        global_clustering.out.swarms,
        global_dereplication.out.fasta,
        list_all_cluster_seeds_of_size_greater_than_2.out.stats,
        basename
    )

    // [S36]/[S37]: post-cleave taxonomy + chimera annotations.
    fake_taxonomic_assignment2(cleaving.out.representatives, basename)
    chimera_detection_post_cleave(
        global_clustering.out.representatives,
        cleaving.out.representatives,
        chimera_detection.out.log,
        basename
    )

    // [S35]: assemble the filtered occurrence table.
    build_occurrence_table(
        global_clustering.out.representatives,
        cleaving.out.representatives,
        global_clustering.out.stats,
        cleaving.out.stats,
        global_clustering.out.swarms,
        cleaving.out.swarms,
        chimera_detection.out.uchime,
        chimera_detection_post_cleave.out.uchime,
        build_expected_error_file.out.qual,
        fake_taxonomic_assignment.out.results,
        fake_taxonomic_assignment2.out.results,
        build_distribution_file.out.distr,
        basename,
        sample_ids,                          // [S09]
    )

    // [S38]/[S39]: collapse sub- and super-string OTUs. The combined
    // vsearch + python step emits <basename>_superstring_clustering.log
    // ([S45]).
    search_for_terminal_gaps(build_occurrence_table.out.table)
    merge_substring_otus(
        build_occurrence_table.out.table,
        search_for_terminal_gaps.out.uc,
        search_for_terminal_gaps.out.log,
        basename,
    )

    // [S40]–[S44]: mumu (ex-lulu) post-clustering filter pass.
    extract_otu_fasta(merge_substring_otus.out.table)
    trim_metadata_for_mumu(merge_substring_otus.out.table)
    find_similar_sequences(extract_otu_fasta.out.fasta)
    run_mumu(trim_metadata_for_mumu.out.table, find_similar_sequences.out.matches, basename)
    rebuild_post_mumu_table(run_mumu.out.table, merge_substring_otus.out.table, basename)
    extract_mumu_fasta(rebuild_post_mumu_table.out.table)

    // [S102]–[S105]/D20: optional, terminal coarse re-clustering. Gated
    // on the master switch --recluster_id (a compile-time param check,
    // D-b). When on, the reclustered table replaces the post-mumu table
    // as the emitted deliverable (D-a) and Part C runs on the reduced
    // set; when off (default) no recluster process is scheduled and the
    // output is byte-identical to before.
    def final_table = rebuild_post_mumu_table.out.table
    if ( params.recluster_id ) {
        recluster_search(extract_mumu_fasta.out.fasta)
        recluster_merge(
            rebuild_post_mumu_table.out.table,
            recluster_search.out.uc,
            recluster_search.out.log,
            basename,
        )
        extract_recluster_fasta(recluster_merge.out.table)
        final_table = recluster_merge.out.table
    }

    emit:
    // [S46]/[S105] final occurrence table — the reclustered table when
    // --recluster_id is set, otherwise the post-mumu table. Consumed by
    // part_C when Part B and Part C run end-to-end.
    table = final_table
}


workflow part_B_shadow {
    // [S56]: shadow Part B. Runs the **same processes** as part_B on
    // the Part A shadow outputs (<sampleId>_notmerged.{fas,qual,stats}
    // — see [S04]). The A-padding emitted by Part A is composed
    // entirely of A/C/G/T (see [S63]), so swarm accepts the sequences
    // as-is and no mask/restore wrapping is needed around
    // global_clustering. The basename carries a trailing _notmerged
    // token so every published artefact is distinguishable from the
    // regular Part B output.

    take:
    fasta_list   // value channel: List<Path>, one _notmerged.fas per sample
    qual_list    // value channel: List<Path>, one _notmerged.qual per sample
    stats_list   // value channel: List<Path>, one _notmerged.stats per sample

    main:
    def basename = fasta_list.map { files ->
        "${params.project_name}_${files.size()}_samples_notmerged"
    }

    // [S09]: same column-list contract as the regular Part B path —
    // empty shadow samples must still appear as zero columns.
    def sample_ids = fasta_list.map { files ->
        files.collect { it.baseName }.toSorted().join(",")
    }

    build_expected_error_file(qual_list, basename, hash_id_length())
    build_distribution_file(fasta_list, basename)
    list_all_cluster_seeds_of_size_greater_than_2(stats_list, basename)
    global_dereplication(fasta_list, basename)

    global_clustering(global_dereplication.out.fasta, basename)
    def representatives = global_clustering.out.representatives

    fake_taxonomic_assignment(representatives, basename)
    chimera_detection(representatives, basename)

    cleaving(
        global_clustering.out.stats,
        global_clustering.out.struct,
        global_clustering.out.swarms,
        global_dereplication.out.fasta,
        list_all_cluster_seeds_of_size_greater_than_2.out.stats,
        basename
    )

    fake_taxonomic_assignment2(cleaving.out.representatives, basename)
    chimera_detection_post_cleave(
        representatives,
        cleaving.out.representatives,
        chimera_detection.out.log,
        basename
    )

    build_occurrence_table(
        representatives,
        cleaving.out.representatives,
        global_clustering.out.stats,
        cleaving.out.stats,
        global_clustering.out.swarms,
        cleaving.out.swarms,
        chimera_detection.out.uchime,
        chimera_detection_post_cleave.out.uchime,
        build_expected_error_file.out.qual,
        fake_taxonomic_assignment.out.results,
        fake_taxonomic_assignment2.out.results,
        build_distribution_file.out.distr,
        basename,
        sample_ids,                          // [S09]
    )

    search_for_terminal_gaps(build_occurrence_table.out.table)
    merge_substring_otus(
        build_occurrence_table.out.table,
        search_for_terminal_gaps.out.uc,
        search_for_terminal_gaps.out.log,
        basename,
    )

    extract_otu_fasta(merge_substring_otus.out.table)
    trim_metadata_for_mumu(merge_substring_otus.out.table)
    find_similar_sequences(extract_otu_fasta.out.fasta)
    run_mumu(trim_metadata_for_mumu.out.table, find_similar_sequences.out.matches, basename)
    rebuild_post_mumu_table(run_mumu.out.table, merge_substring_otus.out.table, basename)
    extract_mumu_fasta(rebuild_post_mumu_table.out.table)

    // [S105]/D20: the re-clustering gate is symmetric with the regular
    // part_B path (the _notmerged basename is already threaded through).
    def final_table = rebuild_post_mumu_table.out.table
    if ( params.recluster_id ) {
        recluster_search(extract_mumu_fasta.out.fasta)
        recluster_merge(
            rebuild_post_mumu_table.out.table,
            recluster_search.out.uc,
            recluster_search.out.log,
            basename,
        )
        extract_recluster_fasta(recluster_merge.out.table)
        final_table = recluster_merge.out.table
    }

    emit:
    // [S46]/[S105] final shadow occurrence table — exposed so callers
    // can splice taxonomy onto the shadow path via a standalone Part C
    // invocation. End-to-end wiring runs Part C on the regular path
    // only; calling it twice from the same scope is forbidden in
    // DSL2.
    table = final_table
}
