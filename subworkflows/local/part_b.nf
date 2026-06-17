
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
    global_clustering(global_dereplication.out[0], basename)

    // Placeholder taxonomy + chimera flag run in parallel — both
    // consume only the swarm representatives ([S33]/[S34]).
    fake_taxonomic_assignment(global_clustering.out[3], basename)
    chimera_detection(global_clustering.out[3], basename)

    // Re-cleave clusters using the per-sample stats ([S22]).
    cleaving(
        global_clustering.out[1],            // .stats
        global_clustering.out[2],            // .struct
        global_clustering.out[0],            // .swarms
        global_dereplication.out[0],         // global .fas
        list_all_cluster_seeds_of_size_greater_than_2.out[0],
        basename
    )

    // [S36]/[S37]: post-cleave taxonomy + chimera annotations.
    fake_taxonomic_assignment2(cleaving.out[2], basename)
    chimera_detection_post_cleave(
        global_clustering.out[3],            // pre-cleave reps
        cleaving.out[2],                     // cleaved reps (fas2)
        chimera_detection.out[1],            // pre-cleave stderr ([S45])
        basename
    )

    // [S35]: assemble the filtered occurrence table.
    build_occurrence_table(
        global_clustering.out[3],            // _1f_representatives.fas
        cleaving.out[2],                     // _1f_representatives.fas2
        global_clustering.out[1],            // _1f.stats
        cleaving.out[0],                     // _1f.stats2
        global_clustering.out[0],            // _1f.swarms
        cleaving.out[1],                     // _1f.swarms2
        chimera_detection.out[0],            // _1f_representatives.uchime
        chimera_detection_post_cleave.out[0],           // _1f_representatives.uchime2
        build_expected_error_file.out[0],    // .qual
        fake_taxonomic_assignment.out[0],    // _1f_representatives.results
        fake_taxonomic_assignment2.out[0],   // _1f_representatives.results2
        build_distribution_file.out[0],      // .distr
        basename,
        sample_ids,                          // [S09]
    )

    // [S38]/[S39]: collapse sub- and super-string OTUs. The combined
    // vsearch + python step emits <basename>_superstring_clustering.log
    // ([S45]).
    search_for_terminal_gaps(build_occurrence_table.out[0])
    merge_substring_otus(
        build_occurrence_table.out[0],
        search_for_terminal_gaps.out[0],   // .uc hits
        search_for_terminal_gaps.out[1],   // vsearch search.log
        basename,
    )

    // [S40]–[S44]: mumu (ex-lulu) post-clustering filter pass.
    extract_otu_fasta(merge_substring_otus.out[0])
    trim_metadata_for_mumu(merge_substring_otus.out[0])
    find_similar_sequences(extract_otu_fasta.out[0])
    run_mumu(trim_metadata_for_mumu.out[0], find_similar_sequences.out[0], basename)
    rebuild_post_mumu_table(run_mumu.out[0], merge_substring_otus.out[0], basename)
    extract_mumu_fasta(rebuild_post_mumu_table.out[0])

    emit:
    // [S46] final occurrence table — consumed by part_C when Part B
    // and Part C run end-to-end.
    table = rebuild_post_mumu_table.out[0]
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

    global_clustering(global_dereplication.out[0], basename)
    def representatives = global_clustering.out[3]

    fake_taxonomic_assignment(representatives, basename)
    chimera_detection(representatives, basename)

    cleaving(
        global_clustering.out[1],            // .stats
        global_clustering.out[2],            // .struct
        global_clustering.out[0],            // .swarms
        global_dereplication.out[0],         // global .fas
        list_all_cluster_seeds_of_size_greater_than_2.out[0],
        basename
    )

    fake_taxonomic_assignment2(cleaving.out[2], basename)
    chimera_detection_post_cleave(
        representatives,
        cleaving.out[2],
        chimera_detection.out[1],
        basename
    )

    build_occurrence_table(
        representatives,
        cleaving.out[2],
        global_clustering.out[1],
        cleaving.out[0],
        global_clustering.out[0],
        cleaving.out[1],
        chimera_detection.out[0],
        chimera_detection_post_cleave.out[0],
        build_expected_error_file.out[0],
        fake_taxonomic_assignment.out[0],
        fake_taxonomic_assignment2.out[0],
        build_distribution_file.out[0],
        basename,
        sample_ids,                          // [S09]
    )

    search_for_terminal_gaps(build_occurrence_table.out[0])
    merge_substring_otus(
        build_occurrence_table.out[0],
        search_for_terminal_gaps.out[0],
        search_for_terminal_gaps.out[1],
        basename,
    )

    extract_otu_fasta(merge_substring_otus.out[0])
    trim_metadata_for_mumu(merge_substring_otus.out[0])
    find_similar_sequences(extract_otu_fasta.out[0])
    run_mumu(trim_metadata_for_mumu.out[0], find_similar_sequences.out[0], basename)
    rebuild_post_mumu_table(run_mumu.out[0], merge_substring_otus.out[0], basename)
    extract_mumu_fasta(rebuild_post_mumu_table.out[0])

    emit:
    // [S46] final shadow occurrence table — exposed so callers can
    // splice taxonomy onto the shadow path via a standalone Part C
    // invocation. End-to-end wiring runs Part C on the regular path
    // only; calling it twice from the same scope is forbidden in
    // DSL2.
    table = rebuild_post_mumu_table.out[0]
}
