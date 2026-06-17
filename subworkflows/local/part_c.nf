// ============================================================================
// Part C — taxonomic assignment (stampa re-implementation)
// ============================================================================
// Skeleton phase: the processes are wired but the test coverage is
// tagged `pending` because D04 still needs to land before the
// stand-alone CLI / output policy is finalised.

include { extract_fasta_sequences_from_occurrence_table } from '../../modules/local/part_c/extract_fasta_sequences_from_occurrence_table.nf'
include { assign_taxonomy_stampa }                         from '../../modules/local/part_c/assign_taxonomy_stampa.nf'
include { assign_taxonomy_sintax }                         from '../../modules/local/part_c/assign_taxonomy_sintax.nf'
include { update_occurrence_table }                        from '../../modules/local/part_c/update_occurrence_table.nf'
include { compute_majority_assignment }                    from '../../modules/local/part_c/compute_majority_assignment.nf'
include { normalize_path }                                 from '../../modules/local/functions.nf'


workflow part_C {
    // Re-usable Part C wiring shared by the standalone workflow and
    // the Part B → Part C end-to-end path.

    take:
    occurrence_table   // value channel: Path to <basename>_table.tsv
    basename           // value channel: <project>_<N>_samples

    main:
    // [S47]/[S64]: pick the reference matching the assignment method.
    // The two formats are not interchangeable; stampa expects a
    // space-separated lineage in the header and sintax expects a
    // `;tax=...;` annotation. The relevant assert in the entry point
    // already guarantees the matching flag is set before we land here.
    def reference = ( params.taxonomy_method == 'sintax' )
        ? file(normalize_path(params.reference_dataset_sintax))
        : file(normalize_path(params.reference_dataset))

    // [S48]: extract a representatives FASTA from the occurrence
    // table. (The fasta-input branch promised by [S48] is blocked
    // on D04 and is not wired in this skeleton.)
    extract_fasta_sequences_from_occurrence_table(occurrence_table)

    if ( params.taxonomy_method == 'sintax' ) {
        // [S50]: shadow path.
        assign_taxonomy_sintax(
            extract_fasta_sequences_from_occurrence_table.out[0],
            reference,
            basename,
        )
        update_occurrence_table(
            occurrence_table,
            assign_taxonomy_sintax.out[0],
            basename,
        )
    } else {
        // [S49]: stampa scatter-gather. Split the representatives
        // fasta into chunks (or pass it through unchanged when
        // params.stampa_chunk_size == 0), fan each chunk through
        // assign_taxonomy_stampa in parallel, then concatenate +
        // sort the slices via collectFile.
        //
        // Each chunk's TSV is paired with the basename value
        // channel via `combine` so the collectFile closure can
        // build "${basename}_taxonomy_stampa.tsv" per-item
        // (basename can't be interpolated into collectFile's
        // `name:` String at workflow-build time).
        //
        // Note on the sort closure: collectFile sorts in JVM heap,
        // which is fine for nf-metabarcoding's representative
        // counts (typically < a few million OTUs). If a real run
        // ever hits the OOM line, swap this for an external
        // `sort_taxonomy` process running `LC_ALL=C sort -k2,2nr
        // -k1,1d` on an unsorted collectFile output. See Plan B
        // (2026-05-19) for the fallback writeup.
        def reps_ch = extract_fasta_sequences_from_occurrence_table.out[0]
        def chunks = (params.stampa_chunk_size > 0)
            ? reps_ch.splitFasta(by: params.stampa_chunk_size, file: true)
            : reps_ch

        def merged = assign_taxonomy_stampa(chunks, reference)
            .combine(basename)
            .collectFile(
                storeDir: normalize_path(params.results_folder),
                sort: { line ->
                    // -k2,2nr -k1,1d == abundance desc, amplicon asc.
                    // collectFile's sort closure takes one argument
                    // (a line) and must return a Comparable key.
                    // ArrayLists aren't Comparable on the JVM side,
                    // so we build a single String key: width-13 inverted
                    // abundance (so descending becomes ascending under
                    // lexicographic ordering), then a tab, then the
                    // amplicon name (already ascending).
                    def f = line.tokenize('\t')
                    String.format(
                        "%013d\t%s",
                        9999999999999L - (f[1] as Long),
                        f[0],
                    )
                },
            ) { chunk_tsv, bn ->
                ["${bn}_taxonomy_stampa.tsv", chunk_tsv.text]
            }

        // [S51] empty-input contract: when the occurrence table has no
        // rows, the extracted representatives fasta is empty,
        // splitFasta emits zero chunks, and `merged` never fires —
        // which would leave update_occurrence_table dangling. Fall
        // back to an empty assignments file so update_occurrence_table
        // still runs and publishes a header-only
        // <basename>_table_assigned.tsv.
        def empty_assignments = file("${workflow.workDir}/empty_assignments.tsv")
        empty_assignments.text = ""

        update_occurrence_table(
            occurrence_table,
            merged.ifEmpty(empty_assignments),
            basename,
        )

        // [S66]: opt-in majority-rule assignment, run on the freshly
        // assigned table. Stampa branch only — the startup assert
        // ([S66]) guarantees we never reach here under sintax, and the
        // shadow path never calls part_C. `reference` is the
        // stampa-formatted --reference_dataset resolved above.
        if ( params.majority_assignment ) {
            compute_majority_assignment(
                update_occurrence_table.out[0],
                reference,
                basename,
            )
        }
    }

    emit:
    // [S51] annotated occurrence table — sibling of Part B's
    // <basename>_table.tsv, named <basename>_table_assigned.tsv.
    table = update_occurrence_table.out[0]
}


workflow part_C_shadow {
    // [S50]: shadow Part C. Runs `vsearch --sintax` (always, no
    // toggle) on the shadow occurrence table produced by shadow
    // Part B ([S56]) or supplied via [S62]'s standalone-mode probe.
    // The basename carries the `_notmerged` token so every
    // published artefact is distinguishable from the regular Part C
    // output. Shares `extract_fasta_sequences_from_occurrence_table`
    // and `update_occurrence_table` with `part_C` — DSL2 allows
    // these processes to be invoked from both sub-workflows because
    // each invocation lives in its own workflow scope.
    //
    // [S50] short-circuit: when the upstream produced no `_notmerged`
    // samples, `part_B_shadow` still emits a header-only table
    // (`_0_samples_notmerged_table.tsv`). The filter on `populated`
    // drops that case so no shadow artefacts are published when there
    // is nothing to assign. Standalone mode reaches this workflow only
    // when the sibling file exists ([S62]), but the filter is also a
    // safety net there (empty sibling → no published shadow artefact).

    take:
    occurrence_table   // value channel: Path to <basename>_notmerged_table.tsv
    basename           // value channel: <project>_<N>_samples_notmerged

    main:
    // [S64]: shadow Part C consumes the sintax-formatted reference.
    // Entry points gate this workflow on params.reference_dataset_sintax
    // being non-null, so the file() call below is safe by the time the
    // workflow is invoked.
    def reference = file(normalize_path(params.reference_dataset_sintax))

    def populated = occurrence_table.filter { tbl ->
        // any non-empty line after the header → there is at least one
        // amplicon to assign.
        tbl.toFile().readLines().drop(1).any { it.trim() }
    }

    extract_fasta_sequences_from_occurrence_table(populated)

    assign_taxonomy_sintax(
        extract_fasta_sequences_from_occurrence_table.out[0],
        reference,
        basename,
    )

    update_occurrence_table(
        populated,
        assign_taxonomy_sintax.out[0],
        basename,
    )

    emit:
    table = update_occurrence_table.out[0]
}
