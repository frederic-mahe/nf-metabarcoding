// ============================================================================
// Part C — taxonomic assignment (stampa re-implementation)
// ============================================================================
// `part_C_assign` is the shared assignment core: given a representatives
// FASTA it runs the method selected by --taxonomy_method (stampa
// scatter-gather or sintax) and publishes the standalone
// <basename>_taxonomy_<method>.tsv. `part_C` wraps it with the
// occurrence-table extraction and join for table input ([S46]→[S51]);
// `part_C_fasta` runs the assignment alone for fasta input ([S48]).
// `part_C_shadow` is the sintax-only shadow path ([S50]) and keeps its
// own wiring (it ignores --taxonomy_method).

include { extract_fasta_sequences_from_occurrence_table } from '../../modules/local/part_c/extract_fasta_sequences_from_occurrence_table.nf'
include { assign_taxonomy_stampa }                         from '../../modules/local/part_c/assign_taxonomy_stampa.nf'
include { sort_taxonomy }                                  from '../../modules/local/part_c/sort_taxonomy.nf'
include { assign_taxonomy_sintax }                         from '../../modules/local/part_c/assign_taxonomy_sintax.nf'
include { update_occurrence_table }                        from '../../modules/local/part_c/update_occurrence_table.nf'
include { compute_majority_assignment }                    from '../../modules/local/part_c/compute_majority_assignment.nf'
include { normalize_path; log_dir; coerce_bool }             from '../../modules/local/functions.nf'


workflow part_C_assign {
    // Shared taxonomic-assignment core for the regular Part C paths:
    // given a representatives FASTA, run the method selected by
    // --taxonomy_method ([S61]) and publish the standalone
    // <basename>_taxonomy_<method>.tsv. Called by `part_C` (table
    // input, which then splices the result back onto the occurrence
    // table) and by `part_C_fasta` ([S48] fasta input, which keeps the
    // standalone table as its only deliverable). The shadow path does
    // NOT use this helper: it is sintax-only regardless of
    // --taxonomy_method, so it keeps its own wiring.

    take:
    representatives   // channel: representatives FASTA (one file)
    basename          // value channel: <basename>

    main:
    // [S47]/[S64]: pick the reference matching the assignment method.
    // The two formats are not interchangeable; stampa expects a
    // space-separated lineage in the header and sintax expects a
    // `;tax=...;` annotation. The relevant assert in the entry point
    // already guarantees the matching flag is set before we land here.
    def reference = ( params.taxonomy_method == 'sintax' )
        ? file(normalize_path(params.reference_dataset_sintax))
        : file(normalize_path(params.reference_dataset))

    def assignments
    if ( params.taxonomy_method == 'sintax' ) {
        // [S50]/[S61]: vsearch --sintax. Emits the canonical 5-column
        // assignments intermediate (consumed by the table-input join)
        // and publishes the standalone 4-column
        // <basename>_taxonomy_sintax.tsv.
        assign_taxonomy_sintax(representatives, reference, basename)
        assignments = assign_taxonomy_sintax.out.taxonomy
    } else {
        // [S49]: stampa scatter-gather. Split the representatives
        // fasta into chunks (or pass it through unchanged when
        // params.stampa_chunk_size == 0), fan each chunk through
        // assign_taxonomy_stampa in parallel, then concatenate +
        // sort the slices via collectFile.
        //
        def chunks = (params.stampa_chunk_size > 0)
            ? representatives.splitFasta(by: params.stampa_chunk_size, file: true)
            : representatives

        assign_taxonomy_stampa(chunks, reference)

        // [S49]/[S59]/D16: gather each chunk's vsearch.log into the
        // single published <basename>_taxonomy.log under logs/part_c/,
        // the stampa counterpart of the sintax path's --log artefact.
        // Empty input (no chunks) → nothing gathered, no log published,
        // matching the "only stages that run produce output" contract.
        assign_taxonomy_stampa.out.log
            .combine(basename)
            .collectFile(storeDir: log_dir('part_c')) { chunk_log, bn ->
                ["${bn}_taxonomy.log", chunk_log.text]
            }

        // Gather the per-chunk slices into a single (unsorted) file,
        // then stabilise the order in the sort_taxonomy process.
        //
        // collectFile cannot do the sorting itself: its `sort:` closure
        // orders whole entries (chunks), not the lines inside a chunk,
        // so it only sorted correctly in the degenerate
        // one-record-per-chunk case and left multi-record chunks (the
        // default and the `local`/`demo` profiles) in vsearch's
        // thread-dependent order. Here collectFile only concatenates;
        // sort_taxonomy applies `LC_ALL=C sort -k2,2nr -k1,1d`
        // (abundance desc, amplicon asc) and publishes the result.
        def gathered = assign_taxonomy_stampa.out.taxonomy
            .collectFile(name: 'taxonomy_stampa.unsorted.tsv')

        sort_taxonomy(gathered, basename)
        assignments = sort_taxonomy.out.taxonomy
    }

    emit:
    // Canonical per-amplicon assignments channel: the stampa
    // <basename>_taxonomy_stampa.tsv (5-column, header) or the sintax
    // <basename>_assignments_sintax.tsv intermediate. Empty (stampa,
    // empty input → zero chunks) when there is nothing to assign;
    // table-input callers fall back to an empty file ([S51]).
    taxonomy = assignments
}


workflow part_C {
    // Table-input Part C ([S46]→[S51]): extract a representatives FASTA
    // from the occurrence table, assign taxonomy via the shared
    // `part_C_assign` helper, then splice the assignment back onto the
    // table. Shared by the standalone Part C workflow and the
    // Part B → Part C end-to-end path.

    take:
    occurrence_table   // value channel: Path to <basename>_table.tsv
    basename           // value channel: <project>_<N>_samples

    main:
    // [S48]: extract a representatives FASTA from the occurrence table.
    extract_fasta_sequences_from_occurrence_table(occurrence_table)

    part_C_assign(
        extract_fasta_sequences_from_occurrence_table.out.fasta,
        basename,
    )

    // [S51] empty-input contract: when the occurrence table has no
    // rows, the extracted representatives fasta is empty, the stampa
    // scatter emits zero chunks, and `part_C_assign` emits nothing —
    // which would leave update_occurrence_table dangling. Fall back to
    // an empty assignments file so update_occurrence_table still runs
    // and publishes a header-only <basename>_table_assigned.tsv. (The
    // sintax branch always emits a file, so ifEmpty is a no-op there.)
    def empty_assignments = file("${workflow.workDir}/empty_assignments.tsv")
    empty_assignments.text = ""

    update_occurrence_table(
        occurrence_table,
        part_C_assign.out.taxonomy.ifEmpty(empty_assignments),
        basename,
    )

    // [S66]: opt-in majority-rule assignment, run on the freshly
    // assigned table. Stampa branch only — the startup assert ([S66])
    // guarantees majority is never combined with sintax, so the
    // stampa-formatted --reference_dataset is the right reference.
    if ( coerce_bool(params.majority_assignment) ) {
        def reference = file(normalize_path(params.reference_dataset))
        compute_majority_assignment(
            update_occurrence_table.out.table,
            reference,
            basename,
        )
    }

    emit:
    // [S51] annotated occurrence table — sibling of Part B's
    // <basename>_table.tsv, named <basename>_table_assigned.tsv.
    table = update_occurrence_table.out.table
}


workflow part_C_fasta {
    // [S48] fasta-input Part C: the user supplies a representatives
    // FASTA via --representatives_fasta instead of an occurrence table.
    // The extraction and the occurrence-table join are skipped; the
    // sole deliverable is the standalone <basename>_taxonomy_<method>.tsv
    // published by the assignment process (sort_taxonomy for stampa,
    // assign_taxonomy_sintax for sintax). No <basename>_table_assigned.tsv
    // is produced — there is no table to splice onto.

    take:
    representatives   // value channel: Path to the representatives FASTA
    basename          // value channel: <basename>

    main:
    part_C_assign(representatives, basename)

    emit:
    taxonomy = part_C_assign.out.taxonomy
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
        extract_fasta_sequences_from_occurrence_table.out.fasta,
        reference,
        basename,
    )

    update_occurrence_table(
        populated,
        assign_taxonomy_sintax.out.taxonomy,
        basename,
    )

    emit:
    table = update_occurrence_table.out.table
}
