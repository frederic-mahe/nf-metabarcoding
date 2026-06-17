#!/usr/bin/env nextflow

include { normalize_path; usage } from './modules/local/functions.nf'
include { part_A } from './subworkflows/local/part_a.nf'
include { part_B; part_B_shadow } from './subworkflows/local/part_b.nf'
include { discover_part_b_fasta } from './modules/local/part_b/discover_part_b_fasta.nf'


// ============================================================================
// Part C — taxonomic assignment (stampa re-implementation)
// ============================================================================
// Skeleton phase: the processes are wired but the test coverage is
// tagged `pending` because D04 still needs to land before the
// stand-alone CLI / output policy is finalised.

process extract_fasta_sequences_from_occurrence_table {
    // [S48]: extract a representatives FASTA from the occurrence
    // table. Column 4 = amplicon ID, column 2 = abundance, column
    // 10 = sequence (same layout as [S40]'s extract_otu_fasta).

    input:
    path occurrence_table

    output:
    path "${occurrence_table.baseName}_representatives.fas"

    shell:
    '''
    awk 'NR > 1 {printf ">"$4";size="$2";\\n"$10"\\n"}' !{occurrence_table} \
        > !{occurrence_table.baseName}_representatives.fas
    '''
}


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
    path "stampa_chunk.tsv"

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
}


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


process update_occurrence_table {
    // [S51]: splice the taxonomy assignments back onto Part B's
    // <basename>_table.tsv. The output is published as a sibling
    // file `<basename>_table_assigned.tsv` so Part B's unannotated
    // table is preserved alongside the annotated one (D04 sub-q2).
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode,
        enabled: params.results_folder != null

    input:
    path occurrence_table
    path assignments
    val basename

    output:
    path "${basename}_table_assigned.tsv"

    shell:
    '''
    update_occurrence_table.py \
        --occurrence_table !{occurrence_table} \
        --assignments !{assignments} \
        > !{basename}_table_assigned.tsv
    '''
}


process compute_majority_assignment {
    // [S66]: opt-in final step of the regular Part C path. Recompute a
    // majority-rule taxonomy per OTU from the reference accessions in
    // the `references` column of the assigned table
    // (<basename>_table_assigned.tsv, [S51]) and the stampa-formatted
    // reference dataset ([S47]). Emits an independent three-column
    // table <basename>_table_assigned_majority.tsv
    // (OTU\tamplicon\ttaxonomy_majority). Never runs on the shadow
    // path; gated on params.majority_assignment and (by the startup
    // assert) on --taxonomy_method=stampa.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode,
        enabled: params.results_folder != null

    input:
    path assigned_table
    path reference_dataset
    val basename

    output:
    path "${basename}_table_assigned_majority.tsv"

    shell:
    '''
    majority_assignment.py \
        --input_table !{assigned_table} \
        --reference_db !{reference_dataset} \
        > !{basename}_table_assigned_majority.tsv
    '''
}


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


workflow part_c {
    // Standalone Part C ([S47]/[S48]): the user provides
    // --occurrence_table (Part B's <basename>_table.tsv) and a
    // reference whose format matches --taxonomy_method ([S47]/[S64]).
    // The fasta-input branch is blocked on D04 and is not exposed by
    // this skeleton.
    if ( params.taxonomy_method == 'sintax' ) {
        assert params.reference_dataset_sintax :
            "--reference_dataset_sintax must be set when --taxonomy_method=sintax"
    } else {
        assert params.reference_dataset :
            "--reference_dataset must be set when --taxonomy_method=stampa"
    }
    assert params.occurrence_table :
        "--occurrence_table must be set when running Part C standalone " +
        "(fasta input is blocked on DECISIONS.md D04)"
    assert params.results_folder :
        "--results_folder must be set when running Part C"

    def results_dir = new File(normalize_path(params.results_folder).toString())
    results_dir.mkdirs()

    def table_path = file(normalize_path(params.occurrence_table))
    def derived_basename = table_path.baseName.replaceFirst(/_table$/, '')
    def basename_ch = Channel.value(derived_basename)
    def table_ch    = Channel.fromPath(normalize_path(params.occurrence_table))

    part_C(table_ch, basename_ch)

    // [S62]/[S64]: probe for a shadow sibling next to the input
    // table. The file lookup is done at workflow-build time (not
    // inside a channel operator) so the toggle is "the file exists on
    // disk now". When the sibling is present **and**
    // --reference_dataset_sintax is set, route the sibling through
    // part_C_shadow alongside the regular part_C invocation — DSL2
    // allows two distinct sub-workflows in the same scope. When
    // either condition is unmet, the shadow workflow is simply not
    // invoked.
    def shadow_table_path = file(
        "${table_path.parent}/${derived_basename}_notmerged_table.tsv"
    )
    if ( shadow_table_path.exists() && params.reference_dataset_sintax ) {
        def shadow_basename_ch = Channel.value("${derived_basename}_notmerged")
        def shadow_table_ch    = Channel.fromPath(shadow_table_path.toString())
        part_C_shadow(shadow_table_ch, shadow_basename_ch)
    }
}


workflow part_b {
    // Standalone Part B ([S25]/[S26]/[S27]/[S56]): discover the
    // per-sample .fas files under params.fasta_folder and route them
    // to two parallel workflows:
    //   - part_B        for samples without _notmerged suffix
    //   - part_B_shadow for samples with the _notmerged suffix
    // .qual and .stats are derived as sibling files of each .fas, so
    // the same per-sample file-naming convention applies to both paths.
    //
    // When --reference_dataset is also set, Part C is wired onto the
    // regular path so the table_assigned.tsv lands alongside the Part B
    // outputs. Shadow taxonomy ([S50]) is wired via the sibling
    // workflow `part_C_shadow` on `part_B_shadow.out.table` — DSL2
    // allows two distinct sub-workflows in the same scope.
    assert params.project_name :
        "--project_name must be set when --fasta_folder is set"
    assert params.results_folder :
        "--results_folder must be set when --fasta_folder is set"

    // [S26]: create the results folder (and any missing parents)
    // before anything publishes into it.
    def results_dir = new File(normalize_path(params.results_folder).toString())
    results_dir.mkdirs()

    discover_part_b_fasta()

    def regular_samples_ch = discover_part_b_fasta.out[0]
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], file(row[1])) }
    def shadow_samples_ch = discover_part_b_fasta.out[1]
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], file(row[1])) }

    def fasta_list = regular_samples_ch.map { _id, f -> f }.collect()
    def qual_list  = regular_samples_ch
        .map { id, f -> file("${f.parent}/${id}.qual") }
        .collect()
    def stats_list = regular_samples_ch
        .map { id, f -> file("${f.parent}/${id}.stats") }
        .collect()
    part_B(fasta_list, qual_list, stats_list)

    def s_fasta_list = shadow_samples_ch.map { _id, f -> f }.collect()
    def s_qual_list  = shadow_samples_ch
        .map { id, f -> file("${f.parent}/${id}.qual") }
        .collect()
    def s_stats_list = shadow_samples_ch
        .map { id, f -> file("${f.parent}/${id}.stats") }
        .collect()
    part_B_shadow(s_fasta_list, s_qual_list, s_stats_list)

    // [S47]/[S64]: invoke Part C when the user supplied the reference
    // matching the selected method. Shadow Part C ([S50]) additionally
    // requires --reference_dataset_sintax — see [S64].
    def regular_ref_present = ( params.taxonomy_method == 'sintax' )
        ? params.reference_dataset_sintax
        : params.reference_dataset
    if ( regular_ref_present ) {
        def basename_ch = fasta_list.map { files ->
            "${params.project_name}_${files.size()}_samples"
        }
        part_C(part_B.out.table, basename_ch)
    }
    if ( params.reference_dataset_sintax ) {
        // [S50]/[S64]: shadow taxonomy on part_B_shadow's occurrence
        // table. The header-only gate inside part_C_shadow suppresses
        // output when the shadow side had no samples upstream.
        def shadow_basename_ch = s_fasta_list.map { files ->
            "${params.project_name}_${files.size()}_samples_notmerged"
        }
        part_C_shadow(part_B_shadow.out.table, shadow_basename_ch)
    }
}


workflow {
    // [S57]: --help short-circuits before any required-param assert.
    // print() (not log.info) so the usage block lands on stdout —
    // the conventional channel for help output and what nf-test
    // captures in workflow.stdout.
    if ( params.help ) {
        print usage()
        return
    }

    // [S58]: validate --publish_mode before any process is wired so a
    // typo aborts the run immediately with a clear message instead of
    // failing on the first PublishDir attempt.
    def allowed_modes = ['copy', 'copyNoFollow', 'link', 'move', 'rellink', 'symlink']
    assert params.publish_mode in allowed_modes :
        "--publish_mode must be one of ${allowed_modes}, got '${params.publish_mode}'"

    // [S61]: validate --taxonomy_method up-front. Only the regular
    // Part C path consults this flag; shadow Part C always uses sintax
    // ([S50]).
    def allowed_taxonomy_methods = ['stampa', 'sintax']
    assert params.taxonomy_method in allowed_taxonomy_methods :
        "--taxonomy_method must be one of ${allowed_taxonomy_methods}, got '${params.taxonomy_method}'"

    // [S66]: majority assignment recomputes a per-OTU taxonomy from the
    // reference accessions listed in the `references` column. Only the
    // stampa method populates that column ([S50] leaves it at the NA
    // placeholder), so --majority_assignment is incompatible with
    // --taxonomy_method=sintax. Fail fast before any process is wired.
    assert !(params.majority_assignment && params.taxonomy_method == 'sintax') :
        "--majority_assignment requires --taxonomy_method=stampa " +
        "(sintax leaves the references column unpopulated)"

    // [S65]: validate --hash_function up-front so an unsupported hash
    // aborts before any process is wired rather than surfacing as an
    // obscure vsearch --relabel error mid-pipeline.
    def allowed_hash_functions = ['sha1', 'md5']
    assert params.hash_function in allowed_hash_functions :
        "--hash_function must be one of ${allowed_hash_functions}, got '${params.hash_function}'"

    // [S63]: validate --join_padding_length as a positive integer
    // before any process is scheduled. Non-positive values, non-integers,
    // and other invalid inputs would otherwise surface as a confusing
    // vsearch error mid-pipeline.
    def jpl = params.join_padding_length
    assert (jpl instanceof Number) && (jpl as int) == jpl && (jpl as int) >= 1 :
        "--join_padding_length must be a positive integer, got '${jpl}'"

    // [S60]: every path-typed param is read through `normalize_path()`
    // at its use site (file(), publishDir, etc.) — see the helper at the
    // top of this file. Nextflow 25's `params` map is read-only from
    // inside a workflow body, so a one-shot mutation here is silently
    // dropped; the wrap-at-use-site pattern is the workaround.

    // [S47]/[S48]: --occurrence_table switches the pipeline into
    // Part C standalone mode (Parts A and B do not run). The
    // mode also requires --reference_dataset and --results_folder.
    if ( params.occurrence_table ) {
        part_c()
        return
    }

    // [S27]: --fasta_folder switches the pipeline into Part B
    // standalone mode. Part A's --fastq_folder requirement is lifted
    // in that mode (Part A does not run). When --reference_dataset is
    // also set, `workflow part_b` chains Part C onto Part B's regular
    // output.
    if ( params.fasta_folder ) {
        part_b()
        return
    }

    // [S18]: forward_primer, reverse_primer and fastq_folder are required
    // and have no default (supply via CLI or project config); the workflow
    // asserts them at startup.
    assert params.fastq_folder : "--fastq_folder must be set (no default)"

    // [S18]/[S20]: primers and --no_trimming are mutually exclusive — when
    // params.no_trimming is true, forward_primer and reverse_primer must be
    // empty and the trim_primers step is skipped; otherwise both are required.
    if ( params.no_trimming ) {
        assert !params.forward_primer : "--forward_primer must be empty when --no_trimming is set"
        assert !params.reverse_primer : "--reverse_primer must be empty when --no_trimming is set"
    } else {
        assert params.forward_primer : "--forward_primer must be set (no default)"
        assert params.reverse_primer : "--reverse_primer must be set (no default)"
    }

    part_A()

    // ----- Part A → Part B handoff ([S27]/[S56]) -----
    // When --project_name is set, Part B runs on the Part A outputs
    // (without going through disk discovery). Regular and shadow
    // samples ([S04]) feed two parallel workflows:
    //   - part_B        on samples whose ID does NOT end in _notmerged
    //   - part_B_shadow on samples whose ID DOES end in _notmerged
    // The shadow workflow's published artefacts carry a _notmerged
    // token in their basename so they sit alongside the regular
    // ones in --results_folder without colliding. When
    // --reference_dataset is also set, Part C runs on the regular
    // Part B table and `part_C_shadow` ([S50]) runs on the shadow
    // table — DSL2 allows two distinct sub-workflows in the same
    // scope.
    if ( params.project_name ) {
        assert params.results_folder :
            "--results_folder must be set when --project_name is set"
        def results_dir = new File(normalize_path(params.results_folder).toString())
        results_dir.mkdirs()

        def regular_fasta = part_A.out.fasta.filter { id, _f -> !id.endsWith("_notmerged") }
        def regular_qual  = part_A.out.qual .filter { id, _f -> !id.endsWith("_notmerged") }
        def regular_stats = part_A.out.stats.filter { id, _f -> !id.endsWith("_notmerged") }

        def shadow_fasta = part_A.out.fasta.filter { id, _f -> id.endsWith("_notmerged") }
        def shadow_qual  = part_A.out.qual .filter { id, _f -> id.endsWith("_notmerged") }
        def shadow_stats = part_A.out.stats.filter { id, _f -> id.endsWith("_notmerged") }

        // join the three streams on sample ID so the lists stay
        // aligned even when processes complete out of order.
        def joined_b = regular_fasta
            .join(regular_qual)
            .join(regular_stats)
        def b_fasta = joined_b.map { _id, f, _q, _s -> f }.collect()
        def b_qual  = joined_b.map { _id, _f, q, _s -> q }.collect()
        def b_stats = joined_b.map { _id, _f, _q, s -> s }.collect()
        part_B(b_fasta, b_qual, b_stats)

        def joined_shadow = shadow_fasta
            .join(shadow_qual)
            .join(shadow_stats)
        def s_fasta = joined_shadow.map { _id, f, _q, _s -> f }.collect()
        def s_qual  = joined_shadow.map { _id, _f, q, _s -> q }.collect()
        def s_stats = joined_shadow.map { _id, _f, _q, s -> s }.collect()
        part_B_shadow(s_fasta, s_qual, s_stats)

        // [S47]/[S64]: same routing as `workflow part_b`. Run regular
        // Part C when the reference matching --taxonomy_method is set;
        // run shadow Part C ([S50]) when --reference_dataset_sintax is
        // set. The two gates are independent: a stampa-only run still
        // produces regular taxonomy; supplying just the sintax flag
        // skips the regular part_C and produces shadow taxonomy only.
        def regular_ref_present = ( params.taxonomy_method == 'sintax' )
            ? params.reference_dataset_sintax
            : params.reference_dataset
        if ( regular_ref_present ) {
            def basename_ch = b_fasta.map { files ->
                "${params.project_name}_${files.size()}_samples"
            }
            part_C(part_B.out.table, basename_ch)
        }
        if ( params.reference_dataset_sintax ) {
            // [S50]/[S64]: shadow taxonomy on part_B_shadow's table.
            // The header-only gate inside part_C_shadow suppresses
            // output when the shadow side had zero _notmerged samples.
            def shadow_basename_ch = s_fasta.map { files ->
                "${params.project_name}_${files.size()}_samples_notmerged"
            }
            part_C_shadow(part_B_shadow.out.table, shadow_basename_ch)
        }
    }
}
