#!/usr/bin/env nextflow

include { normalize_path; validate_params; samplesheet_profile; check_primer_format; resource_size_warnings } from './modules/local/functions.nf'
include { part_A } from './subworkflows/local/part_a.nf'
include { part_B; part_B_shadow } from './subworkflows/local/part_b.nf'
include { discover_part_b_fasta } from './modules/local/part_b/discover_part_b_fasta.nf'
include { validate_samplesheet } from './modules/local/validate_samplesheet.nf'
include { part_C; part_C_fasta; part_C_shadow } from './subworkflows/local/part_c.nf'
include { dump_software_versions } from './modules/local/dump_software_versions.nf'
include { validateParameters; paramsHelp } from 'plugin/nf-schema'


workflow part_c {
    // Standalone table-input Part C ([S47]/[S48]): the user provides
    // --occurrence_table (Part B's <basename>_table.tsv) and a
    // reference whose format matches --taxonomy_method ([S47]/[S64]).
    // The fasta-input variant is the sibling `part_c_fasta` entry,
    // reached via --representatives_fasta.
    if ( params.taxonomy_method == 'sintax' ) {
        assert params.reference_dataset_sintax :
            "--reference_dataset_sintax must be set when --taxonomy_method=sintax"
    } else {
        assert params.reference_dataset :
            "--reference_dataset must be set when --taxonomy_method=stampa"
    }
    assert params.occurrence_table :
        "--occurrence_table must be set when running table-input Part C standalone"

    // [S71]/D08: no startup mkdirs — publishDir materialises
    // <outdir>/<subdir> on first publish (no parse-time filesystem I/O).
    // --outdir always resolves (default 'results'), so no output-folder
    // param is required ([S26] superseded).

    def table_path = file(normalize_path(params.occurrence_table))
    def derived_basename = table_path.baseName.replaceFirst(/_table$/, '')
    def basename_ch = Channel.value(derived_basename)
    def table_ch    = Channel.fromPath(normalize_path(params.occurrence_table))

    part_C(table_ch, basename_ch)

    // [S62]/[S64]/D07: route the shadow sibling through a staged
    // channel instead of a parse-time disk probe. part_C_shadow is
    // wired unconditionally (within the gate), so the DAG shape no
    // longer depends on head-node disk state at parse time; a runtime
    // `.filter { it.exists() }` empties the channel when the
    // <basename>_notmerged_table.tsv sibling is absent, so the branch
    // self-suppresses (no work, no output). The --recover_unmerged
    // ([S78]) and --reference_dataset_sintax gates are param checks
    // (config, not disk state), so they stay an `if`.
    if ( params.recover_unmerged && params.reference_dataset_sintax ) {
        def shadow_table_ch = Channel
            .fromPath(
                "${table_path.parent}/${derived_basename}_notmerged_table.tsv",
                checkIfExists: false
            )
            .filter { it.exists() }
        def shadow_basename_ch = Channel.value("${derived_basename}_notmerged")
        part_C_shadow(shadow_table_ch, shadow_basename_ch)
    }
}


workflow part_c_fasta {
    // Standalone fasta-input Part C ([S48]): the user provides
    // --representatives_fasta (a representatives FASTA, header
    // `<amplicon>;size=<abundance>;`) and a reference whose format
    // matches --taxonomy_method ([S47]/[S64]). The occurrence-table
    // extraction and join are skipped; the deliverable is the
    // standalone <basename>_taxonomy_<method>.tsv. The basename is
    // derived from the fasta filename. No shadow path runs in this
    // mode (there is no _notmerged sibling concept for a bare fasta).
    if ( params.taxonomy_method == 'sintax' ) {
        assert params.reference_dataset_sintax :
            "--reference_dataset_sintax must be set when --taxonomy_method=sintax"
    } else {
        assert params.reference_dataset :
            "--reference_dataset must be set when --taxonomy_method=stampa"
    }

    def fasta_path = file(normalize_path(params.representatives_fasta))
    def basename_ch = Channel.value(fasta_path.baseName)
    def fasta_ch    = Channel.fromPath(normalize_path(params.representatives_fasta))

    part_C_fasta(fasta_ch, basename_ch)
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

    // [S71]/D08: no startup mkdirs — publishDir materialises
    // <outdir>/<subdir> on first publish (no parse-time filesystem I/O).
    // --outdir always resolves (default 'results'), so no output-folder
    // param is required ([S26] superseded).

    // Build the six per-sample lists (regular + shadow × fasta/qual/
    // stats) from either the validated --input samplesheet (fasta
    // profile, [S70]) or the folder scan ([S27]/[S56]). The folder scan
    // derives qual/stats as siblings of each .fas (unchanged); the
    // samplesheet carries them explicitly (defaulting to the same
    // siblings in bin/parse_samplesheet.py). Both paths split into the
    // regular and shadow ([S56]) streams by the _notmerged suffix.
    def fasta_list
    def qual_list
    def stats_list
    def s_fasta_list
    def s_qual_list
    def s_stats_list
    if ( params.input ) {
        validate_samplesheet(file(normalize_path(params.input)))
        def split_ch = validate_samplesheet.out
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(
                row.sample,
                file(row.fasta, checkIfExists: true),
                file(row.qual,  checkIfExists: true),
                file(row.stats, checkIfExists: true)
            ) }
            .branch { id, _f, _q, _s ->
                shadow:  id.endsWith('_notmerged')
                regular: !id.endsWith('_notmerged')
            }
        fasta_list   = split_ch.regular.map { _id, f, _q, _s -> f }.collect()
        qual_list    = split_ch.regular.map { _id, _f, q, _s -> q }.collect()
        stats_list   = split_ch.regular.map { _id, _f, _q, s -> s }.collect()
        s_fasta_list = split_ch.shadow.map { _id, f, _q, _s -> f }.collect()
        s_qual_list  = split_ch.shadow.map { _id, _f, q, _s -> q }.collect()
        s_stats_list = split_ch.shadow.map { _id, _f, _q, s -> s }.collect()
    } else {
        discover_part_b_fasta()

        def regular_samples_ch = discover_part_b_fasta.out.regular
            .splitCsv(sep: '\t')
            .map { row -> tuple(row[0], file(row[1])) }
        def shadow_samples_ch = discover_part_b_fasta.out.shadow
            .splitCsv(sep: '\t')
            .map { row -> tuple(row[0], file(row[1])) }

        fasta_list = regular_samples_ch.map { _id, f -> f }.collect()
        qual_list  = regular_samples_ch
            .map { id, f -> file("${f.parent}/${id}.qual") }
            .collect()
        stats_list = regular_samples_ch
            .map { id, f -> file("${f.parent}/${id}.stats") }
            .collect()
        s_fasta_list = shadow_samples_ch.map { _id, f -> f }.collect()
        s_qual_list  = shadow_samples_ch
            .map { id, f -> file("${f.parent}/${id}.qual") }
            .collect()
        s_stats_list = shadow_samples_ch
            .map { id, f -> file("${f.parent}/${id}.stats") }
            .collect()
    }
    part_B(fasta_list, qual_list, stats_list)

    // [S47]/[S64]: invoke Part C when the user supplied the reference
    // matching the selected method. Shadow Part C ([S50]) additionally
    // requires --recover_unmerged ([S78]) and --reference_dataset_sintax
    // — see [S64].
    def regular_ref_present = ( params.taxonomy_method == 'sintax' )
        ? params.reference_dataset_sintax
        : params.reference_dataset
    if ( regular_ref_present ) {
        def basename_ch = fasta_list.map { files ->
            "${params.project_name}_${files.size()}_samples"
        }
        part_C(part_B.out.table, basename_ch)
    }

    // [S78]: the shadow Part B / Part C path is opt-in. When
    // --recover_unmerged is unset (the default), part_B_shadow is never
    // invoked, so no `_notmerged` Part B table is produced even if
    // `_notmerged.fas` files were discovered in --fasta_folder.
    if ( params.recover_unmerged ) {
        part_B_shadow(s_fasta_list, s_qual_list, s_stats_list)
        if ( params.reference_dataset_sintax ) {
            // [S50]/[S64]: shadow taxonomy on part_B_shadow's occurrence
            // table. The header-only gate inside part_C_shadow
            // suppresses output when the shadow side had no samples.
            def shadow_basename_ch = s_fasta_list.map { files ->
                "${params.project_name}_${files.size()}_samples_notmerged"
            }
            part_C_shadow(part_B_shadow.out.table, shadow_basename_ch)
        }
    }
}


workflow {
    // [S57]: --help short-circuits before any required-param assert.
    // print() (not log.info) so the usage block lands on stdout —
    // the conventional channel for help output and what nf-test
    // captures in workflow.stdout.
    if ( params.help ) {
        // [S57]: a small hand-written banner carries the two things the
        // parameter schema cannot express — the three entry points and
        // the README/SPECIFICATIONS pointer — wrapped around the
        // schema-generated parameter listing from paramsHelp() (every
        // params.* flag grouped by part, with type and default).
        print """\
            nf-metabarcoding — swarm-based metabarcoding pipeline

            Usage:
              nextflow run main.nf --fastq_folder PATH         [Part A → B [→ C]]
              nextflow run main.nf --fasta_folder PATH         [Part B standalone]
              nextflow run main.nf --occurrence_table PATH     [Part C standalone]
              nextflow run main.nf --representatives_fasta PATH [Part C, fasta input]

            """.stripIndent()
        print paramsHelp(
            command: "nextflow run main.nf --fastq_folder PATH " +
                     "--forward_primer SEQ --reverse_primer SEQ"
        )
        print """\

            See README.md for examples and SPECIFICATIONS.md for the behaviour
            contract ([Sxx] IDs).
            """.stripIndent()
        return
    }

    // [S58]/[S61]/[S65]/[S63]/[S72]: per-value type / enum / range checks,
    // declared in nextflow_schema.json and enforced here. Runs first so a
    // bad scalar aborts naming the parameter before validate_params()'
    // file-content check ([S73]) — which reads a reference header — gets a
    // chance to fire on an otherwise-doomed run.
    validateParameters()

    // [S02]/[S66]/[S71]: cross-parameter and file-content guards the
    // schema cannot express (input-mode exclusivity, the majority/sintax
    // incompatibility, reference-format sniffing, the --results_folder
    // deprecation warning). Mode-specific requirements (fastq_folder /
    // primers / reference) stay inline below because they depend on the
    // selected run mode.
    validate_params()

    // [S79]: under -profile slurm, warn (don't abort) when the dataset /
    // reference size hints are unset so a forgotten --dataset_size_gb
    // surfaces here rather than as an OOM kill deep in a large run. The
    // resource tiers only apply under slurm, so the helper is silent
    // otherwise; `workflow.profile` is the active comma-joined profile
    // list.
    resource_size_warnings(
        workflow.profile,
        params.dataset_size_gb,
        params.reference_size_gb,
        (params.reference_dataset || params.reference_dataset_sintax) as boolean
    ).each { System.err.println("WARNING: ${it}") }

    // [S68]/[S71]: record tool versions under <outdir>/pipeline_info/ on
    // every entry point — including a Part A-only run. --outdir always
    // resolves (default 'results'), so the report is published regardless
    // of the selected mode. Invoked once here (not per router/branch) so
    // the single-host process is never scheduled twice in one run.
    dump_software_versions()

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

    // [S48]: --representatives_fasta switches the pipeline into
    // fasta-input Part C mode (Parts A and B do not run; no
    // occurrence-table join). Requires the reference matching
    // --taxonomy_method; the deliverable is <basename>_taxonomy_<method>.tsv.
    if ( params.representatives_fasta ) {
        part_c_fasta()
        return
    }

    // [S27]/[S70]: --fasta_folder, or an --input samplesheet in the
    // fasta profile, switches the pipeline into Part B standalone mode.
    // Part A's --fastq_folder requirement is lifted in that mode (Part A
    // does not run). When --reference_dataset is also set, `workflow
    // part_b` chains Part C onto Part B's regular output.
    def input_profile = samplesheet_profile()
    if ( params.fasta_folder || input_profile == 'fasta' ) {
        part_b()
        return
    }

    // [S18]/[S70]: Part A end-to-end runs from --fastq_folder or from an
    // --input samplesheet (fastq profile). forward_primer / reverse_primer
    // are required (no default); fastq_folder is required only when
    // --input is not given.
    if ( !params.input ) {
        assert params.fastq_folder : "--fastq_folder must be set (no default)"
    }

    // [S18]/[S20]: primers and --no_trimming are mutually exclusive — when
    // params.no_trimming is true, forward_primer and reverse_primer must be
    // empty and the trim_primers step is skipped; otherwise both are required.
    if ( params.no_trimming ) {
        assert !params.forward_primer : "--forward_primer must be empty when --no_trimming is set"
        assert !params.reverse_primer : "--reverse_primer must be empty when --no_trimming is set"
    } else {
        assert params.forward_primer : "--forward_primer must be set (no default)"
        assert params.reverse_primer : "--reverse_primer must be set (no default)"
        // [S74]: both primers must be IUPAC nucleotide strings (the
        // trim_primers / reverse_complement.sh alphabet) — fail fast on a
        // typo rather than mid-run in cutadapt.
        check_primer_format('forward_primer', params.forward_primer)
        check_primer_format('reverse_primer', params.reverse_primer)
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
        // [S71]/D08: no startup mkdirs — publishDir materialises
        // <outdir>/<subdir> on first publish (no parse-time filesystem
        // I/O). --outdir always resolves (default 'results'), so no
        // output-folder param is required ([S26] superseded).

        def regular_fasta = part_A.out.fasta.filter { id, _f -> !id.endsWith("_notmerged") }
        def regular_qual  = part_A.out.qual .filter { id, _f -> !id.endsWith("_notmerged") }
        def regular_stats = part_A.out.stats.filter { id, _f -> !id.endsWith("_notmerged") }

        // join the three streams on sample ID so the lists stay
        // aligned even when processes complete out of order.
        def joined_b = regular_fasta
            .join(regular_qual)
            .join(regular_stats)
        def b_fasta = joined_b.map { _id, f, _q, _s -> f }.collect()
        def b_qual  = joined_b.map { _id, _f, q, _s -> q }.collect()
        def b_stats = joined_b.map { _id, _f, _q, s -> s }.collect()
        part_B(b_fasta, b_qual, b_stats)

        // [S47]/[S64]: run regular Part C when the reference matching
        // --taxonomy_method is set.
        def regular_ref_present = ( params.taxonomy_method == 'sintax' )
            ? params.reference_dataset_sintax
            : params.reference_dataset
        if ( regular_ref_present ) {
            def basename_ch = b_fasta.map { files ->
                "${params.project_name}_${files.size()}_samples"
            }
            part_C(part_B.out.table, basename_ch)
        }

        // [S78]: the shadow Part B / Part C path is opt-in. When
        // --recover_unmerged is unset (the default), Part A emits no
        // `_notmerged` samples and the shadow sub-workflows are never
        // invoked. Shadow Part C ([S50]) additionally requires
        // --reference_dataset_sintax ([S64]).
        if ( params.recover_unmerged ) {
            def shadow_fasta = part_A.out.fasta.filter { id, _f -> id.endsWith("_notmerged") }
            def shadow_qual  = part_A.out.qual .filter { id, _f -> id.endsWith("_notmerged") }
            def shadow_stats = part_A.out.stats.filter { id, _f -> id.endsWith("_notmerged") }

            def joined_shadow = shadow_fasta
                .join(shadow_qual)
                .join(shadow_stats)
            def s_fasta = joined_shadow.map { _id, f, _q, _s -> f }.collect()
            def s_qual  = joined_shadow.map { _id, _f, q, _s -> q }.collect()
            def s_stats = joined_shadow.map { _id, _f, _q, s -> s }.collect()
            part_B_shadow(s_fasta, s_qual, s_stats)

            if ( params.reference_dataset_sintax ) {
                // [S50]/[S64]: shadow taxonomy on part_B_shadow's table.
                def shadow_basename_ch = s_fasta.map { files ->
                    "${params.project_name}_${files.size()}_samples_notmerged"
                }
                part_C_shadow(part_B_shadow.out.table, shadow_basename_ch)
            }
        }
    }
}
