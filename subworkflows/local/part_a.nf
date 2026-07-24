include { discover_inputs }                 from '../../modules/local/part_a/discover_inputs.nf'
include { merge_fastq_pairs }               from '../../modules/local/part_a/merge_fastq_pairs.nf'
include { strip_reads }                     from '../../modules/local/part_a/strip_reads.nf'
include { join_notmerged }                  from '../../modules/local/part_a/join_notmerged.nf'
include { trim_primers }                    from '../../modules/local/part_a/trim_primers.nf'
include { filter_and_convert_to_fasta }     from '../../modules/local/part_a/filter_and_convert_to_fasta.nf'
include { extract_expected_error_values }   from '../../modules/local/part_a/extract_expected_error_values.nf'
include { dereplicate_fasta }               from '../../modules/local/part_a/dereplicate_fasta.nf'
include { list_local_clusters }             from '../../modules/local/part_a/list_local_clusters.nf'
include { summarize_read_counts }           from '../../modules/local/part_a/summarize_read_counts.nf'
include { validate_samplesheet }            from '../../modules/local/validate_samplesheet.nf'
include { normalize_path; hash_relabel_flag; hash_id_length; coerce_bool } from '../../modules/local/functions.nf'


workflow part_A {
    // Part A end-to-end: fastq → per-sample dereplicated fasta + qual
    // + stats. Wraps the regular and shadow ([S04]) pipelines so every
    // process tagged with this part shows up as `part_A:<process>` in
    // the Nextflow tty output.
    //
    // Emits three (sampleId, path) channels — one per per-sample
    // artefact — so the main workflow can split them into regular vs
    // shadow streams when feeding part_B / part_B_shadow.

    main:
    // Build the (sample, r1, r2) stream from either the validated
    // --input samplesheet (fastq profile, [S70]) or the pattern-driven
    // folder scan ([S10]/[S11]/[S12]/[S21]). Either way, pairs go
    // through merge_fastq_pairs and single-end files (null r2) skip the
    // merging step.
    def samples
    if ( params.input ) {
        validate_samplesheet(file(normalize_path(params.input)))
        samples = validate_samplesheet.out
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                def r2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null
                tuple(row.sample, file(row.fastq_1, checkIfExists: true), r2)
            }
    } else {
        discover_inputs()
        samples = discover_inputs.out
            .splitCsv(sep: '\t')
            .map { row ->
                def r2 = (row.size() > 2 && row[2]) ? file(row[2]) : null
                tuple(row[0], file(row[1]), r2)
            }
    }

    def branched = samples
        .branch { _sample, _r1, r2 ->
            paired:   r2 != null
            unpaired: r2 == null
        }

    def paired_ch = branched.paired.map { sample, r1, r2 ->
        tuple(sample, [r1, r2])
    }
    def unpaired_ch = branched.unpaired.map { sample, r1, _r2 ->
        tuple(sample, r1)
    }

    // discover pairs and merge
    merge_fastq_pairs(paired_ch)

    // Build a unified (id, fastq) stream for the rest of Part A.
    def regular_ch = merge_fastq_pairs.out.merged.mix(unpaired_ch)

    // [S04]/[S78] shadow pipeline, opt-in via --recover_unmerged
    // (default false). When enabled, optionally strip the low-quality
    // 3' tails ([S24]) and then join unmerged R1/R2 with A-padding
    // (length params.join_padding_length, default 8 — see [S63]). The
    // shadow sampleId is `<sampleId>_notmerged`, so all subsequent
    // processes naturally publish artefacts at `<sampleId>_notmerged.*`
    // without touching the regular pipeline. Because the padding uses
    // A/C/G/T only, the shadow and regular paths share the same
    // downstream filters (--fastq_maxns 0) and the same swarm
    // invocation — no per-path branching is needed after this point.
    // When disabled (the default), the not-merged reads are simply not
    // consumed: no shadow processes run and no `_notmerged` artefacts
    // are produced.
    def to_process
    if ( coerce_bool(params.recover_unmerged) ) {
        strip_reads(merge_fastq_pairs.out.notmerged)
        join_notmerged(
            strip_reads.out.stripped.map { id, fwd, rev ->
                tuple("${id}_notmerged", fwd, rev)
            }
        )
        to_process = regular_ch.mix(join_notmerged.out.joined)
    } else {
        to_process = regular_ch
    }

    // trim primers (skipped when --no_trimming is set)
    def trimmed_ch
    def trim_fwd_log_ch
    def trim_rev_log_ch
    if ( coerce_bool(params.no_trimming) ) {
        trimmed_ch = to_process
        // no trimming logs in this mode — the read-count summary ([S86])
        // records the F / R columns as 0.
        trim_fwd_log_ch = Channel.empty()
        trim_rev_log_ch = Channel.empty()
    } else {
        trim_primers(to_process)
        trimmed_ch = trim_primers.out.trimmed
        trim_fwd_log_ch = trim_primers.out.log_forward
        trim_rev_log_ch = trim_primers.out.log_reverse
    }

    // convert to fasta with the hash digest + ee ([S65]), apply
    // min-length / max-N=0 filters. The A-padded shadow reads carry no
    // Ns, so they pass the same max-N=0 threshold as the regular path.
    filter_and_convert_to_fasta(trimmed_ch, hash_relabel_flag())

    // set aside EE values
    extract_expected_error_values(
        filter_and_convert_to_fasta.out.fasta,
        hash_id_length()
    )

    // dereplicate
    dereplicate_fasta(filter_and_convert_to_fasta.out.fasta)

    list_local_clusters(dereplicate_fasta.out.fasta)

    // [S86]: per-sample read-count summary. Runs on every fastq path
    // (Part A-only and end-to-end); non-fastq modes never invoke part_A,
    // so it only ever runs on a fastq path. Called here, inside part_A,
    // so it shows up as `part_A:summarize_read_counts` alongside the
    // other Part A steps. Shadow (_notmerged) samples are excluded; the
    // basename uses the Part B construct when --project_name is set,
    // otherwise <N>_samples.
    def rc_ids = dereplicate_fasta.out.fasta
        .filter { id, _f -> !id.endsWith("_notmerged") }
        .map { id, _f -> id }
    def rc_ids_file = rc_ids.collectFile(name: 'sample_ids.txt', newLine: true, sort: true)
    def rc_logs = merge_fastq_pairs.out.log
        .mix(trim_fwd_log_ch)
        .mix(trim_rev_log_ch)
        .filter { logfile -> !logfile.name.contains("_notmerged") }
        .collect()
        .ifEmpty([])
    def rc_basename = rc_ids.count().map { n ->
        params.project_name ? "${params.project_name}_${n}_samples" : "${n}_samples"
    }
    summarize_read_counts(rc_ids_file, rc_logs, rc_basename)

    emit:
    // (sampleId, path) tuples — shadow samples carry a _notmerged
    // suffix so the main workflow can filter them into the shadow
    // stream before calling part_B / part_B_shadow.
    fasta = dereplicate_fasta.out.fasta
    qual  = extract_expected_error_values.out.qual
    stats = list_local_clusters.out.stats
    // per-step log streams feed the read-count summary ([S86]). These
    // are bare path channels (filenames carry the sampleId and the
    // _notmerged token); merge_fastq_pairs runs for every paired sample
    // so merge_log is empty only when every sample is single-end.
    merge_log    = merge_fastq_pairs.out.log
    trim_fwd_log = trim_fwd_log_ch
    trim_rev_log = trim_rev_log_ch
}
