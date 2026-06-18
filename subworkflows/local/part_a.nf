include { discover_inputs }                 from '../../modules/local/part_a/discover_inputs.nf'
include { merge_fastq_pairs }               from '../../modules/local/part_a/merge_fastq_pairs.nf'
include { strip_reads }                     from '../../modules/local/part_a/strip_reads.nf'
include { join_notmerged }                  from '../../modules/local/part_a/join_notmerged.nf'
include { trim_primers }                    from '../../modules/local/part_a/trim_primers.nf'
include { filter_and_convert_to_fasta }     from '../../modules/local/part_a/filter_and_convert_to_fasta.nf'
include { extract_expected_error_values }   from '../../modules/local/part_a/extract_expected_error_values.nf'
include { dereplicate_fasta }               from '../../modules/local/part_a/dereplicate_fasta.nf'
include { list_local_clusters }             from '../../modules/local/part_a/list_local_clusters.nf'
include { hash_relabel_flag; hash_id_length } from '../../modules/local/functions.nf'


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
    // [S10]/[S11]/[S12]/[S21]: pattern-driven discovery. Pairs go
    // through merge_fastq_pairs; single-end files skip the merging
    // step.
    discover_inputs()

    def branched = discover_inputs.out
        .splitCsv(sep: '\t')
        .map { row ->
            def sample = row[0]
            def r1 = file(row[1])
            def r2 = (row.size() > 2 && row[2]) ? file(row[2]) : null
            tuple(sample, r1, r2)
        }
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

    // [S04] shadow pipeline: optionally strip the low-quality 3' tails
    // ([S24]) and then join unmerged R1/R2 with A-padding (length
    // params.join_padding_length, default 8 — see [S63]). The shadow
    // sampleId is `<sampleId>_notmerged`, so all subsequent processes
    // naturally publish artefacts at `<sampleId>_notmerged.*` without
    // touching the regular pipeline. Because the padding uses A/C/G/T
    // only, the shadow and regular paths share the same downstream
    // filters (--fastq_maxns 0) and the same swarm invocation —
    // no per-path branching is needed after this point.
    strip_reads(merge_fastq_pairs.out.notmerged)
    join_notmerged(
        strip_reads.out.stripped.map { id, fwd, rev ->
            tuple("${id}_notmerged", fwd, rev)
        }
    )

    // Build a unified (id, fastq) stream for the rest of Part A.
    def regular_ch = merge_fastq_pairs.out.merged.mix(unpaired_ch)
    def shadow_ch = join_notmerged.out.joined

    def to_process = regular_ch.mix(shadow_ch)

    // trim primers (skipped when --no_trimming is set)
    def trimmed_ch
    if ( params.no_trimming ) {
        trimmed_ch = to_process
    } else {
        trim_primers(to_process)
        trimmed_ch = trim_primers.out.trimmed
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

    emit:
    // (sampleId, path) tuples — shadow samples carry a _notmerged
    // suffix so the main workflow can filter them into the shadow
    // stream before calling part_B / part_B_shadow.
    fasta = dereplicate_fasta.out.fasta
    qual  = extract_expected_error_values.out.qual
    stats = list_local_clusters.out.stats
}
