include { normalize_path } from '../functions.nf'


process discover_inputs {
    // Walk every folder listed in params.fastq_folder (a single path
    // or a comma-separated list) and emit a TSV:
    //     sample_id<TAB>r1<TAB>r2   (r2 empty for single-end samples)
    // bin/discover_fastq.py is the single source of truth for the
    // canonical paired-end pattern table ([S11]); params.fastq_pattern,
    // when set, is the user override checked before that table.

    output:
    path "samples.tsv",          emit: samples
    path "discovery_stderr.txt", emit: warnings

    script:
    // [S10]: resolve each folder against launchDir so a relative
    // fastq_folder (from the CLI or a nextflow.config) points at the
    // launch directory — standard Nextflow file() semantics — rather
    // than at the ephemeral task work dir this script runs in.
    // normalize_path expands a leading `~`; file() then anchors any
    // still-relative remainder to launchDir (absolute paths pass
    // through unchanged).
    def folders = (params.fastq_folder instanceof List)
        ? params.fastq_folder
            .collect { file(normalize_path(it)).toAbsolutePath().toString() }
        : params.fastq_folder
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it }
            .collect { file(normalize_path(it)).toAbsolutePath().toString() }
    def folder_args = folders.collect { "'${it}'" }.join(' ')
    def extra_arg = (params.fastq_pattern
        && !params.fastq_pattern.toString().isEmpty())
            ? "--extra-pattern '${params.fastq_pattern}'"
            : ""
    // [S126]: keep the helper's stderr in a file so part_A can relay
    // its warning: lines (a successful task's stderr is not shown);
    // on failure, replay it on stderr so the error reaches the user.
    """
    if ! discover_fastq.py ${extra_arg} ${folder_args} > samples.tsv 2> discovery_stderr.txt ; then
        cat discovery_stderr.txt >&2
        exit 1
    fi
    """
}
