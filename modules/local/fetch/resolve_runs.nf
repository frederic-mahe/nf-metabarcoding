process resolve_runs {
    // [S99]: resolve one ENA/SRA accession (bioproject or study) to its
    // constituent run accessions before any download, so the download
    // stage can fan out one task per run ([S101]). fastq-dl's
    // --only-download-metadata writes a `fastq-run-info.tsv` whose
    // `run_accession` column lists the runs; we extract that column
    // (minus the header) into <accession>_runs.tsv.
    //
    // [S101]: fastq-dl is pinned as a per-process conda directive (never
    // in environment.yml) so a run that never touches the fetch entry
    // never resolves it, and Wave builds a fetch-only image on demand.
    conda 'bioconda::fastq-dl=4.0.1'

    tag "${accession}"

    input:
    val accession

    output:
    tuple val(accession), path("${accession}_runs.tsv"), emit: runs

    script:
    """
    set -euo pipefail

    fastq-dl \\
        --accession ${accession} \\
        --provider ena \\
        --only-download-metadata

    # extract the run_accession column (skip the header row); tolerate a
    # metadata file whose columns are in any order by keying on the header.
    awk -F '\\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) if (\$i == "run_accession") col = i
            next
        }
        col { print \$col }
    ' fastq-run-info.tsv > ${accession}_runs.tsv
    """

    stub:
    """
    printf '%s\\n' ERR13334567 ERR13334568 > ${accession}_runs.tsv
    """
}
