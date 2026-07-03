include { publish_dir } from '../functions.nf'


process download_run {
    // [S101]: download a single run accession's fastq file(s) via
    // fastq-dl, one task per run so the failure contract is per-run: a
    // run that fails to download fails only its own task; the runs that
    // succeeded are published and cached, so a re-run with -resume
    // retries only the failed runs. --provider ena requests the fastq
    // directly from ENA (no SRA-Lite conversion).
    //
    // [S101]: fastq-dl is pinned as a per-process conda directive (never
    // in environment.yml) so a run that never touches the fetch entry
    // never resolves it, and Wave builds a fetch-only image on demand.
    conda 'bioconda::fastq-dl=4.0.1'

    // [S100]: publish under a per-accession subfolder so runs from
    // distinct accessions in a comma-separated --accession list never
    // collide in a shared directory.
    publishDir path: { publish_dir(accession) }, mode: params.publish_mode, pattern: "*.fastq.gz"

    tag "${accession}/${run}"

    input:
    tuple val(accession), val(run)

    output:
    tuple val(accession), path("*.fastq.gz"), emit: reads

    // [S101]: raise fastq-dl's retry ceiling from its default of 3 to 5
    // attempts, so a transient ENA / network hiccup is retried more times
    // before the run-task (and only that task) fails.
    script:
    """
    set -euo pipefail

    fastq-dl \\
        --accession ${run} \\
        --provider ena \\
        --max-attempts 5
    """

    stub:
    """
    touch ${run}_1.fastq.gz ${run}_2.fastq.gz
    """
}
