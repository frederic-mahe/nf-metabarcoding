include { publish_dir } from './functions.nf'


process dump_software_versions {
    // [S68]: record the versions of the external tools (and the Python
    // interpreter) the run relied on, for reproducibility. Each tool's
    // raw `--version` output is piped through bin/collect_versions.py,
    // which extracts the version token; a tool missing from the active
    // environment (PATH / conda / module / container) is recorded as
    // `n/a` rather than dropped, so the gap is visible in the report.
    // [S71]: published to <outdir>/pipeline_info/ (a sibling of
    // occurrence_table/, not inside it).
    publishDir path: { publish_dir('pipeline_info') }, mode: params.publish_mode

    output:
    path "software_versions.yml"

    script:
    """
    set -euo pipefail

    {
        printf 'vsearch\\t%s\\n'  "\$(vsearch  --version 2>&1 | head -n 1)"
        printf 'swarm\\t%s\\n'    "\$(swarm    --version 2>&1 | head -n 1)"
        printf 'cutadapt\\t%s\\n' "\$(cutadapt --version 2>&1 | head -n 1)"
        printf 'mumu\\t%s\\n'     "\$(mumu     --version 2>&1 | head -n 1)"
        printf 'python\\t%s\\n'   "\$(python3  --version 2>&1 | head -n 1)"
    } | collect_versions.py > software_versions.yml
    """
}
