process validate_samplesheet {
    // [S70]: structurally validate the --input samplesheet and emit the
    // normalized TSV (one row per sample, profile implicit in the
    // header columns). bin/parse_samplesheet.py does the structural
    // checks (columns, duplicates, reserved suffix, single-end
    // inference, sibling qual/stats defaulting); file existence is
    // enforced by the caller's file(checkIfExists: true) when each row
    // is staged. The samplesheet is a declared path input, so Nextflow
    // stages it and -resume tracks changes to it.

    input:
    path samplesheet

    output:
    path "samplesheet.normalized.tsv"

    script:
    """
    parse_samplesheet.py ${samplesheet} > samplesheet.normalized.tsv
    """
}
