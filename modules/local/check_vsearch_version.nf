process check_vsearch_version {
    // [S136] (D25): report the first line of `vsearch --version` (the
    // version is printed on stderr, so both streams are merged). An
    // empty line means vsearch is missing. The comparison to the
    // minimum is done by vsearch_version_error() in functions.nf, and
    // the calling workflow stops the run on an error. Runs as a task
    // (not at parse time) because vsearch may only be visible on the
    // compute nodes, through conda, modules or containers.

    output:
    stdout emit: version

    script:
    """
    { vsearch --version 2>&1 || true ; } | head -n 1
    """

    stub:
    """
    echo "vsearch v2.32.0_stub"
    """
}
