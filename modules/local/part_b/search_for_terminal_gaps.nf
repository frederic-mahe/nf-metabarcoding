process search_for_terminal_gaps {
    // [S38]: self-cluster the OTU table at id=1.0 to find OTU pairs
    // that are identical modulo terminal gaps (sub-strings or
    // super-strings). The output is the `H` lines of vsearch's
    // `--uc` stream — `merge_substring_otus` consumes those to
    // collapse pupil OTUs onto their masters.
    //
    // grep "^H" returns exit 1 when there are no hits; `|| true`
    // keeps the process green (an empty .uc is a legitimate outcome).
    //
    // [S45]: vsearch's --log captures the search step's stderr; the
    // merge_substring_otus process cats it together with the merge
    // step's stderr to produce <basename>_superstring_clustering.log.

    input:
    path otu_table

    output:
    path "${otu_table.baseName}.uc"
    path "search.log"

    shell:
    '''
    awk 'NR > 1 {printf ">"$1"\\n"$10"\\n"}' !{otu_table} | \
        vsearch \
            --threads !{task.cpus} \
            --cluster_smallmem - \
            --id 1.0 \
            --qmask none \
            --usersort \
            --log search.log \
            --uc - | \
        grep "^H" > !{otu_table.baseName}.uc || true
    '''
}
