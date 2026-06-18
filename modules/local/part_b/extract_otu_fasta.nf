process extract_otu_fasta {
    // [S40]: emit a FASTA from an OTU table (every data row).
    // Header `<amplicon>;size=<total>;`; column 4 is the amplicon,
    // column 2 the total abundance, column 10 the sequence.
    //
    // [S59]: pre-mumu fasta is an internal intermediate consumed by
    // find_similar_sequences (mumu match list) — not published.
    // The post-mumu sibling `extract_mumu_fasta` produces the one
    // FASTA that reaches the results folder.

    input:
    path table

    output:
    path "${table.baseName}.fas", emit: fasta

    shell:
    '''
    awk 'NR > 1 {printf ">"$4";size="$2";\\n"$10"\\n"}' !{table} \
        > !{table.baseName}.fas
    '''
}
