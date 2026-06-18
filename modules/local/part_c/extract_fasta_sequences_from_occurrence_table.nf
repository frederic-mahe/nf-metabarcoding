process extract_fasta_sequences_from_occurrence_table {
    // [S48]: extract a representatives FASTA from the occurrence
    // table. Column 4 = amplicon ID, column 2 = abundance, column
    // 10 = sequence (same layout as [S40]'s extract_otu_fasta).

    input:
    path occurrence_table

    output:
    path "${occurrence_table.baseName}_representatives.fas", emit: fasta

    shell:
    '''
    awk 'NR > 1 {printf ">"$4";size="$2";\\n"$10"\\n"}' !{occurrence_table} \
        > !{occurrence_table.baseName}_representatives.fas
    '''
}
