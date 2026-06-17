process trim_metadata_for_mumu {
    // [S41]: keep the amplicon column (col 4) and every sample
    // column (cols 14+). mumu's --otu_table consumes this shape.

    input:
    path table

    output:
    path "${table.baseName}_reduced.table"

    shell:
    '''
    cut -f 4,14- !{table} > !{table.baseName}_reduced.table
    '''
}
