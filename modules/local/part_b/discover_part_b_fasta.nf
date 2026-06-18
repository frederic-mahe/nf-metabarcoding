include { normalize_path } from '../functions.nf'


process discover_part_b_fasta {
    // [S27]/[S56]: walk params.fasta_folder, assert unique sample IDs,
    // and emit two TSVs:
    //   - fastas.tsv         — regular samples (excludes *_notmerged.fas)
    //   - shadow_fastas.tsv  — only *_notmerged.fas samples (shadow Part B)
    // Each row is `sample_id<TAB>fasta_path`. Either file may be empty
    // (e.g. shadow_fastas.tsv is empty when no _notmerged.fas exists).

    output:
    path "fastas.tsv",        emit: regular
    path "shadow_fastas.tsv", emit: shadow

    script:
    def folders = (params.fasta_folder instanceof List)
        ? params.fasta_folder.collect { normalize_path(it) }
        : params.fasta_folder
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it }
            .collect { normalize_path(it) }
    def folder_args = folders.collect { "'${it}'" }.join(' ')
    """
    discover_fasta.py          ${folder_args} > fastas.tsv
    discover_fasta.py --shadow ${folder_args} > shadow_fastas.tsv
    """
}
