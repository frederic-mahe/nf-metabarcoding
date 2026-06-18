include { publish_dir } from '../functions.nf'


process global_dereplication {
    // [S31]: cat every input .fas and pass it through
    // vsearch --fastx_uniques. --sizein/--sizeout preserve the
    // per-sample size annotations so vsearch sums abundances across
    // samples. The vsearch log lands at <basename>_dereplication.log
    // ([S45]).
    //
    // [S59]: only the log reaches the results folder; the
    // dereplicated .fas is an internal intermediate.
    publishDir path: { publish_dir('occurrence_table') }, mode: params.publish_mode, pattern: "*.log"

    input:
    path fastas
    val basename

    output:
    path "${basename}.fas",                 emit: fasta
    path "${basename}_dereplication.log",   emit: log

    shell:
    '''
    cat !{fastas} | \
        vsearch \
            --fastx_uniques - \
            --sizein \
            --sizeout \
            --fasta_width 0 \
            --quiet \
            --log !{basename}_dereplication.log \
            --fastaout !{basename}.fas
    '''
}
