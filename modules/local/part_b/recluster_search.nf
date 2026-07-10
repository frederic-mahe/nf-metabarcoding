process recluster_search {
    // [S103]/D20: optional, terminal coarse-clustering pass for
    // divergent markers. Runs `vsearch --cluster_size` (abundance-based
    // greedy clustering, AGC) on the post-mumu FASTA
    // (`extract_mumu_fasta`'s output, [S40], header
    // `>amplicon;size=total;`) so `--sizein` reads the abundances
    // directly. The output is the `H` lines of vsearch's `--uc` stream —
    // `recluster_merge` consumes those to fold each member OTU onto its
    // centroid.
    //
    // Only `--id` (--recluster_id) and `--iddef` (--recluster_iddef) are
    // user-facing ([S102]); the rest are fixed: `--cluster_size` is the
    // feature, `--sizein`/`--sizeout` carry the abundances,
    // `--qmask none` matches [S38], and `--maxaccepts 0 --maxrejects 0`
    // (0 = unlimited) buys an exhaustive centroid search — accuracy over
    // speed, like [S42]'s self-search.
    //
    // grep "^H" returns exit 1 when there are no hits; `|| true` keeps
    // the process green (an empty .uc is a legitimate outcome — nothing
    // to collapse).
    //
    // [S105]: vsearch's --log captures the search step's stderr;
    // recluster_merge cats it together with the merge step's stderr to
    // produce <basename>_reclustering.log.

    input:
    path fasta

    output:
    path "${fasta.baseName}.uc", emit: uc
    path "recluster_search.log", emit: log

    shell:
    '''
    set -euo pipefail

    vsearch \
        --cluster_size !{fasta} \
        --sizein \
        --sizeout \
        --id !{params.recluster_id} \
        --iddef !{params.recluster_iddef} \
        --qmask none \
        --fasta_width 0 \
        --strand plus \
        --maxaccepts 0 \
        --maxrejects 0 \
        --threads !{task.cpus} \
        --log recluster_search.log \
        --uc - | \
        grep "^H" > !{fasta.baseName}.uc || true
    '''

    stub:
    """
    touch ${fasta.baseName}.uc recluster_search.log
    """
}
