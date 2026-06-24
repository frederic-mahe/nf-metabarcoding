#!/usr/bin/env bats
#
# Unit tests for bin/filter_swarm_stats.awk — the per-sample swarm
# statistics filter ([S17]). It keeps only rows whose column 2 (total
# reads) is strictly greater than `min_cluster_size`, defaulting to 2
# (the legacy "> 2 reads" rule) when the variable is not supplied.
#
# Input is swarm --statistics-file format (TSV, column 2 = total reads);
# the caller (list_local_clusters) passes the threshold with
# `-v min_cluster_size=`.
#
# COVERAGE: [S17]
#
# Run from the repository root:
#   bats tests/bin/

setup() {
    REPO_ROOT="$(git rev-parse --show-toplevel)"
    SCRIPT="${REPO_ROOT}/bin/filter_swarm_stats.awk"
}

# Three clusters with totals 1, 3, 5 in column 2 (swarm stats layout).

@test "default threshold keeps only clusters with more than 2 reads" {
    run bash -c "printf 'c1\t1\ts\nc2\t3\ts\nc3\t5\ts\n' | '${SCRIPT}'"
    [ "${status}" -eq 0 ]
    [ "${output}" = "$(printf 'c2\t3\ts\nc3\t5\ts')" ]
}

@test "min_cluster_size=0 keeps every cluster (> 0 reads)" {
    run bash -c "printf 'c1\t1\ts\nc2\t3\ts\nc3\t5\ts\n' | '${SCRIPT}' -v min_cluster_size=0"
    [ "${status}" -eq 0 ]
    [ "${output}" = "$(printf 'c1\t1\ts\nc2\t3\ts\nc3\t5\ts')" ]
}

@test "min_cluster_size=4 drops clusters at or below 4 reads" {
    run bash -c "printf 'c1\t1\ts\nc2\t3\ts\nc3\t5\ts\n' | '${SCRIPT}' -v min_cluster_size=4"
    [ "${status}" -eq 0 ]
    [ "${output}" = "$(printf 'c3\t5\ts')" ]
}

@test "the TSV layout is preserved verbatim for surviving rows" {
    run bash -c "printf 'c2\t3\tseed2\t2\textra\n' | '${SCRIPT}' -v min_cluster_size=2"
    [ "${status}" -eq 0 ]
    [ "${output}" = "$(printf 'c2\t3\tseed2\t2\textra')" ]
}
