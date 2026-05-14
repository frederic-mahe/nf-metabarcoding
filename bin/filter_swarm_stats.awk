#!/usr/bin/awk -f
#
# filter_swarm_stats.awk
#
# Keep only swarm clusters with more than 2 reads (spec [S17]).
#
# Input: swarm --statistics-file format (TSV, column 2 = total reads).
# Output: same rows, same TSV layout, filtered.

BEGIN { FS = OFS = "\t" }
$2 > 2
