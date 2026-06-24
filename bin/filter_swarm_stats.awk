#!/usr/bin/awk -f
#
# filter_swarm_stats.awk
#
# Keep only swarm clusters with strictly more than `min_cluster_size`
# reads (spec [S17]). `min_cluster_size` is supplied with `-v` by the
# caller (list_local_clusters); when unset it defaults to 2, the legacy
# value, so the script stays usable standalone.
#
# Input: swarm --statistics-file format (TSV, column 2 = total reads).
# Output: same rows, same TSV layout, filtered.

BEGIN { FS = OFS = "\t"; if (min_cluster_size == "") min_cluster_size = 2 }
$2 > min_cluster_size
