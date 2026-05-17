#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
   break swarm clusters, using sample distribution data
"""

from __future__ import annotations

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2021/03/19"
__version__ = "$Revision: 1.2"

import argparse
import copy
import operator
import os
import re
import sys
from typing import Optional


# *************************************************************************** #
#                                                                             #
#                                  Functions                                  #
#                                                                             #
# *************************************************************************** #

def per_sample_stats_parse(
    per_sample_stats_file: str, percentage: float,
) -> tuple[float, list[dict[str, int]]]:
    """
    Map samples, OTU seeds and stats.
    """
    separator = "\t"
    per_sample_stats = [dict() for i in range(0, 256)]
    number_of_samples = 0
    previous_sample = None

    with open(per_sample_stats_file, "r") as stats_file:
        print("PROGRESS: parsing per-sample stats", file=sys.stderr)
        for line in stats_file:
            line = line.strip().split(separator)
            sample, cloud, mass, seed, seed_abundance = line[0:5]
            index = int(seed[0:2], 16)
            if seed in per_sample_stats[index]:
                per_sample_stats[index][seed] += 1
            else:
                per_sample_stats[index][seed] = 1
            if sample != previous_sample:
                previous_sample = sample
                number_of_samples += 1

    # keep only local seeds present in at least "percentage" of samples
    threshold = percentage * number_of_samples
    seeds = [{k: v for k, v in per_sample_stats[i].items() if v >= threshold}
             for i in range(0, 256)]

    return threshold, seeds


def stats_parse(
    global_stats_file: str,
    threshold: float,
    seeds: list[dict[str, int]],
) -> None:
    """
    Find and eliminate global seeds.
    """
    separator = "\t"
    # 1 - a secondary seed cannot be more abundant than the global
    # seed, otherwise it would have been selected to be the global
    # seed.
    # 2 - a secondary seed cannot have less reads than the threshold
    # value, otherwise it would appear in "percentage" samples with
    # less than 1 read per sample which is not possible.
    # 3 - consequently, a global seed cannot have an abundance smaller
    # than the threshold value.
    with open(global_stats_file, "r") as stats_file:
        print("PROGRESS: parsing stats", file=sys.stderr)
        for line in stats_file:
            line = line.strip().split(separator)
            cloud, mass, seed, seed_abundance, singletons = line[0:5]

            # stats file is sorted by decreasing seed_abundance
            if int(seed_abundance) < threshold:
                break  # no need to read more lines (see point 3)

            # only check clusters with at least two unique sequences,
            # and with enough reads to be hosting a secondary seed
            # that passes our threshold
            if int(cloud) > 1 and int(mass) >= 2 * threshold + int(singletons):
                index = int(seed[0:2], 16)
                # eliminate global seeds from the list of seeds
                if seed in seeds[index]:
                    del seeds[index][seed]

    return None


def swarms_parse(
    swarms_file: str,
    seeds: list[dict[str, int]],
) -> tuple[list[dict[str, int]], list[dict[str, str]]]:
    """
    Map amplicons and abundance values. Only keep clusters with local seeds.
    """
    # Can a global seed be absent from the local seed list? Yes, it is
    # possible: an amplicon can be very abundant in one sample, and a
    # nearly-identical amplicon can be present in many samples but
    # with a lower total abundance. Hence, I can only rely on a search
    # through the "swarms" file to get a full list of clusters that
    # can be cleaved.
    separator = "_|;size=|;? "  # parsing of abundance annotations
    swarms = [dict() for i in range(0, 256)]
    global_seeds = [dict() for i in range(0, 256)]
    seeds_set = {k for d in seeds for k in d.keys()}
    number_of_seeds = len(seeds_set)

    with open(swarms_file, "r") as swarms_file:
        print("PROGRESS: parsing swarms", file=sys.stderr)
        for line in swarms_file:
            line = line.strip()
            amplicons = re.split(separator, line)[0::2]
            abundances = re.split(separator, line)[1::2]
            seed = amplicons[0]
            common = seeds_set & set(amplicons)
            if common:
                number_of_seeds -= len(common)
                for amplicon, abundance in zip(amplicons, abundances):
                    index = int(amplicon[0:2], 16)
                    swarms[index][amplicon] = int(abundance)
                for amplicon in common:
                    index = int(amplicon[0:2], 16)
                    global_seeds[index][amplicon] = seed
            if number_of_seeds == 0:
                break

    return swarms, global_seeds


def struct_parse(
    struct_file: str,
    seeds: list[dict[str, int]],
    global_seeds: list[dict[str, str]],
) -> list[dict[str, set[str]]]:
    """
    carve out sub-clusters.
    """
    separator = "\t"
    global_seeds_set = {v for d in global_seeds for v in d.values()}
    number_of_seeds = sum([len(d) for d in seeds])
    previous_cluster_id = 0
    clusters = dict()
    new_clusters = list()

    with open(struct_file, "r") as struct_file:
        print("PROGRESS: parsing struct", file=sys.stderr)
        for line in struct_file:
            line = line.strip()
            father, son, diffs, cluster_id, steps = line.split(separator)

            # initialize per-cluster parameters
            if int(cluster_id) != previous_cluster_id:
                if clusters:  # save previous results, if any
                    new_clusters.append(clusters)
                previous_cluster_id = int(cluster_id)
                global_seed = father
                has_a_local_seed = (True if global_seed in global_seeds_set
                                    else False)
                clusters = dict()
                if has_a_local_seed:
                    clusters[global_seed] = set([global_seed])

            # stop parsing the file as soon as possible
            if number_of_seeds == 0:
                break

            # skip clusters that do not contain local seeds
            if has_a_local_seed is False:
                continue

            # detect local seeds
            index = int(son[0:2], 16)
            if son in seeds[index]:
                number_of_seeds -= 1
                clusters[son] = set([son])
                continue

            # populate cluster (assuming a father-son link)
            for d in clusters:
                if father in clusters[d]:
                    clusters[d].add(son)
                    break
            else:
                # father is not in the sub-clusters (i.e. an "orphan"
                # created by the grafting process)
                for d in clusters:
                    if son in clusters[d]:
                        clusters[d].add(father)

    return new_clusters


def add_abundance_values(
    new_clusters: list[dict[str, set[str]]],
    swarms: list[dict[str, int]],
) -> list[dict[str, list[tuple[str, int]]]]:
    """
    Add abundance values and sort (deal with a rare case).
    """
    print("PROGRESS: sorting each cluster", file=sys.stderr)
    new_clusters_with_abundance = copy.deepcopy(new_clusters)  # suboptimal
    for i, super_cluster in enumerate(new_clusters):
        for cluster in super_cluster:
            swarm = list()
            for amplicon in super_cluster[cluster]:
                index = int(amplicon[0:2], 16)
                abundance = swarms[index][amplicon]
                swarm.append((amplicon, abundance))
            # sort amplicons by decreasing abundance value and by
            # name (fix a rare bug: a tie leading to wrong seed
            # selection)
            swarm.sort(key=lambda x: (-x[1], x[0]))

            # # make sure original seed stays first if its a tie (does not work) [2025-09-08 lun.]
            # new_seed = swarm[0][0]
            # new_abundance = swarm[0][1]
            # index = int(new_seed[0:2], 16)
            # is_a_local_seed = new_seed in swarms[index]
            # # if not a local_seed, search for one in the entries with the same abundance
            # if not is_a_local_seed:
            #     # find first that is a local seed
            #     for a_tuple in swarm[1:]:
            #         if a_tuple[1] == new_abundance:
            #             candidate = a_tuple[0]
            #             if candidate in swarms[int(candidate[0:2], 16)]:
            #                 j = [a_seed[0] for a_seed in swarm].index(candidate)
            #                 swarm[0], swarm[j] = swarm[j], swarm[0]
            #                 break

            # check if the seed changed after sorting (rare case)
            new_seed = swarm[0][0]
            if new_seed != cluster:
                del new_clusters_with_abundance[i][cluster]
                cluster = new_seed
            new_clusters_with_abundance[i][cluster] = swarm

    return new_clusters_with_abundance


def per_cluster_stats(
    out_stats_file: str,
    new_clusters_with_abundance: list[dict[str, list[tuple[str, int]]]],
) -> list[tuple[int, int, str, int, int, str, str]]:
    """
    Compute per-cluster stats.
    """
    print("PROGRESS: computing per-cluster stats", file=sys.stderr)
    new_stats = list()
    for super_cluster in new_clusters_with_abundance:
        for cluster in super_cluster:
            total_abundance = 0
            singletons = 0
            seed_abundance = super_cluster[cluster][0][1]
            for amplicon, abundance in super_cluster[cluster]:
                total_abundance += abundance
                if abundance == 1:
                    singletons += 1
            # number_of_uniques = 1  # 1. number of unique amplicons
            # total_abundance = 0    # 2. total abundance of amplicons
            # seed_label = 0         # 3. label of the initial seed
            # seed_abundance = 0     # 4. initial seed abundance
            # singletons = 0         # 5. number of amps with an abundance of 1
            # number_of_steps = 0    # 6. number of steps in the cluster
            # number of layers       # 7. columns 6 and 7 are not updated
            new_stats.append((len(super_cluster[cluster]),
                              total_abundance,
                              cluster,
                              seed_abundance,
                              singletons,
                              "0",
                              "0"))
    # sort clusters by increasing number of reads first, then in a
    # second step sort by decreasing number of unique amplicons, and
    # amplicon name (sort by amplicon name first, stable sorting will
    # preserve the order for ties)
    new_stats.sort(key=operator.itemgetter(2))
    new_stats.sort(key=operator.itemgetter(1, 0), reverse=True)
    with open(out_stats_file, "w") as new_stats_file:
        for t in new_stats:
            print(*t, sep="\t", file=new_stats_file)

    return new_stats


def per_cluster_swarms(
    out_swarms_file: str,
    new_clusters_with_abundance: list[dict[str, list[tuple[str, int]]]],
) -> None:
    """
    Compute per-cluster swarms.
    """
    print("PROGRESS: computing per-cluster swarms", file=sys.stderr)
    with open(out_swarms_file, "w") as new_swarms_file:
        for super_cluster in new_clusters_with_abundance:
            for cluster in super_cluster:
                print(*[t[0] + ";size=" + str(t[1]) for t in super_cluster[cluster]],
                      sep=" ", file=new_swarms_file)

    return None


def fasta_parse(
    fasta_file: str,
    out_fasta_file: str,
    new_stats: list[tuple[int, int, str, int, int, str, str]],
) -> None:
    """
    Get seed sequences, update abundances.
    """
    # create a dict of target amplicons and abundances
    fasta = [dict() for i in range(0, 256)]
    for t in new_stats:
        index = int(t[2][0:2], 16)
        fasta[index][t[2]] = [t[1]]
    min_abundance = min([t[3] for t in new_stats])
    # filter the fasta file
    index = None
    separator = ";size="
    with open(fasta_file, "r") as fasta_file:
        print("PROGRESS: parsing fasta file", file=sys.stderr)
        for line in fasta_file:
            if line.startswith(">"):
                amplicon, abundance = line.strip(">;\n").split(separator)
                index = int(amplicon[0:2], 16)
                if int(abundance) < min_abundance:
                    break  # no need to read more lines
            else:
                if amplicon in fasta[index]:
                    fasta[index][amplicon].append(line.strip())

    with open(out_fasta_file, "w") as new_representatives_file:
        for t in new_stats:
            amplicon = t[2]
            index = int(amplicon[0:2], 16)
            try:
                abundance, sequence = fasta[index][amplicon]
            except ValueError:
                print(amplicon, fasta[index][amplicon], min_abundance)
                sys.exit(-1)
            print(">" + amplicon + separator + str(abundance)
                  + "\n" + sequence, sep="", file=new_representatives_file)

    return None


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """
    Parse arguments from command line.
    """
    parser = argparse.ArgumentParser(
        description="break swarm clusters, using sample distribution data.")

    parser.add_argument("--global_stats",
                        dest="global_stats_file",
                        required=True,
                        help="cluster statistics")

    parser.add_argument("--per_sample_stats",
                        dest="per_sample_stats_file",
                        required=True,
                        help="per-sample cluster statistics")

    parser.add_argument("--fasta",
                        dest="fasta_file",
                        required=True,
                        help="amplicon sequences")

    parser.add_argument("--struct",
                        dest="struct_file",
                        required=True,
                        help="internal structure of clusters")

    parser.add_argument("--swarms",
                        dest="swarms_file",
                        required=True,
                        help="list of amplicons per cluster")

    parser.add_argument("--percentage",
                        dest="percentage",
                        type=float,
                        default=0.05,
                        help="cleaving threshold: a sub-seed candidate "
                             "must appear as a per-sample cluster seed "
                             "in at least this fraction of samples "
                             "(default: 0.05, the legacy keystone value)")

    parser.add_argument("--fastidious",
                        dest="fastidious",
                        action=argparse.BooleanOptionalAction,
                        default=True,
                        help="upstream swarm run used --fastidious; sets "
                             "the representatives output filename suffix "
                             "to '_1f_representatives.fas2' "
                             "(default: True, matches the legacy "
                             "pipeline). Use --no-fastidious for '_1_'.")

    parser.add_argument("--out-stats",
                        dest="out_stats_file",
                        default=None,
                        help="output path for the cleaved cluster stats "
                             "(default: <--global_stats>2)")

    parser.add_argument("--out-swarms",
                        dest="out_swarms_file",
                        default=None,
                        help="output path for the cleaved cluster swarms "
                             "(default: <--swarms>2)")

    parser.add_argument("--out-fasta",
                        dest="out_fasta_file",
                        default=None,
                        help="output path for the cleaved cluster "
                             "representatives FASTA (default: "
                             "<basename of --fasta>_<1|1f>_representatives.fas2)")

    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    """
    break clusters and output updated rep, stats and swarms files.
    """
    # capture arguments
    args = parse_args(argv)
    global_stats_file = args.global_stats_file
    per_sample_stats_file = args.per_sample_stats_file
    struct_file = args.struct_file
    swarms_file = args.swarms_file
    fasta_file = args.fasta_file

    # fastidious or not? (explicit flag, default True — see [S22])
    swarm_parameters = "1f" if args.fastidious else "1"

    # resolve output paths (legacy `<input>2` defaults if --out-* unset)
    out_stats_file = args.out_stats_file or global_stats_file + "2"
    out_swarms_file = args.out_swarms_file or swarms_file + "2"
    out_fasta_file = args.out_fasta_file or (
        os.path.splitext(fasta_file)[0]
        + "_"
        + swarm_parameters
        + "_representatives.fas2"
    )

    # cleaving threshold (keystone parameter)
    percentage = args.percentage

    # Parse input files
    threshold, seeds = per_sample_stats_parse(per_sample_stats_file,
                                              percentage)
    stats_parse(global_stats_file, threshold, seeds)
    swarms, global_seeds = swarms_parse(swarms_file, seeds)
    new_clusters = struct_parse(struct_file, seeds, global_seeds)

    # Add abundance values
    new_clusters_with_abundance = add_abundance_values(new_clusters,
                                                       swarms)

    # Create output files (stats2, swarms2, fas2)
    new_stats = per_cluster_stats(out_stats_file,
                                  new_clusters_with_abundance)
    per_cluster_swarms(out_swarms_file, new_clusters_with_abundance)
    fasta_parse(fasta_file, out_fasta_file, new_stats)

    return 0


# *************************************************************************** #
#                                                                             #
#                                     Body                                    #
#                                                                             #
# *************************************************************************** #

if __name__ == '__main__':
    sys.exit(main())
