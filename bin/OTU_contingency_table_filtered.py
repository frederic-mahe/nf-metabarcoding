#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge different results into a sorted and filtered OTU contingency table.
"""

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2020/01/14"
__version__ = "$Revision: 4.0"

import re
import sys
from argparse import ArgumentParser
import operator


# *************************************************************************** #
#                                                                             #
#                                  Functions                                  #
#                                                                             #
# *************************************************************************** #
def arg_parse():
    """
    Parse arguments from command line.
    """

    parser = ArgumentParser()

    parser.add_argument("-r", "--representatives",
                        action="store",
                        dest="representatives",
                        required=True)

    parser.add_argument("-s", "--stats",
                        action="store",
                        dest="stats",
                        required=True)

    parser.add_argument("-sw", "--swarms",
                        action="store",
                        dest="swarms",
                        required=True)

    parser.add_argument("-c", "--chimera",
                        action="store",
                        dest="chimera",
                        required=True)

    parser.add_argument("-q", "--quality",
                        action="store",
                        dest="quality",
                        required=True)

    parser.add_argument("-a", "--assignments",
                        action="store",
                        dest="assignments",
                        required=True)

    parser.add_argument("-d", "--distribution",
                        action="store",
                        dest="distribution",
                        required=True)
    args = parser.parse_args()

    return args.representatives, args.stats, args.swarms, \
        args.chimera, args.quality, args.assignments, args.distribution


def representatives_parse(stampa, repre):
    # "${REPRESENTATIVES}" \
    """
    Get seed sequences.
    """
    separator = ";size="
    representatives_file = repre
    representatives = [dict() for i in range(0, 256)]
    with open(representatives_file, "r") as representatives_file:
        print("PROGRESS: parsing fasta representatives", file=sys.stderr)
        for line in representatives_file:
            if line.startswith(">"):
                amplicon = line.strip(">;\n").split(separator)[0]
                index = int(amplicon[0:2], 16)
            else:
                # discard small OTUs not needed in the final table
                if amplicon in stampa[index]:
                    representatives[index][amplicon] = line.strip()

    return representatives


def stats_parse(representatives, stat):
    # "${STATS}" \
    """
    Map OTU seeds and stats.
    """
    separator = "\t"
    stats_file = stat
    stats = dict()
    seeds = dict()
    with open(stats_file, "r") as stats_file:
        print("PROGRESS: parsing stats", file=sys.stderr)
        for line in stats_file:
            line = line.strip().split(separator)
            cloud, mass, seed, seed_abundance = line[0:4]
            index = int(seed[0:2], 16)
            if seed in representatives[index]:
                stats[seed] = int(mass)
                seeds[seed] = (int(seed_abundance), int(cloud))
    # Sort OTUs by decreasing mass
    sorted_stats = sorted(iter(stats.items()),
                          key=operator.itemgetter(1, 0))
    sorted_stats.reverse()

    return stats, sorted_stats, seeds


def swarms_parse(representatives, swarm):
    # "${SWARMS}" \
    """
    Map OTUs.
    """
    separator = "_[0-9]+|;size=[0-9]+;?| "  # parsing of abundance annotations
    swarms_file = swarm
    swarms = [dict() for i in range(0, 256)]
    valid_OTUs = [dict() for i in range(0, 256)]
    with open(swarms_file, "r") as swarms_file:
        print("PROGRESS: parsing swarms", file=sys.stderr)
        for line in swarms_file:
            line = line.strip()
            amplicons = re.split(separator, line)[0::2]
            seed = amplicons[0]
            index = int(seed[0:2], 16)
            if seed in representatives[index]:
                swarms[index][seed] = [amplicons]
                for amplicon in amplicons:
                    index = int(amplicon[0:2], 16)
                    valid_OTUs[index][amplicon] = seed

    return swarms, valid_OTUs


def uchime_parse(representatives, chime):
    # "${UCHIME}" \
    """
    Map OTU's chimera status.
    """
    separator = "\t"
    uchime_file = chime
    uchime = dict()
    with open(uchime_file, "r") as uchime_file:
        print("PROGRESS: parsing uchime", file=sys.stderr)
        for line in uchime_file:
            OTU = line.strip().split(separator)
            try:
                seed = OTU[1].split(";")[0]
                index = int(seed[0:2], 16)
            except IndexError:  # deal with partial line (missing seed)
                continue
            try:
                status = OTU[17]
            except IndexError:  # deal with unfinished chimera detection runs
                status = "NA"
            if seed in representatives[index]:
                uchime[seed] = status

    return uchime


def quality_parse(representatives, qual):
    # "${QUALITY}" \
    """
    List good amplicons.
    """
    quality_file = qual
    quality = dict()
    with open(quality_file, "r") as quality_file:
        print("PROGRESS: parsing amplicon quality (EE)", file=sys.stderr)
        for line in quality_file:
            sha1, qual, length = line.strip().split()
            index = int(sha1[0:2], 16)
            if sha1 in representatives[index]:
                quality[sha1] = float(qual) / int(length)

    return quality


def stampa_parse(assign):
    # "${ASSIGNMENTS}" \
    """
    Map amplicon ids and taxonomic assignments.
    """
    separator = "\t"
    stampa_file = assign
    stampa = [dict() for i in range(0, 256)]
    
    with open(stampa_file, "r") as stampa_file:
        print("PROGRESS: parsing taxonomic assignments", file=sys.stderr)
        for line in stampa_file:
            line = line.strip().split(separator)
            amplicon, abundance, identity, taxonomy, references = line
            # remove rare but annoying character
            taxonomy = taxonomy.replace("#", "")
            index = int(amplicon[0:2], 16)
            stampa[index][amplicon] = (identity, taxonomy, references)

    return stampa


def distribution_parse(valid_OTUs, distr):
    # "${DISTRIBUTION}"
    """
    Map amplicon ids, abundances and samples.
    """
    distr_file = distr
    samples = dict()
    # initialize the seed to sample dict
    seeds2samples = [dict() for i in range(0, 256)]
    seeds = list()
    for i in range(0, 256):
        seeds += list(valid_OTUs[i].values())
    for seed in set(seeds):
        index = int(seed[0:2], 16)
        seeds2samples[index][seed] = dict()

    with open(distr_file, "r") as distr_file:
        print("PROGRESS: parsing distribution file", file=sys.stderr)
        for line in distr_file:
            amplicon, sample, abundance = line.strip().split("\t")
            # deal with duplicated samples
            samples[sample] = samples.get(sample, 0) + 1
            index = int(amplicon[0:2], 16)
            if amplicon in valid_OTUs[index]:
                abundance = int(abundance)
                # update the seed distribution directly
                seed = valid_OTUs[index][amplicon]
                index = int(seed[0:2], 16)
                seeds2samples[index][seed][sample] = (
                    seeds2samples[index][seed].get(sample, 0)
                    + abundance)
        samples = sorted(samples.keys())

    return seeds2samples, samples


def print_table(representatives, stats, sorted_stats,
                swarms, uchime, seeds2samples,
                samples, quality, seeds, stampa):
    """
    Export results.
    """
    print("PROGRESS: filtering and writing OTUs", file=sys.stderr)
    # Print table header
    print("OTU", "total", "cloud",
          "amplicon", "length", "abundance",
          "chimera", "spread", "quality",
          "sequence", "identity", "taxonomy", "references",
          "\t".join(samples),
          sep="\t", file=sys.stdout)

    # Print table content
    i = 1
    for seed, abundance in sorted_stats:
        index = int(seed[0:2], 16)
        sequence = representatives[index][seed]
        occurrences = dict([(sample, 0) for sample in samples])
        ## temporary hotfix! (not good, produces a zero-abundance OTU)
        try:
            occurrences.update(seeds2samples[index][seed])
        except KeyError:
            # In rare cases, the cleaving step can change the seed of
            # a cluster (for example, when one or more amplicons have
            # the same abundance than the seed, but are lower in the
            # alphabetical order). The original cluster still exists
            # in the first "swarm", but appears with another seed in
            # "swarm2". When these files are parsed (in order), a list
            # of valid_otus is populated with the original seed, and
            # then updated to point to the new seed. The stats and
            # sorted_stats list also contain both the original and new
            # seeds. However, downstream mappings such as
            # seeds2samples, do not contain the original seed, only
            # the new one.  When the script tries to print the
            # occurrence table, sorted by decreasing abundance, the
            # original seed is not present in seeds2samples and
            # triggers an error. The solution is to skip the original
            # seed, as the cluster is already represented by the new
            # seed. Not skipping creates an empty cluster (zero reads)
            # with the original seed. Empty clusters should be
            # discarded by downstream analyses, but skipping is much
            # cleaner.
            continue
            # print("KeyError:", seed, "not in seeds2samples", sep=" ", file=sys.stderr)
        spread = len([occurrences[sample] for sample in samples
                      if occurrences[sample] > 0])
        sequence_abundance, cloud = seeds[seed]

        # Quality (note: more digits with python 3)
        if seed in quality:
            high_quality = quality[seed]
        else:
            high_quality = "NA"

        # Chimera checking (deal with incomplete cases. Is it useful?)
        if seed in uchime:
            chimera_status = uchime[seed]
        else:
            chimera_status = "NA"

        # Taxonomic assignment
        if seed in stampa[index]:
            identity, taxonomy, references = stampa[index][seed]
        else:
            identity, taxonomy, references = "NA", "NA", "NA"

        # Apply filters
        max_ee = 0.0002
        if (chimera_status == "N" and
                high_quality <= max_ee
                and (abundance >= 3 or spread >= 2)):
            print(i, abundance, cloud,
                  seed, len(sequence), sequence_abundance,
                  chimera_status, spread, high_quality, sequence,
                  identity, taxonomy, references,
                  "\t".join([str(occurrences[sample]) for sample in samples]),
                  sep="\t", file=sys.stdout)
            i += 1

    return


def main():
    """
    Read swarm files and build a sorted OTU contingency table.
    """
    # Parse arguments from command line
    repre, stat, swarm, chime, qual, assign, distr = arg_parse()

    # Parse taxonomic assignment results (i.e. valid OTUs for the final table)
    stampa = stampa_parse(assign)

    # Parse OTU representatives
    representatives = representatives_parse(stampa, repre)

    # Parse OTU stats
    stats, sorted_stats, seeds = stats_parse(representatives, stat)

    # Parse OTUs (swarms)
    swarms, valid_OTUs = swarms_parse(representatives, swarm)

    # Parse chimera detection results (uchime)
    uchime = uchime_parse(representatives, chime)

    # Parse sequence's best error rates (a.k.a. quality)
    quality = quality_parse(representatives, qual)

    # Parse distribution file
    seeds2samples, samples = distribution_parse(valid_OTUs, distr)

    # Print table header
    print_table(representatives, stats, sorted_stats, swarms,
                uchime, seeds2samples, samples, quality,
                seeds, stampa)

    return


# *************************************************************************** #
#                                                                             #
#                                     Body                                    #
#                                                                             #
# *************************************************************************** #

if __name__ == '__main__':

    main()

sys.exit(0)
