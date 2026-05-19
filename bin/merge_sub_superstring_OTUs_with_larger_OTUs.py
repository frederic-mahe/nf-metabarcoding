#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    merge super or sub-string OTUs with their source-OTUs.
"""

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2019/09/20"
__version__ = "$Revision: 1.2"

import sys
import csv
from optparse import OptionParser


# *************************************************************************** #
#                                                                             #
#                                  Functions                                  #
#                                                                             #
# *************************************************************************** #

def option_parse():
    """
    Parse arguments from command line.
    """
    parser = OptionParser(usage="usage: %prog -t filename -m filename",
                          version="%prog 1.0")

    parser.add_option("-t", "--table",
                      metavar="TABLE",
                      action="store",
                      dest="table",
                      help="set TABLE as input")

    parser.add_option("-m", "--matches",
                      metavar="MATCHES",
                      action="store",
                      dest="matches",
                      help="set MATCHES as input")

    parser.add_option("-o", "--output",
                      metavar="OUTPUT",
                      action="store",
                      dest="output_table",
                      help="set OUTPUT as output")

    (options, args) = parser.parse_args()

    return options.table, options.matches, options.output_table


def main():
    """
    merge super or sub-string OTUs with their source-OTUs.
    """
    input_table, matches, output_table = option_parse()
    # output_table = os.path.splitext(input_table)[0] + "_gene_counts.tsv"

    with open(matches, "r") as matches:
        OTU_connexions = dict()
        for match in matches:
            pupil, master = match.strip().split("\t")[-2:]
            OTU_connexions[pupil] = master

    # Could there be common OTUs? given the clustering process, I
    # don't think so, but I can systematically test for that.
    pupils = set(OTU_connexions.keys())
    masters = set(OTU_connexions.values())
    if masters.intersection(pupils):
        print("WARNING: there are OTUs common to hit and query columns",
              file=sys.stderr)
        print(masters.intersection(pupils), file=sys.stderr)
        sys.exit(-1)

    # Parse the OTU table
    metadata = set(["OTU", "total", "cloud", "amplicon", "length",
                    "abundance", "chimera", "spread", "quality", "sequence",
                    "identity", "taxonomy", "references"])
    with open(input_table, "r") as input_table:
        reader = csv.DictReader(input_table, delimiter="\t")
        fieldnames = reader.fieldnames
        sample_names = set(fieldnames) - metadata
        master_OTUs = dict()
        with open(output_table, "w") as output_table:
            writer = csv.DictWriter(output_table,
                                    fieldnames=reader.fieldnames,
                                    delimiter="\t")
            writer.writeheader()
            for row in reader:
                OTU = row["OTU"]
                if OTU not in masters and OTU not in pupils:
                    writer.writerow(row)
                else:
                    if OTU in masters:
                        # store it
                        master_OTUs[OTU] = row
                        for sample in sample_names:
                            master_OTUs[OTU][sample] = int(row[sample])
                            master_OTUs[OTU]["spread"] = int(row["spread"])
                            master_OTUs[OTU]["total"] = int(row["total"])
                            master_OTUs[OTU]["cloud"] = int(row["cloud"])
                    if OTU in pupils:
                        # who is its master?
                        master = OTU_connexions[OTU]
                        # update master's occurrences
                        for sample in sample_names:
                            master_OTUs[master][sample] += int(row[sample])
                            # update spread
                        new_spread = len([master_OTUs[master][sample]
                                          for sample in sample_names
                                          if master_OTUs[master][sample] > 0])
                        master_OTUs[master]["spread"] = new_spread
                        # update total
                        master_OTUs[master]["total"] += int(row["total"])
                        # update cloud
                        master_OTUs[master]["cloud"] += int(row["cloud"]) + 1
                        # output the updated masters
            for OTU in master_OTUs:
                writer.writerow(master_OTUs[OTU])

    return None


# *************************************************************************** #
#                                                                             #
#                                     Body                                    #
#                                                                             #
# *************************************************************************** #

if __name__ == '__main__':

    main()

sys.exit(0)
