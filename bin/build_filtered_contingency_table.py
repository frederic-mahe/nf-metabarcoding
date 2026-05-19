#!/usr/bin/env python3
"""Merge swarm/uchime/quality/stampa/distribution into a filtered OTU table.

This script consumes the Part B inputs produced by `global_clustering`,
`chimera_detection`, `build_expected_error_file`, `fake_taxonomic_assignment`
(or stampa, once Part C runs), and `build_distribution_file`, and emits
the filtered occurrence table to stdout.

Filter rule: a cluster appears in the output iff

    chimera == "N" AND ee/length <= 0.0002 AND
    (abundance >= 3 OR spread >= 2)

Sample columns are sorted; empty samples contribute a zero-filled
column (see ``[S09]``).

Characterization tests in
``tests/python/test_build_filtered_contingency_table.py`` pin the
byte-exact output — refactors must keep that test green.
"""

from __future__ import annotations

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2020/01/14"
__version__ = "$Revision: 5.0"

import argparse
import operator
import re
import sys
from typing import Iterator, Optional


# Filter constants.
MAX_EE: float = 0.0002


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse arguments from the command line."""
    parser = argparse.ArgumentParser(
        description=(
            "Merge swarm/uchime/quality/stampa/distribution into a "
            "filtered OTU contingency table."
        ),
    )
    parser.add_argument("-r", "--representatives",
                        required=True,
                        help="swarm --seeds FASTA")
    parser.add_argument("-s", "--stats",
                        required=True,
                        help="swarm --statistics-file output")
    parser.add_argument("-sw", "--swarms",
                        required=True,
                        help="swarm --output-file")
    parser.add_argument("-c", "--chimera",
                        required=True,
                        help="vsearch --uchime_denovo --uchimeout")
    parser.add_argument("-q", "--quality",
                        required=True,
                        help="per-amplicon ee/length TSV "
                             "(build_expected_error_file output)")
    parser.add_argument("-a", "--assignments",
                        required=True,
                        help="taxonomic assignments (stampa output)")
    parser.add_argument("-d", "--distribution",
                        required=True,
                        help="sequence-to-sample distribution TSV "
                             "(build_distribution_file output)")
    return parser.parse_args(argv)


def stampa_parse(
    path: str,
) -> dict[str, tuple[str, str, str]]:
    """Map ``amplicon -> (identity, taxonomy, references)``."""
    print("PROGRESS: parsing taxonomic assignments", file=sys.stderr)
    stampa: dict[str, tuple[str, str, str]] = {}
    with open(path) as f:
        for line in f:
            amplicon, _abundance, identity, taxonomy, references = (
                line.strip().split("\t")
            )
            # remove rare but annoying character
            taxonomy = taxonomy.replace("#", "")
            stampa[amplicon] = (identity, taxonomy, references)
    return stampa


def representatives_parse(
    path: str, stampa: dict[str, tuple[str, str, str]],
) -> dict[str, str]:
    """Map ``seed -> sequence`` for every seed present in ``stampa``.

    Seeds absent from ``stampa`` are dropped: every downstream parser
    is keyed on ``representatives``, so the stampa filter cascades to
    them automatically.
    """
    print("PROGRESS: parsing fasta representatives", file=sys.stderr)
    representatives: dict[str, str] = {}
    amplicon: Optional[str] = None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                amplicon = line.strip(">;\n").split(";size=")[0]
            else:
                if amplicon is not None and amplicon in stampa:
                    representatives[amplicon] = line.strip()
    return representatives


def stats_parse(
    path: str, representatives: dict[str, str],
) -> tuple[
    dict[str, int],
    list[tuple[str, int]],
    dict[str, tuple[int, int]],
]:
    """Map seeds -> stats, returning ``(stats, sorted_stats, seeds)``.

    * ``stats`` — seed → total abundance (mass)
    * ``sorted_stats`` — ``[(seed, mass), ...]`` sorted by mass
      descending (and by seed descending as a tiebreaker, matching the
      legacy ``reverse()``-after-ascending-sort behaviour)
    * ``seeds`` — seed → ``(seed_abundance, cloud)``
    """
    print("PROGRESS: parsing stats", file=sys.stderr)
    stats: dict[str, int] = {}
    seeds: dict[str, tuple[int, int]] = {}
    with open(path) as f:
        for line in f:
            cols = line.strip().split("\t")
            cloud, mass, seed, seed_abundance = cols[0:4]
            if seed in representatives:
                stats[seed] = int(mass)
                seeds[seed] = (int(seed_abundance), int(cloud))
    sorted_stats = sorted(stats.items(), key=operator.itemgetter(1, 0))
    sorted_stats.reverse()
    return stats, sorted_stats, seeds


def swarms_parse(
    path: str, representatives: dict[str, str],
) -> tuple[dict[str, list[list[str]]], dict[str, str]]:
    """Parse ``swarm --output-file``.

    Returns ``(swarms, valid_OTUs)`` where ``valid_OTUs`` maps every
    amplicon (seed or not) in a representative-resident cluster to its
    seed, and ``swarms`` keeps the raw amplicon list per seed.
    """
    print("PROGRESS: parsing swarms", file=sys.stderr)
    # Regex matches abundance suffixes and the inter-amplicon space.
    separator = "_[0-9]+|;size=[0-9]+;?| "
    swarms: dict[str, list[list[str]]] = {}
    valid_otus: dict[str, str] = {}
    with open(path) as f:
        for line in f:
            amplicons = re.split(separator, line.strip())[0::2]
            seed = amplicons[0]
            if seed in representatives:
                swarms[seed] = [amplicons]
                for amplicon in amplicons:
                    valid_otus[amplicon] = seed
    return swarms, valid_otus


def uchime_parse(
    path: str, representatives: dict[str, str],
) -> dict[str, str]:
    """Map ``seed -> chimera status``.

    Tolerates two malformed-row shapes the legacy script handled:

    * a row with no second column (``IndexError`` on ``OTU[1]``) — skipped
    * a row with fewer than 18 columns (``IndexError`` on ``OTU[17]``) —
      status falls back to ``"NA"``
    """
    print("PROGRESS: parsing uchime", file=sys.stderr)
    uchime: dict[str, str] = {}
    with open(path) as f:
        for line in f:
            cols = line.strip().split("\t")
            try:
                seed = cols[1].split(";")[0]
            except IndexError:
                continue
            try:
                status = cols[17]
            except IndexError:
                status = "NA"
            if seed in representatives:
                uchime[seed] = status
    return uchime


def quality_parse(
    path: str, representatives: dict[str, str],
) -> dict[str, float]:
    """Map ``seed -> ee/length`` for every seed in ``representatives``."""
    print("PROGRESS: parsing amplicon quality (EE)", file=sys.stderr)
    quality: dict[str, float] = {}
    with open(path) as f:
        for line in f:
            sha1, ee, length = line.strip().split()
            if sha1 in representatives:
                quality[sha1] = float(ee) / int(length)
    return quality


def distribution_parse(
    path: str, valid_otus: dict[str, str],
) -> tuple[dict[str, dict[str, int]], list[str]]:
    """Aggregate amplicon abundances onto their seed, per sample.

    Returns ``(seeds2samples, samples)`` where:

    * ``seeds2samples[seed][sample]`` is the summed abundance of every
      amplicon that resolves to ``seed`` in sample ``sample``; a
      duplicate ``(amplicon, sample)`` row in the distribution simply
      adds to the running total ("deal with duplicated samples" in the
      legacy code).
    * ``samples`` is the sorted list of every sample name seen in the
      distribution file, even for amplicons that are not in
      ``valid_otus`` — so a sample that lives only in filtered-out
      clusters still surfaces as a zero column.
    """
    print("PROGRESS: parsing distribution file", file=sys.stderr)
    seeds2samples: dict[str, dict[str, int]] = {
        seed: {} for seed in set(valid_otus.values())
    }
    sample_set: set[str] = set()
    with open(path) as f:
        for line in f:
            amplicon, sample, abundance = line.strip().split("\t")
            sample_set.add(sample)
            if amplicon in valid_otus:
                seed = valid_otus[amplicon]
                seeds2samples[seed][sample] = (
                    seeds2samples[seed].get(sample, 0) + int(abundance)
                )
    return seeds2samples, sorted(sample_set)


def iter_rows(
    representatives: dict[str, str],
    sorted_stats: list[tuple[str, int]],
    seeds: dict[str, tuple[int, int]],
    seeds2samples: dict[str, dict[str, int]],
    samples: list[str],
    quality: dict[str, float],
    uchime: dict[str, str],
    stampa: dict[str, tuple[str, str, str]],
) -> Iterator[tuple[object, ...]]:
    """Yield one output row per cluster that passes the filter.

    The cleaving step (``[S22]``) can rewrite a cluster's seed; the
    original seed then survives in ``stats`` / ``sorted_stats`` but is
    absent from ``seeds2samples``. Skipping it (rather than emitting a
    zero-abundance row) mirrors the legacy script's
    ``except KeyError: continue`` branch — the rewritten cluster is
    already represented by the new seed.
    """
    i = 1
    for seed, mass in sorted_stats:
        sequence = representatives[seed]
        if seed not in seeds2samples:
            continue
        occurrences = {sample: 0 for sample in samples}
        occurrences.update(seeds2samples[seed])
        spread = sum(1 for sample in samples if occurrences[sample] > 0)
        sequence_abundance, cloud = seeds[seed]

        high_quality: object = quality.get(seed, "NA")
        chimera_status = uchime.get(seed, "NA")
        identity, taxonomy, references = stampa.get(
            seed, ("NA", "NA", "NA"),
        )

        if (
            chimera_status == "N"
            and high_quality <= MAX_EE
            and (mass >= 3 or spread >= 2)
        ):
            yield (
                i, mass, cloud,
                seed, len(sequence), sequence_abundance,
                chimera_status, spread, high_quality, sequence,
                identity, taxonomy, references,
                "\t".join(str(occurrences[sample]) for sample in samples),
            )
            i += 1


def write_table(
    samples: list[str],
    rows: Iterator[tuple[object, ...]],
) -> None:
    """Write header + every yielded row to stdout, tab-separated."""
    print("PROGRESS: filtering and writing OTUs", file=sys.stderr)
    print(
        "OTU", "total", "cloud",
        "amplicon", "length", "abundance",
        "chimera", "spread", "quality",
        "sequence", "identity", "taxonomy", "references",
        "\t".join(samples),
        sep="\t",
    )
    for row in rows:
        print(*row, sep="\t")


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)

    stampa = stampa_parse(args.assignments)
    representatives = representatives_parse(args.representatives, stampa)
    _, sorted_stats, seeds = stats_parse(args.stats, representatives)
    _, valid_otus = swarms_parse(args.swarms, representatives)
    uchime = uchime_parse(args.chimera, representatives)
    quality = quality_parse(args.quality, representatives)
    seeds2samples, samples = distribution_parse(
        args.distribution, valid_otus,
    )

    write_table(
        samples,
        iter_rows(
            representatives, sorted_stats, seeds, seeds2samples,
            samples, quality, uchime, stampa,
        ),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
