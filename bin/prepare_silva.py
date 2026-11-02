#!/usr/bin/env python3
"""Build rank-annotated SILVA SSU references for stampa and sintax.

Reads three files of one SILVA release (plain or gzip-compressed):

- ``tax_slv_ssu_<release>.txt``: one row per taxon, giving its path
  and its rank (``domain``, ``phylum``, ``major_clade``, ...);
- ``taxmap_slv_ssu_ref_nr_<release>.txt``: one row per sequence
  (``primaryAccession``, ``start``, ``stop``, ``path``,
  ``organism_name``, ``taxid``);
- ``SILVA_<release>_SSURef_NR99_tax_silva_trunc.fasta``: the
  sequences, with headers ``>accession.start.stop lineage``. Older
  releases ship it without ``trunc``; for 138.1 the full
  ``SSURef`` fasta can be used with ``--skip_unmapped``, the NR99
  taxmap then selecting the NR99 records (``[S122]``).

And writes up to two references (gzip-compressed if the name ends with
``.gz``):

- stampa (``--stampa``)::

      >HG374866.9.1343 domain__Bacteria|...|Rhizoctonia_solani

- sintax (``--sintax``), for ``vsearch --sintax``::

      >HG374866.9.1343;tax=d:Bacteria,...,t:Rhizoctonia_solani;

Sequences are converted from RNA to DNA (``U`` -> ``T``).

Design notes:

- taxmap rows are keyed on (accession, start, stop), not on the
  accession alone: a genome can contribute several rRNA copies, and
  SILVA may classify them differently;
- a fasta record without a taxmap row (unless ``--skip_unmapped``),
  or a duplicated identifier, stops the script (no silent drop, no
  orphan sequence lines);
- sintax ranks are mapped by exact rank name (``superkingdom`` is not
  ``kingdom``). vsearch only reads the ranks ``d k p c o f g s t``
  and silently ignores any other, so other SILVA ranks
  (``major_clade``, ``subphylum``, ...) are dropped by default;
- in the sintax output only (``[S123]``-``[S125]``): placeholder names
  (``uncultured``, ``Incertae Sedis``, ``X--other``, ...) take their
  nearest named ancestor as a prefix, because vsearch would otherwise
  pool unrelated lineages that share the name; the organism becomes
  ``s:Genus_epithet`` only when it agrees with the lineage's genus
  (SILVA's organism is the source organism of the record: a host
  genome, ``metagenome``...); remaining homonyms are reported;
- characters that break a header grammar are replaced by ``_``:
  whitespace everywhere, ``|`` and ``*`` in stampa lineages (``*`` is
  stampa's wildcard), ``,`` and ``;`` in sintax names;
- outputs are written to a temporary file and renamed on success, so a
  failed run never leaves a truncated reference behind.

The workflow does not call this script (``D21``, ``[S119]``-``[S125]``):
run it once per SILVA release and pass its outputs to
``--reference_dataset`` (stampa) and ``--reference_dataset_sintax``.
"""

import sys

if sys.version_info < (3, 10):
    sys.exit(
        "prepare_silva.py: Python >= 3.10 required, found "
        + sys.version.split()[0]
        + ". Enable a container or conda profile"
        + " (e.g. the pipeline's environment.yml)."
    )

import argparse  # noqa: E402
import csv  # noqa: E402
import gzip  # noqa: E402
import logging  # noqa: E402
import os  # noqa: E402
import re  # noqa: E402
from collections import Counter  # noqa: E402
from contextlib import ExitStack, contextmanager  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import IO, Iterator, NamedTuple, Optional  # noqa: E402

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2026/11/02"
__version__ = "1.0"

LOGGER = logging.getLogger("prepare_silva")

GZIP_MAGIC = b"\x1f\x8b"

# SILVA rank names -> sintax rank letters understood by vsearch
SINTAX_RANKS: dict[str, str] = {
    "domain": "d",
    "kingdom": "k",
    "phylum": "p",
    "class": "c",
    "order": "o",
    "family": "f",
    "genus": "g",
}
# rank letter used for unrecognised ranks with --sintax_other_ranks x
OTHER_RANK = "x"
SPECIES_RANK = "s"
STRAIN_RANK = "t"

# [S123]: SILVA names that stand for "no name at this rank" ...
PLACEHOLDER = re.compile(
    r"^(?:uncultured|Incertae[ -]Sedis\??|Unknown Family)$"
)
# ... and SILVA 144's "other members of X", renamed too, but also an
# anchor for the placeholders below it
OTHER_SUFFIX = "--other"
# [S124]: a species epithet: lowercase letters and hyphens only ...
EPITHET = re.compile(r"^[a-z][a-z-]+$")
# ... but not one of these words, which follow a genus name in SILVA
# organism names without being an epithet
NOT_EPITHETS = frozenset({
    "sp", "spp", "cf", "aff", "nr", "bacterium", "archaeon",
    "endosymbiont", "symbiont", "clone", "isolate", "genomosp", "group",
    "complex", "metagenome", "environmental", "taxon", "primary",
    "secondary", "mycorrhizal", "enrichment",
})
CANDIDATUS = "Candidatus"

TAXMAP_COLUMNS = ("primaryAccession", "start", "stop", "path", "organism_name")

WHITESPACE = re.compile(r"\s+")
STAMPA_FORBIDDEN = re.compile(r"[|*]")
SINTAX_FORBIDDEN = re.compile(r"[,;]")


class PrepareError(Exception):
    """Inconsistent or malformed input; reported without a traceback."""


class Taxon(NamedTuple):
    """One element of a taxonomic path."""

    rank: Optional[str]  # None when the path is absent from the tax file
    name: str


class Lineage(NamedTuple):
    """Ranked path of one sequence, plus its organism name."""

    taxa: tuple[Taxon, ...]
    organism: str


SequenceKey = tuple[str, int, int]


# ------------------------------------------------------------------ io --


def open_text(path: Path) -> IO[str]:
    """Open a plain or gzip-compressed file for reading (magic bytes)."""
    with path.open("rb") as probe:
        compressed = probe.read(2) == GZIP_MAGIC
    if compressed:
        return gzip.open(path, mode="rt", encoding="utf-8", newline="")
    return path.open(mode="rt", encoding="utf-8", newline="")


@contextmanager
def atomic_output(path: Path) -> Iterator[IO[str]]:
    """Write to ``path`` via a temporary sibling, renamed on success."""
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        if path.suffix == ".gz":
            handle: IO[str] = gzip.open(
                tmp, mode="wt", encoding="utf-8", compresslevel=6
            )
        else:
            handle = tmp.open(mode="wt", encoding="utf-8")
        with handle:
            yield handle
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


# ------------------------------------------------------------- parsing --


def split_path(path: str) -> list[str]:
    """``'A;B;C;'`` -> ``['A', 'B', 'C']``."""
    return path.rstrip(";").split(";")


def parse_taxonomy(handle: IO[str]) -> dict[str, str]:
    """Map each taxon path (``'A;B;'``) to its rank (``'phylum'``)."""
    ranks: dict[str, str] = {}
    for number, line in enumerate(handle, start=1):
        if not line.strip():
            continue
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 3 or not fields[0].endswith(";") or not fields[2]:
            raise PrepareError(f"taxonomy line {number}: malformed row")
        ranks[fields[0]] = fields[2]
    if not ranks:
        raise PrepareError("taxonomy file is empty")
    return ranks


def build_lineage(path: str, organism: str, ranks: dict[str, str]) -> Lineage:
    """Attach a rank to every element of ``path``."""
    elements = split_path(path)
    taxa = tuple(
        Taxon(ranks.get(";".join(elements[: i + 1]) + ";"), name)
        for i, name in enumerate(elements)
    )
    return Lineage(taxa, organism)


def parse_taxmap(
    handle: IO[str], ranks: dict[str, str]
) -> dict[SequenceKey, Lineage]:
    """Map (accession, start, stop) to a ranked lineage."""
    reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
    missing = set(TAXMAP_COLUMNS) - set(reader.fieldnames or ())
    if missing:
        raise PrepareError(
            f"taxmap: missing column(s) {', '.join(sorted(missing))}"
        )
    lineages: dict[SequenceKey, Lineage] = {}
    for row in reader:
        try:
            key = (
                row["primaryAccession"], int(row["start"]), int(row["stop"])
            )
        except (TypeError, ValueError) as error:
            raise PrepareError(
                f"taxmap line {reader.line_num}: bad coordinates"
            ) from error
        if key in lineages:
            raise PrepareError(
                f"taxmap line {reader.line_num}: duplicated entry {key}"
            )
        lineages[key] = build_lineage(
            row["path"], row["organism_name"] or "", ranks
        )
    return lineages


def parse_identifier(identifier: str) -> SequenceKey:
    """``'HG374866.9.1343'`` -> ``('HG374866', 9, 1343)``."""
    parts = identifier.rsplit(".", 2)
    try:
        accession, start, stop = parts
        return accession, int(start), int(stop)
    except ValueError as error:
        raise PrepareError(
            f"fasta: identifier {identifier!r} is not accession.start.stop"
        ) from error


def iter_fasta(handle: IO[str]) -> Iterator[tuple[str, list[str]]]:
    """Yield ``(identifier, sequence lines)`` for each fasta record."""
    identifier: Optional[str] = None
    lines: list[str] = []
    for number, raw in enumerate(handle, start=1):
        line = raw.rstrip("\r\n")
        if line.startswith(">"):
            if identifier is not None:
                yield identifier, lines
            identifier = line[1:].split(maxsplit=1)[0] if line[1:] else ""
            if not identifier:
                raise PrepareError(f"fasta line {number}: empty header")
            lines = []
        elif identifier is None:
            if line.strip():
                raise PrepareError(f"fasta line {number}: data before header")
        elif line:
            lines.append(line)
    if identifier is not None:
        yield identifier, lines


# ---------------------------------------------------------- formatting --


def clean(name: str, forbidden: re.Pattern[str]) -> str:
    """Replace whitespace and grammar-breaking characters by ``_``."""
    return forbidden.sub("_", WHITESPACE.sub("_", name.strip()))


def stampa_lineage(lineage: Lineage) -> str:
    """``domain__Bacteria|...|genus__Pseudomonas|Organism_name``."""
    fields = [
        clean(f"{taxon.rank}__{taxon.name}" if taxon.rank else taxon.name,
              STAMPA_FORBIDDEN)
        for taxon in lineage.taxa
    ]
    organism = clean(lineage.organism, STAMPA_FORBIDDEN)
    if organism:
        fields.append(organism)
    return "|".join(fields)


def sintax_fields(
    lineage: Lineage,
    organism: str = "species",
    keep_other_ranks: bool = False,
    rename_placeholders: bool = True,
) -> list[str]:
    """``['d:Bacteria', ..., 'g:Escherichia-Shigella', 's:Escherichia_coli']``.

    ``organism``: ``species`` ([S124]), ``strain`` (the whole organism
    name as ``t:``) or ``none``.
    """
    if rename_placeholders:
        names = anchored_names(lineage)
    else:
        names = [taxon.name for taxon in lineage.taxa]
    fields = []
    for taxon, name in zip(lineage.taxa, names):
        letter = SINTAX_RANKS.get(taxon.rank or "")
        if letter is None:
            if not keep_other_ranks:
                continue
            letter = OTHER_RANK
        fields.append(f"{letter}:{clean(name, SINTAX_FORBIDDEN)}")
    if organism == "species":
        letter, name = SPECIES_RANK, species_name(lineage)
    elif organism == "strain":
        letter, name = STRAIN_RANK, lineage.organism
    else:
        letter, name = "", None
    name = clean(name or "", SINTAX_FORBIDDEN)
    if name:
        fields.append(f"{letter}:{name}")
    return fields


def sintax_lineage(
    lineage: Lineage,
    organism: str = "species",
    keep_other_ranks: bool = False,
    rename_placeholders: bool = True,
) -> str:
    """``d:Bacteria,...,g:Escherichia-Shigella,s:Escherichia_coli``."""
    return ",".join(sintax_fields(
        lineage, organism, keep_other_ranks, rename_placeholders))


class HomonymIndex:
    """Parent lineages of each sintax (rank, name) pair ([S125])."""

    def __init__(self) -> None:
        self.parents: dict[str, set[str]] = {}

    def add(self, fields: list[str]) -> None:
        parent = ""
        for field in fields:
            self.parents.setdefault(field, set()).add(parent)
            parent = f"{parent},{field}" if parent else field

    def homonyms(self) -> list[tuple[str, str, str]]:
        """``(rank, name, parent lineage)`` of every pair under more
        than one parent, sorted."""
        rows = []
        for field, parents in self.parents.items():
            if len(parents) > 1:
                rank, _, name = field.partition(":")
                rows.extend((rank, name, parent) for parent in parents)
        return sorted(rows)

    def write(self, path: Path) -> None:
        with atomic_output(path) as handle:
            print("rank", "name", "parent_lineage", sep="\t", file=handle)
            for row in self.homonyms():
                print(*row, sep="\t", file=handle)


def is_placeholder(name: str) -> bool:
    """``uncultured``, ``Incertae Sedis``, ``Bacillus--other``, ... ."""
    return PLACEHOLDER.match(name) is not None or name.endswith(OTHER_SUFFIX)


def anchored_names(lineage: Lineage) -> list[str]:
    """Taxon names, placeholders prefixed with their nearest named
    ancestor ([S123]): ``Fidelibacterota;Incertae Sedis`` ->
    ``Fidelibacterota``, ``Fidelibacterota Incertae Sedis``.

    An ``X--other`` name is prefixed too, and is then the anchor of
    the placeholders below it, under its own name.
    """
    names = []
    anchor = ""
    for taxon in lineage.taxa:
        name = taxon.name
        if is_placeholder(name):
            names.append(f"{anchor} {name}" if anchor else name)
        else:
            names.append(name)
        if not PLACEHOLDER.match(name):
            anchor = name
    return names


def genus_names(genus: str) -> set[str]:
    """Genus names a SILVA genus stands for ([S124]).

    ``Escherichia-Shigella`` -> {Escherichia, Shigella};
    ``Clostridium sensu stricto 1`` -> {Clostridium};
    ``[Eubacterium] coprostanoligenes group`` -> {Eubacterium};
    ``Candidatus Pelagibacter`` -> {Pelagibacter}.
    """
    words = genus.replace("[", "").replace("]", "").split()
    if words and words[0] == CANDIDATUS:
        words = words[1:]
    return set(filter(None, words[0].split("-"))) if words else set()


def species_name(lineage: Lineage) -> Optional[str]:
    """``Genus epithet`` when the organism agrees with the lineage's
    genus, None otherwise ([S124])."""
    if not lineage.taxa:
        return None
    genus = lineage.taxa[-1]
    if genus.rank != "genus" or is_placeholder(genus.name):
        return None
    words = lineage.organism.split()
    prefix = []
    if words and words[0] == CANDIDATUS:
        prefix = [CANDIDATUS]
        words = words[1:]
    if len(words) < 2 or words[0] == "uncultured":
        return None
    if words[0].replace("[", "").replace("]", "") \
            not in genus_names(genus.name):
        return None
    epithet = words[1]
    if not EPITHET.match(epithet) or epithet in NOT_EPITHETS:
        return None
    return " ".join(prefix + words[:2])


def to_dna(lines: list[str]) -> list[str]:
    """RNA -> DNA alphabet, case preserved."""
    return [line.replace("U", "T").replace("u", "t") for line in lines]


# ----------------------------------------------------------------- run --


def convert(
    fasta: IO[str],
    lineages: dict[SequenceKey, Lineage],
    stampa: Optional[IO[str]],
    sintax: Optional[IO[str]],
    organism: str = "species",
    keep_other_ranks: bool = False,
    skip_unmapped: bool = False,
    rename_placeholders: bool = True,
    index: Optional[HomonymIndex] = None,
) -> Counter[str]:
    """Write one stampa and/or sintax record per fasta record."""
    stats: Counter[str] = Counter()
    seen: set[str] = set()
    for identifier, lines in iter_fasta(fasta):
        if identifier in seen:
            raise PrepareError(f"fasta: duplicated identifier {identifier}")
        seen.add(identifier)
        key = parse_identifier(identifier)
        lineage = lineages.get(key)
        if lineage is None:
            if skip_unmapped:
                stats["skipped"] += 1
                continue
            raise PrepareError(f"fasta: {identifier} absent from the taxmap")
        if not lines:
            raise PrepareError(f"fasta: {identifier} has no sequence")
        if any(taxon.rank is None for taxon in lineage.taxa):
            stats["unranked_taxon"] += 1
        sequence = "\n".join(to_dna(lines))
        if stampa is not None:
            stampa.write(f">{identifier} {stampa_lineage(lineage)}\n")
            stampa.write(f"{sequence}\n")
        if sintax is not None:
            fields = sintax_fields(lineage, organism, keep_other_ranks,
                                   rename_placeholders)
            if index is not None:
                index.add(fields)
            if fields and fields[-1].startswith(f"{SPECIES_RANK}:"):
                stats["species"] += 1
            tax = ",".join(fields)
            sintax.write(f">{identifier};tax={tax};\n{sequence}\n")
        stats["sequences"] += 1
    if not stats["sequences"]:
        raise PrepareError("fasta: no sequence found")
    unused = len(lineages) - stats["sequences"]
    if unused:
        LOGGER.warning("%d taxmap row(s) without a fasta sequence", unused)
    return stats


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Build rank-annotated stampa and/or sintax references from "
            "a SILVA SSU release (tax_slv, taxmap and fasta files)."
        ),
    )
    parser.add_argument(
        "-l", "--taxonomic_levels", "--taxonomy", dest="taxonomy",
        type=Path, required=True, metavar="FILE",
        help="tax_slv_ssu_<release>.txt[.gz]",
    )
    parser.add_argument(
        "-m", "--taxonomic_map", "--taxmap", dest="taxmap",
        type=Path, required=True, metavar="FILE",
        help="taxmap_slv_ssu_ref_nr_<release>.txt[.gz]",
    )
    parser.add_argument(
        "-f", "--fasta_sequences", "--fasta", dest="fasta",
        type=Path, required=True, metavar="FILE",
        help="SILVA_<release>_SSURef_NR99_tax_silva_trunc.fasta[.gz]",
    )
    parser.add_argument(
        "--stampa", type=Path, metavar="FILE",
        help="stampa output (gzip-compressed if FILE ends with .gz)",
    )
    parser.add_argument(
        "--sintax", type=Path, metavar="FILE",
        help="sintax output (gzip-compressed if FILE ends with .gz)",
    )
    parser.add_argument(
        "--sintax_organism", choices=("species", "strain", "none"),
        default="species",
        help=(
            "what sintax does with SILVA's organism name: 'species' "
            "(default) writes s:Genus_epithet when the name agrees with "
            "the lineage's genus, nothing otherwise; 'strain' writes the "
            "whole name as t:; 'none' omits it"
        ),
    )
    parser.add_argument(
        "--sintax_placeholders", choices=("rename", "keep"),
        default="rename",
        help=(
            "sintax placeholder names (uncultured, Incertae Sedis, "
            "X--other, ...): prefix them with their nearest named "
            "ancestor (default) or keep them as is"
        ),
    )
    parser.add_argument(
        "--homonyms", type=Path, metavar="FILE",
        help=(
            "write the sintax (rank, name) pairs found under more than "
            "one parent lineage to FILE (TSV)"
        ),
    )
    parser.add_argument(
        "--sintax_other_ranks", choices=("drop", "x"), default="drop",
        help=(
            "ranks vsearch does not know (major_clade, subphylum, ...): "
            "drop them (default) or keep them as 'x:name' "
            "(vsearch ignores them either way)"
        ),
    )
    parser.add_argument(
        "--skip_unmapped", action="store_true",
        help=(
            "drop fasta records absent from the taxmap instead of "
            "stopping (e.g. SILVA 138.1 full SSURef fasta with the NR99 "
            "taxmap)"
        ),
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version=f"%(prog)s {__version__}",
    )
    args = parser.parse_args(argv)
    if args.stampa is None and args.sintax is None:
        parser.error("at least one of --stampa or --sintax is required")
    return args


def main(argv: Optional[list[str]] = None) -> int:
    logging.basicConfig(format="%(name)s: %(message)s", level=logging.INFO)
    args = parse_args(argv)
    index = HomonymIndex() if args.sintax else None
    try:
        with open_text(args.taxonomy) as handle:
            ranks = parse_taxonomy(handle)
        LOGGER.info("%d taxa read from %s", len(ranks), args.taxonomy)
        with open_text(args.taxmap) as handle:
            lineages = parse_taxmap(handle, ranks)
        LOGGER.info("%d taxmap rows read from %s", len(lineages), args.taxmap)
        with ExitStack() as stack:
            fasta = stack.enter_context(open_text(args.fasta))
            stampa = stack.enter_context(atomic_output(args.stampa)) \
                if args.stampa else None
            sintax = stack.enter_context(atomic_output(args.sintax)) \
                if args.sintax else None
            stats = convert(
                fasta, lineages, stampa, sintax,
                organism=args.sintax_organism,
                keep_other_ranks=args.sintax_other_ranks == "x",
                skip_unmapped=args.skip_unmapped,
                rename_placeholders=args.sintax_placeholders == "rename",
                index=index,
            )
        if index is not None and args.homonyms:
            index.write(args.homonyms)
    except (PrepareError, OSError, UnicodeDecodeError, csv.Error) as error:
        LOGGER.error("%s", error)
        return 1
    LOGGER.info("%d sequences written", stats["sequences"])
    if stats["skipped"]:
        LOGGER.info(
            "%d fasta record(s) absent from the taxmap skipped",
            stats["skipped"],
        )
    if args.sintax and args.sintax_organism == "species":
        LOGGER.info("%d sequence(s) with a sintax species (s:)",
                    stats["species"])
    if index is not None:
        pairs = len({row[:2] for row in index.homonyms()})
        LOGGER.info(
            "%d sintax (rank, name) pair(s) under more than one parent "
            "lineage%s", pairs,
            f" (see {args.homonyms})" if args.homonyms else "",
        )
    if stats["unranked_taxon"]:
        LOGGER.warning(
            "%d sequence(s) with a taxon absent from the taxonomy file",
            stats["unranked_taxon"],
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
