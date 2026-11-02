"""Tests for bin/prepare_silva.py, the SILVA reference builder (D21)."""

# COVERAGE: [S119], [S120], [S121], [S122], [S123], [S124], [S125]

import gzip
import io
import random
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

import pytest

import prepare_silva as ps
from check_reference import check_reference, status

BIN = Path(__file__).resolve().parents[2] / "bin" / "prepare_silva.py"

TAXMAP_HEADER = "primaryAccession\tstart\tstop\tpath\torganism_name\ttaxid\n"

# SILVA 144 layout: kingdom under Bacteria, major clades, super-ranks
TAX_144 = """\
Bacteria;\t3\tdomain\t\t
Bacteria;Pseudomonadati;\t1\tkingdom\t\t144
Bacteria;Pseudomonadati;Pseudomonadota;\t4\tphylum\t\t144
Eukaryota;\t2\tdomain\t\t
Eukaryota;Amorphea;\t5\tmajor_clade\t\t144
Eukaryota;Amorphea;Nucletmycea;\t6\tsuperkingdom\t\t144
Eukaryota;Amorphea;Nucletmycea;Fungi;\t7\tkingdom\t\t144
Eukaryota;Amorphea;Nucletmycea;Fungi;Mag;\t8\tsuperorder\t\t144
"""

TAXMAP_144 = TAXMAP_HEADER + """\
A1\t1\t100\tBacteria;Pseudomonadati;Pseudomonadota;\tRhizoctonia solani\t4
A1\t200\t300\tBacteria;Pseudomonadati;\tother copy\t1
E1\t5\t50\tEukaryota;Amorphea;Nucletmycea;Fungi;Mag;\tName, with comma*\t8
E2\t1\t9\tEukaryota;Amorphea;\t\t5
"""

FASTA_144 = """\
>A1.1.100 Bacteria;Pseudomonadati;Pseudomonadota;Rhizoctonia solani
ACGU
ACGU
>A1.200.300 Bacteria;Pseudomonadati;other copy
uuuu
>E1.5.50 Eukaryota;Amorphea;Nucletmycea;Fungi;Mag;Name, with comma*
GGGG
>E2.1.9 Eukaryota;Amorphea;
CCCC
"""

# SILVA 138.1 / 138.2 layout: no kingdom, phylum right under the domain
TAX_138 = """\
Bacteria;\t3\tdomain\t\t
Bacteria;Proteobacteria;\t2375\tphylum\t\t123
Bacteria;Proteobacteria;Gammaproteobacteria;\t2505\tclass\t\t123
Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;\t\
2580\torder\t\t138
"""

TAXMAP_138 = TAXMAP_HEADER + """\
A16379\t1\t1485\t\
Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;\t\
[Haemophilus] ducreyi\t2580
"""

FASTA_138 = """\
>A16379.1.1485 Bacteria;Proteobacteria;Gammaproteobacteria;\
Enterobacterales;[Haemophilus] ducreyi
ACGU
"""

# 138.1 full SSURef fasta: records the NR99 taxmap does not list
FASTA_138_SSUREF = """\
>BD359735.3.2145 Eukaryota;SAR;Alveolata;Plasmodium malariae
UUUU
""" + FASTA_138 + """\
>X00001.1.10 Bacteria;Proteobacteria;not in NR99
GGGG
"""


def run(fasta: str = FASTA_144, taxmap: str = TAXMAP_144,
        tax: str = TAX_144, **options) -> tuple[list[str], list[str]]:
    """Convert in memory; return the stampa and sintax lines."""
    ranks = ps.parse_taxonomy(io.StringIO(tax))
    lineages = ps.parse_taxmap(io.StringIO(taxmap), ranks)
    stampa, sintax = io.StringIO(), io.StringIO()
    ps.convert(io.StringIO(fasta), lineages, stampa, sintax, **options)
    return stampa.getvalue().splitlines(), sintax.getvalue().splitlines()


def headers(lines: list[str]) -> list[str]:
    return [line for line in lines if line.startswith(">")]


def write_inputs(tmp_path: Path, fasta: str = FASTA_144,
                 taxmap: str = TAXMAP_144, tax: str = TAX_144) -> list[str]:
    """Write the inputs (taxonomy and fasta gzipped) and return CLI args."""
    files = {"tax.txt.gz": tax, "taxmap.txt": taxmap, "seq.fasta.gz": fasta}
    for name, content in files.items():
        path = tmp_path / name
        if name.endswith(".gz"):
            with gzip.open(path, "wt") as handle:
                handle.write(content)
        else:
            path.write_text(content)
    return ["--taxonomy", str(tmp_path / "tax.txt.gz"),
            "--taxmap", str(tmp_path / "taxmap.txt"),
            "--fasta", str(tmp_path / "seq.fasta.gz")]


def run_cli(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(BIN), *args],
                          capture_output=True, text=True)


# ---------------------------------------------------------------------------
# [S119] inputs, outputs, interface


def test_cli_writes_both_references(tmp_path: Path) -> None:
    stampa = tmp_path / "out_stampa.fasta"
    sintax = tmp_path / "out_sintax.fasta.gz"
    result = run_cli(write_inputs(tmp_path) + [
        "--stampa", str(stampa), "--sintax", str(sintax)])
    assert result.returncode == 0, result.stderr
    assert stampa.read_text().startswith(">A1.1.100 domain__Bacteria|")
    with gzip.open(sintax, "rt") as handle:
        assert handle.readline().startswith(">A1.1.100;tax=d:Bacteria,")
    leftovers = [p for p in tmp_path.iterdir() if p.name.startswith(".")]
    assert leftovers == []


def test_cli_requires_an_output(tmp_path: Path) -> None:
    result = run_cli(write_inputs(tmp_path))
    assert result.returncode == 2
    assert "--stampa or --sintax" in result.stderr


def test_outputs_pass_check_reference(tmp_path: Path) -> None:
    stampa = tmp_path / "stampa.fasta.gz"
    sintax = tmp_path / "sintax.fasta"
    assert ps.main(write_inputs(tmp_path) + [
        "--stampa", str(stampa), "--sintax", str(sintax)]) == 0
    assert status(check_reference(str(stampa), "stampa", True)) == "clean"
    assert status(check_reference(str(sintax), "sintax")) == "clean"


# ---------------------------------------------------------------------------
# [S120] lineages and header formats


def test_stampa_headers() -> None:
    stampa, _ = run()
    assert headers(stampa) == [
        ">A1.1.100 domain__Bacteria|kingdom__Pseudomonadati"
        "|phylum__Pseudomonadota|Rhizoctonia_solani",
        ">A1.200.300 domain__Bacteria|kingdom__Pseudomonadati|other_copy",
        ">E1.5.50 domain__Eukaryota|major_clade__Amorphea"
        "|superkingdom__Nucletmycea|kingdom__Fungi|superorder__Mag"
        "|Name,_with_comma_",
        ">E2.1.9 domain__Eukaryota|major_clade__Amorphea",
    ]


def test_sintax_headers() -> None:
    _, sintax = run(organism="strain")
    assert headers(sintax) == [
        ">A1.1.100;tax=d:Bacteria,k:Pseudomonadati,p:Pseudomonadota"
        ",t:Rhizoctonia_solani;",
        ">A1.200.300;tax=d:Bacteria,k:Pseudomonadati,t:other_copy;",
        # superkingdom is not kingdom, superorder is not order
        ">E1.5.50;tax=d:Eukaryota,k:Fungi,t:Name__with_comma*;",
        # empty organism: no trailing separator, no empty rank
        ">E2.1.9;tax=d:Eukaryota;",
    ]


def test_sintax_organism_and_other_ranks_options() -> None:
    _, sintax = run(organism="none", keep_other_ranks=True)
    assert headers(sintax)[2] == (
        ">E1.5.50;tax=d:Eukaryota,x:Amorphea,x:Nucletmycea,k:Fungi,x:Mag;"
    )
    # default "species": this lineage stops at the phylum ([S124])
    _, sintax = run()
    assert headers(sintax)[0].endswith(",p:Pseudomonadota;")


def test_cli_sintax_options(tmp_path: Path) -> None:
    sintax = tmp_path / "sintax.fasta"
    assert ps.main(write_inputs(tmp_path) + [
        "--sintax", str(sintax), "--sintax_organism", "none",
        "--sintax_other_ranks", "x"]) == 0
    assert sintax.read_text().splitlines()[0] == (
        ">A1.1.100;tax=d:Bacteria,k:Pseudomonadati,p:Pseudomonadota;"
    )


def test_each_rrna_copy_keeps_its_own_lineage() -> None:
    stampa, _ = run()
    assert headers(stampa)[0].endswith(
        "|phylum__Pseudomonadota|Rhizoctonia_solani")
    assert headers(stampa)[1].endswith(
        "|kingdom__Pseudomonadati|other_copy")


def test_element_absent_from_taxonomy_has_no_rank() -> None:
    lineage = ps.build_lineage("Archaea;Foo;", "", {"Archaea;": "domain"})
    assert ps.stampa_lineage(lineage) == "domain__Archaea|Foo"
    assert ps.sintax_lineage(lineage, "strain", False) == "d:Archaea"
    assert ps.sintax_lineage(lineage, "strain", True) == "d:Archaea,x:Foo"


def test_rna_to_dna_keeps_line_layout() -> None:
    stampa, sintax = run()
    assert stampa[1:3] == ["ACGT", "ACGT"]
    assert stampa[4] == "tttt"
    assert sintax[1:3] == ["ACGT", "ACGT"]


@pytest.mark.skipif(shutil.which("vsearch") is None, reason="needs vsearch")
def test_vsearch_reads_every_rank(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fasta"
    ps.main(write_inputs(tmp_path) + ["--sintax", str(reference),
                                      "--sintax_organism", "strain"])
    # replace the first toy sequence by a long, non-repetitive one
    sequence = "".join(random.Random(1).choices("ACGT", k=300))
    reference.write_text(
        reference.read_text().replace("ACGT\nACGT", sequence))
    query = tmp_path / "query.fasta"
    query.write_text(f">q\n{sequence}\n")
    result = subprocess.run(
        ["vsearch", "--sintax", str(query), "--db", str(reference),
         "--tabbedout", "-", "--sintax_cutoff", "0.8", "--randseed", "1",
         "--minseqlength", "1", "--quiet"],
        capture_output=True, text=True, check=True,
    )
    assert result.stdout.split("\t")[3].strip() == (
        "d:Bacteria,k:Pseudomonadati,p:Pseudomonadota,t:Rhizoctonia_solani"
    )


# ---------------------------------------------------------------------------
# [S121] failures: exit 1, one-line message, no output


@pytest.mark.parametrize("fasta, taxmap, tax, message", [
    (FASTA_144 + ">Z9.1.10 x\nACGT\n", TAXMAP_144, TAX_144,
     "absent from the taxmap"),
    (FASTA_144 + ">E2.1.9 x\nACGT\n", TAXMAP_144, TAX_144,
     "duplicated identifier"),
    (">A1 x\nACGT\n", TAXMAP_144, TAX_144, "accession.start.stop"),
    (">A1.1.100 x\n>E2.1.9 x\nACGT\n", TAXMAP_144, TAX_144,
     "has no sequence"),
    ("", TAXMAP_144, TAX_144, "no sequence found"),
    (FASTA_144, TAXMAP_144 + "E2\t1\t9\tEukaryota;\t\t2\n", TAX_144,
     "duplicated entry"),
    (FASTA_144, "primaryAccession\tstart\tstop\tpath\n", TAX_144,
     "organism_name"),
    (FASTA_144, TAXMAP_144, "Bacteria;\t3\n", "malformed row"),
], ids=["unmapped", "duplicate_id", "bad_id", "no_sequence", "empty_fasta",
        "duplicate_taxmap_row", "missing_column", "bad_taxonomy_row"])
def test_failures(tmp_path: Path, fasta: str, taxmap: str, tax: str,
                  message: str) -> None:
    result = run_cli(write_inputs(tmp_path, fasta, taxmap, tax) + [
        "--stampa", str(tmp_path / "out_stampa.fasta.gz"),
        "--sintax", str(tmp_path / "out_sintax.fasta")])
    assert result.returncode == 1
    assert message in result.stderr
    assert "Traceback" not in result.stderr
    assert len(result.stderr.strip().splitlines()) <= 3  # info lines + error
    assert not [p for p in tmp_path.iterdir() if "out" in p.name]


def test_skip_unmapped_drops_only_unmapped_records() -> None:
    stampa, _ = run(fasta=FASTA_144 + ">Z9.1.10 x\nACGT\n",
                    skip_unmapped=True)
    assert [h.split()[0] for h in headers(stampa)] == [
        ">A1.1.100", ">A1.200.300", ">E1.5.50", ">E2.1.9"]


# ---------------------------------------------------------------------------
# [S122] SILVA release layouts


def test_layout_138() -> None:
    stampa, sintax = run(FASTA_138, TAXMAP_138, TAX_138, organism="strain")
    assert headers(stampa) == [
        ">A16379.1.1485 domain__Bacteria|phylum__Proteobacteria"
        "|class__Gammaproteobacteria|order__Enterobacterales"
        "|[Haemophilus]_ducreyi",
    ]
    assert headers(sintax) == [
        ">A16379.1.1485;tax=d:Bacteria,p:Proteobacteria"
        ",c:Gammaproteobacteria,o:Enterobacterales,t:[Haemophilus]_ducreyi;",
    ]


def test_layout_138_1_full_ssuref_needs_skip_unmapped(tmp_path: Path) -> None:
    args = write_inputs(tmp_path, FASTA_138_SSUREF, TAXMAP_138, TAX_138)
    output = tmp_path / "out.fasta"
    result = run_cli(args + ["--stampa", str(output)])
    assert result.returncode == 1
    assert "BD359735.3.2145 absent from the taxmap" in result.stderr
    result = run_cli(args + ["--stampa", str(output), "--skip_unmapped"])
    assert result.returncode == 0, result.stderr
    assert "2 fasta record(s) absent from the taxmap skipped" in result.stderr
    assert headers(output.read_text().splitlines()) == [
        ">A16379.1.1485 domain__Bacteria|phylum__Proteobacteria"
        "|class__Gammaproteobacteria|order__Enterobacterales"
        "|[Haemophilus]_ducreyi",
    ]


def test_layout_144() -> None:
    _, sintax = run()
    assert headers(sintax)[0].startswith(
        ">A1.1.100;tax=d:Bacteria,k:Pseudomonadati,p:Pseudomonadota")


# ---------------------------------------------------------------------------
# [S123] sintax placeholders take their nearest named ancestor as prefix

TAX_PLACEHOLDERS = """\
Bacteria;\t3\tdomain\t\t
Bacteria;Fidelibacterota;\t10\tphylum\t\t
Bacteria;Fidelibacterota;Incertae Sedis;\t11\tclass\t\t
Bacteria;Fidelibacterota;Incertae Sedis;Incertae Sedis;\t12\torder\t\t
Bacteria;Fidelibacterota;Incertae Sedis;Incertae Sedis;uncultured;\t\
13\tfamily\t\t
Bacteria;Bacilli;\t20\tclass\t\t
Bacteria;Bacilli;Bacillales--other;\t21\torder\t\t
Bacteria;Bacilli;Bacillales--other;Bacillus--other;\t22\tgenus\t\t
Bacteria;Bacilli;Bacillales--other;Bacillus--other;Incertae Sedis;\t\
23\tspecies_group\t\t
Bacteria;Rhizobiales Incertae Sedis;\t30\tphylum\t\t
Eukaryota;\t2\tdomain\t\t
Eukaryota;Amorphea;\t5\tmajor_clade\t\t
Eukaryota;Amorphea;Unknown Family;\t6\tfamily\t\t
Eukaryota;Amorphea;Unknown Family;Incertae-Sedis;\t7\tgenus\t\t
"""

TAXMAP_PLACEHOLDERS = TAXMAP_HEADER + """\
P1\t1\t9\t\
Bacteria;Fidelibacterota;Incertae Sedis;Incertae Sedis;uncultured;\t\t13
P2\t1\t9\t\
Bacteria;Bacilli;Bacillales--other;Bacillus--other;Incertae Sedis;\t\t23
P3\t1\t9\tBacteria;Rhizobiales Incertae Sedis;\t\t30
P4\t1\t9\tEukaryota;Amorphea;Unknown Family;Incertae-Sedis;\t\t7
"""

FASTA_PLACEHOLDERS = "".join(
    f">P{i}.1.9 x\nACGT\n" for i in range(1, 5))


def run_placeholders(**options) -> tuple[list[str], list[str]]:
    stampa, sintax = run(FASTA_PLACEHOLDERS, TAXMAP_PLACEHOLDERS,
                         TAX_PLACEHOLDERS, **options)
    return headers(stampa), headers(sintax)


def test_placeholders_are_renamed_in_sintax() -> None:
    _, sintax = run_placeholders(keep_other_ranks=True)
    assert sintax == [
        # a run of placeholders: same anchor all the way down
        ">P1.1.9;tax=d:Bacteria,p:Fidelibacterota"
        ",c:Fidelibacterota_Incertae_Sedis"
        ",o:Fidelibacterota_Incertae_Sedis"
        ",f:Fidelibacterota_uncultured;",
        # --other names are renamed, and anchor what is below them
        ">P2.1.9;tax=d:Bacteria,c:Bacilli,o:Bacilli_Bacillales--other"
        ",g:Bacillales--other_Bacillus--other"
        ",x:Bacillus--other_Incertae_Sedis;",
        # a name that only contains "Incertae Sedis" is not a placeholder
        ">P3.1.9;tax=d:Bacteria,p:Rhizobiales_Incertae_Sedis;",
        # the anchor may have a rank sintax does not write (major_clade)
        ">P4.1.9;tax=d:Eukaryota,x:Amorphea,f:Amorphea_Unknown_Family"
        ",g:Amorphea_Incertae-Sedis;",
    ]


def test_placeholders_keep_and_stampa_untouched() -> None:
    stampa, sintax = run_placeholders(rename_placeholders=False)
    assert sintax[0] == (
        ">P1.1.9;tax=d:Bacteria,p:Fidelibacterota,c:Incertae_Sedis"
        ",o:Incertae_Sedis,f:uncultured;")
    renamed_stampa, _ = run_placeholders()
    assert renamed_stampa == stampa
    assert stampa[0] == (
        ">P1.1.9 domain__Bacteria|phylum__Fidelibacterota"
        "|class__Incertae_Sedis|order__Incertae_Sedis|family__uncultured")


def test_placeholder_without_named_ancestor_keeps_its_name() -> None:
    lineage = ps.build_lineage("Incertae Sedis;", "", {
        "Incertae Sedis;": "domain"})
    assert ps.sintax_lineage(lineage, "none", False) == "d:Incertae_Sedis"


def test_cli_sintax_placeholders_keep(tmp_path: Path) -> None:
    sintax = tmp_path / "sintax.fasta"
    args = write_inputs(tmp_path, FASTA_PLACEHOLDERS, TAXMAP_PLACEHOLDERS,
                        TAX_PLACEHOLDERS)
    assert ps.main(args + ["--sintax", str(sintax),
                           "--sintax_placeholders", "keep"]) == 0
    assert sintax.read_text().splitlines()[0].endswith(",f:uncultured;")


# ---------------------------------------------------------------------------
# [S124] organism -> s:Genus_epithet, or nothing

SPECIES_RANKS = {
    "Bacteria;": "domain",
    "Bacteria;Enterobacteriaceae;": "family",
}


@pytest.mark.parametrize("genus, organism, expected", [
    ("Escherichia-Shigella", "Escherichia coli", "Escherichia coli"),
    ("Escherichia-Shigella", "Shigella flexneri 2a str. 2457T",
     "Shigella flexneri"),
    ("Clostridium sensu stricto 1", "Clostridium botulinum A str. ATCC 3502",
     "Clostridium botulinum"),
    ("[Eubacterium] coprostanoligenes group",
     "[Eubacterium] coprostanoligenes", "[Eubacterium] coprostanoligenes"),
    ("Haemophilus", "[Haemophilus] ducreyi", "[Haemophilus] ducreyi"),
    ("Candidatus Pelagibacter", "Candidatus Pelagibacter ubique HTCC1062",
     "Candidatus Pelagibacter ubique"),
    ("Fusobacterium", "Fusobacterium nucleatum subsp. vincentii ATCC 49256",
     "Fusobacterium nucleatum"),
    ("Auricularia", "Auricularia auricula-judae (ear fungus)",
     "Auricularia auricula-judae"),
    ("Plasmodium", "Plasmodium vivax (malaria parasite P. vivax)",
     "Plasmodium vivax"),
    # genus disagrees: host organism, renamed genus
    ("Pseudomonas", "Rhizoctonia solani", None),
    ("Lactiplantibacillus", "Lactobacillus plantarum", None),
    ("Pseudomonas", "metagenome", None),
    # uncultured or empty organism
    ("Pseudomonas", "uncultured bacterium", None),
    ("Pseudomonas", "uncultured Pseudomonas fluorescens", None),
    ("Pseudomonas", "", None),
    # no epithet
    ("Pseudomonas", "Pseudomonas", None),
    ("Pseudomonas", "Pseudomonas sp.", None),
    ("Pseudomonas", "Pseudomonas sp", None),
    ("Pseudomonas", "Pseudomonas HTCC2255", None),
    ("Pseudomonas", "Pseudomonas 2a", None),
    ("Phytophthora", "Phytophthora taxon juncus", None),
    ("Rickettsia", "Rickettsia secondary endosymbiont of Archarius", None),
    ("Glomus", "Glomus mycorrhizal symbiont of Marchantia", None),
    ("Nitrospira", "Nitrospira enrichment culture clone GaII", None),
    ("Wolbachia", "Wolbachia endosymbiont of Drosophila", None),
    # the genus itself is a placeholder
    ("uncultured", "uncultured bacterium", None),
    ("Incertae Sedis", "Incertae Sedis bacterium", None),
])
def test_species_name(genus: str, organism: str,
                      expected: Optional[str]) -> None:
    path = f"Bacteria;Enterobacteriaceae;{genus};"
    ranks = dict(SPECIES_RANKS, **{path: "genus"})
    lineage = ps.build_lineage(path, organism, ranks)
    assert ps.species_name(lineage) == expected


def test_species_needs_a_lineage_ending_at_genus() -> None:
    lineage = ps.build_lineage("Bacteria;Enterobacteriaceae;",
                               "Escherichia coli", SPECIES_RANKS)
    assert ps.species_name(lineage) is None


def test_species_in_sintax_header() -> None:
    path = "Bacteria;Enterobacteriaceae;Escherichia-Shigella;"
    ranks = dict(SPECIES_RANKS, **{path: "genus"})
    lineage = ps.build_lineage(
        path, "Escherichia coli O157:H7 str. Sakai", ranks)
    assert ps.sintax_lineage(lineage) == (
        "d:Bacteria,f:Enterobacteriaceae,g:Escherichia-Shigella"
        ",s:Escherichia_coli")
    assert ps.sintax_lineage(lineage, "strain") == (
        "d:Bacteria,f:Enterobacteriaceae,g:Escherichia-Shigella"
        ",t:Escherichia_coli_O157:H7_str._Sakai")
    assert ps.sintax_lineage(lineage, "none") == (
        "d:Bacteria,f:Enterobacteriaceae,g:Escherichia-Shigella")


# ---------------------------------------------------------------------------
# [S125] homonym report

TAX_HOMONYMS = """\
Eukaryota;\t2\tdomain\t\t
Eukaryota;Tricholomataceae;\t3\tfamily\t\t
Eukaryota;Tricholomataceae;Tricholoma;\t4\tgenus\t\t
Eukaryota;Ulvaceae;\t5\tfamily\t\t
Eukaryota;Ulvaceae;Tricholoma;\t6\tgenus\t\t
Eukaryota;Ulvaceae;Incertae Sedis;\t7\tgenus\t\t
Eukaryota;Tricholomataceae;Incertae Sedis;\t8\tgenus\t\t
"""

TAXMAP_HOMONYMS = TAXMAP_HEADER + """\
H1\t1\t9\tEukaryota;Tricholomataceae;Tricholoma;\tTricholoma terreum\t4
H2\t1\t9\tEukaryota;Ulvaceae;Tricholoma;\tTricholoma sp.\t6
H3\t1\t9\tEukaryota;Ulvaceae;Incertae Sedis;\t\t7
H4\t1\t9\tEukaryota;Tricholomataceae;Incertae Sedis;\t\t8
"""

FASTA_HOMONYMS = "".join(f">H{i}.1.9 x\nACGT\n" for i in range(1, 5))


def test_homonym_report(tmp_path: Path) -> None:
    report = tmp_path / "homonyms.tsv"
    args = write_inputs(tmp_path, FASTA_HOMONYMS, TAXMAP_HOMONYMS,
                        TAX_HOMONYMS)
    result = run_cli(args + ["--sintax", str(tmp_path / "s.fasta"),
                             "--homonyms", str(report)])
    assert result.returncode == 0, result.stderr
    assert "1 sintax (rank, name) pair(s) under more than one parent" \
        in result.stderr
    # the renamed placeholders (g:Ulvaceae_Incertae_Sedis, ...) are not
    # homonyms
    assert report.read_text().splitlines() == [
        "rank\tname\tparent_lineage",
        "g\tTricholoma\td:Eukaryota,f:Tricholomataceae",
        "g\tTricholoma\td:Eukaryota,f:Ulvaceae",
    ]


def test_homonym_report_lists_kept_placeholders(tmp_path: Path) -> None:
    report = tmp_path / "homonyms.tsv"
    args = write_inputs(tmp_path, FASTA_HOMONYMS, TAXMAP_HOMONYMS,
                        TAX_HOMONYMS)
    assert ps.main(args + ["--sintax", str(tmp_path / "s.fasta"),
                           "--sintax_placeholders", "keep",
                           "--homonyms", str(report)]) == 0
    assert "g\tIncertae_Sedis\td:Eukaryota,f:Ulvaceae" \
        in report.read_text().splitlines()


def test_homonym_report_header_only(tmp_path: Path) -> None:
    report = tmp_path / "homonyms.tsv"
    assert ps.main(write_inputs(tmp_path) + [
        "--sintax", str(tmp_path / "s.fasta"),
        "--homonyms", str(report)]) == 0
    assert report.read_text() == "rank\tname\tparent_lineage\n"
