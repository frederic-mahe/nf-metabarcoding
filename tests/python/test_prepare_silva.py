"""Tests for bin/prepare_silva.py, the SILVA reference builder (D21)."""

# COVERAGE: [S119], [S120], [S121], [S122]

import gzip
import io
import random
import shutil
import subprocess
import sys
from pathlib import Path

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
    _, sintax = run()
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
    _, sintax = run(organism_rank=None, keep_other_ranks=True)
    assert headers(sintax)[2] == (
        ">E1.5.50;tax=d:Eukaryota,x:Amorphea,x:Nucletmycea,k:Fungi,x:Mag;"
    )
    _, sintax = run(organism_rank="s")
    assert headers(sintax)[0].endswith(",s:Rhizoctonia_solani;")


def test_cli_sintax_options(tmp_path: Path) -> None:
    sintax = tmp_path / "sintax.fasta"
    assert ps.main(write_inputs(tmp_path) + [
        "--sintax", str(sintax), "--sintax_organism_rank", "none",
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
    assert ps.sintax_lineage(lineage, "t", False) == "d:Archaea"
    assert ps.sintax_lineage(lineage, "t", True) == "d:Archaea,x:Foo"


def test_rna_to_dna_keeps_line_layout() -> None:
    stampa, sintax = run()
    assert stampa[1:3] == ["ACGT", "ACGT"]
    assert stampa[4] == "tttt"
    assert sintax[1:3] == ["ACGT", "ACGT"]


@pytest.mark.skipif(shutil.which("vsearch") is None, reason="needs vsearch")
def test_vsearch_reads_every_rank(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fasta"
    ps.main(write_inputs(tmp_path) + ["--sintax", str(reference)])
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
    stampa, sintax = run(FASTA_138, TAXMAP_138, TAX_138)
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
        ">A1.1.100;tax=d:Bacteria,k:Pseudomonadati,p:Pseudomonadota,")
