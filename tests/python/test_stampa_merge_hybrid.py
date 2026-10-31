"""bin/stampa_merge.py on hits against a hybrid reference header."""

# COVERAGE: [S64]

import io
from pathlib import Path

from stampa_merge import emit_results


def test_hybrid_header_keeps_stampa_lineage_and_full_token(
    tmp_path: Path,
) -> None:
    # [S64]: a hybrid header carries both annotations
    # (`>id;tax=...; lineage`). stampa splits the vsearch target at the
    # first space, so the lineage after the space is assigned and the
    # whole `id;tax=...;` token lands in the references column.
    hits = tmp_path / "hits"
    hits.write_text(
        "amp1;size=10;\t100.0\t"
        "refA;tax=d:Bacteria,p:Proteobacteria; Bacteria|Proteobacteria\n"
    )
    out = io.StringIO()
    emit_results(str(hits), out)
    assert out.getvalue().splitlines() == [
        "amp1\t10\t100.0\tBacteria|Proteobacteria\t"
        "refA;tax=d:Bacteria,p:Proteobacteria;"
    ]
