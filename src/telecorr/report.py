"""Small helpers for notebooks: Markdown tables and sentences built from computed values."""
from __future__ import annotations

from typing import Any, Iterable, Sequence

import numpy as np


def fmt(v: Any, nd: int = 4) -> str:
    if v is None:
        return "-"
    if isinstance(v, np.ndarray):
        v = v.tolist()
    if isinstance(v, np.generic):
        v = v.item()
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return f"{v:.{nd}f}" if abs(v) >= 10 ** (-nd) or v == 0 else f"{v:.2e}"
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(fmt(x, nd) for x in v) + "]"
    return str(v)


def md_table(rows: Iterable[Sequence[Any]], headers: Sequence[str], nd: int = 4) -> str:
    """A GitHub-flavoured Markdown table (floats with ``nd`` decimals)."""
    out = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    for r in rows:
        out.append("| " + " | ".join(fmt(v, nd) for v in r) + " |")
    return "\n".join(out)


def show(text: str) -> None:
    """Display Markdown in a notebook (prints it elsewhere)."""
    try:
        from IPython.display import Markdown, display
        display(Markdown(text))
    except ImportError:                                   # pragma: no cover
        print(text)
