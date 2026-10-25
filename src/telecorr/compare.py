"""Compare a recomputed result with a saved one, value by value (used by the regression tests and notebooks)."""
from __future__ import annotations

import json
import math
from typing import Any


def _plain(x: Any) -> Any:
    """JSON round trip: tuples -> lists, int dict keys -> str, numpy scalars -> Python numbers."""
    return json.loads(json.dumps(x, default=float))


def tree_diff(mine: Any, saved: Any, rtol: float = 1e-11, atol: float = 1e-12, path: str = "") -> tuple[int, list]:
    """Walk ``saved`` and compare every leaf with ``mine``. Numbers agree when
    |a - b| <= atol + rtol |b|; everything else must be equal. Keys only in ``mine`` are ignored.
    Returns (number of compared leaves, list of mismatches (path, mine, saved))."""
    mine, saved = _plain(mine), _plain(saved)
    n, bad = 0, []

    def walk(a, b, p):
        nonlocal n
        if isinstance(b, dict):
            if not isinstance(a, dict):
                bad.append((p, type(a).__name__, "dict")); return
            for k, v in b.items():
                if k not in a:
                    bad.append((f"{p}/{k}", "missing", v)); continue
                walk(a[k], v, f"{p}/{k}")
        elif isinstance(b, list):
            if not isinstance(a, list) or len(a) != len(b):
                bad.append((p, f"len {len(a) if isinstance(a, list) else '-'}", f"len {len(b)}")); return
            for i, (x, y) in enumerate(zip(a, b)):
                walk(x, y, f"{p}[{i}]")
        else:
            n += 1
            num = lambda v: isinstance(v, (int, float)) and not isinstance(v, bool)   # noqa: E731
            if num(a) and num(b):
                if not (math.isfinite(a) and math.isfinite(b)) or abs(a - b) > atol + rtol * abs(b):
                    if not (a == b):
                        bad.append((p, a, b))
            elif a != b:
                bad.append((p, a, b))

    walk(mine, saved, path)
    return n, bad


def assert_tree_close(mine: Any, saved: Any, rtol: float = 1e-11, atol: float = 1e-12, min_leaves: int = 1) -> int:
    """Raise AssertionError listing the first mismatches; return the number of compared values."""
    n, bad = tree_diff(mine, saved, rtol, atol)
    assert not bad, f"{len(bad)} of {n} values differ, e.g. " + "; ".join(map(str, bad[:5]))
    assert n >= min_leaves, f"only {n} values compared"
    return n
