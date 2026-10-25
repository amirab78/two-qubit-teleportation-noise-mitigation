"""Execute the notebooks in place (outputs saved), in order, and report the time of each.

usage:  python scripts/run_notebooks.py [notebook names or numbers ...]     (default: all of notebooks/)
"""
from __future__ import annotations

import os
import sys
import time

import nbformat
from nbclient import NotebookClient

HERE = os.path.dirname(os.path.abspath(__file__))
NB_DIR = os.path.join(os.path.dirname(HERE), "notebooks")


def run(path: str, timeout: int = 1800, kernel: str = "python3") -> float:
    nb = nbformat.read(path, as_version=4)
    t0 = time.time()
    NotebookClient(nb, timeout=timeout, kernel_name=kernel, resources={"metadata": {"path": NB_DIR}}).execute()
    nbformat.write(nb, path)
    return time.time() - t0


def main(args: list[str]) -> int:
    names = sorted(f for f in os.listdir(NB_DIR) if f.endswith(".ipynb"))
    if args:
        names = [n for n in names if any(n.startswith(a) or n == a for a in args)]
    kernel = os.environ.get("TELECORR_KERNEL", "python3")
    for n in names:
        dt = run(os.path.join(NB_DIR, n), kernel=kernel)
        print(f"{n}: {dt:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
