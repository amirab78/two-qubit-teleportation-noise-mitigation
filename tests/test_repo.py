"""Repository checks: the README's numbers equal their sources, no secrets anywhere, the submission script is disabled."""
import json
import os
import re
import subprocess
import sys

import pytest

from telecorr.reference import number, printed, simulation

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
README = open(os.path.join(REPO, "README.md"), encoding="utf8").read()


def ref(name):
    with open(os.path.join(REPO, "data", "reference", name), encoding="utf8") as fh:
        return json.load(fh)


def test_readme_numbers_equal_their_sources():
    nominal = {r["label"]: r for r in simulation("6a_nominal.json")["rows"]}
    rec = next(c for c in ref("9_8_checks.json")["claims"] if c["key"] == "recov")
    rows = {r["state"]: r for r in next(ch for ch in rec["checks"] if "rows" in ch)["rows"]}
    n2 = ref("8e_n2_analysis.json")["hardware"]
    expected = [
        f"{printed('n0001')} / {printed('n0002')}",
        f"{number('n0007')} / {number('n0008')}",
        number("n0009"), number("n0010"), number("n0015"), printed("n0017"),
        f"{nominal['best local (6 parameters)']['product_HT']:.4f} / {nominal['best local (6 parameters)']['bell']:.4f}",
        f"{nominal['best two-qubit (15)']['product_HT']:.4f} / {nominal['best two-qubit (15)']['bell']:.4f}",
        f"gap ≈ {number('n0035')} / shots", number("s7gains")["slow"].split(" / ")[0],
        f"{sum(j['usage_s'] for j in ref('8d_jobs_table.json'))} s of QPU time",
        f"F = {n2['R1_N2none_product_HT'][0]:.3f} to {n2['R2_N2corr_product_HT'][0]:.3f}",
    ]
    for st in ("product_HT", "bell_30_70"):
        r = rows[st]
        lo, hi = r["fieller95"]
        expected.append(f"| `{st}` | {r['F']['none']:.3f} | {r['F']['corr']:.3f} | {r['F']['native']:.3f} | "
                        f"{r['f']:.3f} [{lo:.3f}, {hi:.3f}] |")
    missing = [e for e in expected if e not in README]
    assert not missing, missing


SCANNERS = {os.path.join("tests", "test_repo.py"), os.path.join("tests", "regression", "test_data_integrity.py")}


def text_files():
    skip = {".git", "__pycache__", ".pytest_cache"}
    for d, dirs, files in os.walk(REPO):
        dirs[:] = [x for x in dirs if x not in skip]
        for f in files:
            p = os.path.join(d, f)
            if f.endswith((".gz", ".qpy", ".npz", ".png")) or os.path.relpath(p, REPO) in SCANNERS:
                continue            # binary (data are scanned by test_data_integrity); the scanners hold the patterns
            yield p


SECRET = [r"[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+\.[A-Za-z0-9.-]*[A-Za-z]{2,}", r"(?i)api[_-]?key", r"(?i)\btoken\b",
          r"crn:v1", r"(?i)ibm_cloud", r"(?i)\binstance\s*[:=]", r"C:\\\\Users", r"C:/Users", r"/home/[a-z]"]


@pytest.mark.parametrize("pattern", SECRET)
def test_no_credentials_or_personal_paths(pattern):
    hits = []
    for p in text_files():
        txt = open(p, encoding="utf8", errors="replace").read()
        for m in re.finditer(pattern, txt):
            hits.append(f"{os.path.relpath(p, REPO)}: {m.group(0)[:40]}")
    assert not hits, hits[:10]


def test_submission_script_is_disabled_by_default():
    script = os.path.join(REPO, "scripts", "submit_hardware_jobs.py")
    src = open(script, encoding="utf8").read()
    assert not re.search(r"^(from|import) qiskit_ibm_runtime", src, re.M)      # imported only inside send()
    env = {k: v for k, v in os.environ.items() if k != "TELECORR_ENABLE_HARDWARE"}
    for args in (["--job", "R1_main_ibm_kingston", "--account", "x"],
                 ["--job", "R1_main_ibm_kingston", "--account", "x", "--enable-hardware"]):
        r = subprocess.run([sys.executable, script, *args], capture_output=True, text=True, env=env, timeout=60)
        assert r.returncode == 2 and "Disabled" in r.stdout


def test_readme_quick_start_runs(capsys):
    code = re.search(r"```python\n(.*?)```", README, re.S).group(1)
    exec(compile(code, "README quick start", "exec"), {})
    before, after = capsys.readouterr().out.split("->")
    assert min(float(x) for x in re.findall(r"[0-9.]+", after)) > 0.99 > max(float(x) for x in re.findall(r"[0-9.]+", before))
