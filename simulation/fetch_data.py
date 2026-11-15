"""Fetch the two public datasets and derive the frozen inputs of the recorded run.

This script is not run by run_all.py; run it once before run_all.py. It downloads the
raw files from the Open Science Framework into data/raw/ (not committed), checks their
SHA-256 against the values pinned below, writes the derived files in data/ that
run_all.py reads, and checks that each derived file reproduces its pinned SHA-256
byte for byte:

* data/trolley_first_sessions.csv (committed): counts of the eight answer patterns to
  the Switch, Loop and Footbridge dilemmas by country and design stratum, one session
  per participant (Awad et al. 2020, OSF project mxa6z).
* data/consciousness_items.csv (not committed, rebuilt here): answers of the 307
  submissions to the survey of consciousness researchers to eleven position items and
  thirteen entity items (five-point scales, 1 = definitely yes), with the ASSC-attendance
  flag and the two screening answers (consciousness: "Do you study or research
  consciousness?", 1 = yes, 2 = no; before: "Have you filled out this survey before?",
  1 = yes, 2 = no, 3 = not asked on the 2018 forms), without participant numbers
  (Francken et al. 2022, OSF project xz3ua). The analysis applies the source's screening.
* data/MANIFEST.json: sources, retrieval date, raw and derived hashes, derivation rules.

Usage: uv run python fetch_data.py [--raw-dir DIR]
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import urllib.request
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

SOURCES = {
    "trolley": {
        "url": "https://osf.io/download/c2qx5/",
        "file": "Shared_data_allResponses.csv",
        "project": "https://osf.io/mxa6z/",
        "sha256": "2fe26d6382d43f532f48dead25563161229e1a9d4ed5df1989851c3bdfc5f320",
        "citation": "Awad, Dsouza, Shariff, Rahwan and Bonnefon (2020), PNAS 117(5):2332-2337",
    },
    "consciousness": {
        "url": "https://osf.io/download/ze2hx/",
        "file": "survey_jasp.csv",
        "project": "https://osf.io/xz3ua/",
        "sha256": "aa82da0bc0edbd37b7b1e2d4646b1c73ba71efb246776d9d0a64ae5fd4ef8e7f",
        "citation": "Francken, Beerendonk, Molenaar, Fahrenfort, Kiverstein, Seth and van Gaal (2022), Neuroscience of Consciousness 2022(1):niac011",
    },
}

POSITION_ITEMS = ["driving", "light", "phenomenal", "higher-order", "explanatorygap", "biology",
                  "evolution", "machines", "zombie", "mary", "transition"]
ENTITY_ITEMS = ["entities_amoeba", "entities_baby", "entities_bat", "entities_dog", "entities_fish",
                "entities_monkey", "entities_octopus", "entities_other_people", "entities_rock",
                "entities_tree", "entities_thermostat", "entities_worm", "entities_yourself"]
SCENARIOS = ("Switch", "Loop", "Footbridge")
SCREENING = ["assc", "consciousness", "before"]
RETRIEVED = "2026-11-15"
DERIVED_SHA256 = {
    "trolley_first_sessions.csv": "58aad63805e2ec8eb49a3a33670a126d518a75fe0767484916077990dd3ce4b3",
    "consciousness_items.csv": "91a4995d048aef97fb2e2a05e7477ec6b1e813a6244aea965c7d092d40c7d95a",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "research-data-fetch/1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp, open(dest, "wb") as out:
        out.write(resp.read())


def derive_trolley(raw: Path, out: Path) -> dict:
    """One complete session per participant: the first in file order whose three
    dilemmas were each answered exactly once.

    Besides the answer pattern, each participant's presentation order, the screen side
    of the sacrifice option in each dilemma (left when answerLeft equals Outcome, the
    coding of the original analysis), device template and whether the session followed
    the main Moral Machine task (seenOther) are kept so that the analysis can stratify.
    """
    sessions: dict[tuple[str, str, str], dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    order: list[tuple[str, str, str]] = []
    with open(raw, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            key = (row["UserID"], row["Session_id"], row["country_code"])
            if key not in sessions:
                order.append(key)
            sessions[key][row["Scenario"]].append(row)
    seen_users: set[str] = set()
    counts: dict[tuple, int] = defaultdict(int)
    n_complete = n_used = n_mixed = 0
    for key in order:
        answers = sessions[key]
        if not all(len(answers.get(s, [])) == 1 for s in SCENARIOS):
            continue
        n_complete += 1
        user = key[0]
        if user in seen_users:
            continue
        seen_users.add(user)
        rows = [answers[s][0] for s in SCENARIOS]
        templates = {r["Template"] for r in rows}
        seen = {r["seenOther"] for r in rows}
        if len(templates) != 1 or len(seen) != 1:
            n_mixed += 1
            continue
        n_used += 1
        pres = "".join(s[0] for s, r in sorted(zip(SCENARIOS, rows), key=lambda z: int(z[1]["Scenario_order"])))
        sides = "".join("L" if (r["answerLeft"] == "TRUE") == (r["Outcome"] == "1") else "R" for r in rows)
        pattern = tuple(int(r["Outcome"]) for r in rows)
        counts[(key[2], pres, sides, templates.pop(), seen.pop()) + pattern] += 1
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["country_code", "order", "sacrifice_side_SLF", "template", "after_moral_machine",
                    "switch", "loop", "footbridge", "count"])
        for k in sorted(counts):
            w.writerow(list(k) + [counts[k]])
    return {"sessions_complete_single_answer": n_complete, "participants_used": n_used,
            "participants_dropped_inconsistent_session_fields": n_mixed,
            "rule": "first complete session per UserID in file order; sessions with a repeated "
                    "answer to any dilemma excluded; Outcome 1 = the sacrifice was accepted; "
                    "order = presentation order by Scenario_order (S, L, F initials); "
                    "sacrifice side = L when answerLeft equals Outcome"}


def derive_consciousness(raw: Path, out: Path) -> dict:
    """Columns copied verbatim for every submission; participant numbers dropped."""
    with open(raw, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh, delimiter=";"))
    header, body = rows[0], rows[1:]
    cols = SCREENING + POSITION_ITEMS + ENTITY_ITEMS
    idx = [header.index(c) for c in cols]
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in body:
            w.writerow([r[i] for i in idx])
    rec = [dict(zip(header, r)) for r in body]
    repeat = sum(1 for r in rec if r["before"] == "1")
    eligible = sum(1 for r in rec if r["consciousness"] == "1" and r["before"] != "1")
    return {"submissions": len(body), "repeat_submissions": repeat,
            "not_researcher": sum(1 for r in rec if r["consciousness"] == "2" and r["before"] != "1"),
            "researcher_status_missing": sum(1 for r in rec if r["consciousness"] == "" and r["before"] != "1"),
            "eligible": eligible,
            "rule": "columns copied verbatim; participant numbers dropped; assc 1 = attended ASSC; "
                    "eligible = consciousness 1 (studies or researches consciousness) and before not 1 "
                    "(not a repeat submission), the screening of Francken et al. (2022)"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", type=Path, default=DATA / "raw")
    args = ap.parse_args()
    args.raw_dir.mkdir(parents=True, exist_ok=True)
    manifest = {"retrieved": RETRIEVED, "sources": {}, "derived": {}}
    for key, src in SOURCES.items():
        dest = args.raw_dir / src["file"]
        if not dest.exists():
            download(src["url"], dest)
        digest = sha256(dest)
        if digest != src["sha256"]:
            raise SystemExit(f"{dest}: SHA-256 {digest} differs from pinned {src['sha256']}")
        manifest["sources"][key] = {k: v for k, v in src.items()} | {"bytes": dest.stat().st_size}
    t_out = DATA / "trolley_first_sessions.csv"
    c_out = DATA / "consciousness_items.csv"
    manifest["derived"]["trolley_first_sessions.csv"] = derive_trolley(args.raw_dir / SOURCES["trolley"]["file"], t_out)
    manifest["derived"]["consciousness_items.csv"] = derive_consciousness(args.raw_dir / SOURCES["consciousness"]["file"], c_out)
    for name in ("trolley_first_sessions.csv", "consciousness_items.csv"):
        digest = sha256(DATA / name)
        if digest != DERIVED_SHA256[name]:
            raise SystemExit(f"{name}: SHA-256 {digest} differs from pinned {DERIVED_SHA256[name]}")
        manifest["derived"][name]["sha256"] = digest
    manifest["licence_note"] = ("Neither OSF project states a licence. The derived files are kept for "
                                "reproduction with attribution to the original authors.")
    (DATA / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
