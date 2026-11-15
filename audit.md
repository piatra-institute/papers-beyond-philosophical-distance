# Audit

Dated log of editorial passes and verification runs. Newest first.

## 2026-11-15 — initial research, computation, recorded run and build

Scope: first complete version of the paper from the seed material in `chats/` (an ideation conversation and a draft manuscript, used as leads only), written by an authoring agent (Claude Code, claude-opus-5-5, family anthropic). The draft's propositions were not reused as results: the paper answers the question the draft left open (when pairwise maps suffice) and adds a test on public data.

Built:
  - `simulation/atlas/`: exact set-system engine (order by Helly-independent point sets with hereditary pruning; direct enumeration of minimal incompatible subfamilies as cross-check), commitment classes (yes/no clause classes, nominal and ordinal questions, products, grids, planar point sets with exact integer geometry), lower-bound constructions, background-logic and pairwise-hull checks, maximum-entropy models (iterative proportional fitting for three answers; penalised exact fits over all states for 8 and 11 answers), connected information, triplet tests, joint-avoidance counts, margin-preserving power analysis and the analytic design calculation; figures; 28 tests.
  - `simulation/fetch_data.py`: separate fetch script for the two Open Science Framework files (Awad et al. 2020; Francken et al. 2022), with pinned SHA-256 and derivation rules recorded in `simulation/data/MANIFEST.json`. The PhilPapers survey sites serve bot challenges and individual responses are not released; no attempt was made to bypass the challenges.

Verification:
  - `uv lock` before the first recorded run.
  - Recorded run `analysis` (`papers run ... -- uv run python run_all.py`, three data files declared with `--input`): PASS, 28 tests passed, all twelve internal checks true. Three receipts exist: the first two are archived in `verification/`; the second added the design calculation and count fields, the third relabelled Figure 2(b). Scientific values are identical across the three runs (seeded); only timings differ.
  - refs: `papers refs --verify`: 48 references, all cited IDs resolve; identity matches or advisory detail differences for every DOI; Helly (1923) and Harborth (1978) have no DOI (NOT_APPLICABLE); Nielsen and Parsons (2007) is NEEDS_REVIEW because Crossref's record lacks a year and series volume (manual identity record left to the referee).
  - sources: three verification subagents and the author checked every reference against its registry and read the supporting passages where open copies existed. Deepest depth per reference in claims.yaml: 31 full text, 9 abstract, 6 identity, 2 without a source claim (doignon1973, hoffman1979); details in sources.md.
  - claims: 213 claims (112 computation, 63 source, 26 interpretation, 6 assumption, 5 definition, 1 normative); all 151 non-source bindings PASS; the 63 source claims are left `verified: false` with honest depth, locators, quotes where read in full, and notes; coverage left pending.
  - standards: PASS (no self-reference; abstract 246 words).
  - voice: two inline-contrast candidates removed; no em-dashes; no "rather than".
  - build: 21 pages, empty build log (no missing glyphs); every page rendered at 80 dpi and inspected by the author; figures redrawn twice before the final recorded run to remove label collisions.
  - check: no FAIL; NEEDS_REVIEW only for claims coverage, the 63 unverified source claims, the editorial and visual review records and the Nielsen-Parsons identity record.
  - Not done here: independent editorial and visual review records and the claims-coverage stamp, which require a reviewer from another model family.
