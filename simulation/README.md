# Beyond Philosophical Distance: computations and analyses

This directory computes the obstruction orders of commitment classes by exhaustive search, verifies the lower-bound constructions, and measures higher-order structure in two public datasets of positions on philosophical questions.

## Layout

- `atlas/setsys.py`: finite set systems as bit masks; the closure of an atom list; the search for the largest set of positions with empty intersection of the closures of its one-point deletions (the obstruction order, Proposition 3 of the paper); and a direct enumeration of minimal incompatible subfamilies used as an independent cross-check on small cases.
- `atlas/classes.py`: atomic commitments for yes/no issues (categorical answers, implications, exclusions, two-literal clauses, Horn clauses, parity constraints, unrestricted sets), for single nominal or ordinal questions, and for products of questions.
- `atlas/geometry.py`: exact integer planar geometry: convex hulls, the closure "convex hull intersected with the occupied positions", boxes, and the largest empty convex polygon by brute force.
- `atlas/constructions.py`: explicit obstructions attaining the lower bounds (implication chains, parity families, the two cycles of two-literal clauses, exclusion triples, unit-cell families on grids), the standard pair of parity families with equal pairwise Jaccard distances, and arc consistency on implication chains.
- `atlas/theory.py`: the exhaustive verification run: orders of the yes/no classes for up to five issues, bundling invariance, witnesses, single-question and product classes, grids, random planar maps, the background-logic check of Proposition 4 and the pairwise-hull check of Proposition 5.
- `atlas/lpbounds.py`: what numerical pairwise tables exclude (Proposition 6): for each answer pattern, the largest probability compatible with the pairwise tables by linear programming; for each zero, the dual certificate checked in exact integer arithmetic; the same set computed by one linear program from the support; and the exhaustive comparison with the pairwise hull for all supports on three and four items.
- `atlas/maxent.py`: maximum-entropy models of binary answers (iterative proportional fitting and an exact one-dimensional fit for three answers; penalised maximum likelihood over all states for more), connected information, the likelihood-ratio test for a three-way interaction.
- `atlas/empirical.py`: the trolley analysis (information decomposition, bootstrap, cross-validation, strata), the survey batteries (screening of eligible respondents, nested cross-validation with the corrected resampled variance of Nadeau and Bengio, triplet tests against the fitted pairwise model with 9,999 Monte Carlo draws and Benjamini-Hochberg control, simultaneous intervals for three-way interactions, joint avoidances, power for one pre-specified triplet, and the discovery power of the full screening procedure by simulation) and the analytic design calculation for the PhilPapers survey.
- `atlas/figures.py`: the two figures of the paper.
- `fetch_data.py`: separate data-collection script, not run by `run_all.py`; run it first. It downloads the raw files from the Open Science Framework into `data/raw/` (not committed), checks their SHA-256, writes the derived files in `data/` and `data/MANIFEST.json`, and checks that each derived file reproduces the SHA-256 pinned in the script.
- `run_all.py`: the recorded run. It runs the tests, the exhaustive verification and the analyses, and writes `output/results.json`, `../paper/figures/orders.pdf` and `../paper/figures/data.pdf`.
- `tests/test_atlas.py`: known orders, agreement of the two order computations on random set systems, constructions, geometry, iterative proportional fitting (including the exclusive-or example), the exact no-three-way fit, maximum-entropy margins, the hull identity, the linear-programming certificates and the closed form for three items, the corrected variance, the planted models, the screening counts of the survey and the hashes of the derived data.

## Frozen data

- `data/trolley_first_sessions.csv`: counts of the eight answer patterns to the Switch, Loop and Footbridge dilemmas (1 = the sacrifice was chosen), by country, presentation order, screen side of the sacrifice option, device and session type, with one session per participant (64,988 participants). Derived from `Shared_data_allResponses.csv` of OSF project mxa6z (Awad, Dsouza, Shariff, Rahwan and Bonnefon 2020).
- `data/consciousness_items.csv` (not committed; rebuilt by `fetch_data.py`): answers of the 307 submissions to the eleven position questions and thirteen entity questions of the survey of consciousness science (five-point scales, 1 = definitely yes), with the conference-attendance flag, the two screening answers ("Do you study or research consciousness?" and "Have you filled out this survey before?") and without participant numbers. Derived from `survey_jasp.csv` of OSF project xz3ua (Francken, Beerendonk, Molenaar, Fahrenfort, Kiverstein, Seth and van Gaal 2022). The analysis keeps the 232 eligible submissions of the source's screening: self-identified consciousness researchers who had not filled in the survey before.
- `data/MANIFEST.json`: source URLs, retrieval date, raw and derived SHA-256 hashes and the derivation rules.

Neither OSF project states a licence. The aggregated trolley counts are kept for reproduction with attribution; the per-respondent survey file and the raw files are fetched and rebuilt, not redistributed.

## Reproduce

```bash
uv lock
uv run python fetch_data.py      # first: download, check and rebuild data/ (SHA-256 pinned)
uv run python run_all.py
```

Seeds are fixed in `run_all.py` (20261115) and in the analysis modules. The discovery-power simulation runs in parallel worker processes with one fixed seed per simulated dataset, so its results do not depend on the number of workers.
