# Beyond Philosophical Distance

Maps of philosophical positions record pairwise relations: distances between positions, correlations between survey answers, shared citations. This paper asks when checking commitments in pairs, or k at a time, decides whether a family of commitments can be held together. The answer is the obstruction order of the commitment class, its Helly number: two for categorical answers and ranges on ordered scales when every combination of answers is admissible, n+1 for categorical answers together with implication rules and for parity constraints on n yes/no issues, 2n for arbitrary two-literal clauses, d+1 for convex regions of a continuous space, 2^d on a grid of ordinal positions, and at least six for convex regions evaluated on thirty occupied positions of a planar map. The orders are proved and verified by exhaustive enumeration, and the paper separates locality in commitments from locality in issues. For maps computed from data, it shows that the empty cells of pairwise tables exclude exactly the answer patterns outside the majority-closed hull of the observed patterns, while their numerical values exclude more, which linear programming finds and certifies. It then measures higher-order structure in two public datasets: third-order structure is small and stable in the trolley judgments of 64,988 participants, and none was detected among 223 consciousness researchers, a sample too small to detect an effect of the trolley size.

[Read the paper](https://piatra.institute/papers/beyond-philosophical-distance.pdf) · [Manuscript](paper/PAPER.md) · [Computations and analyses](simulation/README.md)

## Build

```bash
uv run build.py          # -> paper/PAPER.pdf  (vendored canonical recipe)
```

Requires `pandoc` and `xelatex` on PATH. From the workspace you can also run `papers build beyond-philosophical-distance`.

## Reproduce the results

```bash
cd simulation
uv lock
uv run python fetch_data.py      # first: downloads the OSF files, checks their SHA-256, rebuilds data/
uv run python run_all.py
```

Run `fetch_data.py` first. The per-respondent survey file `simulation/data/consciousness_items.csv` is not in the repository: the script downloads the raw file of the survey's Open Science Framework deposit into `simulation/data/raw/`, checks its pinned SHA-256, and rebuilds the derived file byte for byte, checking it against a second pinned SHA-256. The aggregated trolley counts `simulation/data/trolley_first_sessions.csv` are committed, and the script re-derives and checks them as well. `run_all.py` then executes the tests, computes the obstruction orders by exhaustive search, verifies the lower-bound constructions, computes the linear-programming bounds, analyses the data in `simulation/data/`, and writes `simulation/output/results.json`, `paper/figures/orders.pdf` and `paper/figures/data.pdf`.
