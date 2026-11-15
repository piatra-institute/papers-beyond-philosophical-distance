# Beyond Philosophical Distance

Maps of philosophical positions record pairwise relations: distances between positions, correlations between survey answers, shared citations. This paper asks when such maps can decide whether a family of commitments can be held together. The answer is the obstruction order of the commitment class, its Helly number: two for categorical answers and ranges on ordered scales, n+1 for implication rules and parity constraints on n yes/no issues, 2n for arbitrary two-literal clauses, d+1 for convex regions of a continuous space, 2^d on a grid of ordinal positions, and at least six for convex regions evaluated on thirty occupied positions of a planar map. The orders are proved and verified by exhaustive enumeration, and the paper separates locality in commitments from locality in issues. It then measures higher-order structure in two public datasets: third-order structure is small and stable in the trolley judgments of 64,988 participants and undetectable at the sample size of a survey of consciousness researchers.

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
uv run python run_all.py
```

The run executes the tests, computes the obstruction orders by exhaustive search, verifies the lower-bound constructions, analyses the frozen data in `simulation/data/`, and writes `simulation/output/results.json`, `paper/figures/orders.pdf` and `paper/figures/data.pdf`. `simulation/fetch_data.py` re-derives the frozen data from the two Open Science Framework deposits.
