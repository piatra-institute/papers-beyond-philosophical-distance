# Audit

Dated log of editorial passes and verification runs. Newest first.

## 2026-11-15 — revision round 4 (fourth closure check)

Scope: the fourth closure check (Codex; `closure-round4.json`, `gate-round4.json`) closed all eight findings, supported all 69 source claims and passed all 25 pages, and raised one major issue, one minor and a coverage gap; handled by the authoring agent (Claude Code, claude-opus-5-5, family anthropic). Responses are in `fix-round4.json`. Nothing was committed; metadata.yaml, the Nielsen and Parsons record and the availability section were left unchanged. No code changed, so the recorded run was not repeated.

Changes:
  - Calibration (closure N1): the methods now say that the harmonic-corrected variant would control the false discovery rate under any dependence if the p-values were valid, that the Monte Carlo p-values are simulated from a penalized pairwise model fitted to the same respondents and their calibration was not established, so that rate is nominal too, and that the variant's lower thresholds make its rejections a subset of the step-up procedure's. The results sentence no longer calls the variant valid under any dependence. The passage of src-by-any was updated to its rewritten sentence; its citation, locator, quote, depth and note are unchanged and it stays `verified: false`.
  - Enumeration (closure N2): the description now says that three issues are enumerated over all 255 supports, while four issues use one representative of each of the 401 classes under permutations of the issues and exchanges of yes and no, weighted by class size, and that the methods and the second random distribution apply to each set computed.
  - Coverage (closure N3): a script (`coverage_pass.py`, kept in the authoring agent's scratch area) lists every numeral and spelled-out quantity in every paragraph, heading, table row and caption and matches each occurrence against the claims whose passage contains it. A computation claim binds its own numeral; other claims bind only quantities their rationale names verbatim; source claims bind what their quote or note names. Exempt and reported separately: structural numbers and labels, identifiers, names containing number words, formula constants, dates, citation locators, and a short explicit list of pronoun, temporal and enumeration uses. On the round-3 text and claims it found 389 quantities, 264 bound and 125 unbound. On the round-4 text the round-3 claims bind 257 of 391, since the rewritten sentences lost 7 bindings; after this round all 391 are bound, 134 of them newly, and 114 occurrences are exempt.

Verification:
  - claims: 443 claims (171 computation, 69 source, 149 definition, 42 interpretation, 11 assumption, 1 normative); 89 added (83 definition, 4 computation, 1 interpretation, 1 assumption), 2 replaced by rewritten ones, 3 rebound to edited sentences; all computation bindings PASS; all 69 source claims `verified: false`.
  - build: 25 pages, empty build log; changed pages 15 to 17 rendered at 80 dpi and inspected.
  - check: no FAIL (432 PASS, 72 NEEDS_REVIEW, 4 NOT_APPLICABLE); NEEDS_REVIEW only for the 69 unverified source claims, claims coverage and the editorial and visual review records.

## 2026-11-15 — revision round 3 (second closure check)

Scope: the second closure check (Codex) closed all eight findings, supported all 69 source claims and reproduced every result in a scratch rerun; three small items remained, handled here by the authoring agent (Claude Code, claude-opus-5-5, family anthropic). Responses are in `fix-round3.json` in the same ledger. Nothing was committed; metadata.yaml, the source claims, the Nielsen and Parsons identity record and the availability text were left unchanged.

Changes:
  - Figures (simulation/atlas/figures.py only). Figure 2b: the two annotations moved into a legend below the panel, clear of every line; plain-number ticks on the logarithmic axis; the smallest text rose to 14 pt in figure units. Figure 1: the exponents of the 2^n and 2^d labels are drawn as raised italic text at 12 pt instead of mathtext superscripts at 70% of 13 pt. Printed sizes measured in the built PDF with PyMuPDF: Figure 1's smallest text is 7.35 pt, the drawn exponents (previously 5.57 pt), with all other text at 7.96 pt or more; Figure 2's smallest text is 7.74 pt (the tick exponents were 5.47 pt).
  - Bindings: definition claims for "2/3" (line 35, recorded pairwise_jaccard_values), "probability zero" (proof of Proposition 5) and "expectation zero" (the duality bound behind Proposition 6).
  - Lines 194 and 196: the fourth planted triplet cannot reach the two largest sizes with a positive three-way term, the only direction simulated; a separate check of the same fitted model (not part of the recorded run) found its positive-direction information levelling off at 0.0678 bits while negative terms reach 0.222 bits.

Verification:
  - Recorded run `analysis` repeated with the same command, the same three `--input` files and the same artifacts: PASS (16:40 to 16:48 UTC, 518 s), 34 tests passed, 16 checks true. Against the round-2 receipt, which the tooling archived as `verification/analysis-20261115T164845721867Z.json`, the only changed input is `simulation/atlas/figures.py`, and the changed artifacts are the two figures and `simulation/output/results.json`, which differs from the round-2 results only in 42 wall-clock timing fields (a field-by-field comparison found no other difference).
  - claims: 356 claims (167 computation, 69 source, 68 definition, 41 interpretation, 10 assumption, 1 normative); four added (cov-jaccard-23, cov-zero-148, cov-zero-152, cov-direction) and nine rebound to the two edited passages; no source claim changed.
  - build: 25 pages, empty build log; pages 10 and 19 rendered at 80 dpi and inspected, and pages 18 and 19 checked for the caveat text.
  - check: no FAIL (345 PASS, 72 NEEDS_REVIEW, 4 NOT_APPLICABLE); NEEDS_REVIEW only for the 69 unverified source claims, claims coverage and the editorial and visual review records.

## 2026-11-15 — revision round 2 (closure check)

Scope: response to the independent closure check (Codex, gpt-6-astra; ledger `revision/2026-11/ledger/beyond-philosophical-distance/closure-round1.json` and `gate-round1.json`), which closed all eight round-1 findings, passed all 25 pages and supported 61 of 71 source claims, by the authoring agent (Claude Code, claude-opus-5-5, family anthropic). Responses are in `fix-round2.json` in the same ledger. Nothing was committed. The coordinator rewrote the v1 commit without the survey file and set `status: published` and a `web_blurb` in metadata.yaml; both were left alone, as was the Nielsen and Parsons identity record, which Codex re-verified.

Changes:
  - False discovery rate (closure N1): the screening is described as the Benjamini-Hochberg step-up procedure at a nominal 5%. The text states that the guarantee holds for independent tests and under positive regression dependence (Benjamini and Yekutieli 2001, Theorem 1.2) and was not established for the dependent triplet tests, and reports the simulated error rate under the fitted model as an estimate only. Benjamini and Yekutieli's variant for any dependence (Theorem 1.3) rejects a subset of what the step-up procedure rejects, so the null result also holds under it. Reference benjamini2001 added and verified against Crossref. The figure caption, the power paragraph and the two results sentences were reworded accordingly.
  - Greenberg: "convex and continuous preferences". Abramsky and Brandenburger: contextuality is the absence of a global distribution with the given marginals; strong contextuality, the absence of any consistent global assignment, is named separately.
  - Sources an independent reviewer could not inspect: src-helly, src-gardenfors, src-hansson and src-rescher removed (originals stay cited as attribution); src-harborth bound to Heule and Scheucher (open arXiv, p. 2); src-baker bound to Feder and Vardi's author preprint (Section 6.1.2, Theorem 24); src-thomson kept with the open Yale Law School repository copy, quote and page locators. Every other source claim is unchanged.
  - Coverage: 60 new definition claims (65 in all) and one interpretation claim bind the quantities written as words, the worked-example numbers and the method constants named by the reviewer, each with a rationale naming the proposition or construction and, where the run holds the value, its pointer; new computation claims bind the digits in the power paragraph and the Figure 2 caption. The numeral matcher of the tooling accepts digits only, so no number word can carry a computation claim.
  - Minors: "at least 0.0269" (rounded down); the power paragraph states the 75 datasets at the two largest sizes; the PhilPapers sentence now gives the effect needed for 80% power.
  - Abstract unchanged (246 words). No code changed, so the recorded run was not repeated; `verification/analysis.json` (15:08 to 15:16 UTC) still binds every computation claim.

Verification:
  - claims: 352 claims (167 computation, 69 source, 65 definition, 40 interpretation, 10 assumption, 1 normative); all computation bindings PASS; the 69 source claims stay `verified: false`; a diff against round 1 shows every source claim not named by the closure check unchanged.
  - refs: `papers refs --verify`: 54 references; benjamini2001 verified; the Nielsen and Parsons manual record carried over unchanged.
  - standards PASS (abstract 246 words); metadata PASS; build: 25 pages, empty log; the changed pages (12 and 16 to 22) rendered at 80 dpi and inspected.
  - check: no FAIL (341 PASS, 72 NEEDS_REVIEW, 4 NOT_APPLICABLE); NEEDS_REVIEW only for the 69 unverified source claims, claims coverage and the editorial and visual review records.

## 2026-11-15 — revision round 1 (independent referee report)

Scope: response to the Codex referee report (gpt-6-astra, ledger `revision/2026-11/ledger/beyond-philosophical-distance/referee-codex.json`; verdict "useful synthesis"; findings F1-F8) by the authoring agent (Claude Code, claude-opus-5-5, family anthropic). Per-finding responses are in `revision/2026-11/ledger/beyond-philosophical-distance/fix-round1.json`. Nothing was committed; `simulation/data/consciousness_items.csv` was removed from the index (`git rm --cached`, staged) and gitignored by owner decision.

Changes:
  - F1: false converse of Proposition 5 removed; Proposition 5 restated at the level of supports with its true converse; new Proposition 6 (facial set of the observed margins) with credit to Fienberg and Rinaldo (2012), Eriksson et al. (2006) and Wainwright and Jordan (2008). New `simulation/atlas/lpbounds.py`: one LP per profile (HiGHS), dual certificates checked in exact integer arithmetic, an independent one-LP facial-set computation, and the exhaustive comparison with the majority hull on all supports over three and four issues (four issues through the 401 orbits of item permutations and flips).
  - F2: completeness restricted to partial answer profiles on the full product of answer sets; new paragraph on what correlation and respondent-distance maps record; subtitle changed to "When Pairwise Checks Decide Joint Compatibility" (metadata.yaml title and CITATION.cff regenerated by `papers metadata` and `papers citation`; status untouched).
  - F3: source screening applied (consciousness = 1, before != 1): 232 eligible submissions, 166 eligible attendees, 223 complete respondents per battery; the v1 description of the transition item's answer formats was reversed and is corrected.
  - F4: corrected resampled variance (Nadeau and Bengio) for all cross-validated differences.
  - F5: single-triplet power labelled as unadjusted and pre-specified; discovery power of the full screening simulated by rerunning the whole procedure (penalty selection, pairwise fit, Monte Carlo p-values, Benjamini-Hochberg) on 550 planted and 100 null datasets in parallel workers with one seed per dataset; Monte Carlo draws raised to 9,999; the exact one-dimensional no-three-way fit replaces iterative scaling inside the tests.
  - F6: exclusion language replaced; simultaneous Bonferroni-adjusted Wald intervals for every triplet's three-way log interaction.
  - F7: Proposition 3 stated for closure systems with the convention for classes without obstructions.
  - F8: abstract, introduction, README and Figure 1(a) label corrected.
  - Data handling: `fetch_data.py` pins the raw and the derived SHA-256 and fails on any mismatch; the retrieval date is fixed so that reruns reproduce `MANIFEST.json`; the survey file was deleted and rebuilt from the raw file, byte-identical (SHA-256 91a4995d...c7d95a). README reproduce sections say to run `fetch_data.py` first.
  - References: five added (wainwright2008, fienberg2012, eriksson2006, nadeau2003, benjamini1995), each checked against Crossref; nielsen2007 gained its LNCS volume (4766) from the OpenLibrary edition record and the Crossref volume record.

Verification:
  - Development run (unrecorded) and recorded run `analysis` with the v1 command and the same three `--input` files; 34 tests; 16 internal checks, all true. Receipt `verification/analysis.json` (15:08 to 15:16 UTC, 491 s, PASS); the v1 receipt was archived by the tooling as `verification/analysis-20261115T151624485578Z.json`. Scientific values are identical to the development run.
  - claims: 285 claims (162 computation, 71 source, 37 interpretation, 9 assumption, 5 definition, 1 normative), regenerated from the recorded results; all computation bindings PASS; the 71 source claims, including 8 new ones, stay `verified: false` with honest depth, locators, quotes and notes; coverage pending.
  - refs: `papers refs --verify`: 53 references; the five new DOIs verified (benjamini1995 advisory: the journal's current name differs). Nielsen and Parsons (2007): manual identity record (`--record-identity`, verifier claude-code:claude-opus-5-5:anthropic, source the Crossref volume record) after the automatic check found no year in the chapter record. A check command typed by mistake first overwrote that record with placeholder text; it was re-recorded at once with the real note and verifier, and the cache now holds only the correct record.
  - metadata and citation regenerated (`papers metadata`, `papers citation`); standards PASS (no self-reference; abstract 246 words).
  - build: 25 pages, empty build log; every page rendered at 80 dpi and inspected by the author; Figure 1(a) label shortened before the recorded run to clear panel (b); one sentence on p. 18 reworded after inspection and the PDF rebuilt.
  - check: no FAIL (271 PASS, 74 NEEDS_REVIEW, 4 NOT_APPLICABLE); NEEDS_REVIEW only for the 71 unverified source claims, claims coverage and the editorial and visual review records.
  - Not done here: independent editorial and visual review records and the claims-coverage stamp, which require a reviewer from another model family; population-level uncertainty bounds for the LP exclusions.

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
