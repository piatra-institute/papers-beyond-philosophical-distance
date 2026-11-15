# Editorial record

Status: revision rounds 1 and 2 (2026-11-15) complete in response to an independent referee report and closure check; author's read of the revised text done; independent editorial and visual review pending.

## Contribution and reader

The paper answers a precise question for anyone who builds or reads maps of philosophical positions (survey correlation networks, citation maps, embeddings): when does pairwise information about commitments decide whether a family of them can be held together? The answer is the obstruction order of the commitment class, the Helly number of the class, and the paper computes it for the classes such maps use. Categorical answers and ranges have order two when every combination of answers is admissible; exclusions on a question with k options have order k; implications alone have order one, while categorical answers with implications, Horn rules and parity constraints have order n+1 on n yes/no issues; all two-literal clauses 2n; unrestricted commitments 2^n; convex regions d+1 in a continuum, 2^d on a grid, and at least six on thirty occupied planar positions. Mixing categorical answers (order 2) with implications (order 1) gives order n+1. The paper separates this commitment locality from issue locality (local consistency, near-unanimity, acyclic schemes), states the statistical counterpart in two parts (the empty cells of pairwise tables exclude exactly the profiles outside the majority-closed hull of the observed profiles, Proposition 5; their numerical values exclude the profiles outside the facial set of the observed margins, which a linear program per profile finds and certifies, Proposition 6), and measures third-order structure in two public datasets: small and stable in 64,988 trolley judgments; not detected among 223 eligible consciousness researchers, a sample at which the full screening would have missed an effect of the trolley size. A design calculation shows that the 2020 PhilPapers survey would allow the test at the trolley size for triplets chosen in advance.

What is new against the external literature: the statement of Helly-number, near-unanimity, acyclicity and judgment-aggregation results as one condition on maps of positions; the closed-form table of orders for commitment classes with short self-contained proofs; the mixture effect (categorical answers and implications); the occupancy reading of planar maps, which combines Hoffman's characterization with the empty-hexagon theorem; the pairwise-hull proposition (Baker-Pixley at the level of supports, with a true converse) and its numerical counterpart, the comparison of hull and facial set, exhaustive on three and four issues and on the data; the background-versus-commitment distinction; the exhaustive verification; and the measurements, including the design calculation for the PhilPapers survey. The individual mathematical facts are known in their home fields and are credited there: Helly, Doignon, Hoffman/Averkov, Harborth, Horton, Heule-Scheucher; Baker-Pixley; the Helly-type reading of near-unanimity is Feder and Vardi's (l-Helly property, Theorem 24 of their preprint); the 2n bound for minimally unsatisfiable two-literal clause sets follows from Abbasizanjani and Kullmann's classification (Lee gives 4n - 2), and Lee's construction gives order n^k for k-literal clauses; the quota thresholds are Dietrich and List's and Greenberg's; facial sets, likelihood zeros and the linear-programming test for the existence of the MLE are Fienberg and Rinaldo's and Eriksson, Fienberg, Rinaldo and Sullivant's, and the marginal polytope is described by Wainwright and Jordan.

## Material worth preserving

- The opening trolley triad, which is the n = 2 case of the implication chain and connects the theory to the data.
- Table 1, with the basis of each order.
- Proposition 3 and its search, which turns the order of any finite atlas into a computation on its positions.
- The occupancy paragraph: Helly's three for possible positions against at least six for occupied ones on a planar map of thirty positions.
- The separation of commitment order from issue locality, with the chain as the example.
- Propositions 5 and 6 with the two three-answer examples (exclusive or: nothing excluded; six profiles with one or two yes answers: 000 and 111 excluded by the pairwise function (T-1)(T-2)/2).
- The trolley stratification (countries, presentation orders, side placements, device, session type) and the power analysis with margin-preserving planted interactions, now including the discovery power of the full screening procedure.

## Reverse outline

1. Pairwise maps and joint commitments: states the problem with the trolley triad, situates maps of philosophy, poses the two questions, summarizes the answers and the contribution against home literatures.
2. Commitments, obstructions and the order of a class: definitions; the names of obstructions in logic, diagnosis, belief revision, aggregation and metaphilosophy; Proposition 1 (k-wise maps), distance maps and the standard parity pair; Proposition 2 (bundling); Proposition 3 (order from positions) and its search.
3. Orders of the commitment classes: Table 1; answer profiles and the product lemma; inferential links with proofs (implications, Horn, exclusions, parity, two-literal clauses, unrestricted) and the mixture effect; regions of a conceptual space (continuum, grid, occupied positions) with the occupancy interpretation; exhaustive verification and Figure 1.
4. Locality in commitments and in issues: background logic (Proposition 4, judgment aggregation thresholds, Greenberg, doctrinal paradox, exhaustive check); local consistency, width, acyclic schemes, contextuality; the chain as a case of high commitment order and low issue width; argumentation; bounded order and finite duality.
5. Pairwise statistics: connected information; Proposition 5 (pairwise hull, support level, with converse); the linear program and its dual certificates; Proposition 6 (facial set: depends only on the support, between the support and the hull); the two three-answer examples and the facets of the 2 x 2 x 2 marginal cone; exhaustive comparison on three and four issues; sampling caveat; mixtures.
6. Data: availability of PhilPapers microdata; the two public datasets, the source's screening of the survey, and methods (corrected resampled variance, 9,999 Monte Carlo draws, Benjamini-Hochberg, simulated screening power); trolley results and strata; expert batteries, triplet tests, LP bounds on the data, joint avoidances, simultaneous interaction intervals, single-triplet and screening power; Figure 2; design for the PhilPapers survey.
7. What a map should record: the practical rule by order; what correlation and respondent-distance maps record; issue-local propagation; occupancy; what the data say; limits (crisp compatibility, representation change, the open PhilPapers test).

## Strongest rival

The facts are known in their home fields (Helly numbers, Baker-Pixley and local consistency, acyclic schemes, the median property and quota rules, contextuality, connected information), so restating them for maps of philosophy adds nothing. Response: the paper credits each fact where it is used and claims only the translation, the table of orders for the commitment classes such maps use (with short proofs), the mixture effect, the occupancy bound, the support-level hull characterization and its comparison with the facial set, the exhaustive verification and the measurements. The facial-set part in particular is log-linear theory (Haberman, Eriksson et al., Fienberg and Rinaldo) applied to maps of positions, and is credited as such. The second rival: maps of philosophy are exploratory and not meant to decide compatibility. Response: they are read for unoccupied combinations and possible syntheses, which are compatibility questions, and the order says when that reading is sound.

## Decisions and remaining checks

- The draft's agent-fission example was dropped; the trolley triad illustrates the same mixture of categorical answers and a conditional principle and ties the theory to the data.
- The draft's KH² motivation and the ten clusters were not used; maps of philosophy built from surveys and citations motivate the question on the paper's own terms.
- The parity pair is kept as the standard lower-bound example and is credited as standard.
- The two-literal bound 2n and the Horn bound n+1 are proved in the text. The verification pass found the 2n bound implied by the classification of minimally unsatisfiable 2-CNFs (Abbasizanjani and Kullmann) and the weaker 4n - 2 in Lee (2009); both are credited, and Lee's Ω(n^k) construction was added as a row of Table 1.
- Feder and Vardi's l-Helly property (equivalent to near-unanimity and strict width) is credited as the prior form of Proposition 4; the paragraph on what counts as background was added.
- Schaefer (1978) could not be read; the proofs no longer rely on it (the majority-closure step is proved inline) and the classification is cited with Chen (2006), read in full.
- Jeavons, Cohen and Cooper, Barto and Kozik, Vorob'ev and Hansson are stated at the depth actually read (abstracts); Nehring and Puppe and Rescher at identity depth with notes.
- Individual PhilPapers responses were not available; two public OSF datasets were used, and the PhilPapers test is given as a design with computed power.
- Binarization of five-point answers at "probably yes" is stated as a choice in the discussion.
- Figures were redrawn before the final recorded run: direct labels in Figure 1 and Figure 2(b) to avoid collisions with curves.
- Remaining: independent claims-coverage review, verification of every source claim (`verified` stays false here), and editorial and visual review records.

## Revision round 1 (2026-11-15)

An independent referee (Codex, gpt-6-astra) rated v1 a useful synthesis and raised eight findings. Decisions:

- F1 (blocker): v1's claim that absences inside the majority hull are invisible to every pairwise statistic was false; numerical marginals can force zeros inside the hull. Proposition 5 now states the support-level result with its true converse (the uniform distribution on the hull has the same empty cells), and Proposition 6 states what numerical tables add: the facial set of the observed margins, computed by one LP per profile with exact dual certificates, independent of the frequencies given the support. "Exactly" was removed from every summary of what pairwise statistics determine. Exhaustive comparison: on three issues the numerical tables exclude profiles inside the hull for 60 of 255 supports, on four issues for 27,936 of 65,535; in the data they exclude nothing.
- F2: the completeness claim is restricted to compatibility of partial answer profiles on the full product of answer sets; a new paragraph says what correlation matrices (pairwise tables: hull and facial set, not occupied combinations) and respondent-distance matrices (how many answers are shared, not which) record. The subtitle changed from "Pairwise Maps Suffice for Answer Profiles but Not for Inferential Links" to "When Pairwise Checks Decide Joint Compatibility", which states the result without the conflation of commitment checks with data maps.
- F3: the survey analysis now applies the source's screening (self-identified consciousness researchers, no repeat submissions), reproducing its 232 eligible submissions and 166 eligible attendees; 223 complete respondents per battery. No broader sensitivity analysis is reported. The v1 sentence on the transition item had the answer formats reversed; it now says three answers coded 1, 3 and 5, with 43 analysed respondents using 2 or 4 of a five-point version.
- F4: cross-validation uncertainty uses the corrected resampled variance of Nadeau and Bengio; fold counts are descriptive.
- F5: 0.24 is labelled as the power of one pre-specified triplet without correction; the power of the screening actually applied was simulated by rerunning the whole procedure on 550 datasets with a planted term and 100 without; Monte Carlo draws were raised to 9,999 so that the first Benjamini-Hochberg threshold is attainable.
- F6: exclusion language replaced by non-detection; simultaneous Bonferroni-adjusted Wald intervals for the three-way log interaction show that no triplet's interval excludes the trolley magnitude.
- F7: Proposition 3 now assumes a closure system (any class on a finite domain after adding the domain and all intersections) and fixes the convention for classes without obstructions; the referee's open-interval example is cited as the reason.
- F8: abstract and introduction say "categorical answers together with implications"; Figure 1(a)'s label was corrected the same way.
- Data handling (owner decision): the per-respondent survey file is gitignored and rebuilt byte-identically by fetch_data.py, which pins both the raw and the derived SHA-256; the aggregated trolley counts stay in git.
- Nielsen and Parsons (2007): LNCS volume 4766 added; identity established from the Crossref chapter and volume records and the OpenLibrary edition record.

## Revision round 2 (2026-11-15)

The closure check closed all eight round-1 findings and raised one major issue and three minor ones. Decisions:

- False discovery rate: the Benjamini-Hochberg screening is kept as the procedure that was run and simulated, and its 5% is called nominal, because the dependence condition behind the guarantee (Benjamini and Yekutieli 2001, Theorem 1.2) was not established for overlapping triplets. Switching to the Benjamini-Yekutieli variant would have needed about twice the Monte Carlo draws to make its first threshold attainable and would have lowered the reported power. Its rejections are a subset of the step-up procedure's, so the null result holds under it as stated.
- Sources: claims were kept only where an independent reviewer can inspect the passage. Helly, Gärdenfors, Hansson and Rescher are attribution only; Harborth's result is bound to Heule and Scheucher and the Baker-Pixley consequence to Feder and Vardi; Thomson keeps an open copy with locators. Greenberg's continuity condition and the distinction between contextuality and strong contextuality are now in the text.
- Coverage: number words and example quantities are bound by definition claims with explicit rationales; the abstract did not change.

## Revision round 3 (2026-11-15)

The second closure check left three small items. Figure 2b's annotations moved into a legend below the panel and its logarithmic axis got plain-number ticks; Figure 1's exponents are drawn as raised text so that every piece of figure text prints at 7 pt or more. The bindings for "2/3", "probability zero" and "expectation zero" were added, and the power paragraph and Figure 2's caption now say that the fourth planted triplet cannot reach the two largest sizes with a positive three-way term, the only direction simulated.

## Revision round 4 (2026-11-15)

The variant of the screening with harmonic-corrected thresholds is described as nominal like the screening itself: its guarantee would hold for valid p-values, and the calibration of the Monte Carlo p-values was not established. The enumeration description now names the symmetry reduction for four issues. A systematic coverage pass bound every quantity in the manuscript; the exemptions are listed in the audit entry.
