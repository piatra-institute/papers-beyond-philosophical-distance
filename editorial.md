# Editorial record

Status: author's full read complete; independent editorial and visual review pending.

## Contribution and reader

The paper answers a precise question for anyone who builds or reads maps of philosophical positions (survey correlation networks, citation maps, embeddings): when does pairwise information about commitments decide whether a family of them can be held together? The answer is the obstruction order of the commitment class, the Helly number of the class, and the paper computes it for the classes such maps use. Categorical answers and ranges have order two; exclusions on a question with k options have order k; implications, Horn rules and parity constraints have order n+1 on n yes/no issues; all two-literal clauses 2n; unrestricted commitments 2^n; convex regions d+1 in a continuum, 2^d on a grid, and at least six on thirty occupied planar positions. Mixing categorical answers (order 2) with implications (order 1) gives order n+1. The paper separates this commitment locality from issue locality (local consistency, near-unanimity, acyclic schemes), states the statistical counterpart (pairwise marginals determine the majority-closed hull of the observed answer patterns), and measures third-order structure in two public datasets: small and stable in 64,988 trolley judgments; undetectable, at quantified power, in a 294-respondent expert survey. A design calculation shows that the 2020 PhilPapers survey would make the test decisive.

What is new against the external literature: the statement of Helly-number, near-unanimity, acyclicity and judgment-aggregation results as one condition on maps of positions; the closed-form table of orders for commitment classes with short self-contained proofs; the mixture effect (categorical answers and implications); the occupancy reading of planar maps, which combines Hoffman's characterization with the empty-hexagon theorem; the pairwise-hull proposition, which uses Baker-Pixley to say exactly what pairwise statistics of observed positions can certify; the background-versus-commitment distinction; the exhaustive verification; and the measurements, including the design calculation for the PhilPapers survey. The individual mathematical facts are known in their home fields and are credited there: Helly, Doignon, Hoffman/Averkov, Harborth, Horton, Heule-Scheucher; Baker-Pixley; the Helly-type reading of near-unanimity is Feder and Vardi's (l-Helly property, Theorem 24 of their preprint); the 2n bound for minimally unsatisfiable two-literal clause sets follows from Abbasizanjani and Kullmann's classification (Lee gives 4n - 2), and Lee's construction gives order n^k for k-literal clauses; the quota thresholds are Dietrich and List's and Greenberg's.

## Material worth preserving

- The opening trolley triad, which is the n = 2 case of the implication chain and connects the theory to the data.
- Table 1, with the basis of each order.
- Proposition 3 and its search, which turns the order of any finite atlas into a computation on its positions.
- The occupancy paragraph: Helly's three for possible positions against at least six for occupied ones on a planar map of thirty positions.
- The separation of commitment order from issue locality, with the chain as the example.
- Proposition 5 and the exclusive-or example.
- The trolley stratification (countries, presentation orders, side placements, device, session type) and the power analysis with margin-preserving planted interactions.

## Reverse outline

1. Pairwise maps and joint commitments: states the problem with the trolley triad, situates maps of philosophy, poses the two questions, summarizes the answers and the contribution against home literatures.
2. Commitments, obstructions and the order of a class: definitions; the names of obstructions in logic, diagnosis, belief revision, aggregation and metaphilosophy; Proposition 1 (k-wise maps), distance maps and the standard parity pair; Proposition 2 (bundling); Proposition 3 (order from positions) and its search.
3. Orders of the commitment classes: Table 1; answer profiles and the product lemma; inferential links with proofs (implications, Horn, exclusions, parity, two-literal clauses, unrestricted) and the mixture effect; regions of a conceptual space (continuum, grid, occupied positions) with the occupancy interpretation; exhaustive verification and Figure 1.
4. Locality in commitments and in issues: background logic (Proposition 4, judgment aggregation thresholds, Greenberg, doctrinal paradox, exhaustive check); local consistency, width, acyclic schemes, contextuality; the chain as a case of high commitment order and low issue width; argumentation; bounded order and finite duality.
5. Pairwise statistics: connected information; Proposition 5 (pairwise hull); exclusive or; mixtures.
6. Data: availability of PhilPapers microdata; the two public datasets and methods; trolley results and strata; expert batteries, triplet tests, joint avoidances, power; Figure 2; design for the PhilPapers survey.
7. What a map should record: the practical rule by order; issue-local propagation; occupancy; what the data say; limits (crisp compatibility, representation change, the open decisive test).

## Strongest rival

The facts are known in their home fields (Helly numbers, Baker-Pixley and local consistency, acyclic schemes, the median property and quota rules, contextuality, connected information), so restating them for maps of philosophy adds nothing. Response: the paper credits each fact where it is used and claims only the translation, the table of orders for the commitment classes such maps use (with short proofs), the mixture effect, the occupancy bound, the pairwise-hull characterization, the exhaustive verification and the measurements. The second rival: maps of philosophy are exploratory and not meant to decide compatibility. Response: they are read for unoccupied combinations and possible syntheses, which are compatibility questions, and the order says when that reading is sound.

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
- Remaining: independent claims-coverage review, verification of every source claim (`verified` stays false here), editorial and visual review records, and a manual identity record for Nielsen and Parsons (2007), whose Crossref record lacks a year.
