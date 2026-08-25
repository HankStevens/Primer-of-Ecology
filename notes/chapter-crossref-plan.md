# Chapter cross-reference & appendix plan

*Drafted 2026-08-25. Non-content-changing analysis, plus a first pass of Quarto section anchors. Nothing has been committed -- review with `git diff` and a local `quarto render` before committing on `dev`.*

## 1. What's already been done

The book already gives every chapter a header id (e.g. `# Direct competition and mutualism {#comp}`), but almost none of them followed Quarto's `sec-` naming convention, which is required for the live `@sec-id` cross-reference shortcode (auto-numbered, auto-linked, stays correct if chapters get reordered). Only chapter 11 was using any of the existing ids, and only as plain markdown links with hand-written text (`[Chapter 7](#comp)`), which don't renumber themselves.

I renamed the 9 existing chapter-level ids to the `sec-` convention (`{#comp}` -> `{#sec-comp}`, etc.), renamed the one existing subsection id (`{#fr}` -> `{#sec-fr}` in ch. 8, functional response), added **26 new subsection-level anchors** at the equations/models most commonly leaned on by other chapters (list below), and updated chapter 11's 4 existing links to use the new ids -- 3 of them now use real `@sec-` crossrefs instead of hardcoded chapter numbers.

New anchors added, by chapter:

- **01-theory**: `sec-metabolic-theory-ecology` (renamed from `metabolic_scaling`)
- **02-optimal_foraging**: chapter id only (`sec-oft`)
- **03-Expo-growth**: `sec-lambda-r-relation` (lambda <-> r conversion, "Properties of geometric and exponential growth")
- **04-DIDemography**: `sec-transition-matrix`, `sec-eigenanalysis`, `sec-lambda1`, `sec-stable-stage-distribution`, `sec-reproductive-value`, `sec-sensitivity-elasticity` -- this chapter is the single most-leaned-on chapter in the book (eigenanalysis/stability method reused in 6, 7, 8, 9, 13)
- **05-DDgrowth**: `sec-logistic-growth`, `sec-stability-analysis`, `sec-msy-k2`, `sec-discrete-logistic-growth` -- logistic growth is the second most-reused chapter (6, 7, 8, 9, 07b all build on it)
- **06-space**: `sec-levins-model`, `sec-source-sink-model`
- **07-direct_comp_mut**: `sec-lv-competition`, `sec-return-time`, `sec-jacobian-stability`
- **07b-compcolstorage**: `sec-storage-effect`, `sec-comp-col-tradeoff`
- **08-consumerResource**: `sec-rosenzweig-macarthur` (plus the renamed `sec-fr`)
- **09-host_parasitoid**: `sec-nicholson-bailey-model`
- **10-disease**: `sec-sir-closed-model`
- **11-consumer-resource-comp-mut**: `sec-tilman-resource-ratio-model`
- **13-simplefoodwebs**: `sec-pimm-lawton-methods`
- **14-diversity**: `sec-neutral-theory`

Not yet done, on purpose (kept this pass mechanical/low-risk): a handful of equally important candidates sit mid-paragraph with no header to attach to (R*, Shannon/Simpson diversity, R0, interaction strength). Adding those needs a Pandoc span (`[...]{#sec-id}`) around a phrase rather than tagging a header line, which is a bit more delicate -- happy to do a second pass there. Also untouched: the other ~49 places (of the 53 I found) where a chapter still says "Chapter 3" or "recall that..." in plain prose rather than linking to the new anchors -- mechanical to convert, but touches body text in every chapter, so I held off rather than assume you want that yet.

One structural note: **chapter 1 has no chapter-level heading at all** -- its title comes only from YAML front matter (`title: Theory in Ecology`), with no `# ... {#id}` line in the body. I didn't add one, since it's not clear whether adding a heading would duplicate the title Quarto already renders from YAML. Worth a quick check in a local render if you want chapter 1 to be a link target too.

## 2. Prerequisite map

What each chapter assumes from earlier chapters, without re-deriving it:

| Chapter | Leans on (chapter -> concept) |
|---|---|
| 1. Theory | -- (opening chapter) |
| 2. Optimal Foraging | -- (no hard dependency; echoes ch. 1's "efficient theory" framing loosely) |
| 3. Density-independent growth | -- (written to stand alone) |
| 4. Demography | 3 -> discrete growth / lambda notation, for-loop projection |
| 5. Density-dependent growth | 3 -> exponential growth, r_d=b-d-bd, sparrow dataset |
| 6. Space | 3 -> r=b-d; 4 -> eigenanalysis/stable stage distribution; 5 -> logistic growth, stability-by-partial-derivative |
| 7. Competition & mutualism | 3 -> logistic growth, K=1/alpha; 4/5 -> eigenanalysis, stability method |
| 7b. Coexistence/storage | 3 -> density dependence, doubling time; 6 -> Levins metapopulation model |
| 8. Consumer-resource | 2 -> optimal foraging; 3 -> exponential growth; 5 -> logistic growth; 7 -> Jacobian/eigenanalysis, Routh-Hurwitz |
| 9. Host-parasitoid | general predator-prey framing, discrete growth/lambda, continuous-equilibrium convention |
| 10. Disease | 8 -> type I functional response/mass action; stability via Jacobian eigenanalysis |
| 11. Consumer-resource + comp/mut | 7 -> LV competition/mutualism, invasion criterion; 8 -> consumer-resource formalism, functional response types |
| 13. Food webs | 7 -> competition coefficients (alpha_ij); general -> Jacobian/community matrix, eigenanalysis, return time |
| 14. Diversity | contrasts with earlier stable-coexistence criteria; otherwise self-contained (probability/info-theory, not stability analysis) |

The pattern is clear: **chapters 4 and 5 are the book's load-bearing walls** (eigenanalysis/stability method, and logistic growth, respectively) -- almost everything downstream leans on one or both. That's exactly the material an appendix would take the most pressure off of.

## 3. Appendix sketch

Quarto's book format supports an `appendices:` list in `_quarto.yml`, separate from `chapters:`, which auto-labels sections A, B, C... None currently exists (the `0X-chapters/` folder is a dead leftover, not wired into the build). Proposed:

**Appendix A -- Local stability analysis (the Jacobian/eigenanalysis toolkit).** The single most-repeated "recall that..." in the book. Would cover: linearizing at an equilibrium via partial derivatives, building the Jacobian/community matrix, the eigenvalue stability criterion, return time, oscillation period from complex eigenvalues, and -- since it's on your to_do list already -- the damping ratio. Chapters 4-9 and 13 would point here instead of re-explaining the method each time.

**Appendix B -- Notation and mathematical conventions.** A glossary of the recurring symbols (N, r, K, lambda, alpha, beta...) and basic matrix algebra used throughout, matched to the notation you're already planning to align with Caswell for the demography chapter. This is also the natural home for the `Conj(solve(A))` vs. `eigen(t(A))` question on your to_do list.

**Appendix C -- R toolkit / quick model reference.** Two things in one place: (a) the tidyverse/map conventions you want to standardize on (also already on your to_do list), and (b) a one-page-per-model quick-reference -- equation plus a link to the chapter section where it's fully developed, using the new `@sec-` anchors. This is the piece that most directly answers "point back instead of repeating."

I'd suggest treating this as a sketch to react to rather than something to build yet -- want me to scaffold the three appendix files (headers, `appendices:` entry in `_quarto.yml`, empty section stubs under the new anchor-worthy topics) so you can see the shape before writing into them, or would you rather adjust the plan first?

## 4. Second pass: informal "Chapter N" mentions converted to `@sec-` links

Of the 53 informal cross-chapter mentions originally found, most (~34) were self-referential ("as we saw above," "recall that...") pointing back within the *same* chapter, or general math/citation recalls with no chapter attached -- converting those would mean guessing at a target the text doesn't actually specify, so they were left as plain prose.

**16 explicit mentions were converted** to live `@sec-` crossrefs (auto-numbered, auto-linked, stay correct if chapters are reordered):

- `03-Expo-growth.qmd` -- "Ch. 6" (source population) -> `@sec-source-sink-model`
- `05-DDgrowth.qmd` (4x) -- "Chapter 3" mentions (data source, $r=b-d$, $r_d$, sparrow unboundedness) -> `@sec-expo`
- `06-space.qmd` (5x) -- "Chapter 4" (eigenanalysis) -> `@sec-eigenanalysis`; "Chapter 5" (logistic model) -> `@sec-logistic-growth`; "Chapter 5" (stability rules) -> `@sec-stability-analysis`; "Ch. 5" (Darrtown sparrows) -> `@sec-ddgrowth`; "Chapter 3" (r as N->0) -> `@sec-expo`
- `08-consumerResource.qmd` (3x) -- "chapter 3" -> `@sec-expo`; "Chapter 2" (optimal foraging) -> `@sec-oft`; "(see Competition chapter)" -> `(see @sec-comp)`
- `11-consumer-resource-comp-mut.qmd` (2x, beyond the 4 already fixed in the first pass) -- "chapter 8" (the bolded "Please read chapter 8..." instruction) and "Chapter 8" (functional response description) -> `@sec-cr`
- `13-simplefoodwebs.qmd` -- "our earlier chapters on competition" -> `our earlier @sec-lv-competition`

**3 mentions were left untouched and flagged** -- the chapter number given doesn't match the topic being discussed, which looks like a leftover from an earlier chapter reordering rather than something safe for me to silently "fix" by picking a different target:

1. `07-direct_comp_mut.qmd:125` -- "*This is precisely what we saw in Chapter 3 (logistic growth)...*" Logistic growth is chapter 5's topic (`sec-logistic-growth`), not chapter 3's (density-*independent* growth). Probably should say "Chapter 5."
2. `07b-compcolstorage.qmd:282` -- "*(see Chapter 10 for more detail)*" about species-abundance/rank-abundance distributions. That's diversity-chapter (14) material, not disease-chapter (10) material. Probably should say "Chapter 14."
3. `07b-compcolstorage.qmd:958` -- "*recall that in Chapter 3, we began with an examination of negative density-dependence.*" Negative density-dependence is chapter 5's core topic; chapter 3 is specifically density-*independent* growth. Probably should say "Chapter 5" -- though this one is a bit more ambiguous than the other two, worth a second look rather than assuming.

Once you confirm the intended chapter for each, converting them to `@sec-` links is a one-line fix each.
