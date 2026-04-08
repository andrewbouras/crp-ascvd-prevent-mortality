# Complete Audit — crp-ascvd-mortality project (2026-04-07)

Covers both workstreams: the main hsCRP × PREVENT manuscript (Bouras, Jaiswal) and the IRI side-abstract (Jaiswal et al., 9 authors). Audit performed after IRI rebuild was completed earlier the same day.

Severity scale: **CRITICAL** (must fix before submission) / **MAJOR** (should address, reviewer will flag) / **MINOR** (documentation or polish) / **OK** (verified correct).

---

## A. IRI rebuild (this morning's work)

### A1. [OK] Rubin's-rules pooling math
Replicated independently outside the script. Per-imputation Cox betas for the all-cause model: {−0.3528, −0.3532, −0.3531, −0.3529, −0.3585}; ubar (within-imp var) = 0.00219; B (between-imp var) = 6.03e−06; T = 0.00220; pooled SE = 0.0469; FMI = 0.003. Within/between variance split is correct. Pooled HR exactly reproduces the script's reported 0.7018 (0.6402–0.7694). CV mortality HR 0.7359 (0.6222–0.8704) also reproduces. FMI is legitimately tiny because DXA MI only affects C+D subjects and the lean-mass imputation variance is small.

### A2. [OK] No row inflation in v2 outputs
- `iri_cohort_v2_long.csv`: 42,305 rows, 8,461 unique SEQN, exactly 5 _MULT_ rows per SEQN (verified across all 8,461).
- `iri_cohort_v2_perSEQN.csv`: 8,461 rows, 758 all-cause deaths, 202 CV deaths, 170 HD deaths.
- Cross-check `sum(long_cohort events) / 5 == perSEQN events`: exact match (758 = 758; 202 = 202). No duplication survived.

### A3. [OK] Constant-ALM check for non-MI cycles
For cycles I and J (no DXA MI), 100% of subjects have `ALM_kg` constant across _MULT_ rows — confirming the 5× replication is doing what it should. For cycles C and D, 81% and 73% of subjects have constant ALM, meaning ~19–27% were imputed. This matches published NHANES DXA documentation (imputation fraction for valid-DXA subjects is ~20%).

### A4. [MAJOR] Event counts in I/J cycles are essentially zero
| Cycle | N | All-cause deaths | CV deaths |
|---|---|---|---|
| 2003–2004 | 2,666 | 522 | 142 |
| 2005–2006 | 2,416 | 215 | 58 |
| 2015–2016 | 1,889 | **13** | **2** |
| 2017–2018 | 1,490 | **8** | **0** |

The "4-cycle rebuild" is, for CV mortality, effectively a 2-cycle analysis. Cycles I and J together contribute 2 of 202 CV deaths (1%) and 21 of 758 all-cause deaths (3%). This is because NDI follow-up ends Dec 31, 2019 — cycle J subjects have only ~1–2 years of potential follow-up.

**Impact:** The CV mortality story is almost entirely driven by C+D. Adding I/J adjusts the population (helpful for representativeness) but doesn't meaningfully add CV events. The script's pooled HR is *not* inflated by inclusion of I/J, but the framing "CV finding survives in a larger, multi-cycle cohort" should be honest about this.

**How to address:** Either (a) explicitly report cycle-stratified HRs so reviewers see where events come from, or (b) report both the 4-cycle primary and a "2003–2006 only" sensitivity analysis.

### A5. [MAJOR] No proportional-hazards test across cycles
Median follow-up varies dramatically: cycle C/D ≈ 14–16 yrs, cycle I ≈ 4 yrs, cycle J ≈ 1–2 yrs. The rebuild never tests PH by cycle. With the J cycle contributing 0 CV events and very few AC events, the combined-cohort Cox is effectively fitting cycles C/D while treating I/J as adjustment-only. A Schoenfeld residual test or `strata=CYCLE_SFX` sensitivity should be run before any publication.

### A6. [MAJOR] lifelines robust SE ≠ full survey design
Both the IRI rebuild and the main manuscript use `CoxPHFitter(robust=True, weights_col=...)`. This is a sandwich estimator, not a full Taylor-series survey design with PSU + stratum. For NHANES this is acceptable practice and is what's used in many published papers, but reviewers at EJPC, Circulation, EHJ, etc. will sometimes ask for `survey::svycoxph` in R as a confirmatory sensitivity. Neither workstream currently has an R svycoxph comparison.

### A7. [MINOR] IRI z-scores are unweighted
`compute_iri()` computes sex-specific z-scores per imputation from the unweighted sample mean/SD. A "population-level resilience score" would arguably use survey-weighted moments. But the original (missing) script's behavior is unknown, so there's no ground truth. Defensible either way — note in methods.

### A8. [MINOR] Logistic models produce implausibly narrow CIs
`fit_logit_per_imp` uses statsmodels `GLM(..., freq_weights=...)`. statsmodels treats `freq_weights` as if each observation represents `w` identical copies, inflating effective N. The resulting OR CIs for FAIR_POOR_HEALTH (0.796–0.804) and WALK_DIFF (0.690–0.697) are unrealistically narrow, and the pooled FMI of ~0.95–0.98 is a symptom of the between-imputation variance being large relative to the fake-tiny within-imputation variance. **These outcomes should be dropped from the IRI abstract** — they are also conceptually circular for a composite built from albumin + ALMI at the same cross-section. Already flagged in the rebuild summary.

### A9. [OK] Survey weight divisor
`WGTMEC_COMBINED = WTMEC2YR / N_CYCLES` with N_CYCLES = 4. Correct for 4 pooled cycles per CDC analytic guidance. Previous broken version inherited /5 from the 5-cycle main project, which was a real bug.

---

## B. Main hsCRP × PREVENT manuscript (v2)

### B1. [OK] Cohort count matches
`analytic_cohort.csv`: 14,124 unique SEQN, no row inflation, no _MULT_ column. No MI duplication bug here (DXA is not used in the main project). Cycles: 2003–04 (2,257), 2005–06 (2,355), 2007–08 (3,109), 2009–10 (3,325), 2015–16 (3,078). Sum = 14,124 ✓. Manuscript matches.

### B2. [OK] Mortality counts match
Cohort: 1,364 all-cause deaths, 335 CV deaths (UCOD 1 or 5), 271 HD deaths. Manuscript: 1,364 and 335 ✓. Tables 2 and 3 break down correctly by tier × CRP category and sum to the totals.

### B3. [OK] PREVENT tier distribution
- Low: 8,341 (Manuscript: 8,341 ✓)
- Borderline (5–7.5%): 1,420 + Intermediate (7.5–20%): 3,262 = 4,682 B-I (Manuscript: 4,682 ✓)
- High ≥20%: 1,101 (Manuscript: 1,101 ✓)

### B4. [CRITICAL] PREVENT "validation against preventr v0.11.0 with three reference patients" claim is not supported by the code
Manuscript methods (line 45): *"The implementation was validated against the preventr R package (version 0.11.0) [12], which reproduces the published PREVENT coefficients, using three reference patients with exact numerical agreement."*

Script reality: there is no validation block. The only comment on validation is line 181 of `nhanes_crp_analysis.py`: *"Validated against the AHA PREVENT online calculator."* No three reference patients, no expected outputs, no pass/fail assertion. The preventr R package is never invoked.

This is a direct factual misstatement in the manuscript. Either (a) add a validation block to the script with three test patients and demonstrate numerical agreement with preventr (or with the Khan 2024 Table 2 examples), or (b) rewrite the methods sentence to describe what was actually done (validation against the AHA online calculator, presumably by hand-check).

This is the same class of issue as the scaphoid-manuscript citation hallucinations (memory entry 2026-02-25) — claims that look specific and verifiable are fabricated from thin air.

### B5. [CRITICAL] Reference 8 page range is wrong
Manuscript ref 8: *Mach F et al. 2025 Focused Update of the 2019 ESC/EAS Guidelines for the Management of Dyslipidaemias. Eur Heart J. 2025;46(42):**4359-4475**. doi:10.1093/eurheartj/ehaf190.*

Verified (PMID 40878289, CrossRef, Oxford Academic landing page): page range is **4359–4378**, not 4359–4475. The DOI resolves and the rest of the citation is correct, but the end-page 4475 appears to be a transposition or fabrication. Change `4359-4475` → `4359-4378`.

This is exactly the pattern the memory file warns about: correct DOI, wrong page numbers. The citation-accuracy checklist should have caught this; it did not.

### B6. [MINOR] Reference 12 author initial
preventr R package CRAN listing is "Martin Mayer" (no middle initial visible). Manuscript writes "Mayer MG". Minor — may need to drop to "Mayer M" or verify a middle initial exists on the package DESCRIPTION file.

### B7. [MINOR] Reference 13 cannot be verified
*NCHS 2019 Linked Mortality Files Matching Methodology. CDC; 2022.* — no DOI or PMID, and the cached CDC data-linkage page did not resolve during the audit. The citation is plausible (NCHS does publish such methodology documents), but hand-check on the NCHS data-linkage portal before submission.

### B8. [MAJOR] Manuscript claims "CRP prognostic value varies across PREVENT strata" but never formally tests interaction
- Methods and abstract frame the study as testing whether "hsCRP's prognostic association with mortality differs across PREVENT risk tiers" (paraphrased from the original idea_validation).
- The analysis fits **separate Cox models per tier** and compares the point estimates informally.
- No Cox model with a `CRP × tier` interaction term, no Wald test on the interaction coefficient, no formal test of effect modification.
- The manuscript therefore does not quantitatively support its own central claim. A competent reviewer will say: *"you've shown stratum-specific HRs but not that they differ. Please fit a pooled model with an interaction term and report the interaction p-value."*
- **How to fix**: fit a single Cox model on the whole cohort with CRP category, PREVENT tier, and their interaction (plus base covariates). Report the interaction p-value. If null, reframe the abstract as "in a prespecified stratified analysis, the association was strongest in the borderline-intermediate tier" rather than implying formal effect modification.
- Note: The original `idea_critique.md` flagged this exact issue as Flaw 2 ("No Formal Power Calculation for the Interaction Test") and asked for a pilot of events by tier. The events-by-tier counts exist (Table 2/3), but the interaction test was never added.

### B9. [MAJOR] Bootstrap C-statistic CIs do not use survey weights
In `nhanes_crp_analysis.py::compute_c_statistic` (lines 825–909), the base and CRP Cox models fit in the bootstrap loop (lines 876, 881) do NOT pass `weights_col` or `robust=True`. This means:
- The primary HRs are survey-weighted (Table 2/3)
- But the C-stat confidence intervals (Table 5) are computed from unweighted bootstrap fits
- These are inconsistent. The reported ΔC CIs reflect non-survey-weighted sampling variability.

Either apply weights inside the bootstrap fits, or document explicitly in methods that C-statistic CIs are unweighted while HR CIs are weighted. Most reviewers would prefer the former.

### B10. [MAJOR] Survey-weighted Cox HRs are dramatically larger than adjusted Cox HRs
Table 2 row 13 (Adjusted Weighted, low tier, CRP >3 mg/L): **HR 2.38 (1.42–3.98), p=0.001**
Row 7 (Adjusted, same stratum): HR 1.39 (1.06–1.83), p=0.018
Row 1 (Unadjusted, same stratum): HR 1.25 (0.96–1.62), p=0.096

Manuscript line 85 says "Survey-weighted models yielded consistent but somewhat stronger associations" — a 2.38 vs 1.39 HR is not "somewhat stronger", it's a 70% relative difference on the log-HR scale. This pattern (weighting dramatically inflates the high-CRP effect in the low-risk tier) needs a mechanistic explanation. Likely candidates:

- NHANES MEC subsampling differentially oversampled certain subgroups, and the inverse-probability weighting is pulling rare high-CRP individuals who had high mortality.
- Effective N under weighting is much smaller than raw N, so CI widens dramatically — but here the weighted CI is 1.42–3.98 which is *wider* than the unweighted 1.06–1.83, consistent with smaller effective sample.
- Model instability with small weighted event counts in the low-tier/high-CRP cell (10 CV deaths on 2,282 subjects).

Reviewers will ask. Either pick one model (weighted OR unweighted) as primary and justify it, or add a sensitivity discussion explaining the discrepancy.

### B11. [MINOR] Methods text overclaims Python/lifelines versions
Manuscript (line 65): *"Python 3.12 with the lifelines package (version 0.30.1)"*
Audit environment: Python 3.11.3, lifelines 0.30.0 (same interpreter used to regenerate the Cox models).

Not a substantive issue — lifelines 0.30.0 → 0.30.1 is a patch release with no analytical differences — but the exact version string is verifiable and should match.

### B12. [MINOR] Statin sensitivity framed as "not feasible"
Manuscript limitation 4: *"statin sensitivity analyses were not feasible because prescription medication data were not available for all survey cycles."*

NHANES does have prescription data in the RXQ_RX files for most cycles — it's not completely absent. The honest framing is that *harmonizing across cycles is tedious* or that statin coverage in NHANES RXQ cycles is incomplete, not that it is "not available". The `statin` references in the script (8 total per grep) are likely hardcoded `STATIN = 0.0` placeholders in the PREVENT equation — not actual statin ascertainment. This is a real gap: the PREVENT equation was designed for populations where statin status is known, and setting everyone to 0 will systematically bias PREVENT risk upward in the ~30% of adults actually on statins. Worth acknowledging more explicitly as a limitation.

### B13. [MINOR] 2017–2018 (J) cycle never included in main project
The main manuscript stops at 2015–2016. All J-cycle files are already downloaded (I downloaded them today for the IRI rebuild) and would add ~1,500 additional adults and 1–2 years of additional mortality follow-up for the other cycles. Not a bug — the manuscript's scope statement is internally consistent — but a sensitivity extension for the rebuttal would be straightforward now that the files exist locally.

### B14. [OK] seed=42 for bootstrap reproducibility
`np.random.seed(42)` at line 833 inside `compute_c_statistic`. Bootstrap is deterministic given this seed. Manuscript claim is accurate.

### B15. [OK] Cause-specific CV Cox handling
UCOD 001 and 005 are correctly used for the CV outcome. Non-CV deaths are censored (not dropped), which is the correct cause-specific hazards setup. The manuscript's description on line 59 matches the code.

### B16. [OK] CKD-EPI 2021 race-free equation
The eGFR function in `nhanes_crp_analysis.py` uses the 2021 race-free formula with correct kappa/alpha and sex factors. Matches manuscript claim.

### B17. [OK] Sensitivity without CRP > 10 exclusion
The `build_cohort_all_crp` function exists and runs. `sensitivity_cox_all_crp.csv` contains the results. Manuscript's sensitivity analysis numbers (N = 15,660 after dropping that filter) are consistent with the cohort file.

---

## C. IRI abstract (already drafted, 9 authors)

### C1. [CRITICAL, carried over from morning] Every number in `manuscript/NHANES_IRI_abstract_revised.docx` is wrong
Headline N, death counts, HRs, CIs, p-values, ΔC-statistic — all need to be replaced with the v2 rebuild values:
- N = 8,461 (was 11,854)
- All-cause deaths 758 (was 806), CV deaths 202 (was 221)
- All-cause HR 0.702 (0.640–0.769), p < 1e-13 (was 0.71)
- CV HR 0.736 (0.622–0.870), p = 3.4e−4 (was 0.65, p < 0.001)
- HD HR 0.728 (0.605–0.876), p = 8e−4
- ΔC = +0.008 (was +0.002, framed as improved discrimination — originally an overstatement; the new number is also modest but honest)

### C2. [CRITICAL, new] Methods text in IRI abstract says "2003–2006 and 2015–2016"
Must update to "2003–2006, 2015–2018" (4 cycles: C, D, I, J). Must also add: "Cycles C and D ship DXA-derived lean mass as five multiply-imputed copies; cycles I and J ship single observed values. Models were fit independently on each imputation and pooled using Rubin's rules. Cycles I and J observed values were treated as identical across the 5 imputation slots so the pooling formula applied uniformly."

### C3. [CRITICAL, new] Abstract claims survey PSU/stratum-adjusted design; the code does not do this
Same issue as B6 for the main manuscript. The IRI abstract's methods statement about "PSU/stratum-adjusted survey design" is an overclaim for a lifelines robust-SE model. Rewrite or port to R.

### C4. [CRITICAL, new] Drop the functional outcomes
FAIR_POOR_HEALTH and WALK_DIFF should be removed from the abstract. (a) They are circular for a composite built from albumin + ALMI. (b) The statsmodels GLM `freq_weights` inference is broken for survey data (FMI ≈ 0.95–0.98, CIs implausibly narrow). The original abstract used these as secondary endpoints; the rebuild reproduces them for completeness but they should not be published.

### C5. [MAJOR, new] Cycle-stratified event counts expose a power concern
Per Section A4, 98.8% of CV deaths come from cycles C and D. If any reviewer asks for a "cycles 2015–2018 only" sensitivity, it will be impossible (2 CV events). The abstract should pre-empt this by reporting cycle-specific counts in a supplementary table, or by acknowledging the follow-up imbalance explicitly.

### C6. [OK] IRI workstream now logged
`analysis_log.md` and `decisions_log.md` were updated earlier today (2026-04-07) with the rebuild finding and decision. Memory entry `memory/project_iri_abstract.md` is current. This was one of the original audit gaps ("IRI workstream not logged") — now resolved.

---

## D. Other files checked

### D1. [OK] `idea_critique.md` still accurate
The three SERIOUS flaws identified in the idea-critique stage (mediation circularity, interaction power, PREVENT novelty vulnerability) were all partially addressed: mediation was removed (word does not appear in manuscript v2), novelty framing was strengthened. **But** the interaction power issue (Flaw 2) was never formally resolved — see B8.

### D2. [OK] `analysis/results/iri_v2_*` files all present and consistent
Listed earlier. Long cohort, per-SEQN, pooled results, c-statistic, run log, and summary doc. `iri_analysis.py` is self-contained and reruns cleanly.

### D3. [MINOR] Stale files from the broken IRI workstream remain
- `analysis/results/iri_cohort.csv` (broken, 11,854 rows with MI inflation)
- `analysis/results/iri_cohort_full.csv` (broken)
- `analysis/results/iri_model_results.md` (contains the fabricated numbers + the later duplication-inflated numbers)
- `analysis/results/iri_distribution_summary.md` (not verified, likely broken)

These are left in place per "do not delete data without user confirmation". Consider renaming them with a `_DEPRECATED_` prefix or moving to `analysis/results/_old_iri/` so that future runs don't accidentally read them.

### D4. [MINOR] Two manuscript versions, two DOCX exports
`manuscript_v1.md`, `manuscript_v1.docx`, `manuscript_v2.md`, `Bouras_CRP_ASCVD_PREVENT_Manuscript_v2.docx`, and `FINAL_COMPLETE.docx`. Worth confirming which is canonical and removing stale exports before submission to avoid accidentally submitting v1 numbers or the pre-revision DOCX.

---

## Severity summary

| # | Finding | Workstream | Severity |
|---|---|---|---|
| B4 | PREVENT validation claim is fabricated in methods | Main | **CRITICAL** |
| B5 | Ref 8 page range wrong (4475 → 4378) | Main | **CRITICAL** |
| C1 | Every IRI abstract number is wrong | IRI | **CRITICAL** |
| C2 | IRI abstract methods cycles wrong | IRI | **CRITICAL** |
| C3 | IRI abstract survey-design overclaim | IRI | **CRITICAL** |
| C4 | Drop IRI functional outcomes | IRI | **CRITICAL** |
| A4 | I/J cycles contribute ~0 CV events (power transparency) | IRI | MAJOR |
| A5 | No PH test across cycles | IRI | MAJOR |
| A6 | No full survey-design Cox (R svycoxph) | Both | MAJOR |
| B8 | No formal CRP × tier interaction test | Main | MAJOR |
| B9 | Bootstrap C-stat CIs not weight-adjusted | Main | MAJOR |
| B10 | Weighted vs adjusted HR discrepancy unexplained | Main | MAJOR |
| C5 | Cycle-stratified events should be disclosed | IRI | MAJOR |
| A7 | IRI z-scores unweighted | IRI | MINOR |
| A8 | Logistic GLM freq_weights broken CIs | IRI | MINOR |
| B6 | Ref 12 author initial | Main | MINOR |
| B7 | Ref 13 not independently verified | Main | MINOR |
| B11 | Python/lifelines version mismatch | Main | MINOR |
| B12 | Statin "not feasible" overstatement + hardcoded 0 in PREVENT | Main | MINOR |
| B13 | J cycle could be added | Main | MINOR |
| D3 | Stale broken IRI files still on disk | Both | MINOR |
| D4 | Multiple manuscript versions on disk | Main | MINOR |
| ~20 items | Various verified correct | Both | OK |

**Pre-submission blockers (CRITICAL):**
- Main manuscript: B4 (validation claim), B5 (ref 8 pages)
- IRI abstract: C1–C4

**Reviewer pre-empt (MAJOR):**
- Main: B8 (interaction test), B9 (bootstrap weighting), B10 (weighted vs unweighted HR explanation)
- IRI: A4/A5/C5 (cycle-stratified event disclosure + PH test), A6 (R svycoxph sensitivity)

---

## Recommended next actions (priority ordered)

1. **Main manuscript B4**: Add a `validate_prevent()` function to `nhanes_crp_analysis.py` with three fixed test patients and expected PREVENT risk outputs from `preventr` v0.11.0. Assert exact numerical agreement. Rerun the manuscript only after the validation block passes.
2. **Main manuscript B5**: Fix ref 8 page range `4359-4475` → `4359-4378`.
3. **Main manuscript B8**: Add a pooled Cox with `CRP_CAT × PREVENT_TIER` interaction terms. Report the interaction p-values. Reframe abstract if non-significant.
4. **IRI abstract C1–C4**: Rewrite the abstract with the v2 numbers, remove functional outcomes, update methods to 4 cycles + Rubin pooling, soften survey-design claims.
5. **Both A6**: Port one headline Cox from each workstream to R `survey::svycoxph` and report as a confirmatory sensitivity.
6. **Both**: Add a cycle-stratified sensitivity analysis to both workstreams and report in supplementary.
7. **Citation receipts process**: The ref 8 page-number error should have been caught by the existing citation-receipts checklist. Before re-adding the Mach 2025 update ref, confirm the checklist's PubMed lookup step is actually being executed and not skipped.

Audit completed 2026-04-07. Auditor: Claude Code (Opus 4.6), running `python3 /analysis/scripts/iri_analysis.py` and independent re-derivation of Cox/Rubin math against saved cohort files.
