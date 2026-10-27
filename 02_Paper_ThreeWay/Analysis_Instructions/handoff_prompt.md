# Handoff prompt for Claude Code

Copy everything below the line into a fresh Claude Code session started in the `wheathesis/` root.

---

I'm continuing a greenhouse-experiment manuscript you've been helping me analyse. Two things
changed since we last worked together: **I reorganized the folder**, and **I reframed the paper**.
Please read the current state before running anything.

## Read first, in this order

1. `README.md` — folder structure, the current claim, open items
2. `02_Paper_ThreeWay/Manuscript/story_arc_and_figures.md` — the reframing, story arc, rigor
   checklist, and the four-figure plan (this is the most current document)
3. `02_Paper_ThreeWay/Manuscript/what_this_paper_says.md` — one-page plain-language summary

Older context, read only if needed: `Analysis_Instructions/AMF_intercropping_analysis_plan.md`
and `addendum_after_step2.md`.

## The folder moved

- Live analysis pipeline: `02_Paper_ThreeWay/Analyses/` (was `C_Analyses/`)
- Data: `Data/` at the repo root, shared
- Knit from inside `Analyses/` — the `.Rmd` files read `wheat_audited.rds` from the working
  directory. Relative paths were updated for the new depth and verified.
- `02_Paper_ThreeWay/Archive_SignificanceGrids/` is superseded; don't build on it.

## The paper has been reframed — this matters

**Old framing:** "intercropping buffers an AMF-induced yield penalty under drought."

**New framing:** the AMF effect *reverses* across the drought gradient, and only in monoculture.
Per-pot, inoculated vs sterile:

| | Control | Moderate | Severe |
|---|---|---|---|
| Mono — NPP | +26% | +4% | **+38%** |
| Mono — grain | **+51%** | +7% | **−64%** |
| Inter — grain | −7% | +14% | −1% |

So under severe drought in monoculture, AMF made plants ~38% **bigger** and gave ~64% **less
grain** (allocation 0.243 → 0.059). In mixed stands nothing happens at any water level. The claim
is that the mutualism-to-parasitism transition is contingent on the plant's neighborhood.

**Important consequence:** every number above is an Inoc-vs-Sterile contrast computed *within* a
stand type and *within* a water level, so planting density is constant and the density confound
does not affect them. The confound only bites when comparing mono's absolute level against
inter's — which this framing never needs to do. Preserve that property in all new analyses.

**Note:** `Analyses/08.6_Research_Direction.Rmd` is written to the OLD framing. It will need
updating, but not until the tasks below are done.

## What's already established — don't redo

- Design: 2 stand × 2 inoculation × 3 water, n = 5, 60 pots, balanced, no losses
- **Allometry gate passed** — the effect strengthens after size adjustment (three-way
  p = 0.022 → 0.001), because inoculated plants were the *larger* ones
- **Fragility checks passed** — leave-one-out, permutation, exact Wilcoxon
- **Confirmatory #3 is null** — `Inoc:Water` on total NPP. This is a feature: productivity is
  unaffected, only allocation shifts
- **Mechanism test is null** — arbuscule:vesicle ratio shows no shift. Report it, don't rerun it
- **Dose–response is null** — colonization within inoculated pots is range-restricted
  (0.46–1.00, mean 0.92), so this test is uninformative rather than contradictory

## Tasks, in order

### 1. Three rigor items (these decide what the abstract can claim)

**(a) Highest priority — test the +51% at control in monoculture.** It is currently a raw cell
mean at n = 5, and it is what makes this a *reversal* rather than a one-sided penalty. Estimate
the simple effect of inoculation within Mono × Control from the **full model**, with a 95% CI.
- If it clears zero → "AMF swings from benefit to cost across the drought gradient"
- If it doesn't → "no benefit under adequate water, large cost under drought"

Tell me which before I write anything.

**(b) Test the reversal as a linear trend**, not as three separate comparisons — one contrast on
the inoculation effect across water, within monoculture.

**(c) Confirm the reversal holds in all three currencies** (NPP, grain mass, grain allocation)
using models rather than cell means.

Also: check variance homogeneity across all 12 cells; use `varIdent(~1|Poly)` or
`dispformula = ~Poly` if mono and inter differ.

### 2. Build Figure 1

Three-panel effect-size figure — AMF's effect on NPP, grain mass, and grain allocation across
water, by stand type, with 95% CIs and a reference line at zero. **Code is already written** in
`Manuscript/story_arc_and_figures.md` §4. Adapt rather than rewrite.

### 3. Then (not before)

Update `08.6_Research_Direction.Rmd` to the new framing, and draft a 150-word abstract from
`Manuscript/what_this_paper_says.md`.

## Working conventions — please hold to these

- **Never fit subset models for inference.** Fit the full model once and condition within it
  (`joint_tests(m, by = "Poly")`, `emmeans(..., by = c("Poly","Water"))`). Comparing a mono-only
  model to a full-dataset model is not a valid comparison.
- **Don't add confirmatory tests after the fact.** The reversal-trend contrast (1b) is new —
  label it exploratory. The pre-specified confirmatory set stays as it is.
- **Report effect sizes with confidence intervals**, not just p-values. n = 5 per cell.
- **Keep the confirmatory/exploratory split explicit** in every output, in the section headers
  and not only in prose.
- Flag anything that contradicts the framing above rather than working around it.

Start by reading the three documents, then tell me your plan for task 1 before running it.
