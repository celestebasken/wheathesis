# wheathesis

AMF × drought × neighborhood in hard red spring wheat. Greenhouse factorial, 60 pots:
stand composition (mono: 12 wheat / inter: 6 wheat + 1 fava) × AMF inoculation × water
(80 / 60 / 40% field capacity), n = 5 per cell.

The folder holds **two distinct projects** off the same experiment. They share `Data/` and
nothing else.

---

## `Data/` — shared, do not duplicate

| File | Contents |
|---|---|
| `full_dataset.csv` | 60 pots, biomass and allocation |
| `AMF_colonization.csv` | structure-level grid-intersect counts (arbuscules, vesicles, both, hyphae-only, none) |
| `Data_key_column_definitions.csv` | column dictionary |

---

## `01_Thesis_MonoOnly/` — complete, frozen

The undergraduate thesis (Spring 2026). Used **monoculture pots only**, n = 30.
Finding: AMF inoculation increased total productivity but reduced reproductive allocation
under severe drought.

Jan–Aug 2026. Kept for reference and for text that can be reused. Not under active revision.

Note (Nov 2026): the original fava CSVs had the sterile intercrop pots' (31–45) fava weights wrong
(they were copies of the wheat berry weights); corrected in `Data/` and in
`Initial_biomass_analysis/Biomass/`. The thesis results (monoculture only) are unaffected.

```
Initial_biomass_analysis/   early exploratory analysis, figures, per-variable stats
Cross_compare/              subset comparisons, combined NPP/AMF data
AMF_analysis/               colonization overview
Landing.Rmd                 original entry point
```

---

## `02_Paper_ThreeWay/` — active

The manuscript in progress. Uses **all 60 pots**. The finding is the three-way interaction:

> Under severe drought in monoculture, AMF increased wheat productivity 38% (95% CI 11–73%) while
> cutting grain mass 69% (49–80%; allocation 0.24 → 0.06). Under adequate water the same
> inoculation increased grain mass 66% (3–167%), with no supported shift in allocation. In mixed
> stands neither effect occurred at any water level — the mutualism-to-parasitism transition is
> contingent on the plant's neighborhood.

Effect sizes are model estimates from `Analyses/10_Rigor_Reversal/`; 0.24 and 0.06 are raw cell means.

```
Analyses/                   the numbered pipeline, 01 -> 20 (the live analysis)
Analysis_Instructions/      specs handed to Claude Code for execution
Manuscript/                 narrative, story arc, figure plan (what to write, not what to run)
```

### `Analyses/` running order

| File | Purpose |
|---|---|
| `01_data_prep.R` | builds `wheat_audited.rds` from `Data/` |
| `01_Data_Audit.Rmd` | balance, variable construction, sanity checks |
| `02_Allometry_Gate.Rmd` | **gate** — is the effect a size artifact? (passed) |
| `03_Primary_Model.Rmd` | beta regression, planned water contrasts |
| `04_Fragility_Gate.Rmd` | **gate** — leave-one-out, influence |
| `04b_Addendum_Followups.Rmd` | common-support sensitivity, permutation, exact tests |
| `05_Mechanism_Arb_Vesc.Rmd` | arbuscule:vesicle mechanism test (null) |
| `06_Exploratory.Rmd` | colonization associations; fava biomass vs colonization and inoculation (null with corrected data) |
| `08_Figures_and_Supplement.Rmd` | full figure set + supplementary test table |
| `08.6_Research_Direction.Rmd` | current narrative document |
| `10_Rigor_Reversal/` | rigor checks for the reversal framing (simple effects, trend tests) |
| `11_Figure1_Reversal/` | Figure 1: the reversal |
| `12_Narrative_Update/` | abstract draft (markdown, nothing to knit) |
| `13_Figure2_Allometry/` | Figure 2: allometry |
| `15_Figure4_PathDiagram/` | conceptual path diagram — **Figure 3 in the paper** (folder name predates renumbering) |
| `16_Fava_Neighbor/` | **supplementary** (Fig S-fava), exploratory: fava biomass, competitive balance, fava as a covariate |
| `17_Fava_As_Treatment/` | stand type (fava present) as a treatment for wheat berry allocation |
| `20_AMF_Gram_Scale/` | AMF effects on berry, vegetative and total mass in grams per pot |

Numbers are not contiguous (07, 09, 14, 18, 19 were retired); files keep their original numbers.
The numbered subfolders are knit from inside their own folder.

`wheat_audited.rds` and `loo_perm_null_dist.rds` live here; most `.Rmd` files read them from
the working directory, so knit from inside `Analyses/`.

### `Manuscript/` — read in this order

1. `what_this_paper_says.md` — one page, plain language, no statistics
2. `story_arc_and_figures.md` — outline, rigor checklist, the three main figures + supplement
3. `path_decomposition_and_diagram.md` — direct/indirect effect decomposition, conceptual diagram

---

## Open items

- ~~Test the +51% at control in monoculture~~ — done in `10_Rigor_Reversal/`: model estimate
  +66% grain mass (95% CI 3–167%); the allocation shift at control is not supported.
- ~~Test the reversal as a linear trend~~ — done in `10_Rigor_Reversal/` (trend tests in
  `outputs/reversal_trend_tests.csv`).
- ~~Build Figure 1~~ — done in `Analyses/11_Figure1_Reversal/`.
- Senescence timing per pot was never recorded — blocks the leading mechanistic hypothesis
  (phenological delay).
