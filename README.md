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

```
Initial_biomass_analysis/   early exploratory analysis, figures, per-variable stats
Cross_compare/              subset comparisons, combined NPP/AMF data
AMF_analysis/               colonization overview
Landing.Rmd                 original entry point
```

---

## `02_Paper_ThreeWay/` — active

The manuscript in progress. Uses **all 60 pots**. The finding is the three-way interaction:

> Under severe drought in monoculture, AMF increased wheat productivity ~38% while cutting
> grain ~64% (allocation 0.24 → 0.06). Under adequate water the same inoculation increased
> grain ~51%. In mixed stands neither effect occurred at any water level — the
> mutualism-to-parasitism transition is contingent on the plant's neighborhood.

```
Analyses/                   the numbered pipeline, 01 -> 08.6 (the live analysis)
Analysis_Instructions/      specs handed to Claude Code for execution
Manuscript/                 narrative, story arc, figure plan (what to write, not what to run)
Archive_SignificanceGrids/  Aug-Sep 2026 exploratory test batteries; superseded by Analyses/
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
| `06_Exploratory.Rmd` | colonization associations, fava suppression |
| `08_Figures_and_Supplement.Rmd` | full figure set + supplementary test table |
| `08.6_Research_Direction.Rmd` | current narrative document |

`wheat_audited.rds` and `loo_perm_null_dist.rds` live here; most `.Rmd` files read them from
the working directory, so knit from inside `Analyses/`.

### `Manuscript/` — read in this order

1. `what_this_paper_says.md` — one page, plain language, no statistics
2. `story_arc_and_figures.md` — outline, rigor checklist, the four key figures
3. `path_decomposition_and_diagram.md` — direct/indirect effect decomposition, conceptual diagram

---

## Open items

- Test the **+51% at control in monoculture** — load-bearing for the reversal framing, currently
  a raw cell mean at n = 5.
- Test the reversal as a **linear trend** across water within monoculture (one contrast, not three
  comparisons).
- Build **Figure 1** (three-panel effect-size figure; spec in `Manuscript/story_arc_and_figures.md`).
- Senescence timing per pot was never recorded — blocks the leading mechanistic hypothesis
  (phenological delay).

## Known broken path

`Archive_SignificanceGrids/PerPot/scratch/testing.Rmd` points at a `Data/` path that did not
resolve before the reorganization either. It is scratch; left as found.
