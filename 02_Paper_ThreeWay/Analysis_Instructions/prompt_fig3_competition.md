# Claude Code prompt — Figure 3B, competitive balance in intercrop pots

Copy everything below the line.

---

I want to extend Figure 3 to look at what AMF inoculation does to the **competitive balance**
between wheat and fava in the 30 intercrop pots, and how it changes across drought. Work in
`02_Paper_ThreeWay/Analyses/`, following the conventions in `fig_style.R` and the existing
`14_Figure3_FavaSuppression/`.

## Context — what I already know from the raw means

Inoculation effects in intercrop pots (Inoc minus Sterile, per pot):

| Water | Δ fava | Δ wheat | Δ total pot |
|---|---|---|---|
| Control | −6.12 g | −1.38 g | **−7.50 g (−20%)** |
| Moderate | −5.36 g | +4.90 g | −0.46 g |
| Severe | −0.44 g | +1.38 g | +0.94 g |

So wheat does **not** reliably gain when fava is suppressed. Under adequate water the pot simply
loses biomass. Please verify these numbers first and tell me if I've misread anything — the whole
framing depends on them.

## What to build

**A new sub-document `16_Figure3B_CompetitiveBalance/`**, producing a two-panel figure that
replaces the current single-panel Figure 3:

- **Panel A** — the existing fava biomass figure (reuse `14_.../fig3_redesign.R`).
- **Panel B** — pot composition: wheat and fava biomass components by inoculation × water, so the
  total pot productivity change is visible alongside which component moved. Stacked bars or
  side-by-side components, whichever reads more cleanly — show me both if it's ambiguous.

## Models

Three responses, all on the 30 intercrop pots, all fitted **once** with `Inoc * Water` and
conditioned within (never subset models):

1. `log(fava_npp) ~ Inoc * Water` — already exists in `14_`, reuse it
2. `log(W_NPP) ~ Inoc * Water` — does wheat gain?
3. `log(Pot_NPP_g) ~ Inoc * Water` — does the system gain or lose?
4. `log(W_NPP / fava_npp) ~ Inoc * Water` — competitive balance

For (4), use the **log wheat-to-fava ratio**, not `Pot_pcent_fava`. The percentage conflates fava
shrinking with wheat growing and can't distinguish them; the log ratio is symmetric and additive.
Report the ratio effect as a fold-change with 95% CI per water level.

Check Levene's across the six cells for each response; use a variance structure if it's violated.

## Conventions to hold to

- **Everything here is EXPLORATORY.** Label it as such in section headers, not just in prose. Do
  not add it to the confirmatory set.
- **Use the randomized `Inoc` factor as the predictor, not hyphal colonization.** Colonization is a
  plant response, not a treatment, and the within-inoculated dose–response is null. (If you want a
  supplementary check, colonization coverage in intercrop pots is nearly continuous — sterile spans
  0–0.47, inoculated 0.46–1.00, largest gap 0.20 — so a `log(W_NPP/fava_npp) ~ hyphal * Water`
  model within intercrop pots is defensible as an appendix item only, with association language.)
- Report effect sizes with CIs, not just p-values. n = 5 per cell.
- Follow the encoding rule in `fig_style.R`: colour = stand type, fill = inoculation.
- Anchor any annotation to the data it describes — no labels floating in reserved headroom (see the
  note at the top of `14_.../fig3_redesign.R` for what went wrong the first time).

## Two caveats that must appear in the document

1. **Inoculum confound.** Inoculated pots received live INVAM material (spores, hyphae, root
   fragments, carrier medium, co-introduced microbiota); sterile pots received nothing. A 20%
   productivity loss concentrated entirely in the legume has an obvious alternative explanation —
   something in the inoculum harmed fava directly. There is no killed-inoculum-plus-filtrate
   control. State this plainly.
2. **AMF–rhizobia carbon competition** is the more interesting candidate mechanism: fava was
   nodulating, and AMF and rhizobia both draw photosynthate, so a tripartite symbiosis can be net
   costly to the legume. Nodulation was observed at harvest but never quantified, so flag this as a
   hypothesis, not a finding.

## Deliverables

1. Verification of the table above
2. The four models with effect sizes and CIs per water level
3. The two-panel figure, saved to `figures/` as PDF and PNG
4. A proposed caption
5. A short "notes for the Results text" section, in the style of `14_`

Tell me what you find before building the figure — particularly whether the total-pot-productivity
loss at control holds up with a model behind it.
