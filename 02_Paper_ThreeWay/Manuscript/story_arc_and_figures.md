# Story arc, rigor checklist, and the four figures

Decision taken: the paper is the **three-way interaction**, told from the finding and then
dissected. Structure below, gaps to close before drafting, and the figure set.

---

## 1. The claim, in one paragraph

> Under severe drought in monoculture, AMF inoculation increased wheat productivity by 38% while
> reducing grain by 64%; reproductive allocation fell from 0.24 to 0.06. Under adequate water the
> same inoculation increased grain by 51%. In mixed stands with a legume neighbor, neither effect
> occurred at any water level. The mutualism-to-parasitism transition is therefore contingent on
> the plant's neighborhood, not on water alone.

Every number above is an Inoc-vs-Sterile contrast computed **within** a stand type and **within** a
water level. Density is constant within stand type, so the density confound does not affect any of
them — it only affects comparisons of mono's absolute level against inter's. This framing never
makes that comparison. State this explicitly in the methods; it defuses the paper's oldest
weakness.

---

## 2. Story arc

**Introduction** — builds up (conventional).

1. Mycorrhizal outcomes are context-dependent, spanning mutualism to parasitism.
2. That context is almost always studied as *abiotic* — water, phosphorus, light.
3. Plants rarely grow alone. The biotic neighborhood is an unexamined axis of context.
4. Question: does having a legume neighbor change where wheat falls on that continuum under
   drought?

**Methods** — design first, briefly, with the plant counts stated plainly. The reader must know
mono = 12 wheat and inter = 6 wheat + 1 fava *before* seeing any per-pot number.

**Results** — leads with the finding, then dissects. Each section answers a question the reader is
already asking.

| § | Question | Answer | Figure |
|---|---|---|---|
| R1 | What happened? | AMF reverses from benefit to cost across drought — in monoculture only | **Fig 1** |
| R2 | Productivity or allocation? | Allocation. Productivity was *up* 38% where grain fell 64% | Fig 1C |
| R3 | Is it just plant size? | No — size-adjusted, the effect strengthens | **Fig 2** |
| R4 | Is it fragile? | Robustness: LOO, permutation, distribution-free | text + appendix |
| R5 | What did it cost the neighbor? | Fava suppressed, except under severe drought | **Fig 3** |
| R6 | Why? | Unresolved. Structural mechanism tested and null | **Fig 4** |

**Discussion** — states the finding, then: what it means for the continuum literature; the stress
gradient framing; what the null mechanism rules out; the surviving hypotheses (phenological delay
leading); agronomic implications as hypotheses, not recommendations; limitations.

---

## 3. Rigor gaps to close before drafting

Ordered by how load-bearing they are.

**(a) The +51% at control is now a headline number.** It makes the reversal a reversal rather than
a one-sided penalty. It is currently a raw cell mean at n = 5. Test the simple effect of
inoculation within Mono × Control from the full model, with a CI. If it does not clear zero, the
claim softens to "no benefit under adequate water, large cost under drought" — still publishable,
but you must know which sentence you are entitled to.

**(b) Test the reversal as a trend, not as three comparisons.** The claim "AMF swings from
beneficial to costly as drought intensifies" is one linear contrast on the inoculation effect
across water, within monoculture — not three eyeballed numbers.

```r
emm <- emmeans(m_berry, ~ Inoc | Poly * Water)
ct  <- contrast(emm, "revpairwise", by = c("Poly", "Water"))
# linear trend in the Inoc effect across water, within Mono
contrast(ct, "poly", by = "Poly")[1, ]   # linear component
```

**(c) Confirm the reversal holds in all three currencies.** It appears to (NPP, grain mass, grain
fraction), which means you do not have to nominate one currency and defend it. Verify with models,
not cell means, and report all three in Fig 1.

**(d) Nominate a primary response explicitly.** Suggested: grain allocation (proportion) for the
ecological claim, grain mass per pot as the agronomic translation. Say which is primary in the
methods so the other reads as corroboration rather than as a second bite.

**(e) Re-run the confirmatory set under the new framing.** The pre-specified tests were written for
"intercropping buffers a penalty." They still apply, but the emphasis moved. Do not add new
confirmatory tests after the fact — relabel honestly: the three-way interaction remains
confirmatory, the reversal-trend contrast is new and should be labelled exploratory unless your
advisor agrees it was implicit.

**(f) Variance homogeneity across the 12 cells.** Levene's; use `varIdent(~1|Poly)` or
`dispformula = ~Poly` if mono and inter differ. This recovers real power.

---

## 4. The four figures

Collapsing, not adding. The current allocation figure and productivity figure show one result in
two currencies across two figures; they merge into Fig 1.

### Figure 1 — The reversal (three panels)

AMF's effect on (A) total NPP, (B) grain mass, (C) grain allocation, across water, by stand type,
as model-estimated effect sizes with 95% CIs and a reference line at zero.

This is the whole finding in one figure. The reader sees panel A flat-to-positive everywhere and
panel B reversing in monoculture only — which *is* the allocation story, made visible without
needing ratios.

```r
library(tidyverse); library(emmeans); library(glmmTMB); library(nlme)

tidy_contrast <- function(model, panel_label, transform = c("log", "identity")) {
  transform <- match.arg(transform)
  emm <- emmeans(model, ~ Inoc | Poly * Water)
  if (transform == "identity") emm <- regrid(emm)
  ct  <- contrast(emm, "revpairwise", by = c("Poly", "Water"))
  d   <- as.data.frame(confint(ct))
  lo  <- intersect(c("asymp.LCL", "lower.CL"), names(d))[1]
  hi  <- intersect(c("asymp.UCL", "upper.CL"), names(d))[1]
  out <- tibble(Poly = d$Poly, Water = d$Water,
                est = d$estimate, lo = d[[lo]], hi = d[[hi]])
  # log models -> % change; identity (beta) -> percentage points
  if (transform == "log") {
    out <- out %>% mutate(across(c(est, lo, hi), ~ (exp(.x) - 1) * 100))
  } else {
    out <- out %>% mutate(across(c(est, lo, hi), ~ .x * 100))
  }
  out %>% mutate(panel = panel_label)
}

m_npp   <- lm(log_NPP   ~ Poly * Inoc * Water, data = wheat_c)
m_berry <- lm(log_berry ~ Poly * Inoc * Water, data = wheat_c)
m_alloc <- glmmTMB(berry_alloc ~ Poly * Inoc * Water, dispformula = ~ Poly,
                   family = beta_family(), data = wheat_c)

eff <- bind_rows(
  tidy_contrast(m_npp,   "A. Total productivity\n(% change)",        "log"),
  tidy_contrast(m_berry, "B. Grain mass\n(% change)",                "log"),
  tidy_contrast(m_alloc, "C. Grain allocation\n(percentage points)", "identity")
) %>% mutate(panel = factor(panel, levels = unique(panel)))

ggplot(eff, aes(x = Water, y = est, color = Poly, group = Poly)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey55") +
  geom_line(position = position_dodge(0.35), linewidth = 0.8) +
  geom_point(position = position_dodge(0.35), size = 3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15,
                position = position_dodge(0.35)) +
  facet_wrap(~panel, scales = "free_y") +
  scale_color_manual(values = c(Mono = "#b2182b", Inter = "#2166ac")) +
  labs(x = NULL, y = "Effect of AMF inoculation (Inoc - Sterile)",
       title = "AMF reverses from benefit to cost across drought - in monoculture only",
       subtitle = "Model-estimated contrasts with 95% CI. Dashed line = no effect.") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.text = element_text(face = "bold", size = 10))
```

**Reading note for the caption:** every point is an inoculated-vs-sterile contrast computed within
a stand type and water level, so planting density is held constant and does not contribute to any
estimate shown.

### Figure 2 — Allometry

The existing log(berry) vs log(vegetative mass) panel grid. Keep as is. Carries R3.
Keep the non-overlapping-support caveat in the caption.

### Figure 3 — Fava suppression

The existing boxplot of fava biomass by inoculation across water, all 30 mixed pots. Carries R5.
Labelled exploratory.

### Figure 4 — Conceptual path diagram

See `path_decomposition_and_diagram.md`. Carries R6 and the limitations. Shows the supported
paths, the tested-and-null arbuscule:vesicle path, and the untested phenology path.

### To the appendix

Design schematic; colonization intensity and composition; the arbuscule:vesicle models; the
dose-response figure; raw productivity means; the full supplementary test table; all robustness
checks (LOO, permutation, exact Wilcoxon, sensitivity to the allometric slope).

---

## 5. Suggested next three actions

1. Run the rigor items (a)-(c). They decide which sentences the abstract is entitled to.
2. Build Figure 1 and look at it. If panel A is flat-to-positive everywhere while panel B reverses
   in monoculture, the paper writes itself from there.
3. Draft the abstract — 150 words, no hedging, from the paragraph in §1. If it holds together, the
   outline above is the paper.
