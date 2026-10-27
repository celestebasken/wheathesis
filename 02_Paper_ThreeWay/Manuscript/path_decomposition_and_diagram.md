# Path decomposition and the conceptual diagram

Two deliverables. Neither is a full SEM — see the note at the end for why.

---

## 1. Decomposing the total AMF effect into size-mediated and size-independent parts

**What this replaces.** Gate 2 currently reports "p improved from 0.022 to 0.001 after size
adjustment." That's the same information, but it buries the interesting structure: the two
components have *opposite signs* and partially cancel.

**What it says.** The total effect of inoculation on grain mass is the sum of:

- a **positive** indirect component — AMF increases vegetative mass, and larger plants make more
  grain, and
- a **negative** direct component — at a given vegetative mass, inoculated plants make less grain.

The negative component is larger. Reporting them separately makes the mechanism legible.

### Implementation

```r
library(boot)

d <- wheat %>% filter(Poly == "Mono")     # or the full dataset with Poly in both models
contrasts(d$Water) <- cbind(severe = c(-0.5, -0.5, 1), mild = c(-1, 1, 0))

# Path a: treatment -> mediator (conditional on water)
# Path b: mediator -> outcome;  Path c': treatment -> outcome, size held constant
decompose <- function(dat, idx = seq_len(nrow(dat))) {
  s  <- dat[idx, ]
  ma <- lm(log_veg   ~ Inoc * Water, data = s)
  mb <- lm(log_berry ~ log_veg + Inoc * Water, data = s)

  b  <- coef(mb)["log_veg"]
  out <- c()
  for (w in levels(s$Water)) {
    nd0 <- data.frame(Inoc = factor("Strl", levels = levels(s$Inoc)),
                      Water = factor(w, levels = levels(s$Water)))
    nd1 <- nd0; nd1$Inoc <- factor("Inoc", levels = levels(s$Inoc))

    a_w  <- predict(ma, nd1) - predict(ma, nd0)        # AMF effect on size at water w
    veg0 <- predict(ma, nd0)
    c_w  <- predict(mb, cbind(nd1, log_veg = veg0)) -
            predict(mb, cbind(nd0, log_veg = veg0))    # direct, size held fixed

    out <- c(out, setNames(c(a_w * b, c_w, a_w * b + c_w),
                           paste0(w, c("_indirect", "_direct", "_total"))))
  }
  out
}

set.seed(1)
bt <- boot(d, decompose, R = 2000,
           strata = interaction(d$Inoc, d$Water))      # preserve the design

est <- decompose(d)
ci  <- t(sapply(seq_along(est), function(i)
          boot.ci(bt, index = i, type = "perc")$percent[4:5]))
results <- data.frame(effect = names(est), estimate = est,
                      lo = ci[, 1], hi = ci[, 2])
```

### Reporting

Report the three numbers per water level on the log scale, or back-transform to fold-differences.
The sentence worth writing is the severe-drought one: *the indirect (size-mediated) effect of
inoculation on grain was positive, the direct effect was negative and larger, and the two partly
cancelled — which is why unadjusted allocation ratios understate the penalty.*

### Required caveats

- **This is decomposition, not causal mediation.** Vegetative and berry mass share unmeasured
  determinants (germination vigour, pot microenvironment), so the mediator path is not causally
  identified. Use "size-mediated" and "size-independent," not "mediates." Note that the existing
  Gate 2 regression makes exactly the same assumption — this only makes it explicit.
- **`veg` and `berry` are separate components**, not nested (`NPP = veg + berry`), so there is no
  mechanical dependency inflating the size path. Worth stating; a reviewer will check.
- At n = 30 with three water levels, the conditional intervals will be wide. Report them anyway;
  the severe-drought direct effect should be clearly separated from zero.

---

## 2. The conceptual path diagram (a figure, not a model)

The single most useful thing SEM thinking offers this paper. It carries the whole argument,
including the honest gaps, and commits to no statistics.

### What to draw

```
  Water ───────────────────────────────┐
    │                                  │
    ▼                                  ▼
  Inoculation ──[?]──> Colonization ──[X]──> Fungal C sink
    │                   (unchanged)          (NOT SUPPORTED)
    │
    ├──[+]──> Vegetative mass ──[+]──> Grain mass
    │          (size-mediated, positive)
    │
    ├──[−−]──> Grain mass  (direct, size-independent, LARGE)
    │
    └──[?]──> Phenology / senescence timing ──> Grain filling
               (UNTESTED - not recorded)

  Neighborhood (fava present) ──[moderates]──> the direct path
  Inoculation ──[−]──> Fava biomass ──[moderated by Water]
```

### Annotation scheme

Use line style and a legend to mark each arrow's evidentiary status:

| Status | Style | Paths |
|---|---|---|
| Supported | solid, bold | Inoc → veg mass; Inoc → grain (direct); Inoc → fava biomass |
| Tested, null | dashed, greyed, with ✕ | Colonization → functional shift (arbuscule:vesicle) |
| Untested | dotted, question mark | Inoc → phenology → grain filling |
| Moderation | arrow onto an arrow | Water and Neighborhood onto the direct path |

The dead arbuscule:vesicle path belongs in the figure. A diagram that shows what you looked for
and did not find reads as careful; one that shows only confirmed paths reads as selective.

### Code

```r
library(DiagrammeR)
grViz("
digraph paths {
  graph [rankdir = LR, fontname = Helvetica]
  node  [shape = box, style = rounded, fontname = Helvetica]

  Inoc [label = 'AMF\ninoculation']
  Veg  [label = 'Vegetative\nmass']
  Grain[label = 'Grain\nmass']
  Col  [label = 'Colonization\n(intensity + composition)']
  Phen [label = 'Phenology /\nsenescence timing', style = 'rounded,dashed']

  Inoc -> Veg   [label = '+',  penwidth = 2]
  Veg  -> Grain [label = '+',  penwidth = 2]
  Inoc -> Grain [label = '−−', penwidth = 3]
  Inoc -> Col   [label = 'no change', style = dashed, color = grey50]
  Col  -> Grain [label = '✕ not supported', style = dashed, color = grey50]
  Inoc -> Phen  [label = '?', style = dotted]
  Phen -> Grain [label = '?', style = dotted]
}
")
```

Consider using this as **Figure 1** of the manuscript, before any data. A conceptual framework
figure is conventional in ecology and would orient a reader far faster than prose.

---

## 3. Why not an actual SEM

Recorded so the question doesn't come back.

- **Covariance-based SEM (lavaan) needs n ≫ 100.** At n = 60 (or 30 subset), fit indices such as
  RMSEA and CFI are unreliable and the model is effectively saturated.
- **Piecewise SEM handles small n**, but its global test (d-separation / Fisher's C) does not
  accommodate interaction terms coherently — and the paper's central finding *is* a three-way
  interaction. Moderation would require multigroup fitting at n = 10 per group.
- **The design is randomized and factorial.** SEM's comparative advantage is observational data
  with many correlated variables; a factorial model is the better-specified tool here.

If an advisor or reviewer asks for an SEM, the decomposition above plus the diagram is the
defensible answer, and these three bullets are the reason.
