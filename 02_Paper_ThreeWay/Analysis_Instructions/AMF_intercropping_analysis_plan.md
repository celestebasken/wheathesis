# Analysis Plan — AMF × Drought × Neighborhood in Hard Red Spring Wheat

Working title: *Legume neighbors modulate the mycorrhizal yield-allocation tradeoff in wheat under drought*

This document is written to be handed to Claude Code as a specification. It states what the
paper claims, what it cannot claim, the exact models to fit, and the diagnostics that decide
whether the central result survives.

Make all files for this new set of analyses in the folder C_Analyses

Don't commit anything to github. 

---

## 1. The design, stated honestly

| | Mono | Inter |
|---|---|---|
| Wheat | 24 planted → **12 thinned** | 12 planted → **6 thinned** |
| Fava | 0 | 3 planted → **1 thinned** |
| Pot | 3 gal, 13 kg soil | 3 gal, 13 kg soil |

Crossed with: Inoculation (inoc / sterile) × Water (control 80% FC, moderate 60%, severe 40%).
n = 5 pots per cell. **12 cells, 60 pots, residual df = 48.**

Senescence timing did not differ by treatment.

### The confound, precisely

The `Poly` factor changes two things at once: wheat density (12 → 6 per pot) and the presence
of a legume neighbor. There is no 6-wheat-no-fava control. Therefore:

- **The main effect of `Poly` is not interpretable as an "intercropping effect."** Per-plant wheat
  biomass being higher in inter pots is the expected consequence of halving intraspecific
  competition. Do not present it as facilitation.
- **The `Poly × Inoc × Water` interaction is substantially more defensible.** For density alone to
  explain it, density would have to *modulate the AMF × drought interaction*. That is a real and
  unexcludable alternative — but it is a specific, non-trivial mechanism, not a triviality, and it
  is itself interesting ("the carbon cost of AMF under drought depends on per-plant resource
  availability, irrespective of neighbor identity").

### Language rules for the manuscript

- Name the factor **"stand composition"** or **"neighborhood treatment,"** defined as
  *wheat monoculture (12 plants) vs. mixed stand (6 wheat + 1 fava)*. Define it in the abstract.
- Never write "the effect of intercropping" for a main effect.
- The claim the paper makes is: *the direction and magnitude of the AMF effect on wheat
  reproductive allocation under severe drought depends on stand composition.*
- The claim the paper does **not** make: *fava specifically, via N fixation, rescues wheat yield.*
  That requires the density control and/or tissue N.

### Quantify the density offset

Compute a **wheat-equivalent density** for inter pots: express fava biomass per pot in units of
mean per-plant wheat biomass within the same water × inoc cell. From the existing figures fava is
roughly 25% of pot biomass in inter pots, implying ~8 wheat-equivalents vs. 12 in mono, i.e. inter
pots sat at roughly two-thirds of mono density in biomass terms. Report this number. It converts a
vague objection into a stated, bounded quantity.

### The agronomic escape hatch

Substitutive intercropping *is* what farmers do — same land, different composition. Framed as an
agronomic systems comparison, 6 wheat + 1 fava vs. 12 wheat is exactly the practical contrast.
The confound only blocks *mechanistic* claims about facilitation. Consider foregrounding the
agronomic framing and treating mechanism as explicitly speculative.

---

## 2. The threat that matters more than the density confound: allometry

This applies to the undergraduate thesis result too, and it should be checked first.

Allocation measured as a ratio is confounded with plant size, because allocation is size-dependent
(ontogenetic drift; Weiner 2004, *Perspect. Plant Ecol.*; Müller et al. 2000). In this experiment
both treatments that matter also changed plant size:

- AMF inoculation increased NPP ~19.8%.
- Inter pots have larger individual wheat plants (density release).

So "AMF reduced berry allocation" and "inter buffers the reduction" could both be size artifacts.
**This must be ruled out before anything else.**

### Allometric analysis

Fit reproductive output against vegetative output on log-log axes and test whether treatments
shift the *elevation* of the allometric line, not just the ratio:

```r
# veg = ANPP + BNPP per plant; berry = berry mass per plant
m_allo <- lm(log(berry) ~ log(veg) * Poly * Inoc * Water, data = wheat)
# 1. test slope homogeneity (the log(veg):treatment terms)
# 2. if slopes are common, refit with common slope and test elevation shifts
m_allo2 <- lm(log(berry) ~ log(veg) + Poly * Inoc * Water, data = wheat)
emmeans(m_allo2, ~ Inoc | Poly * Water)   # size-adjusted allocation
```

Also run standardized major axis regression as a cross-check, since both variables carry
measurement error:

```r
smatr::sma(log(berry) ~ log(veg) * grp, data = wheat)   # grp = Poly:Inoc:Water
```

**Decision rule:** if the Inoc × Water effect on berry allocation disappears once size is
accounted for, the paper's finding is "AMF changes plant size, and allocation follows" — a real but
different and weaker story. If it survives, you have a genuine allocation-strategy shift, and you
have pre-empted the strongest reviewer objection. Report the size-adjusted result as primary
either way.

---

## 3. Primary model

One model, all 60 pots. Pot is the experimental unit (already correct in the thesis).

**Primary response:** wheat berry allocation = berry / (ANPP + BNPP + berry), per plant.
This is the response *least* mechanically distorted by density, because it is a within-plant ratio.

```r
# proportions -> beta regression with logit link, variable dispersion by stand composition
library(glmmTMB)
m1 <- glmmTMB(berry_alloc ~ Poly * Inoc * Water,
              dispformula = ~ Poly,
              family = beta_family(), data = wheat)
```

If any berry_alloc values are exactly 0 or 1, apply the Smithson–Verkuilen transform
`(y*(n-1) + 0.5)/n` before fitting, and say so. If zeros are structural (severe drought pots with
no grain), consider zero-inflated beta and report both.

**Primary test:** the `Poly:Inoc:Water` term.

**Primary contrast — this is the sentence the paper is built on.** Do not use a compact letter
display across 12 cells (66 comparisons, near-zero power). Use one difference-of-differences:

```r
emm <- emmeans(m1, ~ Inoc | Poly * Water)
contrast(emm, interaction = c("revpairwise", "revpairwise"),
         by = "Water")   # (Inoc - Strl | Inter) - (Inoc - Strl | Mono), within each Water level
```

Report the severe-drought value with a 95% CI on the response scale. That single interval is the
result.

### Water: do not assume linearity

Your existing tables repeatedly show **Moderate > Control > Severe** for mass variables. That hump
is real enough to matter, so do not fit water as a linear covariate. Instead keep it as a 3-level
factor and partition its 2 df into two orthogonal planned contrasts:

- **C1 — severe stress:** Severe vs. mean(Control, Moderate). This is where all the action is.
- **C2 — mild gradient:** Control vs. Moderate. This is where the leaching artifact lives (§7).

```r
contrasts(wheat$Water) <- cbind(severe = c(-0.5, -0.5, 1), mild = c(-1, 1, 0))
```

This recovers most of the power of a focused test without assuming monotonicity, and it makes the
interaction terms interpretable one at a time.

### Secondary responses (same model structure)

`log(NPP per plant)`, `log(ANPP)`, `log(BNPP)`, `log(berry mass)`, ANPP%, BNPP%,
pot-level total productivity.

For Gaussian responses, check variance homogeneity across all 12 cells with Levene's. If mono and
inter differ (likely — different systems), do not pool:

```r
nlme::gls(y ~ Poly * Inoc * Water, weights = varIdent(form = ~ 1 | Poly), data = wheat)
```

### Compositional constraint

ANPP% + BNPP% + berry% = 1, so three separate ANOVAs are three views of two degrees of freedom,
and a significant shift in one mechanically forces compensation in the others. Treat berry
allocation as **the** primary response (it is the agronomically meaningful quantity) and report the
other two as descriptive. As a robustness check, run an isometric log-ratio transform of the
three-part composition and a MANOVA, or `DirichletReg::DirichReg`. Mention the constraint
explicitly in the methods; reviewers in this literature notice.

---

## 4. Fragility check — run this before writing a word

The three-way interaction appears to hinge on one cell (Inoc / Mono / Severe ≈ 0.09 while the
other three sit at 0.21–0.25). With n = 5, that could be one pot.

```r
# leave-one-out refit across all 60 pots; store the three-way p and the severe-drought contrast
loo <- lapply(seq_len(nrow(wheat)), function(i) {
  m <- update(m1, data = wheat[-i, ])
  c(p3 = car::Anova(m, type = 3)["Poly:Inoc:Water", "Pr(>Chisq)"],
    est = <severe difference-of-differences estimate>)
})
```

Also: `influence.measures()`, Cook's distance, and a dotplot of the raw five pots in that cell.

**Decision rule:** if removing any single pot moves the three-way p-value above ~0.10 or halves the
contrast estimate, say so in the results and temper the abstract. A fragile result reported as
fragile is publishable; a fragile result reported as robust is a retraction risk. Include the LOO
range in a supplement regardless of outcome.

---

## 5. Colonization analyses — rebuild these

Three problems with the current versions:

1. **Bimodal predictor.** In the moderate and severe panels, colonization piles at ~0 and ~1 with
   little in between. A regression on that is an inoc-vs-sterile group contrast wearing a
   continuous-predictor costume, with R² inflated by leverage.
2. **Panel-wise p-values.** Comparing p = 0.006 / 0.001 / 0.3 across three facets is the same error
   as comparing mono-only to inter-only. The claim "the relationship depends on water" requires an
   interaction term.
3. **Colonization is a response, not a treatment.** It was not randomized. Bigger plants can support
   more hyphae, so reverse causation is live.

### Fix

Split cleanly by inferential status:

- **Inoculation = randomized → causal language permitted.** This is the treatment.
- **Colonization = observed → association language only.** Fit within inoculated pots only, where
  the predictor is continuous and reasonably spread:

```r
col <- subset(wheat, Inoc == "Inoc")
m_col <- lm(berry_alloc ~ hyphal * Water * Poly, data = col)   # beta/glmmTMB if preferred
```

Report the `hyphal:Water` interaction as the test. Label the section "associations" and use
"was associated with" throughout. Note explicitly that background colonization in sterile pots
(0–30%) biases the inoculation contrast toward the null — i.e. your treatment effects are
conservative. That is a defensible framing, not just a caveat.

### On the sterile control (clarifying the earlier comment)

Two separate issues, often conflated:

1. **What you raised — imperfect sterility.** Sterile pots showed some colonization. This shrinks
   the contrast and makes tests conservative. Fine; state it.
2. **What I meant — the inoculum is not a clean AMF-only treatment.** Inoculated pots received 32 g
   of live INVAM material containing spores, hyphae, root fragments, growth medium, and whatever
   non-AMF bacteria and fungi rode along; sterile pots received nothing. So the contrast is
   *live inoculum vs. no inoculum*, not *AMF vs. no AMF*. The standard control adds autoclaved
   inoculum plus a fine-sieved microbial wash filtrate (passes bacteria, excludes AMF propagules),
   isolating AMF as the only difference. You do not have that. Write the limitation as: "the
   inoculation contrast reflects the effect of a live AMF-dominated inoculum, and we cannot fully
   separate AMF effects from those of co-introduced microbiota or the small nutrient addition
   associated with the inoculum carrier."

---

## 6. Fava analyses (inter pots only, n = 30)

Only wheat roots were stained, so fava's own mycorrhizal status is unmeasured. This actually makes
the framing cleaner: **the strength of the wheat–AMF partnership predicts fava suppression.**

```r
fava <- subset(wheat, Poly == "Inter")
m_f <- lm(log(fava_npp) ~ hyphal * Water, data = subset(fava, Inoc == "Inoc"))
m_f2 <- lm(fava_pct_pot ~ hyphal * Water, data = subset(fava, Inoc == "Inoc"))
```

Caveats to state: n ≈ 5 per water level once restricted to inoculated pots — this is exploratory,
full stop; fava was certainly colonized but unmeasured; the flat slope under severe drought may be
range restriction (fava NPP compressed into 3–5 g) rather than a genuine absence of relationship.
Test the `hyphal:Water` interaction rather than comparing facet p-values.

### The mechanism probe worth running

You have continuous variation in fava biomass across inter pots. Use it to ask whether the rescue
is fava-specific or just density release:

```r
m_mech <- glmmTMB(berry_alloc ~ fava_npp * Inoc * Water, family = beta_family(),
                  data = subset(wheat, Poly == "Inter"))
```

- If **larger fava → stronger rescue** of the inoculation penalty under drought, that is evidence
  for facilitation (consistent with N transfer).
- If **larger fava → weaker rescue**, competition dominates and the mono/inter difference is more
  likely density release.

Not causal (fava size is not randomized and is collinear with wheat performance), but it is a real,
directional test and it is currently the only mechanistic leverage you have. Report as exploratory.

---

## 7. Two design issues to surface in the manuscript

### Explain the Moderate > Control hump, and fix a units ambiguity

Your mass variables repeatedly show **Moderate > Control > Severe**. This non-monotonicity needs an
explanation in the discussion, or a reviewer will supply a less flattering one.

The most likely explanation is **mild hypoxia in the control treatment**. 80% of field capacity is
quite wet for a pot with restricted drainage, and wheat roots respond poorly to low soil oxygen.
Moderate drought at 60% FC may simply sit closer to the optimum for root function than the control
does. This is easy to state, plausible, and does not undermine anything — the treatment you care
about is severe drought, and the two planned contrasts in §3 already isolate it from this.

**Reporting fix.** As written, the methods give 7.57 L/min emitters × 4 min × 6 days = 181.7 L/week,
which reads as per-pot volume delivered into an 11.4 L pot. Any reviewer will do that arithmetic
and conclude the pots were massively over-watered. Since you calibrated by weighing and pots held
at target field capacity, the 7.57 L/min figure is almost certainly a line or manifold flow rate
split across many emitters. **Report the delivered volume per pot per event and per week**, state
the number of pots per line, and state whether pots drained. One sentence removes an objection that
would otherwise look serious.

### Irrigation replication — you are fine, say so

Each irrigation line was present on every bench, and pots were re-randomized monthly across
benches and reconnected. This means water treatment is **not** confounded with bench, which is the
usual fatal flaw in greenhouse drought studies. Write this explicitly into the methods in one
sentence. Reviewers look for it.

### No fava monoculture → no LER

You cannot compute a land equivalent ratio. Say so plainly and substitute pot-level total
productivity as the system metric, noting it is itself density-confounded. Do not attempt a
partial or improvised LER.

---

## 8. Multiplicity and pre-specification

The existing significance tables represent roughly 150 tests across three datasets. That
transparency is an asset *if* the paper designates a small confirmatory set up front, and a
liability if it reads as though the interesting cells were found by looking.

**Designate as confirmatory (3 tests):**

1. `Poly:Inoc:Water` on wheat berry allocation (size-adjusted).
2. `Inoc:Water` on wheat berry allocation.
3. `Inoc:Water` on wheat NPP per plant.

Everything else — ANPP%, BNPP%, colonization models, all fava analyses — is **exploratory**, labeled
as such in the section headers, with the full test table in a supplement. Report effect sizes with
confidence intervals throughout; at n = 5 per cell, estimation is both more honest and more
persuasive than p-values.

---

## 9. Figures

1. **Design schematic** — pot composition, thinning, drought timing on the Zadoks scale. Given the
   density confound, showing the design plainly is a credibility move.
2. **The key figure** — berry allocation, Water × Inoc, faceted by stand composition, with 95% CIs
   (the current version has no error bars) and the difference-of-differences annotated for severe
   drought.
3. **Allometry** — log(berry) vs. log(vegetative), colored by inoculation, faceted by
   water × composition, with fitted lines. This is what makes §2 visible.
4. **Productivity** — NPP and berry mass per plant, plus pot-level total productivity.
5. **Colonization associations** — inoculated pots only, with the interaction test in the caption
   rather than three separate panel p-values.
6. **Fava suppression** — fava biomass vs. wheat hyphal colonization, clearly labeled exploratory.

---

## 10. Mechanism without sending samples out

No further pots are possible, so the density confound stands and is handled by framing (§1).
Tissue chemistry would be the textbook mechanistic addition, but you already hold the material for
a stronger and cheaper one.

### Option A — functional composition of colonization (CONFIRMED AVAILABLE, zero cost)

Arbuscules and vesicles were tallied separately during scoring. This is the most valuable
under-used asset in the dataset and it should become the mechanistic core of the paper.

**Why it matters.** Arbuscules are the site of bidirectional plant–fungus exchange: they are the
symbiosis actually *trading*. Vesicles are fungal lipid storage: the fungus banking carbon for
itself. The **arbuscule:vesicle ratio** is therefore a direct index of whether the association is
functioning as a mutualism or as a carbon sink.

**Why it is better than it looks.** Your thesis reported a null — colonization intensity did not
change across water treatments — and interpreted it as "the symbiosis was maintained under stress."
Separating structures lets you test something much sharper:

> Colonization *intensity* was unchanged across drought, but its *functional composition* shifted:
> under severe drought the association moved from exchange toward fungal storage, and this shift
> was attenuated when a legume neighbor was present.

That converts a null result into a positive mechanistic finding, and it supplies the missing causal
middle of the whole paper — why inoculated wheat lost reproductive allocation under drought.

The data for this part are in Data/AMF_colonization.csv
One thing to note is that the column A_V is a grid subsect view containing both arbuscules and vesicles 
(and likely hyphae). This means that the total vesicles observed, for example, would be the value in
the column Vesicles + the value in A+V. Just so you know for the code.

**Models (inoculated pots only, n = 30; 2 Poly × 3 Water × 5).**

```r
inoc <- subset(wheat, Inoc == "Inoc")
inoc$AV <- log((inoc$arb_pct + 0.5) / (inoc$ves_pct + 0.5))   # log arbuscule:vesicle ratio

# 1. Does functional composition shift with drought and stand composition?
m_av <- lm(AV ~ Poly * Water, data = inoc)          # test Poly:Water, plus the severe contrast

# 2. Structures modelled separately, as counts out of 100 intersects
m_arb <- glm(cbind(arb_n, 100 - arb_n) ~ Poly * Water, family = binomial, data = inoc)
m_ves <- glm(cbind(ves_n, 100 - ves_n) ~ Poly * Water, family = binomial, data = inoc)
# check overdispersion; use quasibinomial or betabinomial if phi >> 1

# 3. Confirm intensity is unchanged while composition shifts (this contrast IS the result)
m_hyph <- glm(cbind(hyph_n, 100 - hyph_n) ~ Poly * Water, family = binomial, data = inoc)

# 4. Does functional composition track the allocation penalty?
m_link <- glmmTMB(berry_alloc ~ AV * Water * Poly, family = beta_family(), data = inoc)
```

**Framing caution.** Model 4 is observational — A:V was not randomized — so this is association,
not mediation. Write "consistent with" and do not run a formal mediation analysis.

**Phenology caveat to check.** Arbuscules turn over in days and decline naturally at senescence,
while vesicles accumulate. Everything was harvested on one date once 90% of pots reached Zadoks 92.
If severe drought accelerated senescence, those pots may have been further past Z92 at harvest,
which would depress arbuscules for reasons unrelated to the symbiosis. If you have any record of
when each treatment reached Z92, check it and report it. If treatments senesced at different rates,
say so — it is a limitation, not a disqualification, because the mono-vs-inter comparison within
severe drought is still made at comparable phenology.

### Option B — yield components from stored grain (a few hours)

Your berries are in coin envelopes. Count grains per pot and compute mean individual grain weight.
This decomposes yield into **grain number** (sink establishment) and **grain weight** (sink filling),
which are set at different developmental stages. Your drought began at Zadoks 39 — right at floret
development — so the decomposition is directly interpretable:

- Fewer grains, normal weight → AMF suppressed floret set at booting.
- Normal grain number, lighter grains → AMF competed with wheat for carbon during grain filling.

These are different stories with different implications, and distinguishing them is a genuine
mechanistic contribution. Thousand-grain weight is also a standard agronomic variable, which helps
with the applied framing.

### Option C — foliar/grain N (only if a reviewer demands it)

Grain N would test the N-facilitation pathway directly. It is the strongest option and the most
work. **Do not run it preemptively.** Options A and B are sufficient for a first submission; hold N
in reserve as a revision response if a reviewer asks for it.

---

## 11. Open questions to resolve before modeling

- **Fava nodulation: resolved.** A rhizobial inoculant was applied to fava seed at planting and
  nodules were observed at harvest. N fixation was therefore active and the facilitation mechanism
  is live. Nodulation was not quantified — state this as a limitation and describe it as an
  observation at harvest, not as data. Note also that rhizobia were added only to fava, so
  nodulation is a property of the mixed-stand treatment rather than a confound with inoculation.
- **Arbuscules vs. vesicles: resolved.** Scored separately and available. See §10 Option A, which
  is now the mechanistic core of the paper rather than an optional extra.
- **5.9% vs. ~9%: resolved** — the same quantity computed two ways. Pick one definition
  (berry ÷ total NPP, per the thesis methods), apply it everywhere, and decide once whether
  numbers quoted in the text are raw cell means or model-estimated marginal means. Mixing the two
  across an abstract and a figure is the kind of inconsistency reviewers flag as sloppiness.
- Confirm the full dataset is balanced at n = 5 across all 12 cells, with no pot losses.
- Confirm whether senescence timing differed by water treatment (see phenology caveat in §10).

---

## 12. Framing — lead with ecology, close with agronomy

Your data supports an ecological paper better than an agronomic one, and this is a strength, not a
retreat. An agronomy readership will object to pots, n = 5, no field validation, no LER, and a
non-standard planting density. An ecology readership takes greenhouse mesocosms at face value,
because the question is mechanism and context-dependency rather than yield prediction.

### The ecological arc

Three literatures converge on your design, and you can hit all three.

**1. The mutualism–parasitism continuum.** AMF outcomes range from strongly beneficial to net
costly depending on context (Johnson et al. 1997, *New Phytologist*; Hoeksema et al. 2010,
*Ecology Letters*). Almost all of that work varies *abiotic* context — nutrients, water, light.
Your contribution: **the biotic neighborhood also sets position on the continuum.**

**2. The stress gradient hypothesis.** SGH predicts that plant–plant interactions shift from
competition toward facilitation as abiotic stress increases (Bertness & Callaway 1994; Maestre et
al. 2009). Your result is a textbook SGH pattern — mono and mixed stands are indistinguishable
under control and moderate water, and diverge only under severe drought. **The twist is the
interesting part:** the facilitation appears *only in inoculated plants*. The mutualist itself
imposes the cost that the neighbor then relieves. That is a novel formulation — stress generated by
a belowground mutualism, buffered by an aboveground neighbor — and it is the sentence I would build
the abstract around.

**3. Mutualist-mediated competition.** Shared mycorrhizal partners restructure competitive
hierarchies between plants (Hart et al. 2003, *TREE*; Wagg et al. 2011; Bever et al. 2010,
*TREE*). Your fava analyses speak directly to this: the stronger the wheat–AMF association, the more
fava is suppressed — except under severe drought, where the hierarchy collapses. Exploratory, but
it makes the paper bidirectional (how the neighbor changes the symbiosis, *and* how the symbiosis
changes the neighbor) rather than one-directional.

### How the two framings divide the manuscript

- **Intro and discussion:** ecological. Context-dependency of mutualism, allocation theory, SGH.
- **A dedicated implications section:** agronomic. Inoculation decisions are made at planting,
  before drought is known; under monoculture this carries yield risk under severe drought, and
  legume intercropping may hedge it. Keep it to a few paragraphs and frame it as hypothesis for
  field testing, not recommendation.

This division lets you write the agronomic angle honestly without staking the paper's credibility
on pot-derived yield claims.

### Target journals

Realistic for a single greenhouse experiment at n = 5: *Plant and Soil*, *Mycorrhiza*,
*Applied Soil Ecology*, *Oecologia*, or *Functional Ecology* (the last two only if the allometric
and SGH framing carries). *Agriculture, Ecosystems & Environment* would want the density control
you cannot run. Field-oriented outlets such as *Field Crops Research* will want field data.

Given the framing above, I would target **Plant and Soil** or **Mycorrhiza** first — both routinely
publish greenhouse mycorrhizal context-dependency work at this scale, and both have readerships
comfortable with an ecology-primary, agronomy-secondary structure.

---

## 13. Running order — the handoff brief

Work through these in sequence. Steps 2 and 4 are gates: the answers determine whether the intended
paper exists in this form, so do not build figures or draft text before clearing them.

**Step 1 — Data audit.** Assemble the full 60-pot dataset. Confirm balance at n = 5 across all 12
cells. Construct: per-plant ANPP, BNPP, berry, NPP; the three allocation proportions; hyphal,
arbuscule, and vesicle percentages plus the log A:V ratio; fava biomass and fava % of pot biomass;
pot-level total productivity; wheat-equivalent density for inter pots (§1). Fix one definition of
berry allocation and one convention for quoted means. Verify no allocation values are exactly 0 or 1
before beta regression.

**Step 2 — GATE: allometry (§2).** Fit `log(berry) ~ log(veg) * Poly * Inoc * Water`. Test slope
homogeneity, then elevation shifts under a common slope. Cross-check with SMA.
*If the Inoc × Water effect vanishes after size adjustment, stop and reconsider the paper's claim
before proceeding.* Expected outcome is that it survives and strengthens, because inoculated plants
under severe drought were larger yet produced less grain — the opposite of the allometric
expectation.

**Step 3 — Primary model (§3).** Beta regression on berry allocation, `Poly * Inoc * Water`,
dispersion by Poly, water partitioned into the two planned contrasts. Extract the severe-drought
difference-of-differences with a 95% CI. Then the secondary responses, with `varIdent` where
Levene's indicates heterogeneity.

**Step 4 — GATE: fragility (§4).** Leave-one-out refit across all 60 pots, tracking both the
three-way p-value and the severe contrast estimate. Cook's distance and a raw dotplot of the five
pots in the Inoc/Mono/Severe cell. Record the LOO range for the supplement regardless of outcome,
and temper the abstract if a single pot moves the result materially.

**Step 5 — Mechanism (§10 Option A).** The arbuscule:vesicle models, inoculated pots only. Confirm
that hyphal intensity is flat across water while functional composition shifts. This is the section
that makes the paper mechanistic rather than descriptive.

**Step 6 — Exploratory.** Colonization associations within inoculated pots; fava suppression models;
the fava-biomass moderator probe (§6). Label all of it exploratory in the section headers.

**Step 7 — Optional bench work.** Grain counts and mean grain weight from stored berries
(§10 Option B), if steps 2–5 hold up and you want the yield-component decomposition.

**Step 8 — Figures and supplementary test table (§9, §8).** Build only after steps 2 and 4 clear.

Throughout: report effect sizes with confidence intervals, keep the confirmatory/exploratory split
explicit, and save every model object so the supplementary table is generated from code rather than
transcribed by hand.
