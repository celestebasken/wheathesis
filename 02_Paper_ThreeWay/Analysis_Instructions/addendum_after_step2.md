# Addendum to the analysis plan — after the Step 2 allometry gate

Insert as **Step 4b / 4c**, before Step 5. Three issues surfaced in the allometry figure that should
be handled before the 3.4× estimate goes into a manuscript.

---

## 4b-i. Non-overlapping size ranges in the key cell (important)

In the Mono × Severe panel, sterile pots occupy log(vegetative mass) ≈ 0.05–0.25 and inoculated pots
occupy ≈ 0.55–1.0. **There is no overlap.**

A size-adjusted contrast in that cell is therefore an **extrapolation**, not an interpolation. The
estimated vertical gap between the two groups depends entirely on the assumed common slope, because
neither group has any data in the other's size range. The point estimate is not precisely identified
by the data.

This does not undo the gate — the direction is stable across every plausible slope — but the paper
should report a sensitivity range rather than a single number:

```r
d_berry <- mean(log_berry[Inoc]) - mean(log_berry[Strl])   # ≈ -1.22
d_veg   <- mean(log_veg[Inoc])   - mean(log_veg[Strl])     # ≈ +0.60

for (b in seq(-1, 1, by = 0.25)) {
  cat(b, exp(-(d_berry - b * d_veg)), "\n")   # fold-difference at assumed slope b
}
```

The expected pattern: under a **positive** allometric slope — the one actually observed in the
Inter/Control and Inter/Moderate panels, and the biologically normal one — the effect gets **larger**,
not smaller. Under a zero slope it is the raw difference. Only a strongly negative slope shrinks it,
and even then it remains large.

**Manuscript language:** report the range, and state that the estimate is conservative under the
allometric slope observed in well-watered panels. Do not quote 3.4× as a precisely identified
quantity.

---

## 4b-ii. Distribution-free confirmation (cheap, high value)

In Mono × Severe the two groups appear **completely separated** on berry mass — every inoculated pot
below every sterile pot. Complete separation means the result does not depend on normality, variance
homogeneity, or the allometric model at all.

```r
coin::wilcox_test(berry ~ Inoc,
                  data = subset(wheat, Poly == "Mono" & Water == "Severe"),
                  distribution = "exact")
```

With n = 5 vs. 5 the smallest attainable two-sided p is ≈ 0.008. Reaching it is the point: an
assumption-free test agreeing with the parametric model is worth more here than a smaller p-value.
Report both.

Also run a permutation test of the three-way interaction itself (shuffle treatment labels within
water level, ~10,000 draws) as a supplementary robustness check.

---

## 4b-iii. Variance heterogeneity

The inoculated pots in Mono × Severe are far more scattered than the sterile ones — roughly a
sevenfold spread in berry mass versus a tight cluster. Use Welch, or `varIdent` / `dispformula`,
rather than a pooled-variance comparison for anything at the cell level. Check whether the scatter
tracks colonization intensity: if the least-colonized inoculated pots are the high performers, that
is a dose–response signal supporting the causal story, and it feeds directly into Step 5.

---

## 4c. New secondary result worth testing: drought decouples size from reproduction

The allometric slope is clearly positive in Inter/Control and Inter/Moderate, and flat-to-negative in
all three Severe panels. Test formally:

```r
m_slope <- lm(log(berry) ~ log(veg) * Water * Poly * Inoc, data = wheat)
# focus on the log(veg):Water terms
```

If it holds, this is a publishable secondary finding in its own right — **severe drought breaks the
normal scaling of reproductive output with plant size** — and it reinforces the main story by showing
the plants had stopped converting growth into grain at all. It also explains why the inoculated
plants could be simultaneously larger and less reproductive.

**Caution:** with n = 5 and compressed x-ranges, individual panel slopes are unreliable (the
Mono/Control negative slope looks driven by a single low point). Make the claim from the interaction
term, not from the six fitted lines.

---

## Standing item

Senescence timing per pot is absent from both data files. This does not block Steps 3–4, but it
matters for interpreting arbuscule counts in Step 5, since arbuscules decline naturally at
senescence. Check lab notebooks, harvest photographs, or greenhouse logs before Step 5. If
unrecoverable: the **mono-vs-inter comparison within severe drought remains sound**, because those
pots were harvested on the same date under the same water regime. The vulnerable comparison is
arbuscule abundance *across* water levels — treat that one as descriptive and state the caveat.
