# =============================================================================
# battery_core.R
# =============================================================================
# Shared setup + model-fitting logic used by BOTH:
#   - ANOVA_battery_DRAFT.Rmd        (prints the full stats output)
#   - Significance_Grid_DRAFT.Rmd    (builds the at-a-glance grid + mini-plots)
#
# Pulling this out into one file means both drafts are guaranteed to be
# testing the exact same models with the exact same settings - there's no
# risk of the two Rmds quietly drifting apart if one gets edited and the
# other doesn't. Each Rmd just does `source("battery_core.R")` as its first
# real step and gets: combined, inoc_only, inter_only, mono_only,
# battery_spec, results_list, assumption_summary, term_summary.
#
# Nothing in this file prints/knits anything - it only builds objects.
# =============================================================================

library(tidyverse)   # dplyr/tidyr/purrr for data wrangling + functional loops, ggplot2 for later plotting
library(broom)        # tidy() turns base-R model/ANOVA objects into data frames we can bind together
library(car)            # Anova() for Type II sums of squares, leveneTest() for variance homogeneity
library(emmeans)         # post-hoc pairwise comparisons (estimated marginal means)
library(effectsize)       # partial eta-squared effect sizes

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# Path is relative to whichever Rmd sources this file - both live in CleanCode/,
# so both reach Data/full_dataset.csv the same way.
combined <- read.csv("../Data/full_dataset.csv")

combined <- combined %>%
  mutate(
    Inoculation = factor(Inoculation),
    Polyculture = factor(Polyculture),
    Water       = factor(Water),   # alphabetical sort already gives
                                     # Control, Moderate, Severe - the
                                     # sensible low-to-high stress order.
    Group       = interaction(Inoculation, Polyculture, sep = "_")
  )

# No manual NA-row removal here (unlike the old combined[-4,] approach) -
# lm() below drops rows missing whatever variable IT needs, per model, via
# its default na.action. That means the one pot with a missing BNPP value
# only gets dropped from BNPP-involving models, not from every model.

# -----------------------------------------------------------------------------
# 2. Logit-transform every proportion column, up front
# -----------------------------------------------------------------------------
logit_transform_prop <- function(p) {
  # qlogis() = log(p / (1 - p)); undefined at exactly 0/1, so nudge boundary
  # values slightly inward first. Matters for hyphal/arbuscular colonization
  # and Pot_pcent_fava, which legitimately hit exactly 0 for Strl/Mono pots.
  p_clamped <- pmin(pmax(p, 0.001), 0.999)
  qlogis(p_clamped)
}

combined <- combined %>%
  mutate(
    W_pcent_ANPP_logit          = logit_transform_prop(W_pcent_ANPP),
    W_pcent_BNPP_logit          = logit_transform_prop(W_pcent_BNPP),
    W_pcent_berries_logit       = logit_transform_prop(W_pcent_berries),
    hyphal_colonization_logit   = logit_transform_prop(hyphal_colonization),
    arb_vesc_colonization_logit = logit_transform_prop(arb_vesc_colonization),
    Pot_pcent_fava_logit        = logit_transform_prop(Pot_pcent_fava)
  )

# -----------------------------------------------------------------------------
# 3. Subsets
# -----------------------------------------------------------------------------
# AMF colonization is a structural zero for Strl pots, and fava biomass is a
# structural zero for Mono pots - so those outcomes are modeled on subsets
# that exclude the treatment level that determines the zero, rather than
# including a predictor that's really just re-detecting the subsetting rule.
inoc_only  <- combined %>% filter(Inoculation == "Inoc")
inter_only <- combined %>% filter(Polyculture == "Inter")
mono_only  <- combined %>% filter(Polyculture == "Mono")   # kept in case needed later

# -----------------------------------------------------------------------------
# 4. The reusable ANOVA-battery function
# -----------------------------------------------------------------------------
run_anova_check <- function(data,             # data frame to fit on (already subsetted if needed)
                             response,         # string: name of the response column (raw or _logit)
                             predictors,       # character vector of predictor column names
                             label,            # human-readable name for this variable
                             family,           # correction family: "biomass"/"allocation"/"amf"/"fava"
                             posthoc_rhs = NULL  # optional emmeans right-hand-side spec
                             ) {

  # ---- Step 1: build the model formula dynamically ----
  rhs  <- paste(predictors, collapse = " * ")
  form <- stats::as.formula(paste(response, "~", rhs))

  # ---- Step 2: fit with lm(), not aov() ----
  # car::Anova() (Step 3) needs an lm object; lm()/aov() give identical
  # fitted values/residuals for a factorial design like this one anyway.
  mod <- stats::lm(form, data = data)

  # ---- Step 3: Type II ANOVA table ----
  # Type II doesn't depend on term order, unlike base aov()'s Type I - matters
  # here because one response variable has a missing value (breaks perfect
  # balance for just that model) and because we want one consistent rule
  # applied to every variable rather than deciding per-chunk.
  anova_raw <- car::Anova(mod, type = 2)

  anova_tbl <- anova_raw %>%
    broom::tidy() %>%
    dplyr::filter(term != "Residuals")

  # ---- Step 4a: normality check ----
  shapiro_p <- stats::shapiro.test(stats::residuals(mod))$p.value

  # ---- Step 4b: variance-homogeneity check ----
  # Built dynamically from ALL predictors passed in, not hardcoded to two
  # of three factors.
  group_formula <- stats::as.formula(
    paste(response, "~", paste0("interaction(", paste(predictors, collapse = ", "), ")"))
  )
  levene_p <- car::leveneTest(group_formula, data = data)[1, "Pr(>F)"]

  # ---- Step 5: effect size (partial eta-squared) ----
  # Fed the car::Anova() object directly (not `type = 2` as an argument to
  # eta_squared() itself - that collides with eta_squared()'s own internal
  # `type` argument, which picks eta/omega/epsilon-squared, not SS type).
  eta_tbl <- effectsize::eta_squared(anova_raw, partial = TRUE) %>%
    as.data.frame()

  # ---- Step 6: optional post-hoc pairwise comparisons ----
  posthoc <- NULL
  if (!is.null(posthoc_rhs)) {
    posthoc <- emmeans::emmeans(
      mod,
      specs  = stats::as.formula(paste("pairwise ~", posthoc_rhs)),
      adjust = "tukey"
    )
  }

  # ---- Step 7: return everything ----
  list(
    label     = label,
    family    = family,
    model     = mod,
    anova     = anova_tbl,
    shapiro_p = shapiro_p,
    levene_p  = levene_p,
    eta2      = eta_tbl,
    posthoc   = posthoc
  )
}

# -----------------------------------------------------------------------------
# 5. The battery spec: one row per response variable
# -----------------------------------------------------------------------------
# Single source of truth for "which variable, on which data, with which
# predictors." `display_group` is NEW here (not a stats decision) - it's
# purely for Significance_Grid_DRAFT.Rmd to know which of the two output
# grids ("mass" columns in raw grams, vs "percent" columns) a variable
# belongs in, mirroring how your old PDF table split across two grids.
# It has no effect on the ANOVA/FDR logic, which still keys off `family`.
battery_spec <- list(

  # ---------------- Biomass family: raw grams, per plant, full 3-way --------
  list(label = "Wheat total NPP per plant",       response = "W_NPP_pp",       data = combined,
       predictors = c("Inoculation", "Polyculture", "Water"), family = "biomass",
       posthoc_rhs = "Water | Inoculation * Polyculture", display_group = "mass"),

  list(label = "Wheat aboveground NPP per plant", response = "W_ANPP_g_pp",   data = combined,
       predictors = c("Inoculation", "Polyculture", "Water"), family = "biomass",
       posthoc_rhs = "Water | Inoculation * Polyculture", display_group = "mass"),

  list(label = "Wheat belowground NPP per plant", response = "W_BNPP_g_pp",  data = combined,
       predictors = c("Inoculation", "Polyculture", "Water"), family = "biomass",
       posthoc_rhs = "Water | Inoculation * Polyculture", display_group = "mass"),   # n=59, see Section 1

  list(label = "Wheat berry biomass per plant",   response = "W_Berries_g_pp", data = combined,
       predictors = c("Inoculation", "Polyculture", "Water"), family = "biomass",
       posthoc_rhs = "Water | Inoculation * Polyculture", display_group = "mass"),

  list(label = "Total pot productivity",          response = "Pot_NPP_g_pp",  data = combined,
       predictors = c("Inoculation", "Polyculture", "Water"), family = "biomass",
       posthoc_rhs = "Water | Inoculation * Polyculture", display_group = "mass"),
  # ^ ADDED vs. the first draft - this is the "Total pot productivity" column
  # from your old PDF table, which wasn't in the original battery_spec.

  # ---------------- Allocation family: logit-transformed %, full 3-way ------
  list(label = "Aboveground allocation (%)", response = "W_pcent_ANPP_logit",    data = combined,
       predictors = c("Inoculation", "Polyculture", "Water"), family = "allocation",
       posthoc_rhs = "Water | Inoculation * Polyculture", display_group = "percent"),

  list(label = "Belowground allocation (%)", response = "W_pcent_BNPP_logit",   data = combined,
       predictors = c("Inoculation", "Polyculture", "Water"), family = "allocation",
       posthoc_rhs = "Water | Inoculation * Polyculture", display_group = "percent"),

  list(label = "Berry allocation (%)",       response = "W_pcent_berries_logit", data = combined,
       predictors = c("Inoculation", "Polyculture", "Water"), family = "allocation",
       posthoc_rhs = "Water | Inoculation * Polyculture", display_group = "percent"),
  # W_pcent_vegetative deliberately left out - it's just 1 - W_pcent_berries.

  # ---------------- AMF family: Inoc-only subset, 2-way ----------------------
  list(label = "Hyphal colonization (%)",        response = "hyphal_colonization_logit",   data = inoc_only,
       predictors = c("Polyculture", "Water"), family = "amf",
       posthoc_rhs = "Water | Polyculture", display_group = "percent"),

  list(label = "Arbuscular/vesicular colonization (%)", response = "arb_vesc_colonization_logit", data = inoc_only,
       predictors = c("Polyculture", "Water"), family = "amf",
       posthoc_rhs = "Water | Polyculture", display_group = "percent"),

  # ---------------- Fava family: Inter-only subset, 2-way --------------------
  list(label = "Fava biomass",          response = "Fava_g",             data = inter_only,
       predictors = c("Inoculation", "Water"), family = "fava",
       posthoc_rhs = "Water | Inoculation", display_group = "mass"),

  list(label = "Fava % of pot biomass", response = "Pot_pcent_fava_logit", data = inter_only,
       predictors = c("Inoculation", "Water"), family = "fava",
       posthoc_rhs = "Water | Inoculation", display_group = "percent")

  # To add another variable: copy one list() entry, change the fields, and
  # add it to this list (comma-separated).
)

# A named-by-label version, so downstream code can look up "what data/
# response goes with this label" without re-scanning the whole list.
spec_by_label <- setNames(battery_spec, purrr::map_chr(battery_spec, "label"))

# -----------------------------------------------------------------------------
# 6. Run the battery
# -----------------------------------------------------------------------------
results_list <- list()

for (spec in battery_spec) {
  results_list[[spec$label]] <- run_anova_check(
    data        = spec$data,
    response    = spec$response,
    predictors  = spec$predictors,
    label       = spec$label,
    family      = spec$family,
    posthoc_rhs = spec$posthoc_rhs
  )
}

# -----------------------------------------------------------------------------
# 7. Assumption-check summary table
# -----------------------------------------------------------------------------
assumption_summary <- purrr::map_dfr(results_list, function(r) {
  tibble::tibble(
    family         = r$family,
    label          = r$label,
    n              = stats::nobs(r$model),
    shapiro_p      = round(r$shapiro_p, 4),
    levene_p       = round(r$levene_p, 4),
    normality_flag = ifelse(r$shapiro_p < 0.05, "check", "ok"),
    variance_flag  = ifelse(r$levene_p  < 0.05, "check", "ok")
  )
})

# -----------------------------------------------------------------------------
# 8. Term-level ANOVA results, with effect sizes and FDR correction
# -----------------------------------------------------------------------------
term_summary <- purrr::map_dfr(results_list, function(r) {
  r$anova %>%
    dplyr::mutate(label = r$label, family = r$family) %>%
    dplyr::select(family, label, term, df, statistic, p.value)
})

eta_summary <- purrr::map_dfr(results_list, function(r) {
  r$eta2 %>%
    dplyr::rename(term = Parameter, eta2_partial = Eta2_partial) %>%
    dplyr::mutate(label = r$label) %>%
    dplyr::select(label, term, eta2_partial)
})

term_summary <- term_summary %>%
  dplyr::left_join(eta_summary, by = c("label", "term"))

# Multiple-comparison correction done WITHIN family (biomass / allocation /
# amf / fava), not across the whole table at once - see the design writeup
# in ANOVA_battery_DRAFT.Rmd Section 8 for the full reasoning.
term_summary <- term_summary %>%
  dplyr::group_by(family) %>%
  dplyr::mutate(p_adj = stats::p.adjust(p.value, method = "BH")) %>%
  dplyr::ungroup()
