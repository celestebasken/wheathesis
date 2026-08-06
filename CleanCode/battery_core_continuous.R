# =============================================================================
# battery_core_continuous.R
# =============================================================================
# Sibling to battery_core.R, for the "use hyphal_colonization as a CONTINUOUS
# predictor" version of the significance grids, instead of the binary
# Inoculation (Inoc/Strl) factor. Sourced by:
#   - Significance_Grid_Hyphal_Full_DRAFT.Rmd
#   - Significance_Grid_Hyphal_Inter_DRAFT.Rmd
#   - Significance_Grid_Hyphal_Mono_DRAFT.Rmd
#
# This file assumes battery_core.R has ALREADY been sourced (for tidyverse,
# combined/inoc_only/inter_only/mono_only, and run_battery()) - it only adds
# what's actually DIFFERENT about a continuous predictor:
#
#   1. build_battery_spec_continuous(scope) - same response variables as
#      build_battery_spec(), but every model swaps the binary `Inoculation`
#      factor for continuous `hyphal_colonization`, and is fit ONLY on
#      Inoc-only pots (hyphal_colonization is a structural 0 for Strl pots,
#      not a real continuous observation there - same reasoning
#      build_battery_spec() already applies to the AMF-colonization rows).
#      There's no more "AMF family" here (hyphal_colonization can't predict
#      itself) - just biomass/allocation/fava, same as before.
#
#   2. run_ancova_check() - same overall shape/return value as
#      run_anova_check(), so it can be dropped straight into run_battery()
#      via `model_fn = run_ancova_check`, but with two real statistical
#      differences a continuous predictor requires:
#        - Levene's test can only check variance homogeneity across the
#          CATEGORICAL predictors (Water, and Polyculture for scope="full") -
#          it can't group on a continuous variable. This is flagged, not
#          silently glossed over.
#        - Post-hoc is emtrends() (comparing SLOPES of hyphal_colonization,
#          optionally within levels of a categorical factor), not emmeans()
#          pairwise comparison of group means - there's no "group mean" to
#          compare when the predictor of interest is continuous.
#
# A REAL POWER CAVEAT, not just a formality: hyphal_colonization only varies
# among Inoc pots, so every model here starts from half the data (n=30) that
# build_battery_spec() gets for the same scope, before even splitting by
# Water/Polyculture. For scope="inter"/"mono" that's n=15 total, meaning ~5
# points per Water level for a slope estimate; for scope="full" the 3-way
# interaction term's per-cell slopes are drawn from ~5 points each. Treat any
# interaction-term slope facet here as exploratory, not conclusive - the
# assumption_summary's `n` column (from run_battery(), unchanged) makes this
# visible rather than hidden.
#
# Nothing in this file prints/knits anything - it only builds objects.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. The battery spec: one row per response variable, per scope
# -----------------------------------------------------------------------------
# Mirrors build_battery_spec() closely - same response variables, same
# display_group/plot_response conventions - the only structural differences
# are: (a) every model is fit on an Inoc-only subset, (b) `Inoculation` is
# replaced by continuous `hyphal_colonization` in every predictor list, and
# (c) there's no AMF family (nothing left for hyphal_colonization to predict
# there - it WOULD be predicting itself).
build_battery_spec_continuous <- function(scope = c("full", "inter", "mono")) {
  scope <- match.arg(scope)

  if (scope == "full") {
    base_data       <- inoc_only
    base_predictors <- c("hyphal_colonization", "Polyculture", "Water")
    base_posthoc    <- "Water | Polyculture"

    include_fava      <- TRUE
    fava_data         <- inter_only %>% dplyr::filter(Inoculation == "Inoc")
    fava_predictors   <- c("hyphal_colonization", "Water")
    fava_posthoc       <- "Water"

    include_total_pot <- TRUE

  } else if (scope == "inter") {
    base_data       <- inter_only %>% dplyr::filter(Inoculation == "Inoc")
    base_predictors <- c("hyphal_colonization", "Water")
    base_posthoc    <- "Water"

    include_fava      <- TRUE
    fava_data         <- base_data
    fava_predictors   <- c("hyphal_colonization", "Water")
    fava_posthoc      <- "Water"

    include_total_pot <- TRUE

  } else { # scope == "mono"
    base_data       <- mono_only %>% dplyr::filter(Inoculation == "Inoc")
    base_predictors <- c("hyphal_colonization", "Water")
    base_posthoc    <- "Water"

    include_fava      <- FALSE   # same reasoning as build_battery_spec(): Fava_g is
                                   # exactly 0 for every Mono pot, regardless of
                                   # colonization level - zero variance to test
    include_total_pot <- FALSE   # same reasoning too: identical to wheat-only NPP for Mono
  }

  spec <- list(

    # ---------------- Biomass family: raw grams, per plant ------------------
    list(label = "Wheat total NPP per plant",       response = "W_NPP_pp",       data = base_data,
         predictors = base_predictors, family = "biomass", posthoc_rhs = base_posthoc, display_group = "mass"),

    list(label = "Wheat aboveground NPP per plant", response = "W_ANPP_g_pp",   data = base_data,
         predictors = base_predictors, family = "biomass", posthoc_rhs = base_posthoc, display_group = "mass"),

    list(label = "Wheat belowground NPP per plant", response = "W_BNPP_g_pp",  data = base_data,
         predictors = base_predictors, family = "biomass", posthoc_rhs = base_posthoc, display_group = "mass"),

    list(label = "Wheat berry biomass per plant",   response = "W_Berries_g_pp", data = base_data,
         predictors = base_predictors, family = "biomass", posthoc_rhs = base_posthoc, display_group = "mass"),

    # ---------------- Allocation family: logit-transformed % ----------------
    list(label = "Aboveground allocation (%)", response = "W_pcent_ANPP_logit",    data = base_data,
         predictors = base_predictors, family = "allocation", posthoc_rhs = base_posthoc, display_group = "percent",
         plot_response = "W_pcent_ANPP"),

    list(label = "Belowground allocation (%)", response = "W_pcent_BNPP_logit",   data = base_data,
         predictors = base_predictors, family = "allocation", posthoc_rhs = base_posthoc, display_group = "percent",
         plot_response = "W_pcent_BNPP"),

    list(label = "Berry allocation (%)",       response = "W_pcent_berries_logit", data = base_data,
         predictors = base_predictors, family = "allocation", posthoc_rhs = base_posthoc, display_group = "percent",
         plot_response = "W_pcent_berries")
  )

  if (include_total_pot) {
    spec <- c(spec, list(
      list(label = "Total pot productivity", response = "Pot_NPP_g_pp", data = base_data,
           predictors = base_predictors, family = "biomass", posthoc_rhs = base_posthoc, display_group = "mass")
    ))
  }

  if (include_fava) {
    spec <- c(spec, list(
      list(label = "Fava biomass",          response = "Fava_g",             data = fava_data,
           predictors = fava_predictors, family = "fava", posthoc_rhs = fava_posthoc, display_group = "mass"),

      list(label = "Fava % of pot biomass", response = "Pot_pcent_fava_logit", data = fava_data,
           predictors = fava_predictors, family = "fava", posthoc_rhs = fava_posthoc, display_group = "percent",
           plot_response = "Pot_pcent_fava")
    ))
  }

  spec
}

# -----------------------------------------------------------------------------
# 2. The continuous-predictor model-fitting function
# -----------------------------------------------------------------------------
# Same argument list and same return shape as run_anova_check() - that's
# what lets run_battery() (in battery_core.R) accept this as its `model_fn`
# unchanged. Two steps differ for real statistical reasons (Levene's test,
# post-hoc) - see the file header above.
run_ancova_check <- function(data,             # data frame (already Inoc-only-subsetted)
                               response,         # string: name of the response column (raw or _logit)
                               predictors,       # e.g. c("hyphal_colonization","Water") or with "Polyculture" too
                               label,            # human-readable name for this variable
                               family,           # correction family: "biomass"/"allocation"/"fava"
                               posthoc_rhs = NULL  # emtrends right-hand-side spec, e.g. "Water | Polyculture"
                               ) {

  # ---- Step 1: build the model formula dynamically ----
  rhs  <- paste(predictors, collapse = " * ")
  form <- stats::as.formula(paste(response, "~", rhs))

  # ---- Step 2: fit with lm() ----
  # This is a standard ANCOVA setup (one continuous predictor crossed with
  # one or two categorical factors) - lm() handles it exactly like any other
  # linear model; nothing special is needed to tell it hyphal_colonization
  # is continuous rather than a factor (it just isn't wrapped in factor()).
  mod <- stats::lm(form, data = data)

  # ---- Step 3: Type II ANOVA table ----
  # Type II SS is valid for ANCOVA models too (continuous + categorical
  # terms) - same reasoning as run_anova_check(): doesn't depend on term
  # order, applied consistently across every variable.
  anova_raw <- car::Anova(mod, type = 2)

  anova_tbl <- anova_raw %>%
    broom::tidy() %>%
    dplyr::filter(term != "Residuals")

  # ---- Step 4a: normality check ----
  # Unchanged from run_anova_check() - Shapiro-Wilk on residuals doesn't
  # care whether the predictors were continuous or categorical.
  shapiro_p <- stats::shapiro.test(stats::residuals(mod))$p.value

  # ---- Step 4b: variance-homogeneity check - ADAPTED for a continuous predictor ----
  # car::leveneTest needs a CATEGORICAL grouping variable to compare
  # variances across - it cannot group observations by a continuous
  # variable like hyphal_colonization. So this only checks homogeneity of
  # residual variance across the categorical predictors (Water, and
  # Polyculture when present) - it does NOT check whether residual variance
  # changes across the range of hyphal_colonization itself (a "fan-shaped"
  # pattern in a residuals-vs-fitted plot would be the sign of that, which
  # this numeric test can't catch). Flagged here explicitly rather than
  # silently treating this Levene result as if it covered the same ground
  # it does in run_anova_check().
  categorical_predictors <- setdiff(predictors, "hyphal_colonization")
  if (length(categorical_predictors) > 0) {
    group_formula <- stats::as.formula(
      paste(response, "~", paste0("interaction(", paste(categorical_predictors, collapse = ", "), ")"))
    )
    levene_p <- car::leveneTest(group_formula, data = data)[1, "Pr(>F)"]
  } else {
    # Every spec in build_battery_spec_continuous() keeps at least Water as
    # a categorical predictor, so this branch shouldn't trigger in practice -
    # kept as a safe fallback rather than letting leveneTest() error out.
    levene_p <- NA_real_
  }

  # ---- Step 5: effect size (partial eta-squared) ----
  # Unchanged from run_anova_check() - works the same way for a continuous
  # term as a categorical one (it's still "share of variance attributable
  # to this term's F-test"), including the one-predictor Eta2/Eta2_partial
  # naming quirk (see the comment in run_anova_check() for why).
  eta_tbl <- effectsize::eta_squared(anova_raw, partial = TRUE) %>%
    as.data.frame()
  if ("Eta2" %in% names(eta_tbl) && !("Eta2_partial" %in% names(eta_tbl))) {
    eta_tbl <- eta_tbl %>% dplyr::rename(Eta2_partial = Eta2)
  }

  # ---- Step 6: post-hoc - SLOPE trends, not pairwise means ----
  # There's no "pairwise comparison of group means" for a continuous
  # predictor - the ANCOVA equivalent is emtrends(): the estimated SLOPE of
  # hyphal_colonization on the response, either overall or within each
  # level of a categorical grouping factor (posthoc_rhs), plus a test of
  # whether those slopes differ from each other. Same `pairwise ~ ...`
  # spec syntax as emmeans() in run_anova_check(), just handed to
  # emtrends() with `var` naming which predictor's slope to estimate.
  posthoc <- NULL
  if (!is.null(posthoc_rhs)) {
    posthoc <- emmeans::emtrends(
      mod,
      specs = stats::as.formula(paste("pairwise ~", posthoc_rhs)),
      var   = "hyphal_colonization"
    )
  }

  # ---- Step 7: return everything, same shape as run_anova_check() ----
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
