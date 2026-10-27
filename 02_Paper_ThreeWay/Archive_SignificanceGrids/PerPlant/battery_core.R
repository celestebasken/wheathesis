# =============================================================================
# battery_core.R  (PerPlant version - the "per plant", _pp-variable core)
# =============================================================================
# Shared setup + model-fitting logic used by every doc in CleanCode/PerPlant/:
#   - ANOVA_battery_DRAFT.Rmd           (prints the full stats output, full dataset)
#   - Significance_Grid_DRAFT.Rmd       (at-a-glance grid, full dataset)
#   - Significance_Grid_Inter_DRAFT.Rmd (at-a-glance grid, Inter pots only)
#   - Significance_Grid_Mono_DRAFT.Rmd  (at-a-glance grid, Mono pots only)
# ...plus, one directory down, every doc in CleanCode/PerPlant/AMF/ (the
# hyphal-colonization-as-predictor grids + the intercropping-focused AMF doc),
# via `source("../battery_core.R")` - see this file's find_data_file() helper
# below for why sourcing it from a deeper folder still works.
#
# Pulling this out into one file means every doc above is guaranteed to be
# testing variables the same way with the same settings - there's no risk of
# them quietly drifting apart if one gets edited and the others don't.
#
# This is the "per plant" fork specifically - mass-family variables here are
# density-adjusted per-plant averages (W_NPP_pp, W_ANPP_g_pp, etc.). The
# sibling core at CleanCode/PerPot/battery_core.R uses the pot-total versions
# instead (W_NPP, W_ANPP_g, etc.) - see that file's header comment.
#
# Sourcing this file always gives you the FULL-dataset battery already run,
# under the usual top-level names (battery_spec, spec_by_label, results_list,
# assumption_summary, term_summary) - that's what the first two docs above
# use directly, unchanged from before.
#
# The Inter-only and Mono-only docs want the SAME variables, but fit as 2-way
# Inoculation*Water models (Polyculture is constant once you've subsetted to
# one polyculture, so it drops out as a predictor) instead of the full 3-way.
# Rather than copy/paste the spec + model-running code two more times with
# that one difference, both are exposed as functions - `build_battery_spec()`
# and `run_battery()` - which the top-level code below also uses to build the
# "full" version. The Inter/Mono docs just call them again with a different
# scope. See the design note above `build_battery_spec()`.
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
# Walks UP from the current working directory (knitr's default working dir is
# wherever the knitting .Rmd itself lives) until it finds Data/full_dataset.csv
# - rather than a hardcoded "../Data/..." that silently breaks the moment a
# doc that sources this file moves to a different folder depth (as happened
# when PerPlant/AMF/ was added one level below where this file's docs
# normally live). Depth-independent by construction, so this line never needs
# to change again just because a doc gets reorganized into a subfolder.
find_data_file <- function(rel_path = file.path("Data", "full_dataset.csv"),
                            start = getwd(), max_up = 6) {
  dir <- normalizePath(start, mustWork = FALSE)
  for (i in seq_len(max_up + 1)) {
    candidate <- file.path(dir, rel_path)
    if (file.exists(candidate)) return(candidate)
    parent <- dirname(dir)
    if (identical(parent, dir)) break
    dir <- parent
  }
  stop("Could not find ", rel_path, " by walking up from ", start,
       " - are you knitting from somewhere inside the wheathesis project?")
}

combined <- read.csv(find_data_file())

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

  # For a one-predictor (one-way) model - which only ever happens for the
  # AMF-colonization variables once Polyculture is fixed by an Inter-only or
  # Mono-only scope, leaving just `Water` to test - partial eta-squared is
  # mathematically identical to plain eta-squared, and effectsize() labels
  # the column "Eta2" instead of "Eta2_partial" to reflect that. Normalize
  # both cases to one column name here, so run_battery() (which joins on
  # "Eta2_partial") doesn't need to know which flavor of model produced this
  # particular table.
  if ("Eta2" %in% names(eta_tbl) && !("Eta2_partial" %in% names(eta_tbl))) {
    eta_tbl <- eta_tbl %>% dplyr::rename(Eta2_partial = Eta2)
  }

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
# 5. The battery spec: one row per response variable, per scope
# -----------------------------------------------------------------------------
# `display_group` and `plot_response` are cosmetic (not a stats decision):
# `display_group` tells the grid docs which of the two output grids ("mass"
# columns in raw grams, vs "percent" columns) a variable belongs in.
# `plot_response` tells the thumbnail plots which column to actually draw,
# separate from `response` (what the model is fit on) - every _logit
# variable needs the logit scale to be a valid model, but a bar chart of
# logit values goes negative whenever the underlying proportion is below
# 0.5, which put "0" at the TOP of the panel for some variables and the
# BOTTOM for others, purely depending on which side of 0.5 that variable's
# mean happened to sit. Plotting the original 0-1 proportion instead keeps
# every percent-family thumbnail non-negative, so 0 is always at the bottom.
# Variables without a `plot_response` just fall back to `response` (see the
# grid docs) - that covers the mass-family variables, never transformed and
# already non-negative grams.
#
# `scope` controls which pots are being modeled, and it changes more than
# just "which data frame": once you subset to only Inter or only Mono pots,
# Polyculture is constant - it can't be a predictor anymore (zero variation
# left for it to explain) - so those models are a 2-way Inoculation*Water,
# not the full dataset's 3-way Inoculation*Polyculture*Water. AMF variables
# similarly drop from `Polyculture * Water` to just `Water` once Polyculture
# is fixed by the scope, on top of already dropping Inoculation as a
# predictor (it's fixed at "Inoc" by the Inoc-only subsetting, same as the
# full-dataset version). Fava_g / Pot_pcent_fava are dropped ENTIRELY for
# scope = "mono" - they're exactly 0 for every Mono pot by definition (no
# fava plant present), so there's zero variance for an ANOVA to test; fitting
# one would only ever be able to report "no effect" on a constant.
build_battery_spec <- function(scope = c("full", "inter", "mono")) {
  scope <- match.arg(scope)

  if (scope == "full") {
    base_data       <- combined
    base_predictors <- c("Inoculation", "Polyculture", "Water")
    base_posthoc    <- "Water | Inoculation * Polyculture"

    amf_data       <- inoc_only
    amf_predictors <- c("Polyculture", "Water")
    amf_posthoc    <- "Water | Polyculture"

    include_fava      <- TRUE
    fava_data         <- inter_only
    fava_predictors   <- c("Inoculation", "Water")
    fava_posthoc      <- "Water | Inoculation"

    include_total_pot <- TRUE   # meaningfully different from wheat-only NPP (includes fava biomass)

  } else if (scope == "inter") {
    base_data       <- inter_only
    base_predictors <- c("Inoculation", "Water")
    base_posthoc    <- "Water | Inoculation"

    amf_data       <- inter_only %>% dplyr::filter(Inoculation == "Inoc")
    amf_predictors <- c("Water")
    amf_posthoc    <- "Water"

    include_fava      <- TRUE   # every Inter pot has a fava plant - meaningful here
    fava_data         <- inter_only
    fava_predictors   <- c("Inoculation", "Water")
    fava_posthoc      <- "Water | Inoculation"

    include_total_pot <- TRUE   # still meaningfully different here too (wheat + fava)

  } else { # scope == "mono"
    base_data       <- mono_only
    base_predictors <- c("Inoculation", "Water")
    base_posthoc    <- "Water | Inoculation"

    amf_data       <- mono_only %>% dplyr::filter(Inoculation == "Inoc")
    amf_predictors <- c("Water")
    amf_posthoc    <- "Water"

    include_fava      <- FALSE   # see the note above the function
    include_total_pot <- FALSE   # identical to "Wheat total NPP per plant" for Mono pots
                                  # (no fava biomass to add) - per your call, dropped rather
                                  # than kept as a redundant column.
  }

  spec <- list(

    # ---------------- Biomass family: raw grams, per plant ------------------
    list(label = "Wheat total NPP per plant",       response = "W_NPP_pp",       data = base_data,
         predictors = base_predictors, family = "biomass", posthoc_rhs = base_posthoc, display_group = "mass"),

    list(label = "Wheat aboveground NPP per plant", response = "W_ANPP_g_pp",   data = base_data,
         predictors = base_predictors, family = "biomass", posthoc_rhs = base_posthoc, display_group = "mass"),

    list(label = "Wheat belowground NPP per plant", response = "W_BNPP_g_pp",  data = base_data,
         predictors = base_predictors, family = "biomass", posthoc_rhs = base_posthoc, display_group = "mass"),
    # ^ n may be one less than the rest of this scope's variables - one pot
    # is missing a BNPP measurement (see Section 1) - lm()'s default
    # na.action drops it only from models that actually use this column.

    list(label = "Wheat berry biomass per plant",   response = "W_Berries_g_pp", data = base_data,
         predictors = base_predictors, family = "biomass", posthoc_rhs = base_posthoc, display_group = "mass"),

    # ---------------- Allocation family: logit-transformed % ----------------
    # (plot_response = the un-transformed 0-1 column, so the thumbnail shows
    # the natural proportion instead of the logit scale used for the model.)
    list(label = "Aboveground allocation (%)", response = "W_pcent_ANPP_logit",    data = base_data,
         predictors = base_predictors, family = "allocation", posthoc_rhs = base_posthoc, display_group = "percent",
         plot_response = "W_pcent_ANPP"),

    list(label = "Belowground allocation (%)", response = "W_pcent_BNPP_logit",   data = base_data,
         predictors = base_predictors, family = "allocation", posthoc_rhs = base_posthoc, display_group = "percent",
         plot_response = "W_pcent_BNPP"),

    list(label = "Berry allocation (%)",       response = "W_pcent_berries_logit", data = base_data,
         predictors = base_predictors, family = "allocation", posthoc_rhs = base_posthoc, display_group = "percent",
         plot_response = "W_pcent_berries"),
    # W_pcent_vegetative deliberately left out - it's just 1 - W_pcent_berries.

    # ---------------- AMF family: Inoc-only subset ---------------------------
    list(label = "Hyphal colonization (%)",        response = "hyphal_colonization_logit",   data = amf_data,
         predictors = amf_predictors, family = "amf", posthoc_rhs = amf_posthoc, display_group = "percent",
         plot_response = "hyphal_colonization"),

    list(label = "Arbuscular/vesicular colonization (%)", response = "arb_vesc_colonization_logit", data = amf_data,
         predictors = amf_predictors, family = "amf", posthoc_rhs = amf_posthoc, display_group = "percent",
         plot_response = "arb_vesc_colonization")
  )

  # ---------------- Total pot productivity ------------------------------------
  # Excluded entirely for scope = "mono" - per your call, since it's
  # numerically identical to "Wheat total NPP per plant" there (max
  # difference ~5e-5 g, i.e. rounding): a Mono pot's total biomass IS its
  # wheat biomass, with no fava plant to add. Kept for "full" and "inter",
  # where it's a genuinely different quantity (wheat + fava).
  if (include_total_pot) {
    spec <- c(spec, list(
      list(label = "Total pot productivity", response = "Pot_NPP_g_pp", data = base_data,
           predictors = base_predictors, family = "biomass", posthoc_rhs = base_posthoc, display_group = "mass")
    ))
  }

  # ---------------- Fava family: Inter-only subset ---------------------------
  # Included for scope "full" and "inter" (fava biomass varies meaningfully
  # there); excluded entirely for scope "mono" (see the note above the
  # function for why).
  if (include_fava) {
    spec <- c(spec, list(
      list(label = "Fava biomass",          response = "Fava_g",             data = fava_data,
           predictors = fava_predictors, family = "fava", posthoc_rhs = fava_posthoc, display_group = "mass"),

      list(label = "Fava % of pot biomass", response = "Pot_pcent_fava_logit", data = fava_data,
           predictors = fava_predictors, family = "fava", posthoc_rhs = fava_posthoc, display_group = "percent",
           plot_response = "Pot_pcent_fava")
    ))
  }

  # To add another variable: copy one list() entry above, change the fields,
  # and add it to `spec` (or one of the conditional blocks above, if it
  # should be excluded from a particular scope for a structural reason).
  spec
}

# -----------------------------------------------------------------------------
# 6-8. Run a battery: fit every model, build the assumption-check summary and
#      the term-level (FDR-corrected) summary
# -----------------------------------------------------------------------------
# Takes a battery_spec (from build_battery_spec()) and returns everything
# downstream code needs, bundled in one list - this is what used to be three
# separate top-level sections (6/7/8) run once against a single hardcoded
# battery_spec. Wrapping it in a function is what lets the Inter/Mono grid
# docs re-run the exact same pipeline against their own scope's spec.
#
# `model_fn` defaults to run_anova_check() (categorical Inoculation/Water/
# Polyculture factors) - everything below it (assumption summary, FDR
# correction) only cares that whatever function it's given returns the same
# shape (list(label, family, model, anova, shapiro_p, levene_p, eta2,
# posthoc)), not what kind of model produced it. This is what lets
# battery_core_continuous.R's run_ancova_check() (continuous hyphal
# colonization as the predictor of interest, instead of the binary
# Inoculation factor) reuse this exact function unchanged, just passing
# `model_fn = run_ancova_check` - rather than a fourth near-duplicate
# copy of run/summarize logic.
run_battery <- function(battery_spec, model_fn = run_anova_check) {

  # A named-by-label version, so downstream code can look up "what data/
  # response goes with this label" without re-scanning the whole list.
  spec_by_label <- setNames(battery_spec, purrr::map_chr(battery_spec, "label"))

  results_list <- list()
  for (spec in battery_spec) {
    results_list[[spec$label]] <- model_fn(
      data        = spec$data,
      response    = spec$response,
      predictors  = spec$predictors,
      label       = spec$label,
      family      = spec$family,
      posthoc_rhs = spec$posthoc_rhs
    )
  }

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
  # in ANOVA_battery_DRAFT.Rmd Section 8 for the full reasoning. Grouping by
  # `family` here means this is correcting within THIS scope's battery only
  # (e.g. Mono-only "biomass" p-values get BH-corrected against each other,
  # not against the full dataset's biomass p-values) - the right behavior,
  # since each grid doc is its own self-contained set of tests.
  term_summary <- term_summary %>%
    dplyr::group_by(family) %>%
    dplyr::mutate(p_adj = stats::p.adjust(p.value, method = "BH")) %>%
    dplyr::ungroup()

  list(
    spec_by_label      = spec_by_label,
    results_list       = results_list,
    assumption_summary = assumption_summary,
    term_summary       = term_summary
  )
}

# -----------------------------------------------------------------------------
# Backward-compatible top-level objects (scope = "full")
# -----------------------------------------------------------------------------
# ANOVA_battery_DRAFT.Rmd and Significance_Grid_DRAFT.Rmd both expect these
# exact top-level names, unchanged from before this file had functions in it.
battery_spec <- build_battery_spec("full")
batt         <- run_battery(battery_spec)

spec_by_label      <- batt$spec_by_label
results_list       <- batt$results_list
assumption_summary <- batt$assumption_summary
term_summary       <- batt$term_summary
