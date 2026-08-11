# =============================================================================
# 01_data_prep.R
# =============================================================================
# Step 1 of Instructions/AMF_intercropping_analysis_plan.md ("Data audit").
#
# Builds the single working dataset (`wheat`) used by every later step in this
# folder - joins the wheat-level data (Data/full_dataset.csv) to the newly
# added structure-count data (Data/AMF_colonization.csv), and constructs every
# derived variable Step 1 asks for. Column names here deliberately match the
# short names used throughout the instructions doc's own code snippets
# (Poly, Inoc, Water, berry, veg, berry_alloc, hyphal, arb_pct, ves_pct,
# arb_n, ves_n, hyph_n, fava_npp, fava_pct_pot, pot_productivity) so later
# steps can follow the doc's snippets close to verbatim, rather than needing
# a name-translation layer in every script.
#
# The ORIGINAL full_dataset.csv columns (W_ANPP_g_pp etc.) are kept alongside
# the short-name aliases, not dropped - so nothing here is a one-way,
# irreversible rename, just an added, easier-to-use interface on top.
#
# Nothing in this file prints/knits anything - see 01_Data_Audit.Rmd for the
# actual audit report (balance check, boundary check, etc.) built from this.
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 1. Load both data files
# -----------------------------------------------------------------------------
# Path is relative to C_Analyses/, a sibling of Data/.
full  <- read.csv("../Data/full_dataset.csv")
amf   <- read.csv("../Data/AMF_colonization.csv")

# AMF_colonization.csv has two columns both literally named "Pot" in the raw
# CSV header - read.csv() auto-disambiguates the second one to "Pot.1". The
# first ("Pot") is the compound name string (e.g. "1_Inoc_inter_c_1"); the
# second ("Pot.1") is the numeric pot ID (1-60) - the same numbering used in
# full_dataset.csv's own `Pot` column, and the actual join key.
amf <- amf %>%
  rename(AMF_Name = Pot, Pot = Pot.1)

# -----------------------------------------------------------------------------
# 2. Join by Pot number
# -----------------------------------------------------------------------------
# Only pull in the columns full_dataset.csv doesn't already have (the raw
# structure counts) - hyphal_colonization/arb_vesc_colonization/Inoculation/
# Polyculture/Water/Replicate already exist in `full` and are identical
# across both files (verified during the audit, see 01_Data_Audit.Rmd) - no
# need to duplicate or risk a silent mismatch by re-joining them too.
wheat <- full %>%
  left_join(
    amf %>% select(Pot, Arbuscules, Vesicles, A_V, Hyphae_only, None, Date_assessed),
    by = "Pot"
  )

# -----------------------------------------------------------------------------
# 3. Short-name aliases matching the instructions doc's own code snippets
# -----------------------------------------------------------------------------
wheat <- wheat %>%
  mutate(
    Poly  = factor(Polyculture),
    Inoc  = factor(Inoculation),
    Water = factor(Water),   # alphabetical order already gives Control,
                               # Moderate, Severe - the order the doc's planned
                               # contrasts (set in 03_primary_model, not here)
                               # assume.

    # ---- Per-plant biomass (doc's berry/veg/NPP) ----
    berry = W_Berries_g_pp,
    ANPP  = W_ANPP_g_pp,
    BNPP  = W_BNPP_g_pp,
    veg   = W_ANPP_g_pp + W_BNPP_g_pp,   # NA for the one pot missing BNPP - see Section 5
    NPP   = W_NPP_pp,

    # ---- Berry allocation: THE resolved definition (see Section 6) ----
    # berry_alloc = berry / (ANPP + BNPP + berry) = berry / total NPP.
    # full_dataset.csv's own W_pcent_berries already IS this definition -
    # recomputed independently here from the raw per-plant masses as a
    # cross-check that they agree, rather than trusting the label alone.
    berry_alloc = berry / NPP,

    # ---- Pot-level productivity ----
    pot_productivity = Pot_NPP_g_pp,

    # ---- Fava (doc's fava_npp / fava_pct_pot) ----
    fava_npp     = Fava_g,
    fava_pct_pot = Pot_pcent_fava,

    # ---- AMF colonization: overall presence (unchanged from full_dataset.csv) ----
    hyphal = hyphal_colonization,

    # ---- AMF colonization: TRUE arbuscule/vesicle counts ----
    # A_V is a "both structures visible in this grid intersect" bucket, not a
    # third structure type - per the instructions doc, the true arbuscule
    # count is Arbuscules + A_V, and the true vesicle count is Vesicles + A_V.
    # These are raw counts out of 100 grid intersects (verified in the audit:
    # Arbuscules + Vesicles + A_V + Hyphae_only + None == 100 for all 60 pots).
    arb_n = Arbuscules + A_V,
    ves_n = Vesicles + A_V,
    hyph_n = 100 - None,   # total intersects showing ANY fungal structure -
                            # matches hyphal_colonization * 100 (cross-checked
                            # in the audit)
    arb_pct = arb_n / 100,
    ves_pct = ves_n / 100,

    # ---- log arbuscule:vesicle ratio ----
    # Continuity-corrected on the RAW COUNTS (Haldane-Anscombe style: +0.5 to
    # numerator and denominator), not on the 0-1 proportions - a +0.5 offset
    # only behaves like a sensible small-sample correction when it's added to
    # a count out of 100, not to a value already scaled to [0,1]. This avoids
    # log(0) for any pot with zero arbuscules or zero vesicles while staying
    # a small correction relative to a 100-count denominator.
    AV = log((arb_n + 0.5) / (ves_n + 0.5))
  )

# -----------------------------------------------------------------------------
# 4. Wheat-equivalent density for Inter pots (Section 1 of the plan)
# -----------------------------------------------------------------------------
# "Express fava biomass per pot in units of mean per-plant wheat biomass
# within the same water x inoc cell." Read literally: for each Inter pot,
# convert its OWN fava biomass into "how many wheat plants would have
# produced this much biomass," using the mean per-plant wheat productivity
# of the 5 Inter pots sharing that same Inoculation x Water cell as the
# conversion rate (not a single pot's own W_NPP_pp, to avoid the conversion
# rate being noisy/circular for any one pot).
cell_mean_wheat_pp <- wheat %>%
  filter(Poly == "Inter") %>%
  group_by(Inoc, Water) %>%
  summarise(mean_wheat_pp = mean(NPP, na.rm = TRUE), .groups = "drop")

wheat <- wheat %>%
  left_join(cell_mean_wheat_pp, by = c("Inoc", "Water")) %>%
  mutate(
    # NA (not 0) for Mono pots - "wheat-equivalent density" isn't a
    # meaningful quantity for a pot that has no fava to convert.
    wheat_equiv_fava    = if_else(Poly == "Inter", fava_npp / mean_wheat_pp, NA_real_),
    wheat_equiv_density = if_else(Poly == "Inter", Num_wheat_pot + wheat_equiv_fava, NA_real_),
    # As a fraction of Mono's 12-plant density, for the "roughly two-thirds"
    # framing in the instructions doc.
    wheat_equiv_density_frac_of_mono = wheat_equiv_density / 12
  ) %>%
  select(-mean_wheat_pp)   # only needed as an intermediate join column
