# ============================================================
# 01_basic_figures.R
# Figures without p-values or letter labels
# ============================================================

source("00_setup_helpers.R")

# ------------------------------------------------------------
# A. Biomass/yield by precipitation x inoculation
# ------------------------------------------------------------

mono_biomass_long <- mono_only %>%
  dplyr::select(Water_gal, Water, Inoculation, Berries_g_pp, NPP_g_pp) %>%
  pivot_longer(
    cols = c(Berries_g_pp, NPP_g_pp),
    names_to = "Response",
    values_to = "Value"
  ) %>%
  mutate(
    Response = case_when(
      Response == "NPP_g_pp" ~ "Total Biomass",
      Response == "Berries_g_pp" ~ "Yield"
    ),
    Response = factor(Response, levels = c("Total Biomass", "Yield"))
  )

inoc_ppt_interact_basic <- ggplot(
  mono_biomass_long,
  aes(x = Water, y = Value, color = Inoculation, group = Inoculation)
) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.12, linewidth = 0.8) +
  facet_wrap(~ Response, scales = "free_y") +
  scale_color_manual(values = inoc_colors) +
  labs(
    title = "Interacting Effects of Precipitation and AMF Inoculation",
    x = "Drought treatment",
    y = "Biomass (g per plant)",
    color = "Inoculation"
  ) +
  theme_clean_fig()

ggsave(
  "PresentationFigs/inoc_ppt_interact_basic.png",
  inoc_ppt_interact_basic,
  width = 9,
  height = 4.8,
  dpi = 300
)


# ------------------------------------------------------------
# B. Resource allocation by precipitation x inoculation
# ------------------------------------------------------------

mono_alloc_long <- mono_only %>%
  dplyr::select(
    Water, Water_gal, Inoculation,
    percent_berries_of_wheat,
    percent_ANPP_of_wheat,
    percent_BNPP_of_wheat
  ) %>%
  pivot_longer(
    cols = c(
      percent_berries_of_wheat,
      percent_ANPP_of_wheat,
      percent_BNPP_of_wheat
    ),
    names_to = "Response",
    values_to = "Value"
  ) %>%
  mutate(
    Response = case_when(
      Response == "percent_berries_of_wheat" ~ "Yield",
      Response == "percent_ANPP_of_wheat" ~ "Aboveground Biomass",
      Response == "percent_BNPP_of_wheat" ~ "Belowground Biomass"
    ),
    Response = factor(
      Response,
      levels = c("Yield", "Aboveground Biomass", "Belowground Biomass")
    )
  )

resource_allocation_basic <- ggplot(
  mono_alloc_long,
  aes(x = Water, y = Value, color = Inoculation, group = Inoculation)
) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.12, linewidth = 0.8) +
  facet_wrap(~ Response, scales = "free_y", nrow = 1) +
  scale_color_manual(values = inoc_colors) +
  labs(
    title = "Impacts of Precipitation and Inoculation on Resource Allocation",
    x = "Drought treatment",
    y = "Percent allocation",
    color = "Inoculation"
  ) +
  theme_clean_fig()

ggsave(
  "PresentationFigs/resource_allocation_basic.png",
  resource_allocation_basic,
  width = 12,
  height = 4.8,
  dpi = 300
)