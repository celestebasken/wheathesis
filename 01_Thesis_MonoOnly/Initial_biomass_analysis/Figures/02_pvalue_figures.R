# ============================================================
# 02_pvalue_figures.R
# Figures with p-values, R², or compact letter displays
# ============================================================

source("00_setup_helpers.R")

# ------------------------------------------------------------
# Helper: compact letters comparing inoculation within water
# ------------------------------------------------------------

get_inoc_letters <- function(data) {
  data %>%
    group_by(Response, Water) %>%
    group_modify(~ {
      mod <- lm(Value ~ Inoculation, data = .x)
      em <- emmeans::emmeans(mod, ~ Inoculation)
      cld <- multcomp::cld(em, adjust = "sidak", Letters = letters)
      as.data.frame(cld)
    }) %>%
    mutate(.group = stringr::str_trim(.group)) %>%
    dplyr::select(Response, Water, Inoculation, .group)
}


# ------------------------------------------------------------
# A. Three-panel inoculation figure:
# Total biomass, yield, yield allocation
# ------------------------------------------------------------

inoc_response_long <- mono_only %>%
  dplyr::select(
    Water, Inoculation,
    NPP_g_pp,
    Berries_g_pp,
    percent_berries_of_wheat
  ) %>%
  pivot_longer(
    cols = c(NPP_g_pp, Berries_g_pp, percent_berries_of_wheat),
    names_to = "Response",
    values_to = "Value"
  ) %>%
  mutate(
    Response = case_when(
      Response == "NPP_g_pp" ~ "Total Biomass",
      Response == "Berries_g_pp" ~ "Yield",
      Response == "percent_berries_of_wheat" ~ "Yield Allocation"
    ),
    Response = factor(
      Response,
      levels = c("Total Biomass", "Yield", "Yield Allocation")
    )
  )

summary_inoc <- inoc_response_long %>%
  group_by(Response, Water, Inoculation) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    se = se(Value),
    ymin = mean - se,
    ymax = mean + se,
    .groups = "drop"
  )

letters_inoc <- get_inoc_letters(inoc_response_long)

label_positions_inoc <- summary_inoc %>%
  left_join(letters_inoc, by = c("Response", "Water", "Inoculation")) %>%
  group_by(Response) %>%
  mutate(
    panel_range = max(ymax, na.rm = TRUE) - min(ymin, na.rm = TRUE),
    label_y = ymax + 0.10 * panel_range
  ) %>%
  ungroup()

pd <- position_dodge(width = 0.3)

inoc_ppt_interact_3panel <- ggplot(
  summary_inoc,
  aes(x = Water, y = mean, color = Inoculation, group = Inoculation)
) +
  geom_line(linewidth = 1.2, position = pd) +
  geom_point(size = 3.8, position = pd) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = 0.12,
    linewidth = 0.8,
    position = pd
  ) +
  geom_text(
    data = label_positions_inoc,
    aes(x = Water, y = label_y, label = .group, color = Inoculation),
    position = pd,
    inherit.aes = FALSE,
    fontface = "bold",
    size = 5
  ) +
  facet_wrap(~ Response, scales = "free_y", nrow = 1) +
  scale_color_manual(values = inoc_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.25))) +
  labs(
    title = "Interacting Effects of Precipitation and AMF Inoculation",
    x = "Drought treatment",
    y = NULL,
    color = "Inoculation",
    caption = "Points show means ± SE. Different letters indicate Sidak-adjusted pairwise differences between inoculation treatments within each drought level (p < 0.05)."
  ) +
  theme_clean_fig()

ggsave(
  "PresentationFigs/Plabels/inoc_ppt_interact_3panel.png",
  inoc_ppt_interact_3panel,
  width = 12,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# B. Total Biomass + Yield figure
# ------------------------------------------------------------

biomass_yield_long <- mono_only %>%
  dplyr::select(Water, Inoculation, NPP_g_pp, Berries_g_pp) %>%
  pivot_longer(
    cols = c(NPP_g_pp, Berries_g_pp),
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

summary_biomass <- make_summary(biomass_yield_long)
letters_biomass <- make_inoc_letters(biomass_yield_long)
labels_biomass <- add_label_positions(summary_biomass, letters_biomass)

fig_biomass_yield_p <- plot_inoc_line_letters(
  summary_df = summary_biomass,
  labels_df = labels_biomass,
  title = "Interacting Effects of Precipitation and AMF Inoculation",
  ylab = "Biomass (g per plant)",
  caption = "Points show means ± SE. Different letters indicate Sidak-adjusted differences between inoculation treatments within each drought level.",
  nrow = 1
)

ggsave(
  "PresentationFigs/Plabels/inoc_ppt_interact.png",
  fig_biomass_yield_p,
  width = 11,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# c. ANPP / Berries / BNPP allocation figure
# ------------------------------------------------------------

allocation_long <- mono_only %>%
  dplyr::select(
    Water, Inoculation,
    percent_ANPP_of_wheat,
    percent_berries_of_wheat,
    percent_BNPP_of_wheat
  ) %>%
  pivot_longer(
    cols = c(
      percent_ANPP_of_wheat,
      percent_berries_of_wheat,
      percent_BNPP_of_wheat
    ),
    names_to = "Response",
    values_to = "Value"
  ) %>%
  mutate(
    Response = case_when(
      Response == "percent_ANPP_of_wheat" ~ "ANPP",
      Response == "percent_berries_of_wheat" ~ "Berries",
      Response == "percent_BNPP_of_wheat" ~ "BNPP"
    ),
    Response = factor(Response, levels = c("ANPP", "Berries", "BNPP"))
  )

summary_alloc <- make_summary(allocation_long)
letters_alloc <- make_inoc_letters(allocation_long)
labels_alloc <- add_label_positions(summary_alloc, letters_alloc)

fig_allocation_p <- plot_inoc_line_letters(
  summary_df = summary_alloc,
  labels_df = labels_alloc,
  title = "Impacts of Precipitation and Inoculation on Resource Allocation Within Wheat",
  ylab = "Percent allocation",
  caption = "Points show means ± SE. Different letters indicate Sidak-adjusted differences between inoculation treatments within each drought level.",
  nrow = 1
)

ggsave(
  "PresentationFigs/Plabels/resource_allocation.png",
  fig_allocation_p,
  width = 12,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# D. hyphal colonization regression figure
# ------------------------------------------------------------

hyphal_long <- mono_only %>%
  filter(Water %in% c("Control", "Severe")) %>%
  dplyr::select(
    hyphal_colonization,
    Water,
    percent_berries_of_wheat,
    NPP_g_pp
  ) %>%
  pivot_longer(
    cols = c(percent_berries_of_wheat, NPP_g_pp),
    names_to = "Response",
    values_to = "Value"
  ) %>%
  mutate(
    Water = factor(Water, levels = c("Control", "Severe")),
    Response = case_when(
      Response == "percent_berries_of_wheat" ~ "Allocation to Yield (%)",
      Response == "NPP_g_pp" ~ "Total Biomass (g per plant)"
    )
  )

hyphal_labels <- hyphal_long %>%
  group_by(Response, Water) %>%
  group_modify(~ {
    mod <- lm(Value ~ hyphal_colonization, data = .x)
    tidy_mod <- broom::tidy(mod)
    glance_mod <- broom::glance(mod)
    
    p_val <- tidy_mod %>%
      filter(term == "hyphal_colonization") %>%
      pull(p.value)
    
    tibble(
      label = paste0("p = ", signif(p_val, 2), "\nR² = ", round(glance_mod$r.squared, 2)),
      x = -Inf,
      y = -Inf
    )
  })

fig_hyphal_p <- ggplot(
  hyphal_long,
  aes(x = hyphal_colonization, y = Value)
) +
  geom_point(size = 2.8, alpha = 0.8, color = "black") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#3B6C8E",
    fill = "#3B6C8E",
    alpha = 0.2,
    linewidth = 1.2
  ) +
  geom_text(
    data = hyphal_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.15,
    vjust = -0.4,
    size = 4.2
  ) +
  facet_grid(Response ~ Water, scales = "free_y") +
  labs(
    x = "Root hyphal colonization (%)",
    y = NULL,
    title = "AMF Hyphal Colonization Impact on Yield Allocation and Total Biomass"
  ) +
  coord_cartesian(clip = "off") +
  theme_clean_fig()

ggsave(
  "PresentationFigs/Plabels/hyphal_berries_NPP.png",
  fig_hyphal_p,
  width = 9,
  height = 6.5,
  dpi = 300
)
