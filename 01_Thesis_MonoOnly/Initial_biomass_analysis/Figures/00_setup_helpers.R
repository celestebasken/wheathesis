# 00_setup_helpers.R

library(tidyverse)
library(broom)
library(emmeans)
library(multcomp)
library(stringr)

combined <- read.csv("../../Cross_compare/NPP_AMF_combined.csv")

combined <- combined[-4, ] %>%
  mutate(
    Water = factor(Water, levels = c("Control", "Moderate", "Severe")),
    Inoculation = factor(Inoculation, levels = c("Inoc", "Strl")),
    Polyculture = factor(Polyculture),
    percent_vegetative_of_wheat = percent_ANPP_of_wheat + percent_BNPP_of_wheat
  )

mono_only <- combined %>% filter(Polyculture != "Inter")
inter_only <- combined %>% filter(Polyculture != "Mono")

dir.create("PresentationFigs", recursive = TRUE, showWarnings = FALSE)
dir.create("PresentationFigs/Plabels", recursive = TRUE, showWarnings = FALSE)

se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

inoc_colors <- c("Inoc" = "#3B6C8E", "Strl" = "#C2A46D")

theme_clean_fig <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      strip.text = element_text(face = "bold", size = 13),
      strip.background = element_rect(fill = "#EAF1F5", color = NA),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black"),
      panel.spacing = unit(1.4, "lines"),
      legend.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 10)
    )
}

make_summary <- function(df) {
  df %>%
    group_by(Response, Water, Inoculation) %>%
    summarise(
      mean = mean(Value, na.rm = TRUE),
      se = se(Value),
      ymin = mean - se,
      ymax = mean + se,
      .groups = "drop"
    )
}

make_inoc_letters <- function(df, adjust_method = "sidak") {
  df %>%
    group_by(Response, Water) %>%
    group_modify(~ {
      mod <- lm(Value ~ Inoculation, data = .x)
      em <- emmeans::emmeans(mod, ~ Inoculation)
      cld <- multcomp::cld(em, adjust = adjust_method, Letters = letters)
      as.data.frame(cld)
    }) %>%
    mutate(.group = stringr::str_trim(.group)) %>%
    dplyr::select(Response, Water, Inoculation, .group)
}

add_label_positions <- function(summary_df, letters_df) {
  summary_df %>%
    left_join(letters_df, by = c("Response", "Water", "Inoculation")) %>%
    group_by(Response) %>%
    mutate(
      panel_range = max(ymax, na.rm = TRUE) - min(ymin, na.rm = TRUE),
      label_y = ymax + 0.10 * panel_range
    ) %>%
    ungroup()
}

plot_inoc_line_letters <- function(summary_df, labels_df, title, ylab, caption, nrow = NULL) {
  pd <- position_dodge(width = 0.3)
  
  ggplot(
    summary_df,
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
      data = labels_df,
      aes(x = Water, y = label_y, label = .group, color = Inoculation),
      inherit.aes = FALSE,
      position = pd,
      fontface = "bold",
      size = 5
    ) +
    facet_wrap(~ Response, scales = "free_y", nrow = nrow) +
    scale_color_manual(values = inoc_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.25))) +
    labs(
      title = title,
      x = "Drought treatment",
      y = ylab,
      color = "Inoculation",
      caption = caption
    ) +
    theme_clean_fig()
}