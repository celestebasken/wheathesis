# =============================================================================
# Figure 3 redesign - drop-in replacement for the `figure` chunk in
# 14_Fig3_FavaSuppression.Rmd
#
# WHY THE OLD VERSION LOOKED DETACHED
#   The labels sat at a fixed y = max(fava) * 1.2 with a second header at
#   * 1.35, and coord_cartesian reserved ylim up to * 1.37. That is ~37% of the
#   panel height held empty purely to park text, with nothing visually linking
#   "-72%" to the pair of boxes underneath it. Hence the floating.
#
# WHAT CHANGED
#   1. Labels are ANCHORED - a bracket spans the two groups being compared and
#      sits just above that water level's own data. No reserved dead zone.
#   2. Headroom cut from 37% to ~18%, and the remaining space is used.
#   3. The floating italic header moved into the subtitle, where it belongs.
#   4. CI on ONE line, not two.
#   5. Non-significant effect greyed, so "big / big / nothing" reads instantly.
#   6. Boxplots -> points + model mean with 95% CI. At n = 5 a box displays
#      quartiles computed from five numbers, which implies more distributional
#      information than exists; the crossbar shows exactly what the model
#      estimates. (To keep boxes instead, see the note at the bottom.)
# =============================================================================

library(tidyverse)
library(emmeans)
library(car)
source("../fig_style.R")

inter <- readRDS("../wheat_audited.rds") %>% filter(Poly == "Inter")
m     <- lm(log(fava_npp) ~ Inoc * Water, data = inter)
p_iw  <- car::Anova(m, type = 2)["Inoc:Water", "Pr(>F)"]

# --- model means (geometric, back-transformed) and effect sizes --------------
mm <- emmeans(m, ~ Inoc | Water) %>%
  as.data.frame() %>%
  transmute(Water = factor(Water, levels(inter$Water)),
            Inoc  = factor(Inoc,  levels(inter$Inoc)),
            fit = exp(emmean), lo = exp(lower.CL), hi = exp(upper.CL))

eff <- emmeans(m, ~ Inoc | Water) %>%
  contrast("pairwise", by = "Water", adjust = "none") %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as.data.frame() %>%
  transmute(Water = factor(Water, levels(inter$Water)),
            pct = (exp(estimate) - 1) * 100,
            lo  = (exp(lower.CL) - 1) * 100,
            hi  = (exp(upper.CL) - 1) * 100,
            p   = p.value)

# --- bracket geometry: anchored to each water level's own data ---------------
DODGE <- 0.6
OFF   <- DODGE / 4            # two groups dodge to x +/- width/4

brk <- inter %>%
  group_by(Water) %>%
  summarise(top = max(fava_npp), .groups = "drop") %>%
  mutate(Water = factor(Water, levels(inter$Water)),
         xn    = as.numeric(Water),
         ybar  = top * 1.12,        # horizontal bar
         ytick = top * 1.05) %>%    # little downward ticks
  left_join(eff, by = "Water") %>%
  mutate(
    label = sprintf("%s%.0f%%  (%s%.0f to %s%.0f)",
                    ifelse(pct  < 0, "−", "+"), abs(pct),
                    ifelse(lo   < 0, "−", "+"), abs(lo),
                    ifelse(hi   < 0, "−", "+"), abs(hi)),
    sig   = p < 0.05,
    ylab_ = ybar * 1.045
  )

YMAX <- max(brk$ylab_) * 1.06
col  <- pal_poly[["Inter"]]

# --- figure ------------------------------------------------------------------
fig3 <- ggplot(inter, aes(x = Water, y = fava_npp)) +

  # bracket: horizontal bar + two downward ticks, per water level
  geom_segment(data = brk, aes(x = xn - OFF, xend = xn + OFF, y = ybar, yend = ybar),
               inherit.aes = FALSE, colour = "grey55", linewidth = 0.4) +
  geom_segment(data = brk, aes(x = xn - OFF, xend = xn - OFF, y = ybar, yend = ytick),
               inherit.aes = FALSE, colour = "grey55", linewidth = 0.4) +
  geom_segment(data = brk, aes(x = xn + OFF, xend = xn + OFF, y = ybar, yend = ytick),
               inherit.aes = FALSE, colour = "grey55", linewidth = 0.4) +
  geom_text(data = brk, aes(x = xn, y = ylab_, label = label,
                            colour = sig, fontface = ifelse(sig, "bold", "plain")),
            inherit.aes = FALSE, size = 3.4, vjust = 0, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "grey15", `FALSE` = "grey60"), guide = "none") +

  # raw pots
  geom_point(aes(shape = Inoc), colour = col, size = 2.1, stroke = 0.9, fill = NA,
             alpha = 0.75,
             position = position_jitterdodge(jitter.width = 0.10,
                                             dodge.width = DODGE, seed = 1)) +

  # model geometric mean + 95% CI
  geom_linerange(data = mm, aes(y = fit, ymin = lo, ymax = hi, group = Inoc),
                 colour = col, linewidth = 0.6,
                 position = position_dodge(width = DODGE)) +
  geom_point(data = mm, aes(y = fit, shape = Inoc, fill = Inoc),
             colour = col, size = 3.4, stroke = 1.1,
             position = position_dodge(width = DODGE)) +

  scale_shape_manual(values = c(Inoc = 21, Strl = 21), labels = inoc_lab) +
  scale_fill_manual(values = c(Inoc = col, Strl = "white"), labels = inoc_lab) +
  scale_x_discrete(labels = water_lab) +
  scale_y_continuous(breaks = seq(0, 10, 2), limits = c(0, YMAX),
                     expand = expansion(mult = c(0.01, 0))) +
  guides(shape = guide_legend(override.aes = list(size = 3.4, colour = col)),
         fill  = guide_legend(override.aes = list(size = 3.4, colour = col))) +
  labs(x = NULL, y = "Fava biomass per plant (g)",
       title = "Inoculated wheat suppressed its fava neighbor — except under severe drought",
       subtitle = paste0("Change in fava biomass with inoculated vs sterile wheat (95% CI). ",
                         "Large points = model mean; small points = individual pots (n = 5)."),
       caption = "EXPLORATORY") +
  theme_ms() +
  theme(plot.title    = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(colour = "grey35", size = 9.3, margin = margin(b = 8)),
        plot.caption  = element_text(colour = "grey55", size = 8, hjust = 0),
        panel.grid.major.x = element_blank())

print(fig3)

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Fig3_FavaSuppression.pdf", fig3, width = 7.5, height = 5)
ggsave("figures/Fig3_FavaSuppression.png", fig3, width = 7.5, height = 5,
       dpi = 300, bg = "white")

# -----------------------------------------------------------------------------
# To keep boxplots instead of mean + CI, delete the two mm layers above and use:
#
#   geom_boxplot(aes(fill = Inoc), colour = col, linewidth = 0.7, width = 0.5,
#                position = position_dodge(DODGE), outlier.shape = NA) +
#
# and set scale_fill_manual(values = c(Inoc = scales::alpha(col, 0.45),
#                                      Strl = "white")).
# The bracket geometry needs no change.
# -----------------------------------------------------------------------------
