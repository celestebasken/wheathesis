# =============================================================================
# fig_style.R - shared look for manuscript Figures 2-4 (Figure 1 predates this
# file and defines the same palette inline in 11_Figure1_Reversal/).
#
# Encoding rule used in EVERY figure so a reader never has to relearn it:
#   color  = stand type   (orange = monoculture, blue = wheat/ fava intercrop)
#   fill   = inoculation  (filled symbol / solid line = inoculated,
#                          open symbol / dashed line = sterile)
# Inoculation is deliberately NOT given its own hue: any third hue fails the
# colorblind-separation check against the orange/blue pair.
# =============================================================================
library(ggplot2)

pal_poly <- c(Mono = "#CC5500", Inter = "#2166ac")   # validated with validate_palette.js
poly_lab <- c(Mono = "Monoculture (12 wheat)",
              Inter = "Wheat/ fava intercrop (6 wheat + 1 fava)")
shp_inoc <- c(Inoc = 16, Strl = 1)                   # filled vs open circle
lty_inoc <- c(Inoc = "solid", Strl = "22")           # solid vs dashed
inoc_lab <- c(Inoc = "Inoculated", Strl = "Sterile")
water_lab <- c(Control = "Control (80% FC)", Moderate = "Moderate (60% FC)", Severe = "Severe (40% FC)")

theme_ms <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          strip.text = element_text(face = "bold", size = base_size - 1),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3))
}

# =============================================================================
# CONVENTION (Nov 2026): figures that cross inoculation with stand type
#
#   panels = stand type, ALWAYS Mono on the LEFT, Inter on the RIGHT  -> facet_poly()
#   inside each panel: x = water; the contrast is inoculation:
#       Inoculated = filled point, solid line,  bright orange #CC5500 (same as pal_poly Mono;
#                    was burnt orange #8F350A until Nov 10 2026, brightened at Celeste's request)
#       Sterile    = open point,   dashed line, slate blue   #6B90CC
#   p-values above each water level = inoculated vs sterile within that
#   stand type x water cell (unadjusted, from the figure's own model),
#   bold when p < 0.05, grey otherwise -> p_label(), p_text_colour()
#   No subtitles or captions on the plot itself; explain in the Rmd text.
#
# Colours checked with colorspace: hue difference stays large under deutan /
# protan / tritan simulation (deltaE 86-100 for #CC5500 vs #6B90CC) and both are
# >= 3:1 against white (4.3, 3.2). Greyscale lightness is close (L* 51 vs 59),
# so in greyscale fill (filled vs open) and line type (solid vs dashed) carry the
# contrast; they encode the same information redundantly by design.
#
# This is ADDITIVE: pal_poly / shp_inoc / lty_inoc above are unchanged and
# still used by figures that colour by stand type.
# =============================================================================
poly_levels    <- c("Mono", "Inter")
poly_lab_strip <- c(Mono = "Monoculture", Inter = "Wheat/ fava intercrop")   # as in 13_Figure2_Allometry
pal_inoc       <- c(Inoc = "#CC5500", Strl = "#6B90CC")

scale_inoc_color    <- function(...) scale_colour_manual(values = pal_inoc, labels = inoc_lab, name = NULL, ...)
scale_inoc_fill     <- function(...) scale_fill_manual(values = c(Inoc = pal_inoc[["Inoc"]], Strl = "white"),
                                                       labels = inoc_lab, name = NULL, ...)
scale_inoc_shape    <- function(...) scale_shape_manual(values = c(Inoc = 21, Strl = 21), labels = inoc_lab, name = NULL, ...)
scale_inoc_linetype <- function(...) scale_linetype_manual(values = lty_inoc, labels = inoc_lab, name = NULL, ...)
inoc_scales <- function() list(scale_inoc_color(), scale_inoc_fill(), scale_inoc_shape(), scale_inoc_linetype())

# set Poly to the Mono | Inter order before plotting: mutate(Poly = poly_factor(Poly))
poly_factor <- function(x) factor(as.character(x), poly_levels)
facet_poly  <- function(...) facet_wrap(~ Poly, nrow = 1, labeller = labeller(Poly = poly_lab_strip), ...)

p_label       <- function(p) ifelse(p < 1e-4, "p < 0.0001", paste0("p = ", signif(p, 2)))
p_text_colour <- function(p) ifelse(p < 0.05, "grey15", "grey60")
p_text_face   <- function(p) ifelse(p < 0.05, "bold", "plain")

theme_fig <- function(base_size = 12) {
  theme_ms(base_size) +
    theme(strip.text = element_text(face = "bold", size = base_size - 1, hjust = 0))
}
