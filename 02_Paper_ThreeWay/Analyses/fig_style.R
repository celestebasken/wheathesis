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
