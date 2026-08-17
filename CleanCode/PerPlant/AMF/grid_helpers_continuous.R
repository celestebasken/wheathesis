# =============================================================================
# grid_helpers_continuous.R
# =============================================================================
# Sibling to grid_helpers.R, for the "hyphal_colonization as a CONTINUOUS
# predictor" grids. Sourced by:
#   - Significance_Grid_Hyphal_Full_DRAFT.Rmd
#   - Significance_Grid_Hyphal_Inter_DRAFT.Rmd
#   - Significance_Grid_Hyphal_Mono_DRAFT.Rmd
#
# This file assumes grid_helpers.R has ALREADY been sourced - it reuses
# sig_style(), p_to_stars(), format_p(), apply_significance_subtitle(),
# build_cell_data(), and render_grid() UNCHANGED (none of them care whether
# the underlying predictor is categorical or continuous - build_cell_data()
# was generalized with `direction_fn`/`plot_fn` arguments specifically so
# this file only needs to supply the two things that genuinely differ:
#
#   1. direction_label_continuous() - a "direction" for a continuous slope
#      has to come from the fitted model's actual regression coefficient
#      (via emtrends()), not from averaging raw group means the way
#      direction_label() does - there's no "high group"/"low group" to rank.
#
#   2. plot_cell_continuous() - a scatter + linear regression line instead
#      of a bar/line-of-means chart, since the x-axis is now a continuous
#      measurement, not a handful of discrete treatment levels.
#
# Nothing in this file prints/knits anything - it only builds objects.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Helper: turn a hyphal_colonization term into a slope-direction sentence
# -----------------------------------------------------------------------------
# Same `(spec, model, term_vars)` signature as direction_label() (see the
# comment on that function in grid_helpers.R for why) - this version uses
# `model`, not `spec`, since a slope has to come from the fitted regression,
# not from group means.
direction_label_continuous <- function(spec, model, term_vars) {
  cont_var <- "hyphal_colonization"

  if (!(cont_var %in% term_vars)) {
    # ---- Purely categorical term (e.g. "Polyculture", "Water", or
    # "Polyculture:Water") - hyphal_colonization isn't part of THIS
    # particular term at all (it's an ANCOVA: some terms are still plain
    # categorical main effects/interactions), so there's no slope to
    # describe here - it's a mean/intercept comparison, exactly what
    # direction_label() already computes. Delegate to it directly rather
    # than letting `setdiff(term_vars, cont_var)` silently fall through to
    # the slope logic below, which would compute the wrong thing (an
    # hyphal_colonization slope per group) for a term that has nothing to
    # do with that slope. (This was a real bug in the first version of this
    # function - the "Polyculture" row and the "hyphal_colonization:
    # Polyculture" row rendered identically before this check existed.)
    return(direction_label(spec, model, term_vars))
  }

  group_vars <- setdiff(term_vars, cont_var)

  if (length(group_vars) == 0) {
    # ---- Main effect of hyphal_colonization alone: one overall slope ----
    # `specs = ~1` marginalizes over any other factors in the model (e.g.
    # Water), giving a single average slope rather than one per group.
    tr    <- as.data.frame(emmeans::emtrends(model, specs = ~1, var = cont_var))
    slope <- tr[[paste0(cont_var, ".trend")]][1]
    arrow <- if (slope >= 0) "↑ response" else "↓ response"
    paste0("↑ hyphal colonization → ", arrow,
           " (slope = ", formatC(slope, digits = 3, format = "f"), ")")

  } else {
    # ---- Interaction: one slope per level of the grouping factor(s) ----
    # e.g. group_vars = c("Water") -> one slope per Water level;
    #      group_vars = c("Polyculture","Water") -> one slope per combination
    # (only occurs in scope = "full", where the model still has both).
    specs <- stats::as.formula(paste("~", paste(group_vars, collapse = "+")))
    tr    <- as.data.frame(emmeans::emtrends(model, specs = specs, var = cont_var))

    trend_col <- paste0(cont_var, ".trend")
    lines <- character(nrow(tr))
    for (i in seq_len(nrow(tr))) {
      lvl_txt  <- paste(sapply(group_vars, function(v) as.character(tr[[v]][i])), collapse = ", ")
      slope    <- tr[[trend_col]][i]
      sign_txt <- if (slope >= 0) "+" else "-"
      lines[i] <- paste0(lvl_txt, ": slope ", sign_txt, formatC(abs(slope), digits = 3, format = "f"))
    }
    paste(lines, collapse = "<br>")
  }
}

# -----------------------------------------------------------------------------
# 2. Helper: build + save one cell's thumbnail plot (scatter + regression)
# -----------------------------------------------------------------------------
plot_cell_continuous <- function(data, response, term_vars, out_dir, cell_id, p_val) {

  cont_var <- "hyphal_colonization"

  if (!(cont_var %in% term_vars)) {
    # ---- Purely categorical term (e.g. "Polyculture", "Water", or
    # "Polyculture:Water") - hyphal_colonization isn't part of THIS term, so
    # a scatter-vs-hyphal_colonization plot would be showing the wrong
    # thing (it would look identical to whatever hyphal_colonization
    # interaction row shares the same categorical factor, which is exactly
    # the bug this check fixes - see the matching note in
    # direction_label_continuous()). Fall back to the categorical grid's
    # own bar/line-of-means chart (plot_cell(), from grid_helpers.R) - a
    # mean+SE comparison across this term's factor levels, averaged over
    # hyphal_colonization, is the honest way to visualize a term that has
    # nothing to do with the continuous predictor.
    return(plot_cell(data, response, term_vars, out_dir, cell_id, p_val))
  }

  group_vars <- setdiff(term_vars, cont_var)

  # ---- Build the plot ----
  # theme_minimal() FIRST, same ordering rule as plot_cell() in
  # grid_helpers.R (see that function's comment for why order matters here -
  # theme_minimal() is a complete theme, adding it after legend.position
  # would silently wipe that setting back to the default).
  if (length(group_vars) == 0) {
    # ---- Main effect: one scatter + one regression line ----
    p <- ggplot(data, aes(x = .data[[cont_var]], y = .data[[response]])) +
      theme_minimal(base_size = 6) +
      geom_point(size = 0.8, alpha = 0.6, color = "grey40") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.6, color = "steelblue")

  } else if (length(group_vars) == 1) {
    # ---- 2-way: one regression line per level of the grouping factor ----
    p <- ggplot(data, aes(x = .data[[cont_var]], y = .data[[response]], color = .data[[group_vars[1]]])) +
      theme_minimal(base_size = 6) +
      geom_point(size = 0.8, alpha = 0.6) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.6) +
      theme(legend.position = "bottom",
            legend.title    = element_blank(),
            legend.text     = element_text(size = 5),
            legend.key.size = unit(3, "mm"),
            legend.margin   = margin(t = -4))

  } else {
    # ---- 3-way (scope = "full" only): color by the 2-level factor, facet by
    # the 3-level one ----  same "more levels -> facet, fewer -> color" rule
    # plot_cell() uses, so a 2x3 grid of tiny scatters isn't cluttered
    # further by trying to cram a third dimension into color/shape too.
    n_levels     <- sapply(group_vars, function(v) dplyr::n_distinct(data[[v]]))
    ordered_vars <- group_vars[order(-n_levels)]
    facet_var    <- ordered_vars[1]
    color_var    <- ordered_vars[2]

    p <- ggplot(data, aes(x = .data[[cont_var]], y = .data[[response]], color = .data[[color_var]])) +
      theme_minimal(base_size = 6) +
      geom_point(size = 0.6, alpha = 0.6) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.5) +
      facet_wrap(stats::as.formula(paste("~", facet_var))) +
      theme(legend.position = "bottom",
            legend.title    = element_blank(),
            legend.text     = element_text(size = 5),
            legend.key.size = unit(3, "mm"),
            legend.margin   = margin(t = -4))
  }

  # p-value/stars subtitle + 0-baseline guarantee - identical treatment to
  # plot_cell(), via the function shared between them.
  p <- apply_significance_subtitle(p, p_val)

  out_path <- file.path(out_dir, paste0(cell_id, ".png"))
  ggplot2::ggsave(out_path, p, width = 1.3, height = 1.2, dpi = 150)
  out_path
}

# -----------------------------------------------------------------------------
# 3. Row labels/order per scope
# -----------------------------------------------------------------------------
# Same 3-vs-7-term structure as term_labels_for_scope() in grid_helpers.R -
# "full" keeps Polyculture and everything it interacts with; "inter"/"mono"
# only ever have hyphal_colonization*Water (2-way), since Polyculture is
# constant within either subset.
term_labels_for_scope_continuous <- function(scope = c("full", "inter", "mono")) {
  scope <- match.arg(scope)
  if (scope == "full") {
    c(
      "hyphal_colonization"                  = "Hyphal colonization",
      "Polyculture"                          = "Polyculture",
      "Water"                                = "Water (PPT)",
      "hyphal_colonization:Polyculture"      = "Hyphal:Poly",
      "Polyculture:Water"                    = "Poly:Water",
      "hyphal_colonization:Water"            = "Hyphal:Water",
      "hyphal_colonization:Polyculture:Water" = "Hyphal:Poly:Water (all 3)"
    )
  } else {
    c(
      "hyphal_colonization"        = "Hyphal colonization",
      "Water"                      = "Water (PPT)",
      "hyphal_colonization:Water"  = "Hyphal:Water"
    )
  }
}
