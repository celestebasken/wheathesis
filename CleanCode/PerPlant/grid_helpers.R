# =============================================================================
# grid_helpers.R
# =============================================================================
# Shared plot/table-rendering logic used by ALL THREE significance-grid docs:
#   - Significance_Grid_DRAFT.Rmd        (full dataset)
#   - Significance_Grid_Inter_DRAFT.Rmd  (Inter pots only)
#   - Significance_Grid_Mono_DRAFT.Rmd   (Mono pots only)
#
# This file assumes battery_core.R has ALREADY been sourced (for tidyverse/
# ggplot2, and for run_anova_check()/build_battery_spec()/run_battery()) -
# it only adds the "turn a fitted battery into an HTML grid" layer on top.
#
# Pulling this out (rather than each grid doc defining its own copy of
# plot_cell()/direction_label()/etc., which is how the first draft did it)
# is the same reasoning as battery_core.R itself: three near-identical Rmds
# each hand-editing their own copy of this logic is exactly the kind of
# drift this whole cleanup project is trying to get away from.
#
# Nothing in this file prints/knits anything - it only builds objects.
# =============================================================================

library(kableExtra)   # kable_styling()/column_spec() for the final HTML table

# -----------------------------------------------------------------------------
# 1. Helper: turn "Inoculation:Water" into a direction sentence
# -----------------------------------------------------------------------------
# Signature is `(spec, model, term_vars)` - not just `(data, response,
# term_vars)` - even though this categorical version only ever needs
# `spec$data`/`spec$response` and ignores `model` entirely. The reason: the
# continuous-predictor variant (direction_label_continuous(), in
# grid_helpers_continuous.R) needs the FITTED MODEL instead of raw group
# means (a "direction" for a continuous slope has to come from a regression
# coefficient, not from averaging groups) - giving both versions the same
# signature is what lets build_cell_data() call whichever one it's given
# (via its `direction_fn` argument) without knowing which kind it is.
direction_label <- function(spec, model, term_vars) {
  data     <- spec$data
  response <- spec$response
  # term_vars is the term string split on ":", e.g. c("Water") for a main
  # effect, or c("Inoculation","Water") for a 2-way interaction.

  if (length(term_vars) == 1) {
    # ---- Main effect: just rank the group means, high to low ----
    form  <- stats::as.formula(paste(response, "~", term_vars[1]))
    means <- stats::aggregate(form, data = data, FUN = mean, na.rm = TRUE)
    means <- means[order(-means[[response]]), ]
    paste(means[[term_vars[1]]], collapse = " > ")

  } else if (length(term_vars) == 2) {
    # ---- 2-way interaction: rank factor A's levels SEPARATELY within each
    # level of factor B, e.g. "Inter: Inoc > Strl" / "Mono: Strl > Inoc" -
    # this is exactly the by-hand format you used in the old PDF table.
    fA <- term_vars[1]
    fB <- term_vars[2]
    form  <- stats::as.formula(paste(response, "~", fA, "+", fB))
    means <- stats::aggregate(form, data = data, FUN = mean, na.rm = TRUE)

    splits <- split(means, means[[fB]])
    lines <- vapply(names(splits), function(lvl) {
      sub <- splits[[lvl]]
      sub <- sub[order(-sub[[response]]), ]
      paste0(lvl, ": ", paste(sub[[fA]], collapse = " > "))
    }, character(1))
    paste(lines, collapse = "<br>")

  } else {
    # ---- 3-way interaction: too many combinations for a one-line summary -
    # leave a pointer to the thumbnail instead of guessing at a caption.
    "see facets ->"
  }
}

# -----------------------------------------------------------------------------
# 2. Helper: format a p-value into stars + text
# -----------------------------------------------------------------------------
# Standard significance-star convention. Deliberately keyed off p_adj (the
# same BH-adjusted value that decides the Yes/Kinda(<.1)/No label below the
# plot) - not raw p - so the star on the figure can never disagree with the
# text label next to it.
p_to_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.10  ~ ".",     # conventional "marginal" marker
    TRUE      ~ "ns"
  )
}

format_p <- function(p) {
  if (p < 0.001) "p < .001" else paste0("p = ", formatC(p, digits = 3, format = "f"))
}

# -----------------------------------------------------------------------------
# 2b. Helper: apply the p-value subtitle + significance-based styling to any
#     ggplot thumbnail, plus the 0-baseline guarantee
# -----------------------------------------------------------------------------
# Shared by plot_cell() (categorical factors) and plot_cell_continuous()
# (continuous hyphal colonization, in grid_helpers_continuous.R) - the two
# build completely different geometries (bar/line vs. scatter+regression),
# but this last "stamp the p-value on it and make it pop if significant"
# step is identical either way, so it's factored out once rather than
# copy/pasted between them.
apply_significance_subtitle <- function(p, p_val) {
  subtitle_txt <- paste0(format_p(p_val), "  ", p_to_stars(p_val))

  # `expand_limits(y = 0)` guarantees "0 is always at the bottom": geom_col()
  # already anchors bars at 0 by default, but line/scatter plots don't -
  # their y-axis just ranges over whatever the data happens to span, which
  # does NOT necessarily reach down to 0. This forces every panel's axis to
  # include 0 as its floor. (It only ever EXTENDS the range to reach 0 - it
  # can't clip data - so nothing gets cut off.)
  p +
    ggplot2::expand_limits(y = 0) +
    labs(subtitle = subtitle_txt) +
    theme(
      axis.title      = element_blank(),
      plot.margin     = margin(1, 1, 1, 1),
      # Bold + darker when the BH-adjusted result is significant (p < .05),
      # so a significant panel visually pops out of the grid even before
      # you read the number - grey/plain otherwise so "ns" doesn't compete
      # for attention.
      plot.subtitle   = element_text(
        size   = 5.5,
        face   = ifelse(p_val < 0.05, "bold", "plain"),
        colour = ifelse(p_val < 0.05, "black", "grey40"),
        hjust  = 0.5,
        margin = margin(b = 1)
      )
    )
}

# -----------------------------------------------------------------------------
# 3. Helper: build + save one cell's thumbnail plot
# -----------------------------------------------------------------------------
plot_cell <- function(data, response, term_vars, out_dir, cell_id, p_val) {

  # ---- Step 1: mean + standard error per group, however many factors ----
  summary_df <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(term_vars))) %>%
    dplyr::summarise(
      mean_y = mean(.data[[response]], na.rm = TRUE),
      se_y   = sd(.data[[response]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[response]]))),
      .groups = "drop"
    )

  # ---- Step 2: decide which factor goes on x / color / facet ----
  # Whichever factor has the MOST levels goes on the x-axis (in our design
  # that's always Water, when Water is involved, since it has 3 levels vs.
  # 2 for Inoculation/Polyculture). Second-most goes to color. Third (only
  # for the 3-way term, full-dataset scope only) goes to facets.
  n_levels     <- sapply(term_vars, function(v) dplyr::n_distinct(data[[v]]))
  ordered_vars <- term_vars[order(-n_levels)]

  x_var     <- ordered_vars[1]
  color_var <- if (length(ordered_vars) >= 2) ordered_vars[2] else NA
  facet_var <- if (length(ordered_vars) == 3) ordered_vars[3] else NA

  # ---- Step 3: build the plot ----
  # IMPORTANT ORDERING NOTE: theme_minimal() is a COMPLETE theme, not an
  # incremental tweak - adding it AFTER a theme(legend.position = "bottom")
  # call silently wipes that setting back to theme_minimal()'s own default
  # (legend on the right), which is what originally crammed the legend into
  # the corner of the interaction-plot thumbnails. Fix: apply theme_minimal()
  # FIRST, immediately on creation, and only ever use plain theme(...) calls
  # (which merge into the existing theme rather than replacing it) after
  # that point.
  if (length(term_vars) == 1) {
    p <- ggplot(summary_df, aes(x = .data[[x_var]], y = mean_y)) +
      theme_minimal(base_size = 6) +
      geom_col(fill = "grey30", width = 0.6) +
      geom_errorbar(aes(ymin = mean_y - se_y, ymax = mean_y + se_y), width = 0.15)
  } else {
    p <- ggplot(summary_df, aes(x = .data[[x_var]], y = mean_y,
                                 color = .data[[color_var]], group = .data[[color_var]])) +
      theme_minimal(base_size = 6) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 1) +
      geom_errorbar(aes(ymin = mean_y - se_y, ymax = mean_y + se_y), width = 0.15) +
      theme(legend.position = "bottom",
            legend.title    = element_blank(),
            legend.text     = element_text(size = 5),
            legend.key.size = unit(3, "mm"),
            legend.margin   = margin(t = -4))

    if (!is.na(facet_var)) {
      p <- p + facet_wrap(stats::as.formula(paste("~", facet_var)))
    }
  }

  # ---- Step 4: p-value + stars as a subtitle, + the 0-baseline guarantee ----
  # Baked into the figure itself - so if this thumbnail is ever pulled out
  # and used on its own (e.g. dropped into the thesis as a standalone
  # figure), it's still self-labeled and doesn't depend on the surrounding
  # table. (Paired with plot_response in build_battery_spec() - drawing the
  # natural 0-1 proportion scale instead of the logit scale used for the
  # model - every panel in every grid doc has a consistent, non-negative
  # baseline.)
  p <- apply_significance_subtitle(p, p_val)

  # ---- Step 5: save + return the file path ----
  out_path <- file.path(out_dir, paste0(cell_id, ".png"))
  ggplot2::ggsave(out_path, p, width = 1.3, height = 1.2, dpi = 150)
  out_path
}

# -----------------------------------------------------------------------------
# 4. Helper: cell background/label styling per significance level
# -----------------------------------------------------------------------------
# Yes/Kinda get a tinted cell background + colored bold label so they pop out
# of the grid at a glance; No stays plain/muted so it visually recedes
# instead of competing for attention. Plain if/else (not dplyr::case_when)
# because this runs once per row on a single scalar `sig` value each time -
# case_when is for vectorized use.
sig_style <- function(sig) {
  if (sig == "Yes") {
    list(bg = "#e9f7ec", fg = "#1a7431", weight = "bold")     # green tint
  } else if (sig == "Kinda (&lt;.1)") {
    list(bg = "#fff8e1", fg = "#8a6100", weight = "bold")     # amber tint
  } else {
    list(bg = "transparent", fg = "#999999", weight = "normal")  # No: muted, no tint
  }
}

# -----------------------------------------------------------------------------
# 5. Build one row per (response variable, ANOVA term) cell
# -----------------------------------------------------------------------------
# Takes this scope's term_summary + spec_by_label + results_list (all three
# from run_battery()), and where to save this scope's thumbnails - returns
# the full cell_data data frame, one row per label x term, ready to pivot
# into a grid.
#
# `direction_fn`/`plot_fn` default to the categorical-factor versions
# (direction_label()/plot_cell()) - the continuous-hyphal-colonization grid
# docs call this same function with direction_fn = direction_label_continuous,
# plot_fn = plot_cell_continuous (both in grid_helpers_continuous.R) instead,
# rather than needing their own copy of everything else in this function
# (styling, HTML assembly, the "always generate a figure" rule, etc.).
build_cell_data <- function(term_summary, spec_by_label, results_list, plot_dir,
                             direction_fn = direction_label, plot_fn = plot_cell) {

  cell_data <- term_summary %>%
    dplyr::mutate(
      term_vars = strsplit(term, ":"),
      sig = dplyr::case_when(
        p_adj < 0.05 ~ "Yes",
        p_adj < 0.10 ~ "Kinda (&lt;.1)",
        TRUE         ~ "No"
      )
    )

  # purrr::map2_chr walks label + term_vars together, row by row. `spec`
  # (data/response) comes from spec_by_label; `model` (the actual fitted
  # lm object) comes from results_list - direction_label() only needs the
  # former, direction_label_continuous() only needs the latter, but both are
  # always passed so build_cell_data() doesn't need to know which one this
  # particular direction_fn actually uses.
  cell_data$direction <- purrr::map2_chr(
    cell_data$label, cell_data$term_vars,
    function(lbl, tv) {
      spec  <- spec_by_label[[lbl]]
      model <- results_list[[lbl]]$model
      direction_fn(spec, model, tv)
    }
  )

  # A figure is generated for EVERY cell, "No" included - per your call,
  # you'd rather always have the picture available to glance at than save
  # the render time.
  #
  # spec$plot_response is used instead of spec$response when it's set - the
  # natural 0-1 scale for the logit-transformed percent/AMF variables (see
  # the comment on plot_response in battery_core.R).
  cell_data$plot_path <- purrr::pmap_chr(
    list(cell_data$label, cell_data$term, cell_data$term_vars, cell_data$p_adj),
    function(lbl, trm, tv, p_val) {
      spec          <- spec_by_label[[lbl]]
      plot_response <- if (is.null(spec$plot_response)) spec$response else spec$plot_response
      cell_id       <- paste(gsub("[^A-Za-z0-9]", "_", lbl), gsub(":", "_", trm), sep = "__")
      plot_fn(spec$data, plot_response, tv, out_dir = plot_dir, cell_id = cell_id, p_val = p_val)
    }
  )

  # ---- Assemble the actual HTML that will go inside each table cell ----
  # knitr::image_uri() base64-encodes the PNG straight into the <img> tag, so
  # the knitted HTML file is self-contained (no separate image files it
  # depends on to display correctly if you move/share the .html).
  cell_data$cell_html <- purrr::pmap_chr(
    list(cell_data$sig, cell_data$direction, cell_data$plot_path),
    function(sig, dir, path) {
      sty <- sig_style(sig)
      img_html <- paste0('<img src="', knitr::image_uri(path), '" width="70"><br>')
      dir_html <- paste0("<br><span style='font-size:85%;color:#555'>", dir, "</span>")

      paste0(
        "<div style='background-color:", sty$bg, "; padding:4px; border-radius:4px;'>",
        img_html,
        "<b style='color:", sty$fg, "; font-weight:", sty$weight, ";'>", sig, "</b>",
        dir_html,
        "</div>"
      )
    }
  )

  cell_data
}

# -----------------------------------------------------------------------------
# 6. Row labels/order per scope, and the grid-rendering function
# -----------------------------------------------------------------------------
# The full dataset's 3-way models produce up to 7 distinct ANOVA terms
# (including Polyculture and everything it interacts with); the Inter-only
# and Mono-only docs only ever have a 2-way Inoculation*Water model (3
# possible terms), since Polyculture is constant within either subset and
# can't be a predictor. Passing the right set in avoids ever displaying a
# "Poly:..." row that literally cannot occur in a single-polyculture grid.
term_labels_for_scope <- function(scope = c("full", "inter", "mono")) {
  scope <- match.arg(scope)
  if (scope == "full") {
    c(
      "Polyculture"                   = "Polyculture",
      "Inoculation"                   = "AMF inoculation",
      "Water"                         = "Water (PPT)",
      "Inoculation:Polyculture"       = "Poly:Inoc",
      "Polyculture:Water"             = "Poly:Water",
      "Inoculation:Water"             = "Inoc:Water",
      "Inoculation:Polyculture:Water" = "Poly:Inoc:Water (all 3)"
    )
  } else {
    c(
      "Inoculation"       = "AMF inoculation",
      "Water"             = "Water (PPT)",
      "Inoculation:Water" = "Inoc:Water"
    )
  }
}

# Renders one grid ("mass" or "percent") as a kable HTML table.
#   - cell_data:    from build_cell_data()
#   - battery_spec: from build_battery_spec() (same scope as cell_data) -
#                   only used here to get the column order/membership per
#                   display_group, so this file never has to repeat that
#                   list itself.
#   - term_labels:  from term_labels_for_scope() (same scope as cell_data)
render_grid <- function(group, cell_data, battery_spec, term_labels) {
  term_order <- names(term_labels)

  labels_in_group <- battery_spec %>%
    purrr::keep(~ .x$display_group == group) %>%
    purrr::map_chr("label")

  grid_df <- cell_data %>%
    dplyr::filter(label %in% labels_in_group) %>%
    dplyr::mutate(term = factor(term, levels = term_order)) %>%
    dplyr::select(term, label, cell_html) %>%
    tidyr::pivot_wider(names_from = label, values_from = cell_html) %>%
    dplyr::arrange(term) %>%
    # Blank out any term row this response family never fits (NA from the
    # pivot = "this term isn't in that variable's model") with a muted
    # "NA" tag instead of a literal R "NA" string.
    dplyr::mutate(dplyr::across(
      -term,
      ~ ifelse(is.na(.x), "<span style='color:#aaa'><i>NA</i></span>", .x)
    )) %>%
    dplyr::mutate(term = term_labels[as.character(term)]) %>%
    dplyr::select(term, dplyr::all_of(labels_in_group))

  grid_df %>%
    knitr::kable(format = "html", escape = FALSE, col.names = c("", labels_in_group)) %>%
    kableExtra::kable_styling(full_width = FALSE, font_size = 11) %>%
    kableExtra::column_spec(1, bold = TRUE, background = "#f0f0f0")
}
