# Compatibility support for binary PRO-CTCAE indicator (_IND) items.
# Loaded after the primary implementations so existing score/composite behavior
# remains unchanged.

.toxFigures_score_items <- toxFigures
.toxTables_score_items <- toxTables

.proae_indicator_items <- function(dsn) {
  item_names <- toupper(names(dsn))
  ref <- proctcae_vars
  ref$name <- as.character(ref$name)
  ref$fmt <- as.character(ref$fmt)
  ref$name[ref$name %in% item_names & ref$fmt %in% c("yn_2_fmt", "yn_3_fmt", "yn_4_fmt")]
}

.proae_non_indicator_data <- function(dsn, indicator_items, required) {
  keep <- unique(c(required, names(dsn)[!(toupper(names(dsn)) %in% indicator_items)]))
  dsn[, keep[keep %in% names(dsn)], drop = FALSE]
}

.proae_indicator_summary <- function(dsn, id_var, cycle_var, item, baseline_val, type) {
  post <- dsn[!is.na(dsn[[cycle_var]]) & dsn[[cycle_var]] != baseline_val,
              c(id_var, cycle_var, item), drop = FALSE]
  ids <- unique(dsn[[id_var]])
  out <- data.frame(id = ids, stringsAsFactors = FALSE)
  names(out)[1] <- id_var

  max_post <- stats::aggregate(post[[item]], by = list(post[[id_var]]),
                               FUN = function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  names(max_post) <- c(id_var, "max_post_bl")
  out <- merge(out, max_post, by = id_var, all.x = TRUE)

  base <- dsn[dsn[[cycle_var]] == baseline_val & !is.na(dsn[[cycle_var]]),
              c(id_var, item), drop = FALSE]
  names(base)[2] <- "base_val"
  out <- merge(out, base, by = id_var, all.x = TRUE)
  out$max_val <- pmax(out$base_val, out$max_post_bl, na.rm = TRUE)
  out$max_val[is.infinite(out$max_val)] <- NA_real_
  out$bl_adjusted <- ifelse(is.na(out$max_post_bl), NA_real_,
                            ifelse(!is.na(out$base_val) & out$base_val >= out$max_post_bl,
                                   0, out$max_post_bl))
  out$score <- switch(tolower(type),
                      bl_adjusted = out$bl_adjusted,
                      max_post_bl = out$max_post_bl,
                      max = out$max_val,
                      stop("param type must be one of 'bl_adjusted', 'max_post_bl', or 'max'"))
  out
}

#' @export
toxTables <- function(dsn, id_var, cycle_var, baseline_val,
                      type = "bl_adjusted", test = "c", riskdiff = FALSE,
                      risk_ci = "wald", risk_ci_alpha = 0.05,
                      arm_var = NA, cycle_limit = NA) {
  ind_items <- .proae_indicator_items(dsn)
  if (!length(ind_items)) {
    return(.toxTables_score_items(dsn, id_var, cycle_var, baseline_val, type,
                                  test, riskdiff, risk_ci, risk_ci_alpha,
                                  arm_var, cycle_limit))
  }

  if (!is.na(cycle_limit) && is.numeric(cycle_limit)) {
    dsn <- dsn[dsn[[cycle_var]] <= cycle_limit, , drop = FALSE]
  }

  required <- c(id_var, cycle_var, if (!is.na(arm_var)) arm_var)
  score_dsn <- .proae_non_indicator_data(dsn, ind_items, required)
  score_vars <- toupper(names(score_dsn)) %in% c(as.character(proctcae_vars$name),
    paste0(substr(as.character(proctcae_vars$name), 1,
                  nchar(as.character(proctcae_vars$name)) - 5), "_COMP"))
  score_out <- if (any(score_vars)) {
    .toxTables_score_items(score_dsn, id_var, cycle_var, baseline_val, type,
                           test, riskdiff, risk_ci, risk_ci_alpha, arm_var, NA)
  } else {
    list(individual = data.frame(), composite = data.frame())
  }

  rows <- lapply(ind_items, function(item) {
    summary <- .proae_indicator_summary(dsn, id_var, cycle_var, item,
                                        baseline_val, type)
    if (!is.na(arm_var)) {
      arm_map <- unique(dsn[, c(id_var, arm_var), drop = FALSE])
      summary <- merge(summary, arm_map, by = id_var, all.x = TRUE)
    }
    label <- as.character(proctcae_vars$short_label[match(item, as.character(proctcae_vars$name))])

    if (is.na(arm_var)) {
      valid <- !is.na(summary$score)
      n <- sum(valid)
      freq <- sum(summary$score[valid] >= 1)
      pct <- if (n) 100 * round(freq / n, 2) else NA_real_
      data.frame(item_lab = label, overall = n,
                 overall_pres = if (n) paste0(freq, " (", pct, "%)") else "",
                 overall_sev = "", stringsAsFactors = FALSE)
    } else {
      arms <- unique(as.character(dsn[[arm_var]]))
      arms <- arms[!is.na(arms)]
      vals <- list(item_lab = label)
      for (arm in arms) {
        x <- summary$score[as.character(summary[[arm_var]]) == arm]
        n <- sum(!is.na(x)); freq <- sum(x >= 1, na.rm = TRUE)
        pct <- if (n) 100 * round(freq / n, 2) else NA_real_
        key <- gsub(" ", "", arm, fixed = TRUE)
        vals[[paste0(key, "_n")]] <- n
        vals[[paste0(key, "_pres")]] <- if (n) paste0(freq, " (", pct, "%)") else ""
      }
      vals[[if (riskdiff) "rdiff_pres" else "pv_pres"]] <- NA
      for (arm in arms) vals[[paste0(gsub(" ", "", arm, fixed = TRUE), "_sev")]] <- ""
      vals[[if (riskdiff) "rdiff_sev" else "pv_sev"]] <- NA
      as.data.frame(vals, stringsAsFactors = FALSE, check.names = FALSE)
    }
  })

  ind_tab <- do.call(rbind, rows)
  if (!nrow(score_out$individual)) {
    score_out$individual <- ind_tab
  } else {
    all_cols <- union(names(score_out$individual), names(ind_tab))
    for (nm in setdiff(all_cols, names(score_out$individual))) score_out$individual[[nm]] <- NA
    for (nm in setdiff(all_cols, names(ind_tab))) ind_tab[[nm]] <- NA
    score_out$individual <- rbind(score_out$individual[, all_cols, drop = FALSE],
                                  ind_tab[, all_cols, drop = FALSE])
  }
  score_out
}

.proae_indicator_figure <- function(dsn, id_var, cycle_var, item, baseline_val,
                                    arm_var, plot_limit, colors, bar_label,
                                    summary_only, cycles_only, x_label, y_label,
                                    x_lab_angle, x_lab_vjust, x_lab_hjust,
                                    suppress_legend, add_item_title) {
  dat <- dsn[, c(id_var, cycle_var, if (!is.na(arm_var)) arm_var, item), drop = FALSE]
  dat <- dat[!is.na(dat[[item]]) & !is.na(dat[[cycle_var]]), , drop = FALSE]
  if (!is.na(plot_limit)) dat <- dat[dat[[cycle_var]] <= plot_limit, , drop = FALSE]
  arm_name <- if (is.na(arm_var)) "overall_" else arm_var
  if (is.na(arm_var)) dat$overall_ <- "Overall"

  max_post <- stats::aggregate(dat[dat[[cycle_var]] > baseline_val, item],
    by = list(dat[dat[[cycle_var]] > baseline_val, id_var]), FUN = max)
  names(max_post) <- c(id_var, item)
  max_post[[cycle_var]] <- "Maximum*"
  max_post <- merge(max_post, unique(dat[, c(id_var, arm_name), drop = FALSE]), by = id_var)

  base <- dat[dat[[cycle_var]] == baseline_val, c(id_var, item), drop = FALSE]
  names(base)[2] <- "base_val"
  adj <- merge(dat, base, by = id_var, all.x = TRUE)
  adj <- stats::aggregate(adj[adj[[cycle_var]] > baseline_val, item],
    by = list(adj[adj[[cycle_var]] > baseline_val, id_var]), FUN = max)
  names(adj) <- c(id_var, "max_post_bl")
  adj <- merge(adj, base, by = id_var, all.x = TRUE)
  adj[[item]] <- ifelse(!is.na(adj$base_val) & adj$base_val >= adj$max_post_bl, 0, adj$max_post_bl)
  adj[[cycle_var]] <- "Adjusted**"
  adj <- merge(adj[, c(id_var, cycle_var, item), drop = FALSE],
               unique(dat[, c(id_var, arm_name), drop = FALSE]), by = id_var)

  plot_dat <- rbind(dat[, c(id_var, cycle_var, arm_name, item), drop = FALSE],
                    max_post[, c(id_var, cycle_var, arm_name, item), drop = FALSE],
                    adj[, c(id_var, cycle_var, arm_name, item), drop = FALSE])
  if (summary_only) plot_dat <- plot_dat[plot_dat[[cycle_var]] %in% c("Maximum*", "Adjusted**"), ]
  if (cycles_only && !summary_only) plot_dat <- plot_dat[!plot_dat[[cycle_var]] %in% c("Maximum*", "Adjusted**"), ]
  plot_dat$grade <- factor(plot_dat[[item]], levels = c(0, 1), labels = c("0", "1+"))
  plot_dat$cycle_plot <- ifelse(plot_dat[[cycle_var]] %in% c("Maximum*", "Adjusted**"),
                                plot_dat[[cycle_var]], paste0(cycle_var, " ", plot_dat[[cycle_var]]))
  plot_dat$arm <- plot_dat[[arm_name]]

  palette <- if (colors == 3) c("0" = "white", "1+" = "black") else
    if (colors == 2) c("0" = "white", "1+" = "#0072B2") else
      c("0" = "white", "1+" = "#13478C")

  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = arm, fill = grade, color = grade)) +
    ggplot2::geom_bar(position = "fill", width = .8) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_color_manual(values = c("0" = "grey60", "1+" = palette[["1+"]]), drop = FALSE) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(100 * x, "%"), limits = c(0, 1)) +
    ggplot2::facet_grid(. ~ cycle_plot) +
    ggplot2::xlab(if (is.na(arm_var) && x_label == "Randomized Treatment Assignment") "Overall" else x_label) +
    ggplot2::ylab(y_label) +
    ggplot2::theme(panel.background = ggplot2::element_rect(color = "grey", fill = "white"),
                   strip.background = ggplot2::element_rect(colour = "grey", fill = "#ededed"),
                   axis.text.x = ggplot2::element_text(angle = x_lab_angle,
                     vjust = x_lab_vjust, hjust = x_lab_hjust),
                   legend.position = if (suppress_legend) "none" else "right")
  label <- as.character(proctcae_vars$short_label[match(item, as.character(proctcae_vars$name))])
  if (add_item_title) p <- p + ggplot2::ggtitle(label)
  list(label, p)
}

#' @export
toxFigures <- function(dsn, id_var, cycle_var, baseline_val, arm_var = NA,
                       plot_limit = NA, colors = 1, bar_label = 0,
                       cycle_label = FALSE, cycle_vals = NA, cycle_labs = NA,
                       summary_only = FALSE, summary_highlight = FALSE,
                       cycles_only = TRUE, x_lab_angle = 0, x_lab_vjust = 1,
                       x_lab_hjust = 0,
                       x_label = "Randomized Treatment Assignment",
                       y_label = "Percent of Total Frequency",
                       footnote_break = FALSE, suppress_legend = FALSE,
                       add_item_title = FALSE) {
  ind_items <- .proae_indicator_items(dsn)
  if (!length(ind_items)) {
    return(.toxFigures_score_items(dsn, id_var, cycle_var, baseline_val, arm_var,
      plot_limit, colors, bar_label, cycle_label, cycle_vals, cycle_labs,
      summary_only, summary_highlight, cycles_only, x_lab_angle, x_lab_vjust,
      x_lab_hjust, x_label, y_label, footnote_break, suppress_legend,
      add_item_title))
  }

  required <- c(id_var, cycle_var, if (!is.na(arm_var)) arm_var)
  score_dsn <- .proae_non_indicator_data(dsn, ind_items, required)
  score_vars <- toupper(names(score_dsn)) %in% c(as.character(proctcae_vars$name),
    paste0(substr(as.character(proctcae_vars$name), 1,
                  nchar(as.character(proctcae_vars$name)) - 5), "_COMP"))
  score_figs <- if (any(score_vars)) {
    .toxFigures_score_items(score_dsn, id_var, cycle_var, baseline_val, arm_var,
      plot_limit, colors, bar_label, cycle_label, cycle_vals, cycle_labs,
      summary_only, summary_highlight, cycles_only, x_lab_angle, x_lab_vjust,
      x_lab_hjust, x_label, y_label, footnote_break, suppress_legend,
      add_item_title)
  } else list()

  ind_figs <- lapply(ind_items, function(item) .proae_indicator_figure(
    dsn, id_var, cycle_var, item, baseline_val, arm_var, plot_limit, colors,
    bar_label, summary_only, cycles_only, x_label, y_label, x_lab_angle,
    x_lab_vjust, x_lab_hjust, suppress_legend, add_item_title))
  c(score_figs, ind_figs)
}
