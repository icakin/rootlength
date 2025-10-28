#' Run rootlength analysis pipeline
#'
#' Reads a strict 5-column CSV:
#' 1) treatment (character)
#' 2) concentration (numeric)
#' 3) unit (dose unit, e.g., "g/L")
#' 4) root length (numeric)
#' 5) unit (length unit; e.g., "cm","mm","m","in","inch","inches","ft","feet")
#'
#' Converts length units to centimeters, builds groups, and runs either
#' Dunnett vs user-defined control (with HC3 fallback) or Welch ANOVA + Games-Howell.
#' Optionally bootstraps GI change %.
#'
#' @param file Path to the input CSV (strict 5-column schema).
#' @param out_dir Output directory to write CSVs/plots/summary (created if missing).
#' @param mode "control" (Dunnett vs control) or "explore" (Welch + Games-Howell).
#' @param control_name Character. Required if mode = "control".
#'   Use exact label (e.g., "RB (4 g/L)") or pooled treatment (e.g., "RB").
#' @param boot Logical. If TRUE and mode = "control", bootstrap GI change% CIs.
#' @param boot_B Integer; bootstrap iterations.
#' @param seed Integer; random seed for bootstrap.
#' @param device Plot device passed to ggplot2::ggsave ("png","pdf","svg","jpeg","tiff").
#' @param plot_w,plot_h Numeric plot width/height.
#' @param plot_dpi Numeric DPI for raster devices.
#' @param quiet Logical. If TRUE, minimize console messages.
#'
#' @return Invisibly, a list with primary_key, primary_method, results,
#'   games_howell, dunnett_standard, dunnett_hc3.
#'
#' @import dplyr tidyr stringr tibble ggplot2 emmeans multcomp rstatix sandwich broom
#' @importFrom car leveneTest
#' @importFrom utils read.csv
#' @importFrom stats median reorder var
#' @export

run_pipeline <- function(
    file,
    out_dir   = "rootlength_outputs",
    mode      = "control",                 # "control" or "explore"
    control_name = "Control (0 g/L)",      # required if mode="control"
    boot      = TRUE,
    boot_B    = 2000,
    seed      = 123,
    device    = getOption("rootlength.device", "png"),
    plot_w    = getOption("rootlength.plot_w", 9),
    plot_h    = getOption("rootlength.plot_h", 5),
    plot_dpi  = getOption("rootlength.plot_dpi", 300),
    quiet     = TRUE
) {
  # avoid clashing with base::mode()
  mode_chr <- mode

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # ---- I/O guardrails ----
  if (!file.exists(file)) stop(paste("File not found:", file), call. = FALSE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- Read raw (allow duplicate 'unit' headers) ----
  raw <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(raw) < 5) {
    stop("CSV must have 5 columns in this order: treatment, concentration, unit(dose), root length, unit(length).", call. = FALSE)
  }

  # ---- Strict positional mapping ----
  strict_df <- tibble::tibble(
    treatment         = as.character(raw[[1]]),
    concentration     = suppressWarnings(as.numeric(raw[[2]])),
    dose_unit         = as.character(raw[[3]]),
    root_length_value = suppressWarnings(as.numeric(raw[[4]])),
    length_unit       = as.character(raw[[5]])
  )

  if (any(!is.na(strict_df$concentration) & !is.finite(strict_df$concentration)))
    warning("Some 'concentration' values are non-numeric; set to NA.", call. = FALSE)
  if (any(!is.na(strict_df$root_length_value) & !is.finite(strict_df$root_length_value)))
    warning("Some 'root length' values are non-numeric; set to NA.", call. = FALSE)

  # ---- Unit conversion helpers ----
  .normalize_unit <- function(u) tolower(trimws(gsub("\\.", "", as.character(u))))
  length_multiplier_to_cm <- function(u) {
    u <- .normalize_unit(u)
    ifelse(u %in% c("cm", "centimeter", "centimeters"), 1.0,
           ifelse(u %in% c("mm", "millimeter", "millimeters"), 0.1,
                  ifelse(u %in% c("m", "meter", "meters"), 100.0,
                         ifelse(u %in% c("in", "inch", "inches"), 2.54,
                                ifelse(u %in% c("ft", "foot", "feet"), 30.48, NA_real_)))))
  }

  strict_df <- strict_df %>%
    dplyr::mutate(
      len_mult_cm   = length_multiplier_to_cm(length_unit),
      unknown_unit  = is.na(len_mult_cm) & !is.na(length_unit)
    )

  if (any(strict_df$unknown_unit, na.rm = TRUE)) {
    bad <- strict_df %>% dplyr::filter(unknown_unit) %>% dplyr::distinct(length_unit) %>% dplyr::pull(length_unit)
    stop(paste0(
      "Unknown length unit(s): ", paste(bad, collapse = ", "),
      "\nAccepted: cm, mm, m, in, inch, inches, ft, feet."
    ), call. = FALSE)
  }

  # Convert to centimeters
  strict_df <- strict_df %>%
    dplyr::mutate(root_length_cm = root_length_value * len_mult_cm)

  # ---- Label helpers ----
  dose_lab  <- function(conc, unit) ifelse(is.na(conc) | is.na(unit) | unit == "", "-", paste0(conc, " ", unit))
  group_lab <- function(treat, conc, unit) {
    dl <- dose_lab(conc, unit)
    ifelse(dl == "-", treat, paste0(treat, " (", dl, ")"))
  }

  # ---- Optional banner with detected groups (quiet by default) ----
  if (!isTRUE(quiet)) {
    control_options <- strict_df %>%
      dplyr::transmute(group = group_lab(treatment, concentration, dose_unit)) %>%
      dplyr::distinct() %>% dplyr::arrange(group) %>% dplyr::pull(group)
    banner <- paste(
      "================ CONTROL SELECTION ================",
      "Steps:",
      "  1) Review the detected group labels below.",
      "  2) Choose mode & control_name via function arguments.",
      "     - Exact label (recommended): e.g., RB (4 g/L)",
      "     - Pooled treatment: e.g., RB",
      "  Tip: mode = 'explore' ignores control_name.",
      "Detected options:",
      paste(sprintf("  - %s", control_options), collapse = "\n"),
      "==================================================",
      sep = "\n"
    )
    message("\n", banner, "\n")
  }

  # ---- Build tidy data ----
  prepare_rootlength_data <- function(df_strict) {
    df <- df_strict %>%
      dplyr::transmute(
        treatment      = as.character(treatment),
        concentration  = suppressWarnings(as.numeric(concentration)),
        unit           = dplyr::if_else(is.na(dose_unit) | dose_unit == "", NA_character_, as.character(dose_unit)),
        root_length_cm = suppressWarnings(as.numeric(root_length_cm))
      ) %>%
      dplyr::mutate(group = group_lab(.data$treatment, .data$concentration, .data$unit))
    df
  }

  parse_group_like <- function(x) {
    x  <- as.character(x)
    m1 <- stringr::str_match(x, "^(.*?)\\s*\\(([^\\)]+)\\)\\s*$")
    if (!is.na(m1[1,1])) tibble::tibble(treatment = trimws(m1[,2]), dose_str = trimws(m1[,3]))
    else                  tibble::tibble(treatment = trimws(x),     dose_str = NA_character_)
  }

  format_p <- function(p) {
    if (!is.finite(p)) return(NA_character_)
    if (p <= .Machine$double.xmin) return("<1e-16")
    if (p < 0.001) return("<0.001")
    if (p < 0.01)  return(formatC(p, digits = 3, format = "f"))
    if (p < 0.1)   return(formatC(p, digits = 3, format = "f"))
    formatC(p, digits = 3, format = "f")
  }
  signif_from_p <- function(p) {
    if (!is.finite(p)) return("ns")
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    "ns"
  }
  add_p_columns <- function(tab, p_col_candidates = c("p.adj","adj.p.value","p.value","p")) {
    pcol <- p_col_candidates[p_col_candidates %in% names(tab)][1] %||% NA_character_
    if (is.na(pcol)) return(tab)
    dplyr::mutate(tab,
                  p_display      = vapply(.data[[pcol]], format_p, character(1)),
                  `p.adj.signif` = vapply(.data[[pcol]], signif_from_p, character(1))
    )
  }

  compute_vs_control <- function(df, control_rows){
    groups <- df %>%
      dplyr::group_by(.data$treatment, .data$concentration, .data$unit) %>%
      dplyr::summarise(
        n_total       = dplyr::n(),
        n_germinated  = sum(is.finite(.data$root_length_cm) & .data$root_length_cm > 0),
        germ_pct      = 100 * n_germinated / n_total,
        avg_length_cm = mean(.data$root_length_cm, na.rm = TRUE),
        sum_length_cm = sum(.data$root_length_cm, na.rm = TRUE),
        .groups = "drop"
      )
    if (nrow(control_rows) == 0 || all(is.na(control_rows$root_length_cm)))
      stop("Selected control has no valid measurements.", call. = FALSE)

    ctrl_by_unit <- control_rows %>%
      dplyr::group_by(.data$unit) %>%
      dplyr::summarise(
        ctrl_germ_pct = 100 * sum(is.finite(.data$root_length_cm) & .data$root_length_cm > 0) / dplyr::n(),
        ctrl_mean     = mean(.data$root_length_cm, na.rm = TRUE),
        ctrl_sum      = sum(.data$root_length_cm, na.rm = TRUE),
        .groups = "drop"
      )
    ctrl_global <- control_rows %>%
      dplyr::summarise(
        ctrl_germ_pct_g = 100 * sum(is.finite(.data$root_length_cm) & .data$root_length_cm > 0) / dplyr::n(),
        ctrl_mean_g     = mean(.data$root_length_cm, na.rm = TRUE),
        ctrl_sum_g      = sum(.data$root_length_cm, na.rm = TRUE),
        .groups = "drop"
      )

    groups %>%
      dplyr::left_join(ctrl_by_unit, by = "unit") %>%
      dplyr::mutate(
        ctrl_germ_pct = ifelse(is.na(.data$ctrl_germ_pct), ctrl_global$ctrl_germ_pct_g, .data$ctrl_germ_pct),
        ctrl_mean     = ifelse(is.na(.data$ctrl_mean),     ctrl_global$ctrl_mean_g,     .data$ctrl_mean),
        ctrl_sum      = ifelse(is.na(.data$ctrl_sum),      ctrl_global$ctrl_sum_g,      .data$ctrl_sum),

        relative_seed_germination_percent = ifelse(is.finite(.data$ctrl_germ_pct) & .data$ctrl_germ_pct > 0,
                                                   100 * (.data$germ_pct / .data$ctrl_germ_pct), NA_real_),
        relative_growth_percent           = ifelse(is.finite(.data$ctrl_mean) & .data$ctrl_mean > 0,
                                                   100 * (.data$avg_length_cm / .data$ctrl_mean), NA_real_),

        germination_index_percent = pmax(0, (.data$relative_seed_germination_percent *
                                               .data$relative_growth_percent) / 100),

        gi_change_percent     = pmax(.data$germination_index_percent - 100, -100),
        rsg_change_percent    = .data$relative_seed_germination_percent - 100,
        growth_change_percent = .data$relative_growth_percent - 100
      ) %>%
      dplyr::select(
        .data$treatment, .data$concentration, .data$unit,
        .data$n_total, .data$n_germinated, .data$germ_pct,
        .data$avg_length_cm, .data$sum_length_cm,
        .data$relative_seed_germination_percent, .data$rsg_change_percent,
        .data$relative_growth_percent, .data$growth_change_percent,
        .data$germination_index_percent, .data$gi_change_percent
      )
  }

  compute_groups_only <- function(df){
    df %>%
      dplyr::group_by(.data$treatment, .data$concentration, .data$unit) %>%
      dplyr::summarise(
        n_total       = dplyr::n(),
        n_germinated  = sum(is.finite(.data$root_length_cm) & .data$root_length_cm > 0),
        germ_pct      = 100 * n_germinated / n_total,
        avg_length_cm = mean(.data$root_length_cm, na.rm = TRUE),
        sum_length_cm = sum(.data$root_length_cm, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # --- Pooled-control homogeneity check helper ---
  .check_pooled_control <- function(df_len, ctrl) {
    # Only relevant if pooled control
    if (!isTRUE(ctrl$pooled)) return(NULL)

    mk_dose <- function(conc, unit) ifelse(is.na(conc) | is.na(unit) | unit == "", "(no dose)", paste0(conc, " ", unit))

    ctrl_df <- df_len %>%
      dplyr::filter(.data$treatment == ctrl$meta$treatment) %>%
      dplyr::mutate(dose_lbl = mk_dose(.data$concentration, .data$unit)) %>%
      dplyr::filter(is.finite(.data$root_length_cm))

    ndose <- ctrl_df %>% dplyr::distinct(dose_lbl) %>% nrow()
    if (ndose < 2) {
      return(list(ndose = ndose, p = NA_real_, d = NA_real_, note = "Only one dose present; pooling trivial."))
    }

    anova <- tryCatch(rstatix::welch_anova_test(root_length_cm ~ dose_lbl, data = ctrl_df), error = function(e) NULL)
    pval  <- if (!is.null(anova) && "p" %in% names(anova)) anova$p[1] else NA_real_

    dval <- NA_real_
    if (ndose == 2) {
      dtab <- tryCatch(rstatix::cohens_d(root_length_cm ~ dose_lbl, data = ctrl_df, var.equal = FALSE), error = function(e) NULL)
      if (!is.null(dtab) && "effsize" %in% names(dtab)) dval <- dtab$effsize[1]
    }

    note <- if (is.finite(pval) && pval < 0.05) {
      "WARNING: pooled control doses differ (Welch ANOVA p < 0.05). Consider a single-dose control."
    } else if (is.finite(dval) && abs(dval) >= 0.5) {
      "CAUTION: pooled control shows >= medium effect between doses (|d| >= 0.5)."
    } else {
      "OK: no strong evidence of differences among pooled control doses."
    }

    list(ndose = ndose, p = pval, d = dval, note = note)
  }

  # ================== ANALYSIS ==================
  df  <- prepare_rootlength_data(strict_df)

  if (any(df$root_length_cm < 0, na.rm = TRUE))
    warning("Negative root_length_cm detected; kept but likely a data issue.")

  choose_control_rows <- function(df, control_name) {
    cn <- trimws(control_name)
    rows <- df %>% dplyr::filter(.data$group == cn)
    if (nrow(rows) > 0) {
      meta <- rows %>% dplyr::slice_head(n=1) %>% dplyr::select(treatment, concentration, unit) %>% as.list()
      label <- unique(rows$group)[1]
      return(list(rows = rows, label = label, meta = meta, pooled = FALSE))
    }
    parsed <- parse_group_like(cn)
    trt  <- parsed$treatment[1]
    dose <- parsed$dose_str[1]
    if (!nzchar(trt)) stop("Control name couldn't be parsed.", call. = FALSE)
    if (is.na(dose)) {
      rows <- df %>% dplyr::filter(.data$treatment == trt)  # pooled across doses
      label <- trt
      u1 <- NA_character_; conc1 <- NA_real_; pooled <- TRUE
    } else {
      rows <- df %>% dplyr::filter(group == paste0(trt, " (", dose, ")"))
      if (nrow(rows) > 0) {
        u1    <- rows %>% dplyr::filter(!is.na(.data$unit)) %>% dplyr::distinct(unit) %>% dplyr::slice_head(n=1) %>% dplyr::pull(unit)
        conc1 <- rows %>% dplyr::filter(!is.na(.data$concentration)) %>% dplyr::distinct(concentration) %>% dplyr::slice_head(n=1) %>% dplyr::pull(concentration)
      } else {
        u1 <- NA_character_; conc1 <- NA_real_
      }
      label <- paste0(trt, " (", dose, ")"); pooled <- FALSE
    }
    if (nrow(rows) == 0) stop("Control not found. Use one of the detected options.", call. = FALSE)
    list(rows = rows, label = label, meta = list(treatment = trt, concentration = conc1, unit = u1 %||% NA_character_), pooled = pooled)
  }

  if (!mode_chr %in% c("control","explore")) stop("mode must be 'control' or 'explore'.", call. = FALSE)
  ctrl <- NULL
  if (identical(mode_chr, "control")) {
    if (!nzchar(trimws(control_name))) stop("mode='control' requires control_name.", call. = FALSE)
    ctrl <- choose_control_rows(df, control_name)
  }

  # Results table
  if (identical(mode_chr, "control")) {
    res <- compute_vs_control(df, ctrl$rows)
    control_mask <- if (isTRUE(ctrl$pooled)) {
      res$treatment == ctrl$meta$treatment
    } else {
      (res$treatment == ctrl$meta$treatment) &
        ( (is.na(ctrl$meta$concentration) & is.na(res$concentration)) |
            (!is.na(ctrl$meta$concentration) & !is.na(res$concentration) &
               res$concentration == ctrl$meta$concentration) ) &
        ( is.na(ctrl$meta$unit) | (!is.na(res$unit) & res$unit == ctrl$meta$unit) )
    }
    if (any(control_mask, na.rm = TRUE)) {
      res$gi_change_percent[control_mask]     <- 0
      res$rsg_change_percent[control_mask]    <- 0
      res$growth_change_percent[control_mask] <- 0
    }
  } else {
    res <- compute_groups_only(df)
  }

  # Diagnostics to decide robust vs standard
  lev_obj <- tryCatch(
    car::leveneTest(root_length_cm ~ group,
                    data = transform(df, group = factor(group)),
                    center = median),
    error = function(e) e
  )
  lev_p   <- tryCatch({ as.numeric(lev_obj[["Pr(>F)"]][1]) }, error = function(e) NA_real_)
  v_tbl <- df %>%
    dplyr::filter(is.finite(root_length_cm)) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(var = stats::var(root_length_cm), n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(is.finite(var))
  var_ratio <- tryCatch(max(v_tbl$var, na.rm = TRUE) / max(min(v_tbl$var[v_tbl$var>0], na.rm = TRUE), .Machine$double.eps), error = function(e) NA_real_)
  n_ratio   <- tryCatch(max(v_tbl$n,  na.rm = TRUE) / max(min(v_tbl$n[v_tbl$n>0],     na.rm = TRUE), 1L), error = function(e) NA_real_)
  hetero    <- (is.finite(lev_p) && lev_p < 0.05) || (is.finite(var_ratio) && var_ratio > 3)
  imbalance <- is.finite(n_ratio) && n_ratio > 3

  # Root-length statistics
  df_len <- df %>% dplyr::filter(is.finite(root_length_cm))

  # Pooled-control homogeneity check (only if control mode)
  pool_check <- NULL
  if (identical(mode_chr, "control")) {
    pool_check <- .check_pooled_control(df_len, ctrl)
  }

  # Games-Howell
  stats_gh_raw <- tryCatch(
    rstatix::games_howell_test(df_len, root_length_cm ~ group),
    error = function(e) tibble::tibble(error = as.character(e))
  )
  stats_gh <- add_p_columns(stats_gh_raw, c("p.adj","p.value","p"))

  # Dunnett (if control)
  stats_dun_std <- stats_dun_rob <- NULL
  if (identical(mode_chr, "control")) {
    if (isTRUE(ctrl$pooled)) {
      df_len <- df_len %>% dplyr::mutate(group_dun = dplyr::if_else(.data$treatment == ctrl$meta$treatment, ctrl$meta$treatment, .data$group))
      control_group <- ctrl$meta$treatment
    } else {
      df_len <- df_len %>% dplyr::mutate(group_dun = .data$group)
      control_group <- ctrl$label
    }
    df_len$group_dun <- droplevels(factor(df_len$group_dun))
    if (!(control_group %in% levels(df_len$group_dun))) {
      stop("Control group not found in data levels: ", control_group, call. = FALSE)
    }
    df_len$group_dun <- stats::relevel(df_len$group_dun, ref = control_group)

    fit_len <- stats::lm(root_length_cm ~ group_dun, data = df_len)
    em_len  <- emmeans::emmeans(fit_len, ~ group_dun)
    stats_dun_std_raw <- as.data.frame(
      summary(emmeans::contrast(em_len, method = "trt.vs.ctrl", ref = 1, adjust = "holm"))
    )
    stats_dun_std     <- add_p_columns(stats_dun_std_raw, c("p.value","p"))

    rob_vcov <- sandwich::vcovHC(fit_len, type = "HC3")
    glht_out <- multcomp::glht(fit_len, linfct = multcomp::mcp(group_dun = "Dunnett"), vcov = rob_vcov)
    rob_tbl  <- broom::tidy(glht_out)
    if (!"contrast" %in% names(rob_tbl)) if ("term" %in% names(rob_tbl)) rob_tbl$contrast <- rob_tbl$term
    p_col <- if ("adj.p.value" %in% names(rob_tbl)) "adj.p.value" else if ("p.value" %in% names(rob_tbl)) "p.value" else NA_character_
    keep  <- c("contrast","estimate","std.error","statistic", p_col); keep <- keep[keep %in% names(rob_tbl)]
    stats_dun_rob_raw <- rob_tbl[, keep, drop = FALSE]
    if (!is.na(p_col) && p_col != "p.value") names(stats_dun_rob_raw)[names(stats_dun_rob_raw) == p_col] <- "p.value"
    stats_dun_rob     <- add_p_columns(stats_dun_rob_raw, c("p.value"))
  }

  # Primary selection
  if (identical(mode_chr, "control")) {
    if (hetero || imbalance) {
      primary_table  <- stats_dun_rob
      primary_method <- "Dunnett vs control (HC3 robust)"
      primary_key    <- "PRIMARY_rootlength_Dunnett_HC3.csv"
    } else {
      primary_table  <- stats_dun_std
      primary_method <- "Dunnett vs control (standard)"
      primary_key    <- "PRIMARY_rootlength_Dunnett_standard.csv"
    }
  } else {
    primary_table  <- stats_gh
    primary_method <- "Welch ANOVA + Games-Howell (all-pairs)"
    primary_key    <- "PRIMARY_rootlength_Welch_GH.csv"
  }

  # GI bootstrap
  gi_bootstrap <- NULL
  if (isTRUE(boot) && identical(mode_chr, "control")) {
    set.seed(seed)
    summarize_grp <- function(d) tibble::tibble(
      germ_pct = 100 * mean(is.finite(d$root_length_cm) & d$root_length_cm > 0, na.rm = TRUE),
      mean_len = mean(d$root_length_cm, na.rm = TRUE)
    )
    ctrl_stats <- summarize_grp(ctrl$rows)
    if (!is.finite(ctrl_stats$germ_pct) || ctrl_stats$germ_pct <= 0 ||
        !is.finite(ctrl_stats$mean_len) || ctrl_stats$mean_len <= 0) {
      warning("Control has zero/NA germination or mean length; GI bootstrap skipped.")
    } else {
      groups_df  <- df %>% dplyr::distinct(group, treatment, concentration, unit)
      outB <- vector("list", nrow(groups_df))
      for (i in seq_len(nrow(groups_df))) {
        g  <- groups_df$group[i]
        di <- df %>% dplyr::filter(.data$group == g)
        if (nrow(di) < 2) next
        gi_star <- numeric(boot_B)
        for (b in seq_len(boot_B)) {
          di_b <- di[sample.int(nrow(di), replace = TRUE), , drop = FALSE]
          si   <- summarize_grp(di_b)
          GI   <- (si$germ_pct / ctrl_stats$germ_pct) * (si$mean_len / ctrl_stats$mean_len) * 100
          gi_star[b] <- GI - 100
        }
        ci <- stats::quantile(gi_star, c(0.025, 0.5, 0.975), na.rm = TRUE)
        outB[[i]] <- tibble::tibble(group = g,
                                    gi_change_median = unname(ci[2]),
                                    gi_change_lo95   = unname(ci[1]),
                                    gi_change_hi95   = unname(ci[3]))
      }
      gi_bootstrap <- dplyr::bind_rows(outB)
    }
  }

  # ---- Plot helper ----
  ext <- switch(device, png="png", pdf="pdf", svg="svg", jpeg="jpg", tiff="tiff", "png")
  savep <- function(p, fname) ggplot2::ggsave(file.path(out_dir, paste0(fname, ".", ext)), p,
                                              width = plot_w, height = plot_h, dpi = plot_dpi, device = device)

  # Raw plot
  df_plot <- df %>% dplyr::filter(is.finite(.data$root_length_cm))
  p_raw <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$group, y = .data$root_length_cm)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.35, width = 0.45) +
    ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 1) +
    ggplot2::labs(x = "Treatment (dose)", y = "Root length (cm)",
                  title = "Root length by Treatment x Dose") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (nrow(df_plot)) {
    n_tbl <- df_plot %>% dplyr::count(.data$group, name = "n")
    p_raw <- p_raw + ggplot2::annotate("text",
                                       x = seq_len(nrow(n_tbl)), y = -Inf, vjust = -0.6,
                                       label = paste0("n=", n_tbl$n), size = 3
    )
  }
  savep(p_raw, "plot_primary_raw")

  # GI bars (control mode)
  if (identical(mode_chr, "control") && "gi_change_percent" %in% names(res)) {
    res_lbl <- res %>%
      dplyr::mutate(group = {
        dl <- dose_lab(.data$concentration, .data$unit)
        ifelse(dl == "-", .data$treatment, paste0(.data$treatment, " (", dl, ")"))
      })
    p_gi <- ggplot2::ggplot(res_lbl, ggplot2::aes(x = .data$group, y = .data$gi_change_percent)) +
      ggplot2::geom_col(width = 0.65) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::labs(x = "Treatment (dose)", y = "GI change (%) vs control",
                    title = paste0("Germination Index change - control: ", ctrl$label)) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    savep(p_gi, "plot_primary_gi")
  }

  # Dunnett forest (control mode)
  if (identical(mode_chr, "control")) {
    tab <- if (hetero || imbalance) stats_dun_rob else stats_dun_std
    if (!is.null(tab) && nrow(tab)) {
      if (!"contrast" %in% names(tab)) tab$contrast <- tab$contrast %||% tab$group %||% tab$term
      est <- if ("estimate" %in% names(tab)) "estimate" else "emmean"
      se  <- if ("std.error" %in% names(tab)) "std.error" else "SE"
      lo  <- if ("lower.CL" %in% names(tab)) "lower.CL" else NULL
      hi  <- if ("upper.CL" %in% names(tab)) "upper.CL" else NULL
      if (is.null(lo) || is.null(hi)) {
        tab <- dplyr::mutate(tab,
                             lower.CL = .data[[est]] - 1.96*.data[[se]],
                             upper.CL = .data[[est]] + 1.96*.data[[se]]
        )
        lo <- "lower.CL"; hi <- "upper.CL"
      }
      p_dun <- ggplot2::ggplot(tab, ggplot2::aes(y = reorder(.data$contrast, .data[[est]]), x = .data[[est]])) +
        ggplot2::geom_vline(xintercept = 0, linetype = 2) +
        ggplot2::geom_point() +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[lo]], xmax = .data[[hi]]), height = 0.15) +
        ggplot2::labs(x = "Difference vs control (cm)", y = NULL,
                      title = primary_method) +
        ggplot2::theme_classic(base_size = 12)
      savep(p_dun, "plot_primary_length")
    }
  }

  # ---- WRITE CSVs ----
  if (identical(mode_chr, "control") && "gi_change_percent" %in% names(res)) {
    utils::write.csv(res, file.path(out_dir, "rootlength_results_vs_control.csv"), row.names = FALSE)
  } else {
    utils::write.csv(res, file.path(out_dir, "rootlength_group_summaries.csv"), row.names = FALSE)
  }

  if (primary_key == "PRIMARY_rootlength_Welch_GH.csv") {
    utils::write.csv(stats_gh,      file.path(out_dir, primary_key), row.names = FALSE)
    if (!is.null(stats_dun_std)) utils::write.csv(stats_dun_std, file.path(out_dir, "stats_rootlength_dunnett.csv"), row.names = FALSE)
    if (!is.null(stats_dun_rob)) utils::write.csv(stats_dun_rob, file.path(out_dir, "stats_rootlength_dunnett_robustHC3.csv"), row.names = FALSE)
  } else if (primary_key == "PRIMARY_rootlength_Dunnett_standard.csv") {
    utils::write.csv(stats_dun_std, file.path(out_dir, primary_key), row.names = FALSE)
    utils::write.csv(stats_gh,      file.path(out_dir, "stats_rootlength_games_howell.csv"), row.names = FALSE)
    if (!is.null(stats_dun_rob)) utils::write.csv(stats_dun_rob, file.path(out_dir, "stats_rootlength_dunnett_robustHC3.csv"), row.names = FALSE)
  } else {
    utils::write.csv(stats_dun_rob, file.path(out_dir, primary_key), row.names = FALSE)
    utils::write.csv(stats_gh,      file.path(out_dir, "stats_rootlength_games_howell.csv"), row.names = FALSE)
    if (!is.null(stats_dun_std)) utils::write.csv(stats_dun_std, file.path(out_dir, "stats_rootlength_dunnett.csv"), row.names = FALSE)
  }
  if (!is.null(gi_bootstrap))
    utils::write.csv(gi_bootstrap, file.path(out_dir, "primary_gi_change_bootstrap.csv"), row.names = FALSE)

  # ---- selection_summary.txt ----
  lines_summary <- c(
    "rootlength - Selection Summary",
    paste0("Mode: ", mode_chr),
    paste0("User-defined control: ", if (identical(mode_chr, "control")) ctrl$label else "none"),
    paste0("Control pooled across doses: ", if (identical(mode_chr, "control") && isTRUE(ctrl$pooled)) "YES" else "NO"),
    paste0("Primary root-length test: ", primary_method),
    "",
    "Diagnostics:",
    paste0("- Levene (median-centered) p = ", if (is.na(lev_p)) "NA" else formatC(lev_p, format='f', digits=4)),
    paste0("- Variance ratio (max var / min var>0) = ", if (is.na(var_ratio)) "NA" else formatC(var_ratio, format='f', digits=3)),
    paste0("- n-ratio (max n / min n>0) = ", if (is.na(n_ratio)) "NA" else formatC(n_ratio, format='f', digits=3)),
    paste0("- Heteroscedasticity detected: ", if (hetero) "YES" else "NO"),
    paste0("- Group-size imbalance detected: ", if (imbalance) "YES" else "NO")
  )

  # Append pooled-control homogeneity info if applicable
  if (!is.null(pool_check)) {
    lines_summary <- c(
      lines_summary,
      "",
      "Pooled-control homogeneity check:",
      paste0("- # of doses in pooled control: ", pool_check$ndose),
      paste0("- Welch ANOVA p-value: ", if (is.na(pool_check$p)) "NA" else formatC(pool_check$p, format="f", digits=3)),
      paste0("- Cohen's d (if 2 doses): ", if (is.na(pool_check$d)) "NA" else formatC(pool_check$d, format="f", digits=3)),
      paste0("- Verdict: ", pool_check$note)
    )
  }

  lines_summary <- c(
    lines_summary,
    "",
    "Files produced:",
    paste0(" - ", normalizePath(out_dir)),
    paste0("   - ", if (identical(mode_chr, "control") && "gi_change_percent" %in% names(res)) "rootlength_results_vs_control.csv" else "rootlength_group_summaries.csv"),
    paste0("   - ", primary_key),
    "   - stats_rootlength_games_howell.csv",
    if (!is.null(stats_dun_std)) "   - stats_rootlength_dunnett.csv" else NULL,
    if (!is.null(stats_dun_rob)) "   - stats_rootlength_dunnett_robustHC3.csv" else NULL,
    if (!is.null(gi_bootstrap)) "   - primary_gi_change_bootstrap.csv" else NULL,
    paste0("   - plot_primary_raw.", ext),
    if (identical(mode_chr, "control") && "gi_change_percent" %in% names(res)) paste0("   - plot_primary_gi.", ext) else NULL,
    if (identical(mode_chr, "control")) paste0("   - plot_primary_length.", ext) else NULL,
    "",
    "Notes:",
    " - Length units were converted to centimeters internally (cm).",
    " - Dose units are used for labeling only.",
    " - Control is never inferred; only used if you set mode='control' and provide control_name.",
    " - All stats tables include p_display and p.adj.signif (ns, *, **, ***).",
    " - Cutoffs: 0.05 (*), 0.01 (**), 0.001 (***)."
  )
  writeLines(lines_summary, con = file.path(out_dir, "selection_summary.txt"))

  if (!isTRUE(quiet)) {
    message("\nrootlength - analysis complete")
    message("Mode: ", mode_chr)
    if (identical(mode_chr, "control")) {
      message("User-defined control: ", ctrl$label)
      message("Control pooled across doses: ", if (isTRUE(ctrl$pooled)) "YES" else "NO")
    }
    message("Primary test: ", primary_method)
    message("Outputs in: ", normalizePath(out_dir), "\n")
  }

  invisible(list(
    primary_key       = primary_key,
    primary_method    = primary_method,
    results           = res,
    games_howell      = stats_gh,
    dunnett_standard  = stats_dun_std,
    dunnett_hc3       = stats_dun_rob
  ))
}
