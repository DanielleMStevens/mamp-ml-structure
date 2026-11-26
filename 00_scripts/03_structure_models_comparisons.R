#!/usr/bin/env Rscript
# Compare Boltz and RF3 confidence metrics across receptor groups. pTM is shown
# per receptor group only; all other metrics are split by immunogenicity.
#
# Inputs (paths relative to the mamp-ml-structure repo root):
#   08_output_fold_primary_only/boltz{1,2}/boltz_results_boltz/predictions/*/confidence_*.json
#   08_output_fold_primary_only/rf3/*/<target>_summary_confidences.json (top-ranked model)
#   SERK_coreceptor.tsv (row_id, Receptor, Species_normalized, Immunogenicity)
# Outputs (PDF, written to 09_model_plots/<model>_plots/):
#   <model>_pTM.pdf                  mean pTM per receptor group
#   <model>_<metric>_P_vs_NP.pdf     weakly immunogenic counted as immunogenic
#
# Usage (from repo root or anywhere):
#   Rscript 00_scripts/03_structure_models_comparisons.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(purrr)
  library(readr)
  library(showtext)
  library(sysfonts)
})

# Repo root = parent of 00_scripts/ when run via Rscript --file=...
get_repo_root <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
    if (basename(script_dir) == "00_scripts") return(dirname(script_dir))
    return(script_dir)
  }
  if (file.exists("SERK_coreceptor.tsv")) return(normalizePath("."))
  if (file.exists(file.path("..", "SERK_coreceptor.tsv"))) return(normalizePath(".."))
  getwd()
}

repo_root <- get_repo_root()
tsv_path <- file.path(repo_root, "SERK_coreceptor.tsv")
# RF3 reports no interface pLDDT, and its chain_pair_pae averages over all
# receptor-peptide residue pairs rather than interface contacts only.
models <- list(
  boltz1 = list(
    dir = file.path(repo_root, "08_output_fold_primary_only/boltz1/boltz_results_boltz/predictions"),
    reader = "read_boltz"
  ),
  boltz2 = list(
    dir = file.path(repo_root, "08_output_fold_primary_only/boltz2/boltz_results_boltz/predictions"),
    reader = "read_boltz"
  ),
  rf3 = list(
    dir = file.path(repo_root, "08_output_fold_primary_only/rf3"),
    reader = "read_rf3",
    labels = c(complex_ipae = "Inter-chain PAE (\u00c5)"),
    tags = c(complex_ipae = "interchain_PAE"),
    cutoffs = c(complex_ipae = NA)
  )
)

metrics <- c(
  ptm              = "pTM",
  iptm             = "ipTM",
  complex_plddt    = "Complex pLDDT",
  complex_iplddt   = "Interface pLDDT",
  complex_ipae     = "Interface PAE (\u00c5)"
)
metric_file_tags <- c(
  ptm = "pTM", iptm = "ipTM",
  complex_plddt = "complex_pLDDT", complex_iplddt = "complex_ipLDDT",
  complex_ipae = "complex_iPAE"
)
split_metrics <- setdiff(names(metrics), "ptm")

# High-confidence cutoffs drawn as dashed lines. ipTM and pLDDT follow the
# EMBL-EBI AlphaFold guidance (ipTM > 0.8; pLDDT > 90, i.e. 0.9 on Boltz's 0-1
# scale); interface PAE uses the common < 10 A cutoff for mean interface PAE.
confidence_cutoffs <- c(
  ptm            = 0.8,
  iptm           = 0.8,
  complex_plddt  = 0.9,
  complex_iplddt = 0.9,
  complex_ipae   = 10
)

# Receptors with solved receptor-ligand structures first, then the rest.
solved_structure_groups <- c("AtFLS2", "Non-AtFLS2", "MIK2")
unsolved_fill <- "grey94"
group_order <- c(solved_structure_groups,
                 "FLS3", "EFR", "PEPR", "CORE", "RLP23", "RLP42", "PERU", "NUT",
                 "CuRe1", "INR", "LNR")

# Palette and styling follow Evo 2 (Brixi et al. 2026, Nature, Fig. 2):
# highlighted series in blue, comparison series in light grey, black outlines.
evo_blue   <- "#68AFD7"
evo_grey   <- "#D0D2D3"
evo_points <- "#343434"
# Register a TTF and draw text via showtext so legend / axis glyphs are
# outlined in the PDF (avoids Preview kerning bugs with subsetted ArialMT).
base_font <- "plotfont"
font_candidates <- list(
  c("/System/Library/Fonts/Supplemental/Arial.ttf",
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
    "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf"),
  c("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf")
)
font_registered <- FALSE
for (paths in font_candidates) {
  if (file.exists(paths[[1]])) {
    args <- list(family = base_font, regular = paths[[1]])
    if (file.exists(paths[[2]])) args$bold <- paths[[2]]
    if (file.exists(paths[[3]])) args$italic <- paths[[3]]
    if (file.exists(paths[[4]])) args$bolditalic <- paths[[4]]
    do.call(font_add, args)
    font_registered <- TRUE
    break
  }
}
if (!font_registered) {
  base_font <- "sans"
  message("No Arial/DejaVu TTF found; falling back to device sans")
}
showtext_auto()
# Match ggplot pt sizes under cairo_pdf (96 dpi reference).
showtext_opts(dpi = 96)

classes <- c("Immunogenic", "Non-immunogenic")
class_colors <- c("Immunogenic" = "#7C9885", "Non-immunogenic" = "#B5B682")
# Points and median bars drawn on top of each violin, chosen for contrast.
class_ink <- c("Immunogenic" = "black", "Non-immunogenic" = "black")

theme_evo <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = base_font) +
    theme(
      text = element_text(colour = "black", family = base_font),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.text = element_text(colour = "black", size = base_size - 1, family = base_font),
      # vjust = 1 keeps 45°-rotated labels anchored to the tick.
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title = element_text(size = base_size, family = base_font),
      plot.title = element_text(face = "bold", size = base_size + 1, family = base_font),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(
        size = base_size - 1,
        family = base_font,
        colour = "black",
        margin = margin(l = 2, r = 6, unit = "pt")
      ),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.key.size = unit(0.35, "cm"),
      legend.spacing.x = unit(0.25, "cm"),
      legend.margin = margin(t = 2, b = 0),
      plot.caption = element_text(size = base_size - 2, colour = "black", family = base_font)
    )
}

save_plot_pdf <- function(path, plot, width, height) {
  ggplot2::ggsave(path, plot, width = width, height = height, device = grDevices::cairo_pdf)
}

# ---- Metadata ---------------------------------------------------------------

meta <- read_tsv(tsv_path, col_types = cols(.default = col_character())) %>%
  transmute(
    row_id = as.integer(row_id),
    Receptor,
    Species_normalized,
    Class = case_when(
      tolower(Immunogenicity) %in% c("immunogenic", "weakly immunogenic") ~ "Immunogenic",
      tolower(Immunogenicity) == "non-immunogenic" ~ "Non-immunogenic",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Receptor), tolower(Receptor) != "unknown", !is.na(Class)) %>%
  mutate(Group = case_when(
    Receptor == "FLS2" & Species_normalized == "Arabidopsis thaliana" ~ "AtFLS2",
    Receptor == "FLS2" ~ "Non-AtFLS2",
    Receptor == "INR-like" ~ "INR",
    Receptor == "LNR-like" ~ "LNR",
    TRUE ~ Receptor
  ))

# ---- Confidence JSONs -------------------------------------------------------

row_id_from_dir <- function(f) {
  as.integer(sub("^r0*([0-9]+)_.*$", "\\1", basename(dirname(f))))
}

read_boltz <- function(pred_dir) {
  files <- list.files(pred_dir, pattern = "^confidence_.*\\.json$",
                      recursive = TRUE, full.names = TRUE)
  map_dfr(files, function(f) {
    j <- fromJSON(f)
    tibble(
      row_id = row_id_from_dir(f),
      !!!map(set_names(names(metrics)), ~ as.numeric(j[[.x]]))
    )
  })
}

read_rf3 <- function(pred_dir) {
  target_dirs <- list.dirs(pred_dir, recursive = FALSE)
  files <- file.path(target_dirs, paste0(basename(target_dirs), "_summary_confidences.json"))
  files <- files[file.exists(files)]
  map_dfr(files, function(f) {
    j <- fromJSON(f)
    tibble(
      row_id = row_id_from_dir(f),
      ptm = as.numeric(j$ptm),
      iptm = as.numeric(j$iptm),
      complex_plddt = as.numeric(j$overall_plddt),
      complex_iplddt = NA_real_,
      complex_ipae = as.numeric(j$chain_pair_pae[1, 2])
    )
  })
}

# ---- Plotting ---------------------------------------------------------------

order_groups <- function(groups) {
  present <- unique(groups)
  c(intersect(group_order, present), sort(setdiff(present, group_order)))
}

plot_ptm <- function(df) {
  df <- df %>% filter(!is.na(ptm))
  if (nrow(df) == 0) return(NULL)
  df <- df %>% mutate(Group = factor(Group, levels = order_groups(Group)))

  summary_df <- df %>%
    group_by(Group) %>%
    summarise(mean = mean(ptm), n = n(), .groups = "drop")

  set.seed(1)
  ggplot(df, aes(x = Group, y = ptm)) +
    geom_col(data = summary_df, aes(y = mean), width = 0.75,
             fill = evo_blue, colour = "black", linewidth = 0.3) +
    geom_jitter(width = 0.2, height = 0, size = 0.6, alpha = 0.6,
                colour = evo_points, stroke = 0) +
    geom_hline(yintercept = confidence_cutoffs[["ptm"]], linetype = "dashed", linewidth = 0.4,
               colour = "black") +
    geom_text(data = summary_df, aes(y = 0, label = n), vjust = -0.5,
              size = 2.6, colour = "black", family = base_font) +
    scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "pTM",
         caption = paste("Bar = mean; dot = one model; number = n.",
                         "Dashed line = pTM 0.8 (high confidence).", sep = "\n")) +
    theme_evo() +
    theme(plot.caption.position = "plot")
}

p_to_stars <- function(p) {
  case_when(is.na(p) ~ NA_character_, p < 0.05 ~ "*", TRUE ~ "ns")
}

wilcox_tests <- function(df, metric, classes) {
  pairs <- combn(classes, 2, simplify = FALSE)
  res <- map_dfr(split(df, df$Group), function(g) {
    map_dfr(pairs, function(pr) {
      a <- g[[metric]][g$Class == pr[1]]
      b <- g[[metric]][g$Class == pr[2]]
      a <- a[!is.na(a)]
      b <- b[!is.na(b)]
      if (length(a) < 2 || length(b) < 2) return(NULL)
      tibble(Group = g$Group[1], class1 = pr[1], class2 = pr[2],
             p = wilcox.test(a, b, exact = FALSE)$p.value)
    })
  })
  if (nrow(res) == 0) return(res)
  res %>% mutate(p_adj = p.adjust(p, method = "BH"), label = p_to_stars(p_adj))
}

plot_metric <- function(df, metric, classes, y_label = metrics[[metric]],
                        cutoff = confidence_cutoffs[metric]) {
  df <- df %>%
    filter(!is.na(.data[[metric]]), Class %in% classes) %>%
    mutate(Class = factor(Class, levels = classes))
  if (nrow(df) == 0) return(NULL)

  groups_present <- order_groups(df$Group)
  k <- length(classes)
  dodge <- 0.8
  slot <- dodge / k
  offsets <- set_names((seq_len(k) - (k + 1) / 2) * slot, classes)

  df <- df %>%
    mutate(Group = factor(Group, levels = groups_present),
           x = as.integer(Group) + offsets[as.character(Class)])

  y_min <- min(df[[metric]], cutoff, na.rm = TRUE)
  y_max <- max(df[[metric]], cutoff, na.rm = TRUE)
  y_rng <- max(y_max - y_min, 1e-6)

  counts <- df %>%
    count(Group, Class, x) %>%
    mutate(y = y_min - 0.06 * y_rng)

  stats <- wilcox_tests(df %>% mutate(Class = as.character(Class)), metric, classes)
  if (nrow(stats) > 0) {
    stats <- stats %>%
      mutate(
        gx = match(Group, groups_present),
        xmin = gx + offsets[class1] + 0.08 * slot,
        xmax = gx + offsets[class2] - 0.08 * slot,
        y = y_max + 0.06 * y_rng
      )
  }

  set.seed(1)
  # Ticks are chosen on the tighter range (default 5% padding) so the extra
  # headroom above the brackets doesn't change them.
  tick_lo <- y_min - 0.06 * y_rng
  tick_hi <- if (nrow(stats) > 0) y_max + 0.075 * y_rng else y_max
  tick_pad <- 0.05 * (tick_hi - tick_lo)
  y_breaks <- scales::extended_breaks()(c(tick_lo - tick_pad, tick_hi + tick_pad))

  unsolved_x <- which(!groups_present %in% solved_structure_groups)
  unsolved_bands <- tibble(xmin = unsolved_x - 0.5, xmax = unsolved_x + 0.5)

  p <- ggplot(df, aes(x = x, y = .data[[metric]])) +
    geom_rect(data = unsolved_bands,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = unsolved_fill, colour = NA) +
    geom_hline(yintercept = cutoff, linetype = "dashed", linewidth = 0.4,
               colour = "black", na.rm = TRUE) +
    geom_violin(aes(group = interaction(Group, Class), fill = Class),
                width = slot * 0.95, scale = "width", trim = TRUE,
                colour = NA) +
    geom_jitter(aes(colour = Class), width = slot * 0.18, height = 0, size = 0.6,
                alpha = 0.6, stroke = 0) +
    stat_summary(aes(group = interaction(Group, Class), colour = Class), fun = median,
                 geom = "crossbar", width = slot * 0.5, linewidth = 0.3) +
    geom_text(data = counts, aes(x = x, y = y, label = n),
              size = 8 / .pt, colour = "black", family = base_font) +
    scale_fill_manual(values = class_colors[classes], drop = FALSE) +
    scale_colour_manual(values = class_ink[classes], guide = "none") +
    scale_x_continuous(breaks = seq_along(groups_present), labels = groups_present,
                       expand = expansion(add = 0.5)) +
    scale_y_continuous(breaks = y_breaks, expand = expansion(mult = c(0.05, 0.14))) +
    labs(x = NULL, y = y_label, fill = NULL,
         caption = paste("Bar = median; numbers below violins = n.",
                         "Wilcoxon rank-sum, BH-adjusted:",
                         "* p<0.05, ns = not significant.",
                         sep = if (length(groups_present) < 8) "\n" else " ")) +
    theme_evo(base_size = 10)

  if (nrow(stats) > 0) {
    tick <- 0.015 * y_rng
    p <- p +
      geom_segment(data = stats, aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_segment(data = stats, aes(x = xmin, xend = xmin, y = y, yend = y - tick),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_segment(data = stats, aes(x = xmax, xend = xmax, y = y, yend = y - tick),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_text(data = stats, aes(x = (xmin + xmax) / 2, y = y + tick, label = label),
                inherit.aes = FALSE, size = 8 / .pt, vjust = 0, family = base_font)
  }
  p
}

# Compact panel pooling receptor groups by whether a receptor-ligand structure is
# solved: one violin per immunogenicity class within each structure category.
plot_metric_summary <- function(df, metric, classes) {
  structure_levels <- c("Solved", "Unsolved")
  df <- df %>%
    filter(!is.na(.data[[metric]]), Class %in% classes) %>%
    mutate(Class = factor(Class, levels = classes),
           Structure = factor(ifelse(Group %in% solved_structure_groups,
                                     "Solved", "Unsolved"),
                              levels = structure_levels))
  if (nrow(df) == 0) return(NULL)

  short_labels <- c("Immunogenic" = "P", "Non-immunogenic" = "NP")
  cutoff <- confidence_cutoffs[metric]
  y_min <- min(df[[metric]], cutoff, na.rm = TRUE)
  y_max <- max(df[[metric]], cutoff, na.rm = TRUE)
  y_rng <- max(y_max - y_min, 1e-6)
  tick <- 0.015 * y_rng

  counts <- df %>% count(Structure, Class) %>% mutate(y = y_min - 0.08 * y_rng)

  stats <- map_dfr(split(df, df$Structure, drop = TRUE), function(g) {
    a <- g[[metric]][g$Class == classes[1]]
    b <- g[[metric]][g$Class == classes[2]]
    if (length(a) < 2 || length(b) < 2) return(NULL)
    tibble(Structure = g$Structure[1], p = wilcox.test(a, b, exact = FALSE)$p.value)
  })
  if (nrow(stats) > 0) {
    stats <- stats %>%
      mutate(label = p_to_stars(p.adjust(p, method = "BH")),
             y = y_max + 0.06 * y_rng)
  }

  tick_lo <- y_min - 0.08 * y_rng
  tick_hi <- if (nrow(stats) > 0) y_max + 0.075 * y_rng else y_max
  tick_pad <- 0.05 * (tick_hi - tick_lo)
  y_breaks <- scales::extended_breaks()(c(tick_lo - tick_pad, tick_hi + tick_pad))

  unsolved_band <- tibble(Structure = factor("Unsolved", levels = structure_levels))

  set.seed(1)
  p <- ggplot(df, aes(x = Class, y = .data[[metric]])) +
    geom_rect(data = unsolved_band, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = unsolved_fill, colour = NA) +
    geom_hline(yintercept = cutoff, linetype = "dashed", linewidth = 0.3,
               colour = "black", na.rm = TRUE) +
    geom_violin(aes(fill = Class), width = 0.85, scale = "width", trim = TRUE,
                colour = NA) +
    geom_jitter(aes(colour = Class), width = 0.12, height = 0, size = 0.15,
                alpha = 0.2, stroke = 0) +
    stat_summary(aes(colour = Class), fun = median, geom = "crossbar", width = 0.45,
                 linewidth = 0.25) +
    scale_colour_manual(values = class_ink[classes], guide = "none") +
    geom_text(data = counts, aes(x = Class, y = y, label = n),
              size = 5 / .pt, colour = "black", family = base_font) +
    scale_fill_manual(values = class_colors[classes], guide = "none") +
    scale_x_discrete(labels = short_labels[classes]) +
    scale_y_continuous(breaks = y_breaks, expand = expansion(mult = c(0.05, 0.14))) +
    facet_grid(~ Structure, switch = "x") +
    labs(x = NULL, y = metrics[[metric]]) +
    theme_evo(base_size = 7) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(size = 7, colour = "black",
                                    margin = margin(1, 0, 0, 0)),
          panel.spacing = unit(0, "pt"),
          plot.margin = margin(2, 2, 2, 2))

  if (nrow(stats) > 0) {
    p <- p +
      geom_segment(data = stats, aes(x = 1, xend = 2, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.25) +
      geom_segment(data = stats, aes(x = 1, xend = 1, y = y, yend = y - tick),
                   inherit.aes = FALSE, linewidth = 0.25) +
      geom_segment(data = stats, aes(x = 2, xend = 2, y = y, yend = y - tick),
                   inherit.aes = FALSE, linewidth = 0.25) +
      geom_text(data = stats, aes(x = 1.5, y = y + tick, label = label),
                inherit.aes = FALSE, vjust = 0, size = 6 / .pt, family = base_font)
  }
  p
}

for (model_name in names(models)) {
  model <- models[[model_name]]
  if (!dir.exists(model$dir)) {
    message(sprintf("%s: %s not found, skipping", model_name, model$dir))
    next
  }
  conf <- match.fun(model$reader)(model$dir)
  labels <- replace(metrics, names(model$labels), model$labels)
  tags <- replace(metric_file_tags, names(model$tags), model$tags)
  cutoffs <- replace(confidence_cutoffs, names(model$cutoffs), model$cutoffs)
  dat <- inner_join(conf, meta, by = "row_id")
  message(sprintf("%s: %d confidence files, %d matched to known receptors",
                  model_name, nrow(conf), nrow(dat)))

  out_dir <- file.path(repo_root, "10_model_plots", paste0(model_name, "_plots"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  n_groups <- n_distinct(dat$Group)

  p <- plot_ptm(dat)
  if (!is.null(p)) {
    out <- file.path(out_dir, sprintf("%s_pTM.pdf", model_name))
    save_plot_pdf(out, p, width = max(3, 0.22 * n_groups + 1.2), height = 2.8)
    message("  wrote ", basename(out))
  }

  for (metric in split_metrics) {
    p <- plot_metric(dat, metric, classes, y_label = labels[[metric]],
                     cutoff = cutoffs[metric])
    if (is.null(p)) next
    out <- file.path(out_dir, sprintf("%s_%s_P_vs_NP.pdf", model_name, tags[[metric]]))
    save_plot_pdf(out, p, width = max(4.5, 0.5 * n_groups + 1.3), height = 3.1)
    message("  wrote ", basename(out))
  }

  p <- plot_metric_summary(dat, "iptm", classes)
  if (!is.null(p)) {
    out <- file.path(out_dir, sprintf("%s_ipTM_all_receptors_P_vs_NP.pdf", model_name))
    save_plot_pdf(out, p, width = 1.3, height = 1.7)
    message("  wrote ", basename(out))
  }
}
