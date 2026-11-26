#!/usr/bin/env Rscript
# Compute DunbrackLab IPSAE for top-ranked folds and plot ipSAE (Type == max)
# vs immunogenicity, matching the styling of 03_structure_models_comparisons.R.
#
# Models (08_output_fold_primary_only/):
#   boltz1, boltz2  -> pae_*.npz + *_model_0.cif          (cutoffs 10, 10)
#   colabfold_af2   -> *_scores_rank_001_*.json + *_unrelaxed_rank_001_*.pdb
#                                                    (cutoffs 15, 15)
#   rf3             -> *_confidences.json + *_model.cif (cutoffs 10, 10);
#                     summary is patched with chain_pair_iptm for the AF3 parser
#
# IPSAE: 00_scripts/IPSAE/ipsae.py
#   https://github.com/DunbrackLab/IPSAE
#
# Outputs under 11_ipsae_plots/:
#   ipsae_scores/<model>/<pair_id>/   inputs (symlinks/copies) + IPSAE txt/pml
#   <model>_plots/<model>_ipSAE_P_vs_NP.pdf
#   <model>_plots/<model>_ipSAE_all_receptors_P_vs_NP.pdf
#
# Existing chain-pair score files are skipped. Metadata filter matches script 03.
#
# Usage (from repo root):
#   Rscript 00_scripts/04_model_comparisons_ipsae.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(stringr)
  library(showtext)
  library(sysfonts)
})

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
ipsae_py <- file.path(repo_root, "00_scripts/IPSAE/ipsae.py")
fold_root <- file.path(repo_root, "08_output_fold_primary_only")
out_root <- file.path(repo_root, "11_ipsae_plots")
scores_root <- file.path(out_root, "ipsae_scores")

if (!file.exists(ipsae_py)) {
  stop("Missing IPSAE script: ", ipsae_py)
}

python <- Sys.which("python3")
if (!nzchar(python)) python <- Sys.which("python")
if (!nzchar(python)) stop("python3/python not found on PATH")

# Model-specific PAE / distance cutoffs (document in methods.md).
models <- list(
  boltz1 = list(
    kind = "boltz",
    dir = file.path(fold_root, "boltz1/boltz_results_boltz/predictions"),
    pae_cut = 10, dist_cut = 10
  ),
  boltz2 = list(
    kind = "boltz",
    dir = file.path(fold_root, "boltz2/boltz_results_boltz/predictions"),
    pae_cut = 10, dist_cut = 10
  ),
  colabfold_af2 = list(
    kind = "af2",
    dir = file.path(fold_root, "colabfold_af2"),
    pae_cut = 15, dist_cut = 15
  ),
  rf3 = list(
    kind = "rf3",
    dir = file.path(fold_root, "rf3"),
    pae_cut = 10, dist_cut = 10
  )
)

cut_tag <- function(pae_cut, dist_cut) {
  sprintf("%02d_%02d", as.integer(pae_cut), as.integer(dist_cut))
}

# ---- Styling (same as 03_structure_models_comparisons.R) --------------------

solved_structure_groups <- c("AtFLS2", "Non-AtFLS2", "MIK2")
unsolved_fill <- "grey94"
group_order <- c(solved_structure_groups,
                 "FLS3", "EFR", "PEPR", "CORE", "RLP23", "RLP42", "PERU", "NUT",
                 "CuRe1", "INR", "LNR")

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
showtext_opts(dpi = 96)

classes <- c("Immunogenic", "Non-immunogenic")
class_colors <- c("Immunogenic" = "#7C9885", "Non-immunogenic" = "#B5B682")
class_ink <- c("Immunogenic" = "black", "Non-immunogenic" = "black")

theme_evo <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = base_font) +
    theme(
      text = element_text(colour = "black", family = base_font),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.text = element_text(colour = "black", size = base_size - 1, family = base_font),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title = element_text(size = base_size, family = base_font),
      plot.title = element_text(face = "bold", size = base_size + 1, family = base_font),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(
        size = base_size - 1, family = base_font, colour = "black",
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

# ---- Metadata (same filter as script 03) ------------------------------------

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

row_id_from_pair <- function(pair_id) {
  as.integer(sub("^r0*([0-9]+)_.*$", "\\1", pair_id))
}

# ---- Discover top-ranked structure + PAE inputs -----------------------------

discover_boltz <- function(pred_dir) {
  dirs <- list.dirs(pred_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^r[0-9]+_", basename(dirs))]
  map_dfr(dirs, function(d) {
    pair_id <- basename(d)
    cif <- file.path(d, paste0(pair_id, "_model_0.cif"))
    pae <- file.path(d, paste0("pae_", pair_id, "_model_0.npz"))
    conf <- file.path(d, paste0("confidence_", pair_id, "_model_0.json"))
    if (!file.exists(cif) || !file.exists(pae)) return(NULL)
    tibble(pair_id = pair_id, structure = cif, pae_file = pae,
           extra = if (file.exists(conf)) conf else NA_character_)
  })
}

discover_af2 <- function(af_dir) {
  scores <- list.files(af_dir, pattern = "_scores_rank_001_.*\\.json$", full.names = TRUE)
  map_dfr(scores, function(j) {
    bn <- basename(j)
    pair_id <- sub("_scores_rank_001_.*$", "", bn)
    # scores_rank_001_...json -> unrelaxed_rank_001_...pdb (same model/seed)
    pdb_name <- sub("_scores_rank_001_", "_unrelaxed_rank_001_", bn)
    pdb_name <- sub("\\.json$", ".pdb", pdb_name)
    pdb <- file.path(af_dir, pdb_name)
    if (!file.exists(pdb)) return(NULL)
    tibble(pair_id = pair_id, structure = pdb, pae_file = j, extra = NA_character_)
  })
}

discover_rf3 <- function(rf3_dir) {
  dirs <- list.dirs(rf3_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^r[0-9]+_", basename(dirs))]
  map_dfr(dirs, function(d) {
    pair_id <- basename(d)
    cif <- file.path(d, paste0(pair_id, "_model.cif"))
    conf <- file.path(d, paste0(pair_id, "_confidences.json"))
    summary <- file.path(d, paste0(pair_id, "_summary_confidences.json"))
    if (!file.exists(cif) || !file.exists(conf)) return(NULL)
    tibble(pair_id = pair_id, structure = cif, pae_file = conf,
           extra = if (file.exists(summary)) summary else NA_character_)
  })
}

discover_pairs <- function(model) {
  switch(model$kind,
         boltz = discover_boltz(model$dir),
         af2 = discover_af2(model$dir),
         rf3 = discover_rf3(model$dir),
         stop("unknown kind: ", model$kind))
}

# ---- Stage inputs + run IPSAE ----------------------------------------------

link_input <- function(from, to) {
  from <- normalizePath(from, mustWork = TRUE)
  if (file.exists(to)) {
    # Replace broken/empty placeholders from prior failed runs.
    if (!isTRUE(file.info(to)$isdir) && file.info(to)$size == 0) {
      file.remove(to)
    } else if (identical(normalizePath(to, mustWork = FALSE), from) ||
               identical(Sys.readlink(to), from)) {
      return(invisible(TRUE))
    } else {
      file.remove(to)
    }
  }
  # Prefer symlink; fall back to copy (e.g. if filesystem disallows links).
  status <- system2("ln", c("-sf", from, to))
  if (!identical(as.integer(status), 0L) || !file.exists(to)) {
    if (!isTRUE(file.copy(from, to, overwrite = TRUE)) || !file.exists(to)) {
      stop("Could not link/copy ", from, " -> ", to)
    }
  }
  invisible(TRUE)
}

# RF3 summary lacks chain_pair_iptm; IPSAE's AF3 path requires it. Inject a
# 2x2 matrix from overall iptm (used only for ipTM_af columns; we plot ipSAE).
rf3_patch_py <- file.path(repo_root, "00_scripts/IPSAE/patch_rf3_summary_for_ipsae.py")
write_rf3_summary_for_ipsae <- function(src_summary, dest_summary) {
  status <- system2(python, c(rf3_patch_py, src_summary, dest_summary))
  if (!identical(as.integer(status), 0L)) {
    stop("Failed to patch RF3 summary for IPSAE: ", src_summary)
  }
}

stage_and_run_ipsae <- function(model_name, model, pair) {
  tag <- cut_tag(model$pae_cut, model$dist_cut)
  work <- file.path(scores_root, model_name, pair$pair_id)
  dir.create(work, recursive = TRUE, showWarnings = FALSE)

  struct_bn <- basename(pair$structure)
  pae_bn <- basename(pair$pae_file)
  score_path <- file.path(work, paste0(tools::file_path_sans_ext(struct_bn), "_", tag, ".txt"))

  if (file.exists(score_path) && file.info(score_path)$size > 0) {
    return(list(status = "skipped", score_path = score_path))
  }

  # Drop empty score stubs left by prior failed IPSAE runs.
  for (stub in c(score_path,
                 sub("\\.txt$", "_byres.txt", score_path),
                 sub("\\.txt$", ".pml", score_path))) {
    if (file.exists(stub) && isTRUE(file.info(stub)$size == 0)) file.remove(stub)
  }

  link_input(pair$structure, file.path(work, struct_bn))
  link_input(pair$pae_file, file.path(work, pae_bn))

  if (model$kind == "boltz" && !is.na(pair$extra) && file.exists(pair$extra)) {
    link_input(pair$extra, file.path(work, basename(pair$extra)))
  }
  if (model$kind == "rf3") {
    dest_sum <- file.path(work, paste0(pair$pair_id, "_summary_confidences.json"))
    if (!is.na(pair$extra) && file.exists(pair$extra)) {
      write_rf3_summary_for_ipsae(normalizePath(pair$extra), dest_sum)
    } else {
      writeLines('{"iptm":0,"chain_pair_iptm":[[0,0],[0,0]]}', dest_sum)
    }
  }

  log_file <- file.path(work, "ipsae.log")
  status <- system2(
    python,
    c(ipsae_py, pae_bn, struct_bn,
      as.character(model$pae_cut), as.character(model$dist_cut)),
    stdout = log_file, stderr = log_file,
    wait = TRUE
  )
  exit <- if (is.null(status)) 0L else as.integer(status)
  if (!identical(exit, 0L)) {
    return(list(status = paste0("error:", exit), score_path = score_path,
                log = if (file.exists(log_file)) paste(readLines(log_file, warn = FALSE), collapse = "\n") else ""))
  }
  if (!file.exists(score_path) || file.info(score_path)$size == 0) {
    return(list(status = "error:no_score_file", score_path = score_path,
                log = if (file.exists(log_file)) paste(readLines(log_file, warn = FALSE), collapse = "\n") else ""))
  }
  list(status = "ran", score_path = score_path)
}

# Parse chain-pair table; keep Type == max for A-B (or first max if labels differ).
read_ipsae_max <- function(score_path) {
  lines <- readLines(score_path, warn = FALSE)
  hdr_i <- which(grepl("^Chn1\\s+Chn2", lines))[1]
  if (is.na(hdr_i)) return(NA_real_)
  tab <- tryCatch(
    utils::read.table(text = lines[hdr_i:length(lines)], header = TRUE,
                      stringsAsFactors = FALSE, fill = TRUE),
    error = function(e) NULL
  )
  if (is.null(tab) || !"Type" %in% names(tab) || !"ipSAE" %in% names(tab)) {
    return(NA_real_)
  }
  max_rows <- tab[tolower(tab$Type) == "max", , drop = FALSE]
  if (nrow(max_rows) == 0) return(NA_real_)
  ab <- max_rows[max_rows$Chn1 == "A" & max_rows$Chn2 == "B", , drop = FALSE]
  row <- if (nrow(ab) > 0) ab[1, ] else max_rows[1, ]
  as.numeric(row$ipSAE)
}

# ---- Plotting ---------------------------------------------------------------

order_groups <- function(groups) {
  present <- unique(groups)
  c(intersect(group_order, present), sort(setdiff(present, group_order)))
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

plot_ipsae <- function(df, classes, cutoff = NA_real_) {
  metric <- "ipSAE"
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

  counts <- df %>% count(Group, Class, x) %>% mutate(y = y_min - 0.06 * y_rng)
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

  tick_lo <- y_min - 0.06 * y_rng
  tick_hi <- if (nrow(stats) > 0) y_max + 0.075 * y_rng else y_max
  tick_pad <- 0.05 * (tick_hi - tick_lo)
  y_breaks <- scales::extended_breaks()(c(tick_lo - tick_pad, tick_hi + tick_pad))
  unsolved_x <- which(!groups_present %in% solved_structure_groups)
  unsolved_bands <- tibble(xmin = unsolved_x - 0.5, xmax = unsolved_x + 0.5)

  set.seed(1)
  p <- ggplot(df, aes(x = x, y = .data[[metric]])) +
    geom_rect(data = unsolved_bands,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = unsolved_fill, colour = NA) +
    geom_hline(yintercept = cutoff, linetype = "dashed", linewidth = 0.4,
               colour = "black", na.rm = TRUE) +
    geom_violin(aes(group = interaction(Group, Class), fill = Class),
                width = slot * 0.95, scale = "width", trim = TRUE, colour = NA) +
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
    labs(x = NULL, y = "ipSAE", fill = NULL,
         caption = paste("Bar = median; numbers below violins = n.",
                         "Wilcoxon rank-sum, BH-adjusted:",
                         "* p<0.05, ns = not significant.",
                         "IPSAE Type = max (DunbrackLab/IPSAE).",
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

plot_ipsae_summary <- function(df, classes) {
  structure_levels <- c("Solved", "Unsolved")
  df <- df %>%
    filter(!is.na(ipSAE), Class %in% classes) %>%
    mutate(Class = factor(Class, levels = classes),
           Structure = factor(ifelse(Group %in% solved_structure_groups,
                                     "Solved", "Unsolved"),
                              levels = structure_levels))
  if (nrow(df) == 0) return(NULL)

  short_labels <- c("Immunogenic" = "P", "Non-immunogenic" = "NP")
  y_min <- min(df$ipSAE, na.rm = TRUE)
  y_max <- max(df$ipSAE, na.rm = TRUE)
  y_rng <- max(y_max - y_min, 1e-6)
  tick <- 0.015 * y_rng
  counts <- df %>% count(Structure, Class) %>% mutate(y = y_min - 0.08 * y_rng)

  stats <- map_dfr(split(df, df$Structure, drop = TRUE), function(g) {
    a <- g$ipSAE[g$Class == classes[1]]
    b <- g$ipSAE[g$Class == classes[2]]
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
  p <- ggplot(df, aes(x = Class, y = ipSAE)) +
    geom_rect(data = unsolved_band, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = unsolved_fill, colour = NA) +
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
    labs(x = NULL, y = "ipSAE") +
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

# ---- Main -------------------------------------------------------------------

# Optional: Rscript ... --limit 20  (smoke-test first N pairs per model)
args <- commandArgs(trailingOnly = TRUE)
limit <- NA_integer_
if ("--limit" %in% args) {
  i <- match("--limit", args)
  if (i < length(args)) limit <- as.integer(args[[i + 1]])
}

dir.create(scores_root, recursive = TRUE, showWarnings = FALSE)
# IPSAE writes next to the structure path stem; run each job from its work dir.
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)

for (model_name in names(models)) {
  model <- models[[model_name]]
  if (!dir.exists(model$dir)) {
    message(sprintf("%s: %s not found, skipping", model_name, model$dir))
    next
  }

  pairs <- discover_pairs(model)
  if (is.null(pairs) || nrow(pairs) == 0) {
    message(sprintf("%s: no top-ranked structures found", model_name))
    next
  }
  if (!is.na(limit)) {
    pairs <- pairs %>% slice_head(n = limit)
  }
  message(sprintf("%s: %d top-ranked pairs (PAE/dist cutoffs %d/%d)",
                  model_name, nrow(pairs), model$pae_cut, model$dist_cut))

  n_ran <- n_skip <- n_err <- 0L
  scores <- vector("list", nrow(pairs))

  for (i in seq_len(nrow(pairs))) {
    pair <- pairs[i, ]
    setwd(old_wd)
    work <- file.path(scores_root, model_name, pair$pair_id)
    dir.create(work, recursive = TRUE, showWarnings = FALSE)
    setwd(work)
    res <- stage_and_run_ipsae(model_name, model, pair)
    if (identical(res$status, "skipped")) {
      n_skip <- n_skip + 1L
    } else if (startsWith(res$status, "error")) {
      n_err <- n_err + 1L
      if (n_err <= 5) {
        message(sprintf("  ERROR %s: %s", pair$pair_id, res$status))
        if (!is.null(res$log) && nzchar(res$log)) {
          message("    ", gsub("\n", "\n    ", res$log))
        }
      }
      scores[[i]] <- NULL
      next
    } else {
      n_ran <- n_ran + 1L
    }
    ipsae <- read_ipsae_max(res$score_path)
    scores[[i]] <- tibble(
      pair_id = pair$pair_id,
      row_id = row_id_from_pair(pair$pair_id),
      ipSAE = ipsae,
      score_file = res$score_path
    )
    if (i %% 100 == 0) {
      message(sprintf("  ... %d/%d (ran=%d skip=%d err=%d)",
                      i, nrow(pairs), n_ran, n_skip, n_err))
    }
  }
  setwd(old_wd)

  message(sprintf("  IPSAE done: ran=%d skipped=%d errors=%d", n_ran, n_skip, n_err))
  conf <- bind_rows(scores)
  if (nrow(conf) == 0) {
    message("  no scores to plot")
    next
  }

  dat <- inner_join(conf, meta, by = "row_id")
  message(sprintf("  %d scored, %d matched to known receptors",
                  nrow(conf), nrow(dat)))

  plot_dir <- file.path(out_root, paste0(model_name, "_plots"))
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  n_groups <- n_distinct(dat$Group)

  p <- plot_ipsae(dat, classes)
  if (!is.null(p)) {
    out <- file.path(plot_dir, sprintf("%s_ipSAE_P_vs_NP.pdf", model_name))
    save_plot_pdf(out, p, width = max(4.5, 0.5 * n_groups + 1.3), height = 3.1)
    message("  wrote ", basename(out))
  }

  p <- plot_ipsae_summary(dat, classes)
  if (!is.null(p)) {
    out <- file.path(plot_dir, sprintf("%s_ipSAE_all_receptors_P_vs_NP.pdf", model_name))
    save_plot_pdf(out, p, width = 1.3, height = 1.7)
    message("  wrote ", basename(out))
  }

  # Persist a tidy score table for downstream use.
  write_tsv(
    dat %>% select(pair_id, row_id, Group, Class, ipSAE, score_file),
    file.path(plot_dir, sprintf("%s_ipSAE_scores.tsv", model_name))
  )
}

message("Finished. Outputs under ", out_root)
