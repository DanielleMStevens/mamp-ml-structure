#!/usr/bin/env Rscript
# Classify each structure model by ipTM >= 0.8 as Immunogenic, else
# Non-immunogenic, and score against observed immunogenicity labels.
# By default, only entries that appear in the mamp-ml test set are scored, and
# only the intersection of row_ids folded by every structure model (same n).
# Optional SeqOnly / mamp-ml classification_report.txt scores are plotted beside
# the structure models (report metrics used as published).
#
# Inputs:
#   input_files_fold/boltz{1,2}/boltz_results_boltz/predictions/*/confidence_*.json
#   input_files_fold/rf3/*/<target>_summary_confidences.json (top-ranked model)
#   SERK_coreceptor.tsv (row_id, Receptor_locus_ID, Ligand_Sequence, Immunogenicity)
#   --test_data CSV with columns: Locus ID/Genbank, Sequence, Known Outcome
#   --seqonly / --mamp-ml classification_report.txt paths (optional)
# Outputs (PDF, written to model_comparision_stats/ next to this script):
#   accuracy.pdf, f1_macro.pdf, f1_weighted.pdf, all_metrics.pdf
#
# Usage:
#   Rscript 16_structural_models/01_structure_model_accuracy.R \
#     --seqonly ../model_results/02_immuno_stratify_esm2_with_receptor/classification_report.txt \
#     --mamp-ml ../model_results/06_immuno_stratify_esm2_bfactor_weighted/classification_report.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(showtext)
  library(sysfonts)
})

# Register Arial for correct PDF glyph metrics (same approach as 08_model_analysis plots)
font_add(
  "Arial",
  regular = "/System/Library/Fonts/Supplemental/Arial.ttf",
  bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
  italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
  bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf"
)
showtext_auto()
# High dpi: these figures are ~1–3" wide; showtext rasterizes glyphs into the PDF.
showtext_opts(dpi = 600)

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  if (dir.exists("16_structural_models")) return(normalizePath("16_structural_models"))
  getwd()
}

resolve_path <- function(path, roots) {
  if (is.null(path) || !nzchar(path)) return(NA_character_)
  if (file.exists(path)) return(normalizePath(path))
  for (root in roots) {
    candidate <- file.path(root, path)
    if (file.exists(candidate)) return(normalizePath(candidate))
  }
  # Common alias: ../model_results/... <-> 07_model_results/...
  alt <- sub("^(\\.\\./)?model_results/", "07_model_results/", path)
  if (!identical(alt, path)) {
    for (root in roots) {
      candidate <- file.path(root, alt)
      if (file.exists(candidate)) return(normalizePath(candidate))
    }
    if (file.exists(alt)) return(normalizePath(alt))
  }
  NA_character_
}

parse_args <- function(args) {
  out <- list(
    test_data = "05_datasets/test_data_with_all_test_immuno_stratify.csv",
    use_all = FALSE,
    seqonly = NULL,
    mamp_ml = NULL
  )
  i <- 1L
  while (i <= length(args)) {
    if (args[[i]] %in% c("--test_data", "-t")) {
      if (i == length(args)) stop("--test_data requires a path")
      out$test_data <- args[[i + 1L]]
      i <- i + 2L
    } else if (args[[i]] %in% c("--seqonly", "--seq-only")) {
      if (i == length(args)) stop("--seqonly requires a path")
      out$seqonly <- args[[i + 1L]]
      i <- i + 2L
    } else if (args[[i]] %in% c("--mamp-ml", "--mamp_ml")) {
      if (i == length(args)) stop("--mamp-ml requires a path")
      out$mamp_ml <- args[[i + 1L]]
      i <- i + 2L
    } else if (args[[i]] == "--all") {
      out$use_all <- TRUE
      i <- i + 1L
    } else if (args[[i]] %in% c("--help", "-h")) {
      message(paste(
        "Usage: Rscript 01_structure_model_accuracy.R [options]",
        "  --test_data PATH   mamp-ml test CSV (default:",
        "                     05_datasets/test_data_with_all_test_immuno_stratify.csv)",
        "  --seqonly PATH     SeqOnly classification_report.txt",
        "  --mamp-ml PATH     mamp-ml classification_report.txt",
        "  --all              score all known receptors (ignore --test_data)",
        sep = "\n"
      ))
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", args[[i]], " (try --help)")
    }
  }
  out
}

base_dir <- get_script_dir()
repo_root <- dirname(base_dir)
args <- parse_args(commandArgs(trailingOnly = TRUE))
path_roots <- c(getwd(), repo_root, base_dir, dirname(repo_root))
tsv_path <- file.path(base_dir, "SERK_coreceptor.tsv")
iptm_cutoff <- 0.8
classes <- c("Immunogenic", "Non-immunogenic")
base_font <- "Arial"

structure_models <- list(
  boltz1 = list(
    dir = file.path(base_dir, "input_files_fold/boltz1/boltz_results_boltz/predictions"),
    reader = "read_boltz"
  ),
  boltz2 = list(
    dir = file.path(base_dir, "input_files_fold/boltz2/boltz_results_boltz/predictions"),
    reader = "read_boltz"
  ),
  rf3 = list(
    dir = file.path(base_dir, "input_files_fold/rf3"),
    reader = "read_rf3"
  )
)

metric_specs <- tribble(
  ~key,          ~label,         ~file,
  "accuracy",    "Accuracy",     "accuracy.pdf",
  "f1_macro",    "F1 macro",     "f1_macro.pdf",
  "f1_weighted", "F1 weighted",  "f1_weighted.pdf"
)

# ML models keep the original green palette; structure models are grey ramps.
model_colors <- c(
  "SeqOnly" = "#7C9885",
  "mamp-ml" = "#B5B682",
  boltz1    = "#B0B0B0",
  boltz2    = "#787878",
  rf3       = "#404040"
)
chance_line <- 0.5

theme_evo <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = base_font) +
    theme(
      text = element_text(colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line(colour = "black", linewidth = 0.25),
      axis.ticks.length = unit(0.1, "cm"),
      axis.text = element_text(colour = "black", size = base_size - 1),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = base_size),
      axis.title.x = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 2, colour = "black"),
      legend.key.size = unit(0.2, "cm"),
      legend.key = element_rect(colour = "black", linewidth = 0.2),
      legend.background = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.x = unit(0.1, "cm"),
      plot.caption = element_text(size = base_size - 2, colour = "black", hjust = 0),
      plot.margin = margin(4, 4, 2, 4)
    )
}

normalize_label <- function(x) {
  case_when(
    tolower(x) %in% c("immunogenic", "weakly immunogenic") ~ "Immunogenic",
    tolower(x) == "non-immunogenic" ~ "Non-immunogenic",
    TRUE ~ NA_character_
  )
}

read_classification_report <- function(path, model_name) {
  resolved <- resolve_path(path, path_roots)
  if (is.na(resolved)) stop("classification_report not found: ", path)
  lines <- readLines(resolved, warn = FALSE)
  acc_line <- grep("^\\s*accuracy\\s+", lines, value = TRUE)
  macro_line <- grep("^\\s*macro avg\\s+", lines, value = TRUE)
  weighted_line <- grep("^\\s*weighted avg\\s+", lines, value = TRUE)
  if (!length(acc_line) || !length(macro_line) || !length(weighted_line)) {
    stop("Could not parse accuracy/macro/weighted lines in: ", resolved)
  }
  # Lines look like: accuracy <acc> <support>
  #                  macro avg <p> <r> <f1> <support>
  nums <- function(line) {
    vals <- suppressWarnings(as.numeric(strsplit(trimws(line), "\\s+")[[1]]))
    vals[!is.na(vals)]
  }
  acc <- nums(acc_line[[1]])
  macro <- nums(macro_line[[1]])
  weighted <- nums(weighted_line[[1]])
  tibble(
    model = model_name,
    accuracy = acc[[1]],
    f1_macro = macro[[3]],
    f1_weighted = weighted[[3]],
    n = as.integer(acc[[2]]),
    source = resolved
  )
}

# ---- Metadata ---------------------------------------------------------------

meta_raw <- read_tsv(tsv_path, col_types = cols(.default = col_character())) %>%
  transmute(
    row_id = as.integer(row_id),
    Receptor,
    locus = Receptor_locus_ID,
    lig_seq = Ligand_Sequence,
    Observed = normalize_label(Immunogenicity)
  )

test_filter_label <- NULL
if (args$use_all) {
  meta <- meta_raw %>%
    filter(!is.na(Receptor), tolower(Receptor) != "unknown", !is.na(Observed)) %>%
    select(row_id, Receptor, Observed)
  message("Scoring all known receptors (n_meta=", nrow(meta), ")")
} else {
  test_path <- resolve_path(args$test_data, path_roots)
  if (is.na(test_path)) {
    stop("Test data not found: ", args$test_data,
         "\nPass --test_data PATH or --all to score every known receptor.")
  }
  test <- read_csv(test_path, col_types = cols(.default = col_character()))
  required <- c("Locus ID/Genbank", "Sequence", "Known Outcome")
  missing <- setdiff(required, names(test))
  if (length(missing)) {
    stop("Test data missing columns: ", paste(missing, collapse = ", "))
  }
  test_keys <- test %>%
    transmute(
      locus = `Locus ID/Genbank`,
      lig_seq = Sequence,
      Observed = normalize_label(`Known Outcome`)
    ) %>%
    filter(!is.na(Observed), !is.na(locus), !is.na(lig_seq))
  meta <- meta_raw %>%
    select(row_id, Receptor, locus, lig_seq) %>%
    inner_join(test_keys, by = c("locus", "lig_seq")) %>%
    select(row_id, Receptor, Observed) %>%
    distinct()
  unmatched <- nrow(anti_join(test_keys, meta_raw, by = c("locus", "lig_seq")))
  message(sprintf(
    "Test filter: %s\n  test rows=%d  matched row_ids=%d  unmatched test rows=%d",
    test_path, nrow(test_keys), n_distinct(meta$row_id), unmatched
  ))
  test_filter_label <- basename(test_path)
}

# ---- Confidence JSONs -------------------------------------------------------

row_id_from_dir <- function(f) {
  as.integer(sub("^r0*([0-9]+)_.*$", "\\1", basename(dirname(f))))
}

read_boltz <- function(pred_dir) {
  files <- list.files(pred_dir, pattern = "^confidence_.*\\.json$",
                      recursive = TRUE, full.names = TRUE)
  map_dfr(files, function(f) {
    j <- fromJSON(f)
    tibble(row_id = row_id_from_dir(f), iptm = as.numeric(j$iptm))
  })
}

read_rf3 <- function(pred_dir) {
  target_dirs <- list.dirs(pred_dir, recursive = FALSE)
  files <- file.path(target_dirs, paste0(basename(target_dirs), "_summary_confidences.json"))
  files <- files[file.exists(files)]
  map_dfr(files, function(f) {
    j <- fromJSON(f)
    tibble(row_id = row_id_from_dir(f), iptm = as.numeric(j$iptm))
  })
}

# ---- Metrics ----------------------------------------------------------------

f1_one_class <- function(truth, pred, positive) {
  tp <- sum(truth == positive & pred == positive)
  fp <- sum(truth != positive & pred == positive)
  fn <- sum(truth == positive & pred != positive)
  precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
  if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
}

score_predictions <- function(truth, pred, class_levels = classes) {
  truth <- factor(truth, levels = class_levels)
  pred <- factor(pred, levels = class_levels)
  keep <- !is.na(truth) & !is.na(pred)
  truth <- truth[keep]
  pred <- pred[keep]
  n <- length(truth)
  if (n == 0) {
    return(tibble(accuracy = NA_real_, f1_macro = NA_real_, f1_weighted = NA_real_, n = 0L))
  }
  f1s <- map_dbl(class_levels, ~ f1_one_class(truth, pred, .x))
  support <- as.numeric(table(truth)[class_levels])
  tibble(
    accuracy = mean(truth == pred),
    f1_macro = mean(f1s),
    f1_weighted = sum(f1s * support) / sum(support),
    n = as.integer(n)
  )
}

# ---- Run --------------------------------------------------------------------

per_model <- list()
for (model_name in names(structure_models)) {
  model <- structure_models[[model_name]]
  if (!dir.exists(model$dir)) {
    message(sprintf("%s: %s not found, skipping", model_name, model$dir))
    next
  }
  conf <- match.fun(model$reader)(model$dir) %>%
    filter(!is.na(iptm)) %>%
    mutate(Predicted = if_else(iptm >= iptm_cutoff, "Immunogenic", "Non-immunogenic"))
  dat <- inner_join(conf, meta, by = "row_id")
  message(sprintf("%s: %d predictions matched evaluation set", model_name, nrow(dat)))
  per_model[[model_name]] <- dat
}

if (length(per_model) == 0) stop("No structure models scored.")

shared_ids <- Reduce(intersect, map(per_model, "row_id"))
if (length(shared_ids) == 0) {
  stop("No overlapping row_ids across structure models; cannot compare on a shared set.")
}
message(sprintf(
  "Shared structure evaluation set: n=%d row_ids present in all of %s",
  length(shared_ids), paste(names(per_model), collapse = ", ")
))

structure_results <- map_dfr(names(per_model), function(model_name) {
  dat <- per_model[[model_name]] %>% filter(row_id %in% shared_ids)
  scores <- score_predictions(dat$Observed, dat$Predicted)
  message(sprintf(
    "%s: n=%d  accuracy=%.3f  f1_macro=%.3f  f1_weighted=%.3f",
    model_name, scores$n, scores$accuracy, scores$f1_macro, scores$f1_weighted
  ))
  mutate(scores, model = model_name, .before = 1)
})

report_results <- list()
if (!is.null(args$seqonly)) {
  report_results[["SeqOnly"]] <- read_classification_report(args$seqonly, "SeqOnly")
  message(sprintf(
    "SeqOnly: n=%d  accuracy=%.3f  f1_macro=%.3f  f1_weighted=%.3f\n  %s",
    report_results[["SeqOnly"]]$n, report_results[["SeqOnly"]]$accuracy,
    report_results[["SeqOnly"]]$f1_macro, report_results[["SeqOnly"]]$f1_weighted,
    report_results[["SeqOnly"]]$source
  ))
}
if (!is.null(args$mamp_ml)) {
  report_results[["mamp-ml"]] <- read_classification_report(args$mamp_ml, "mamp-ml")
  message(sprintf(
    "mamp-ml: n=%d  accuracy=%.3f  f1_macro=%.3f  f1_weighted=%.3f\n  %s",
    report_results[["mamp-ml"]]$n, report_results[["mamp-ml"]]$accuracy,
    report_results[["mamp-ml"]]$f1_macro, report_results[["mamp-ml"]]$f1_weighted,
    report_results[["mamp-ml"]]$source
  ))
}

report_df <- if (length(report_results)) {
  bind_rows(report_results) %>% select(model, accuracy, f1_macro, f1_weighted, n)
} else {
  tibble(model = character(), accuracy = double(), f1_macro = double(),
         f1_weighted = double(), n = integer())
}

model_levels <- c(
  intersect(c("SeqOnly", "mamp-ml"), report_df$model),
  intersect(names(structure_models), structure_results$model)
)

results <- bind_rows(report_df, structure_results) %>%
  mutate(model = factor(model, levels = model_levels)) %>%
  arrange(model)

out_dir <- file.path(base_dir, "model_comparision_stats")
dir.create(out_dir, showWarnings = FALSE)

structure_n <- unique(structure_results$n)
report_n <- unique(report_df$n)
plot_caption <- if (nrow(report_df) > 0) {
  sprintf("ipTM\u2265%.1f; struct n=%d; ML n=%s",
          iptm_cutoff, structure_n[[1]], paste(unique(report_n), collapse = "/"))
} else {
  sprintf("ipTM\u2265%.1f; shared n=%d", iptm_cutoff, structure_n[[1]])
}

fill_values <- model_colors[levels(results$model)]

plot_one_metric <- function(key, label) {
  df <- results %>%
    transmute(model, value = .data[[key]])
  bar_w <- 0.55
  edge_pad <- 0.4
  n <- nrow(df)
  # All bars stacked flush; only outer edge padding.
  df$x <- (seq_len(n) - 0.5) * bar_w
  x_min <- 0 - edge_pad
  x_max <- n * bar_w + edge_pad

  ggplot(df, aes(x = x, y = value, fill = model)) +
    geom_hline(yintercept = chance_line, colour = "grey70", linewidth = 0.35) +
    geom_col(colour = "black", linewidth = 0.2, width = bar_w) +
    geom_text(aes(label = sprintf("%.2f", value), y = value + 0.03),
              angle = 90, hjust = 0, vjust = 0.5,
              size = 6 / .pt, family = base_font, colour = "black") +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_x_continuous(breaks = df$x, labels = rep("", n), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(0, 1.12), clip = "off") +
    labs(y = label, caption = plot_caption) +
    theme_evo(base_size = 8) +
    theme(plot.caption.position = "plot") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
}

n_bars <- nlevels(results$model)
plot_width <- max(1.0, 0.13 * n_bars + 0.4)

pwalk(metric_specs, function(key, label, file) {
  out <- file.path(out_dir, file)
  ggsave(out, plot_one_metric(key, label),
         width = plot_width, height = 1.7, device = "pdf", dpi = 300)
  message("wrote ", out)
})

plot_all_metrics <- function() {
  bar_w <- 0.55
  group_gap <- 0.65
  edge_pad <- 0.4
  n_models <- n_bars
  group_width <- n_models * bar_w

  df <- results %>%
    select(model, Accuracy = accuracy, `F1 macro` = f1_macro, `F1 weighted` = f1_weighted) %>%
    pivot_longer(-model, names_to = "metric", values_to = "value") %>%
    mutate(
      metric = factor(metric, levels = c("Accuracy", "F1 macro", "F1 weighted")),
      model = factor(model, levels = levels(results$model)),
      i = as.integer(model),
      group_id = as.integer(metric),
      x = (group_id - 1) * (group_width + group_gap) + (i - 0.5) * bar_w
    )

  metric_labs <- df %>%
    group_by(metric) %>%
    summarise(x = mean(x), .groups = "drop") %>%
    mutate(label = as.character(metric))

  x_min <- min(df$x) - bar_w / 2 - edge_pad
  x_max <- max(df$x) + bar_w / 2 + edge_pad

  ggplot(df, aes(x = x, y = value, fill = model)) +
    geom_hline(yintercept = chance_line, colour = "grey70", linewidth = 0.35) +
    geom_col(colour = "black", linewidth = 0.2, width = bar_w) +
    geom_text(aes(label = sprintf("%.2f", value), y = value + 0.03),
              angle = 90, hjust = 0, vjust = 0.5,
              size = 5.5 / .pt, family = base_font, colour = "black") +
    geom_text(
      data = metric_labs,
      aes(x = x, y = 1.08, label = label),
      inherit.aes = FALSE,
      size = 7 / .pt, family = base_font, colour = "black", fontface = "plain"
    ) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_x_continuous(breaks = df$x, labels = rep("", nrow(df)), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(0, 1.2), clip = "off") +
    labs(y = "Score", caption = plot_caption) +
    theme_evo(base_size = 8) +
    theme(plot.caption.position = "plot") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
}

combined_out <- file.path(out_dir, "all_metrics.pdf")
ggsave(combined_out, plot_all_metrics(),
       width = max(2.3, 0.13 * n_bars * 3 + 0.9), height = 1.7,
       device = "pdf", dpi = 300)
message("wrote ", combined_out)
