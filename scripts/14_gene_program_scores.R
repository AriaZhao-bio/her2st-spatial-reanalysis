# Purpose: Score four biologically relevant spatial gene programs and compare
# invasive-cancer with nonmalignant regions using patient-level replication.

library(Seurat)
library(ggplot2)

theme_set(
  theme_classic(base_size = 8, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      plot.title = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 7),
      strip.text = element_text(size = 7, face = "bold"),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6.5),
      panel.grid = element_blank()
    )
)

read_single_gmt <- function(path) {
  line <- readLines(path, warn = FALSE)
  stopifnot(length(line) == 1)
  fields <- strsplit(line, "\t", fixed = FALSE)[[1]]
  stopifnot(length(fields) >= 3)
  list(
    set_name = fields[1],
    source = fields[2],
    genes = unique(fields[-c(1, 2)])
  )
}

gene_set_files <- c(
  HER2_17q12 = "data/raw/gene_sets/HER2_17Q12_CORE.gmt",
  Proliferation_E2F = "data/raw/gene_sets/HALLMARK_E2F_TARGETS.gmt",
  Hypoxia = "data/raw/gene_sets/HALLMARK_HYPOXIA.gmt",
  EMT_stromal = paste0(
    "data/raw/gene_sets/",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.gmt"
  )
)
stopifnot(all(file.exists(gene_set_files)))

gene_sets <- lapply(gene_set_files, read_single_gmt)
program_labels <- c(
  HER2_17q12 = "HER2/17q12 amplicon",
  Proliferation_E2F = "Proliferation (E2F)",
  Hypoxia = "Hypoxia",
  EMT_stromal = "EMT/stromal"
)

gene_set_manifest <- data.frame(
  program = names(gene_sets),
  display_label = unname(program_labels[names(gene_sets)]),
  gene_set_name = vapply(gene_sets, `[[`, character(1), "set_name"),
  source = vapply(gene_sets, `[[`, character(1), "source"),
  number_of_requested_genes = vapply(
    gene_sets,
    function(x) length(x$genes),
    integer(1)
  ),
  file = unname(gene_set_files),
  md5 = unname(tools::md5sum(gene_set_files)),
  stringsAsFactors = FALSE
)
write.csv(
  gene_set_manifest,
  "results/tables/gene_program_gene_set_manifest.csv",
  row.names = FALSE
)

input_files <- sort(list.files(
  "data/processed",
  pattern = "HER2ST_patient_[A-H]_lognorm_pca\\.rds$",
  full.names = TRUE
))
stopifnot(length(input_files) == 8)

score_patient_program <- function(object, requested_genes) {
  available_genes <- intersect(requested_genes, rownames(object))
  if (length(available_genes) < 5) {
    stop("Fewer than five requested genes are available for a program.")
  }

  expression_matrix <- as.matrix(LayerData(
    object,
    assay = "RNA",
    layer = "data",
    features = available_genes
  ))
  gene_sds <- apply(expression_matrix, 1, stats::sd)
  variable_gene <- is.finite(gene_sds) & gene_sds > 0
  expression_matrix <- expression_matrix[variable_gene, , drop = FALSE]
  gene_sds <- gene_sds[variable_gene]

  if (nrow(expression_matrix) < 5) {
    stop("Fewer than five variable genes remain for a program.")
  }

  gene_means <- rowMeans(expression_matrix)
  z_matrix <- sweep(expression_matrix, 1, gene_means, "-")
  z_matrix <- sweep(z_matrix, 1, gene_sds, "/")

  list(
    score = colMeans(z_matrix),
    available_genes = available_genes,
    variable_genes = rownames(expression_matrix)
  )
}

spot_score_list <- vector("list", length(input_files))
coverage_list <- list()
coverage_index <- 1L

for (i in seq_along(input_files)) {
  object <- readRDS(input_files[i])
  patient_id <- unique(object$patient)
  stopifnot(length(patient_id) == 1)
  metadata <- object[[]]

  patient_scores <- data.frame(
    unique_spot_id = rownames(metadata),
    patient = metadata$patient,
    sample_id = metadata$sample_id,
    section = metadata$section,
    original_spot_id = metadata$original_spot_id,
    pixel_x = metadata$pixel_x,
    pixel_y = metadata$pixel_y,
    stringsAsFactors = FALSE
  )

  for (program in names(gene_sets)) {
    scored <- score_patient_program(object, gene_sets[[program]]$genes)
    patient_scores[[program]] <- scored$score[patient_scores$unique_spot_id]

    coverage_list[[coverage_index]] <- data.frame(
      patient = patient_id,
      program = program,
      requested_genes = length(gene_sets[[program]]$genes),
      available_genes = length(scored$available_genes),
      variable_genes_used = length(scored$variable_genes),
      coverage_fraction = length(scored$variable_genes) /
        length(gene_sets[[program]]$genes),
      genes_used = paste(scored$variable_genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
    coverage_index <- coverage_index + 1L
  }

  spot_score_list[[i]] <- patient_scores
  rm(object, metadata, patient_scores)
  gc()
}

spot_scores <- do.call(rbind, spot_score_list)
row.names(spot_scores) <- NULL
coverage <- do.call(rbind, coverage_list)
row.names(coverage) <- NULL

write.csv(
  spot_scores,
  "results/tables/gene_program_spot_scores.csv",
  row.names = FALSE
)
write.csv(
  coverage,
  "results/tables/gene_program_gene_coverage.csv",
  row.names = FALSE
)

pathology <- read.csv(
  "results/tables/pathology_spot_annotations.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
score_index <- match(pathology$unique_spot_id, spot_scores$unique_spot_id)
stopifnot(!any(is.na(score_index)))

annotated_scores <- pathology
for (program in names(gene_sets)) {
  annotated_scores[[program]] <- spot_scores[[program]][score_index]
}

nonmalignant_labels <- c(
  "adipose tissue",
  "breast glands",
  "connective tissue",
  "immune infiltrate"
)

patient_program_list <- list()
result_index <- 1L

for (patient_id in sort(unique(annotated_scores$patient))) {
  patient_data <- annotated_scores[annotated_scores$patient == patient_id, ]
  for (program in names(gene_sets)) {
    invasive_values <- patient_data[[program]][
      patient_data$pathology_label == "invasive cancer"
    ]
    nonmalignant_values <- patient_data[[program]][
      patient_data$pathology_label %in% nonmalignant_labels
    ]
    patient_program_list[[result_index]] <- data.frame(
      patient = patient_id,
      program = program,
      invasive_spots = length(invasive_values),
      nonmalignant_spots = length(nonmalignant_values),
      invasive_median = median(invasive_values),
      nonmalignant_median = median(nonmalignant_values),
      delta_invasive_minus_nonmalignant =
        median(invasive_values) - median(nonmalignant_values),
      stringsAsFactors = FALSE
    )
    result_index <- result_index + 1L
  }
}

patient_program_values <- do.call(rbind, patient_program_list)
row.names(patient_program_values) <- NULL

test_list <- lapply(names(gene_sets), function(program) {
  values <- patient_program_values[
    patient_program_values$program == program,
    "delta_invasive_minus_nonmalignant"
  ]
  test <- wilcox.test(
    values,
    mu = 0,
    exact = TRUE,
    conf.int = TRUE,
    conf.level = 0.95
  )
  data.frame(
    program = program,
    display_label = unname(program_labels[program]),
    independent_unit = "patient",
    number_of_patients = length(values),
    patients_with_positive_delta = sum(values > 0),
    median_patient_delta = median(values),
    hodges_lehmann_pseudomedian = unname(test$estimate),
    confidence_interval_low = test$conf.int[1],
    confidence_interval_high = test$conf.int[2],
    statistic_V = unname(test$statistic),
    exact_p_value = test$p.value,
    stringsAsFactors = FALSE
  )
})
program_tests <- do.call(rbind, test_list)
program_tests$BH_FDR <- p.adjust(program_tests$exact_p_value, method = "BH")

write.csv(
  patient_program_values,
  "results/tables/gene_program_patient_pathology_values.csv",
  row.names = FALSE
)
write.csv(
  program_tests,
  "results/tables/gene_program_patient_level_tests.csv",
  row.names = FALSE
)

print(coverage[, 1:6])
print(patient_program_values)
print(program_tests)

# Patient-level delta plot: each point is one independent patient.
delta_plot_data <- patient_program_values
delta_plot_data$display_label <- factor(
  program_labels[delta_plot_data$program],
  levels = unname(program_labels)
)
test_labels <- program_tests
test_labels$display_label <- factor(
  test_labels$display_label,
  levels = unname(program_labels)
)
test_labels$label <- paste0(
  "positive: ", test_labels$patients_with_positive_delta, "/8\n",
  "FDR = ", formatC(test_labels$BH_FDR, format = "f", digits = 3)
)
label_y <- tapply(
  delta_plot_data$delta_invasive_minus_nonmalignant,
  delta_plot_data$display_label,
  max
)
test_labels$label_y <- unname(label_y[as.character(test_labels$display_label)]) + 0.12

p_delta <- ggplot(
  delta_plot_data,
  aes(
    x = display_label,
    y = delta_invasive_minus_nonmalignant,
    colour = patient
  )
) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "grey45") +
  geom_point(size = 2.0, position = position_jitter(width = 0.07, height = 0)) +
  stat_summary(
    aes(group = display_label),
    fun = median,
    geom = "crossbar",
    width = 0.55,
    linewidth = 0.45,
    colour = "black"
  ) +
  geom_text(
    data = test_labels,
    aes(x = display_label, y = label_y, label = label),
    inherit.aes = FALSE,
    size = 2.1,
    vjust = 0
  ) +
  scale_colour_brewer(palette = "Dark2", name = "Patient") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
  labs(
    title = "Gene-program shifts in invasive cancer regions",
    subtitle = paste0(
      "Patient-level differences: invasive-cancer median minus nonmalignant median; ",
      "n = 8 patients"
    ),
    x = NULL,
    y = "Within-patient score difference",
    caption = paste0(
      "Scores are means of patient-standardized gene expression. ",
      "Two-sided exact Wilcoxon signed-rank tests; BH correction across four programs.\n",
      "EMT/stromal scores can reflect fibroblast-rich tissue and are not tumor-cell-specific."
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.caption = element_text(size = 6, hjust = 0),
    plot.margin = margin(3, 3, 3, 3, unit = "mm")
  )

ggsave(
  "results/figures/gene_program_invasive_vs_nonmalignant.pdf",
  p_delta,
  width = 145 / 25.4,
  height = 95 / 25.4,
  device = cairo_pdf
)
ggsave(
  "results/figures/gene_program_invasive_vs_nonmalignant.png",
  p_delta,
  width = 145 / 25.4,
  height = 95 / 25.4,
  dpi = 400
)

# Spatial maps for the eight pathology-annotated sections. Values are capped
# only for display; the uncapped scores are retained in the output tables.
spatial_long <- do.call(rbind, lapply(names(gene_sets), function(program) {
  data.frame(
    patient = annotated_scores$patient,
    sample_id = annotated_scores$sample_id,
    pixel_x = annotated_scores$pixel_x,
    pixel_y = annotated_scores$pixel_y,
    program = program,
    display_label = unname(program_labels[program]),
    score = annotated_scores[[program]],
    stringsAsFactors = FALSE
  )
}))
spatial_long$display_label <- factor(
  spatial_long$display_label,
  levels = unname(program_labels)
)
spatial_long$sample_id <- factor(
  spatial_long$sample_id,
  levels = c("A1", "B1", "C1", "D1", "E1", "F1", "G2", "H1")
)
display_limit <- as.numeric(quantile(abs(spatial_long$score), 0.99, na.rm = TRUE))
spatial_long$display_score <- pmax(
  -display_limit,
  pmin(spatial_long$score, display_limit)
)

p_spatial <- ggplot(
  spatial_long,
  aes(x = pixel_x, y = pixel_y, colour = display_score)
) +
  geom_point(size = 0.35, alpha = 0.95) +
  scale_y_reverse() +
  scale_colour_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-display_limit, display_limit),
    name = "Program\nscore"
  ) +
  coord_fixed() +
  facet_grid(display_label ~ sample_id) +
  theme_void(base_size = 7, base_family = "sans") +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 7, hjust = 0),
    strip.text = element_text(size = 6.2, face = "bold"),
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6),
    panel.spacing = grid::unit(0.8, "mm"),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  ) +
  labs(
    title = "Spatial gene-program scores in pathology-annotated sections",
    subtitle = "Patient-standardized scores; colour limits capped at the overall 99th percentile for display"
  )

ggsave(
  "results/figures/gene_program_spatial_annotated_sections.pdf",
  p_spatial,
  width = 183 / 25.4,
  height = 105 / 25.4,
  device = cairo_pdf
)
ggsave(
  "results/figures/gene_program_spatial_annotated_sections.png",
  p_spatial,
  width = 183 / 25.4,
  height = 105 / 25.4,
  dpi = 400
)

cat("\nGene-program analysis completed.\n")
