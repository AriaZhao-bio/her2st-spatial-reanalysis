# Purpose: Match pathology-region annotations to spatial spots and summarize ERBB2

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

annotation_files <- sort(list.files(
  "data/raw/meta",
  pattern = "^[A-H][1-6]_labeled_coordinates\\.tsv$",
  full.names = TRUE
))
stopifnot(length(annotation_files) == 8)

annotation_list <- vector("list", length(annotation_files))
match_report_list <- vector("list", length(annotation_files))

for (i in seq_along(annotation_files)) {
  annotation_file <- annotation_files[i]
  sample_id <- sub("_labeled_coordinates\\.tsv$", "", basename(annotation_file))
  patient_id <- substr(sample_id, 1, 1)

  object_file <- file.path(
    "data", "processed",
    paste0("HER2ST_patient_", patient_id, "_lognorm_pca.rds")
  )
  object <- readRDS(object_file)
  metadata <- object[[]]
  metadata <- metadata[metadata$sample_id == sample_id, , drop = FALSE]

  pathology <- read.delim(annotation_file, check.names = FALSE)
  valid_coordinate <- is.finite(pathology$x) & is.finite(pathology$y)
  pathology_valid <- pathology[valid_coordinate, , drop = FALSE]
  invalid_pathology_rows <- sum(!valid_coordinate)
  pathology_valid$original_spot_id <- paste(
    round(pathology_valid$x),
    round(pathology_valid$y),
    sep = "x"
  )

  if (anyDuplicated(pathology_valid$original_spot_id)) {
    stop("Duplicated rounded pathology spot IDs in ", sample_id)
  }

  match_index <- match(metadata$original_spot_id, pathology_valid$original_spot_id)
  matched <- !is.na(match_index)

  annotation_list[[i]] <- data.frame(
    unique_spot_id = rownames(metadata),
    patient = patient_id,
    sample_id = sample_id,
    original_spot_id = metadata$original_spot_id,
    pixel_x = metadata$pixel_x,
    pixel_y = metadata$pixel_y,
    pathology_label = pathology_valid$label[match_index],
    stringsAsFactors = FALSE
  )

  match_report_list[[i]] <- data.frame(
    patient = patient_id,
    sample_id = sample_id,
    expression_spots = nrow(metadata),
    pathology_rows = nrow(pathology),
    pathology_rows_with_missing_grid_coordinates = invalid_pathology_rows,
    matched_spots = sum(matched),
    unmatched_expression_spots = sum(!matched),
    unused_valid_pathology_rows = sum(
      !pathology_valid$original_spot_id %in% metadata$original_spot_id
    )
  )

  rm(object, metadata, pathology, pathology_valid)
  gc()
}

pathology_annotations <- do.call(rbind, annotation_list)
match_report <- do.call(rbind, match_report_list)
row.names(pathology_annotations) <- NULL
row.names(match_report) <- NULL

if (any(is.na(pathology_annotations$pathology_label))) {
  stop("At least one expression spot did not receive a pathology label.")
}

erbb2_spots <- read.csv(
  "results/tables/ERBB2_spot_level.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

erbb2_index <- match(pathology_annotations$unique_spot_id, erbb2_spots$unique_spot_id)
if (any(is.na(erbb2_index))) {
  stop("At least one annotated spot was absent from ERBB2_spot_level.csv.")
}

annotated_erbb2 <- pathology_annotations
annotated_erbb2$ERBB2_normalized <- erbb2_spots$ERBB2_normalized[erbb2_index]
annotated_erbb2$ERBB2_raw_counts <- erbb2_spots$ERBB2_raw_counts[erbb2_index]

summary_groups <- split(
  annotated_erbb2,
  interaction(
    annotated_erbb2$sample_id,
    annotated_erbb2$pathology_label,
    drop = TRUE
  )
)

pathology_summary <- do.call(rbind, lapply(summary_groups, function(d) {
  data.frame(
    patient = d$patient[1],
    sample_id = d$sample_id[1],
    pathology_label = d$pathology_label[1],
    number_of_spots = nrow(d),
    median_ERBB2_normalized = median(d$ERBB2_normalized),
    mean_ERBB2_normalized = mean(d$ERBB2_normalized),
    fraction_ERBB2_detected = mean(d$ERBB2_raw_counts > 0),
    stringsAsFactors = FALSE
  )
}))
row.names(pathology_summary) <- NULL
pathology_summary <- pathology_summary[
  order(pathology_summary$patient, pathology_summary$pathology_label),
]

write.csv(
  pathology_annotations,
  "results/tables/pathology_spot_annotations.csv",
  row.names = FALSE
)
write.csv(
  match_report,
  "results/tables/pathology_annotation_match_report.csv",
  row.names = FALSE
)
write.csv(
  pathology_summary,
  "results/tables/ERBB2_by_pathology_region.csv",
  row.names = FALSE
)

# Patient-level paired comparison. Spots are first summarized within each
# patient, so the independent unit for inference is the patient (n = 8), not
# the individual spot. "Nonmalignant" excludes cancer in situ and
# undetermined regions.
nonmalignant_labels <- c(
  "adipose tissue",
  "breast glands",
  "connective tissue",
  "immune infiltrate"
)

patient_comparison <- do.call(rbind, lapply(
  split(annotated_erbb2, annotated_erbb2$patient),
  function(d) {
    invasive_values <- d$ERBB2_normalized[
      d$pathology_label == "invasive cancer"
    ]
    nonmalignant_values <- d$ERBB2_normalized[
      d$pathology_label %in% nonmalignant_labels
    ]
    data.frame(
      patient = d$patient[1],
      invasive_spots = length(invasive_values),
      nonmalignant_spots = length(nonmalignant_values),
      invasive_median = median(invasive_values),
      nonmalignant_median = median(nonmalignant_values),
      delta_invasive_minus_nonmalignant =
        median(invasive_values) - median(nonmalignant_values)
    )
  }
))
row.names(patient_comparison) <- NULL

paired_test <- wilcox.test(
  patient_comparison$delta_invasive_minus_nonmalignant,
  mu = 0,
  paired = FALSE,
  exact = TRUE,
  conf.int = TRUE,
  conf.level = 0.95
)

paired_test_summary <- data.frame(
  comparison = "Invasive cancer vs pooled nonmalignant regions",
  independent_unit = "patient",
  number_of_patients = nrow(patient_comparison),
  patients_with_positive_delta = sum(
    patient_comparison$delta_invasive_minus_nonmalignant > 0
  ),
  median_patient_delta = median(
    patient_comparison$delta_invasive_minus_nonmalignant
  ),
  hodges_lehmann_pseudomedian = unname(paired_test$estimate),
  confidence_level = 0.95,
  confidence_interval_low = paired_test$conf.int[1],
  confidence_interval_high = paired_test$conf.int[2],
  test = "Two-sided exact one-sample Wilcoxon signed-rank test on patient deltas",
  statistic_V = unname(paired_test$statistic),
  exact_p_value = paired_test$p.value,
  stringsAsFactors = FALSE
)

write.csv(
  patient_comparison,
  "results/tables/ERBB2_invasive_vs_nonmalignant_patient_values.csv",
  row.names = FALSE
)
write.csv(
  paired_test_summary,
  "results/tables/ERBB2_invasive_vs_nonmalignant_test.csv",
  row.names = FALSE
)

print(match_report)
print(pathology_summary)
print(patient_comparison)
print(paired_test_summary)

label_counts <- sort(table(annotated_erbb2$pathology_label), decreasing = TRUE)
label_levels <- names(label_counts)
annotated_erbb2$pathology_label <- factor(
  annotated_erbb2$pathology_label,
  levels = label_levels
)

pathology_palette <- setNames(
  grDevices::hcl.colors(length(label_levels), palette = "Dark 3"),
  label_levels
)

p_map <- ggplot(
  annotated_erbb2,
  aes(x = pixel_x, y = pixel_y, colour = pathology_label)
) +
  geom_point(size = 0.55, alpha = 0.95) +
  scale_y_reverse() +
  scale_colour_manual(values = pathology_palette, name = "Pathology region") +
  coord_fixed() +
  facet_wrap(~sample_id, ncol = 4) +
  theme_void(base_size = 8, base_family = "sans") +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 7, hjust = 0),
    strip.text = element_text(size = 7, face = "bold"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.2),
    panel.spacing = grid::unit(1.2, "mm"),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  ) +
  labs(
    title = "Pathology regions in the eight annotated HER2ST sections",
    subtitle = "One pathology-annotated section per patient; colours show the published region labels"
  )

ggsave(
  "results/figures/pathology_regions_8_sections.pdf",
  p_map,
  width = 183 / 25.4,
  height = 105 / 25.4,
  device = cairo_pdf
)
ggsave(
  "results/figures/pathology_regions_8_sections.png",
  p_map,
  width = 183 / 25.4,
  height = 105 / 25.4,
  dpi = 400
)

# Spot distributions are shown descriptively. No spot-level p-values are used
# because spots from the same patient/section are not biological replicates.
p_box <- ggplot(
  annotated_erbb2,
  aes(x = pathology_label, y = ERBB2_normalized, fill = pathology_label)
) +
  geom_boxplot(
    width = 0.68,
    outlier.shape = NA,
    linewidth = 0.28
  ) +
  facet_wrap(~sample_id, ncol = 4, scales = "free_x") +
  scale_fill_manual(values = pathology_palette, guide = "none") +
  labs(
    title = "ERBB2 expression across published pathology regions",
    subtitle = "Descriptive spot distributions; each patient contributes one annotated section",
    x = NULL,
    y = "Normalized ERBB2"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(margin = margin(r = 3, unit = "mm")),
    plot.margin = margin(3, 3, 3, 7, unit = "mm")
  )

ggsave(
  "results/figures/ERBB2_by_pathology_region.pdf",
  p_box,
  width = 183 / 25.4,
  height = 120 / 25.4,
  device = cairo_pdf
)
ggsave(
  "results/figures/ERBB2_by_pathology_region.png",
  p_box,
  width = 183 / 25.4,
  height = 120 / 25.4,
  dpi = 400
)

patient_long <- rbind(
  data.frame(
    patient = patient_comparison$patient,
    region = "Nonmalignant",
    median_ERBB2 = patient_comparison$nonmalignant_median
  ),
  data.frame(
    patient = patient_comparison$patient,
    region = "Invasive cancer",
    median_ERBB2 = patient_comparison$invasive_median
  )
)
patient_long$region <- factor(
  patient_long$region,
  levels = c("Nonmalignant", "Invasive cancer")
)

p_paired <- ggplot(
  patient_long,
  aes(x = region, y = median_ERBB2, group = patient)
) +
  geom_line(colour = "#7A7A7A", linewidth = 0.45) +
  geom_point(aes(colour = patient), size = 2.0) +
  scale_colour_brewer(palette = "Dark2", name = "Patient") +
  labs(
    title = "ERBB2 is consistently higher in invasive cancer regions",
    subtitle = paste0(
      "Patient-level medians; two-sided exact Wilcoxon signed-rank test\n",
      "n = 8 patients; exact p = ",
      formatC(paired_test$p.value, format = "f", digits = 4)
    ),
    x = NULL,
    y = "Median normalized ERBB2",
    caption = paste0(
      "Nonmalignant pools adipose tissue, breast glands, connective tissue\n",
      "and immune infiltrate. Patient is the independent unit."
    )
  ) +
  theme(
    legend.position = "right",
    plot.caption = element_text(size = 6, hjust = 0),
    plot.margin = margin(3, 3, 3, 3, unit = "mm")
  )

ggsave(
  "results/figures/ERBB2_invasive_vs_nonmalignant_patient_paired.pdf",
  p_paired,
  width = 120 / 25.4,
  height = 88 / 25.4,
  device = cairo_pdf
)
ggsave(
  "results/figures/ERBB2_invasive_vs_nonmalignant_patient_paired.png",
  p_paired,
  width = 120 / 25.4,
  height = 88 / 25.4,
  dpi = 400
)

cat("\nPathology annotation and ERBB2 regional summaries completed.\n")
