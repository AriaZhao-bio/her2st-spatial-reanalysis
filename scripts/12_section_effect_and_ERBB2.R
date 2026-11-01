# Purpose: Quantify section effects and map normalized ERBB2 expression

library(Seurat)
library(ggplot2)

theme_set(
  theme_classic(base_size = 7, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      plot.title = element_text(size = 8, face = "bold"),
      strip.text = element_text(size = 7, face = "bold"),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      panel.grid = element_blank()
    )
)

input_files <- sort(list.files(
  "data/processed",
  pattern = "HER2ST_patient_[A-H]_lognorm_pca\\.rds$",
  full.names = TRUE
))
stopifnot(length(input_files) == 8)

effect_list <- vector("list", length(input_files))
spot_list <- vector("list", length(input_files))

for (i in seq_along(input_files)) {
  object <- readRDS(input_files[i])
  patient_id <- unique(object$patient)
  stopifnot(length(patient_id) == 1, "ERBB2" %in% rownames(object))

  # Descriptive section effect in the first 10 PCs:
  # between-section sum of squares / total sum of squares.
  pc <- Embeddings(object, reduction = "pca")[, 1:10, drop = FALSE]
  section_group <- factor(object$sample_id)
  grand_mean <- colMeans(pc)
  total_ss <- sum(sweep(pc, 2, grand_mean, "-")^2)
  between_ss <- 0
  for (section_id in levels(section_group)) {
    idx <- section_group == section_id
    section_mean <- colMeans(pc[idx, , drop = FALSE])
    between_ss <- between_ss + sum(idx) * sum((section_mean - grand_mean)^2)
  }
  effect_list[[i]] <- data.frame(
    patient = patient_id,
    number_of_sections = nlevels(section_group),
    number_of_spots = ncol(object),
    section_R2_first10PCs = between_ss / total_ss
  )

  normalized_erbb2 <- FetchData(object, vars = "ERBB2", layer = "data")[[1]]
  raw_erbb2 <- FetchData(object, vars = "ERBB2", layer = "counts")[[1]]
  metadata <- object[[]]

  spot_list[[i]] <- data.frame(
    unique_spot_id = rownames(metadata),
    patient = metadata$patient,
    sample_id = metadata$sample_id,
    section = metadata$section,
    original_spot_id = metadata$original_spot_id,
    pixel_x = metadata$pixel_x,
    pixel_y = metadata$pixel_y,
    ERBB2_normalized = normalized_erbb2,
    ERBB2_raw_counts = raw_erbb2,
    stringsAsFactors = FALSE
  )

  rm(object, pc, metadata)
  gc()
}

section_effect <- do.call(rbind, effect_list)
erbb2_spots <- do.call(rbind, spot_list)
row.names(erbb2_spots) <- NULL

section_groups <- split(erbb2_spots, erbb2_spots$sample_id)
section_summary <- do.call(rbind, lapply(section_groups, function(d) {
  data.frame(
    patient = d$patient[1],
    sample_id = d$sample_id[1],
    section = d$section[1],
    number_of_spots = nrow(d),
    median_ERBB2_normalized = median(d$ERBB2_normalized),
    mean_ERBB2_normalized = mean(d$ERBB2_normalized),
    fraction_ERBB2_detected = mean(d$ERBB2_raw_counts > 0)
  )
}))
row.names(section_summary) <- NULL
section_summary <- section_summary[order(section_summary$patient, section_summary$section), ]

write.csv(section_effect, "results/tables/section_effect_first10PCs.csv", row.names = FALSE)
write.csv(erbb2_spots, "results/tables/ERBB2_spot_level.csv", row.names = FALSE)
write.csv(section_summary, "results/tables/ERBB2_section_summary.csv", row.names = FALSE)

print(section_effect)
print(section_summary)

# Technical diagnostic: descriptive section effect only; no p-values because
# spatial spots are not independent biological replicates.
p_effect <- ggplot(section_effect, aes(x = patient, y = section_R2_first10PCs)) +
  geom_col(fill = "#4C78A8", width = 0.72) +
  geom_text(
    aes(label = sprintf("%.2f", section_R2_first10PCs)),
    vjust = -0.35,
    size = 2.2
  ) +
  coord_cartesian(ylim = c(0, max(section_effect$section_R2_first10PCs) * 1.14)) +
  labs(
    title = "Section-associated variation within each patient",
    subtitle = paste(
      "Descriptive section-associated variation across the first 10 PCs;",
      "spots are not independent biological replicates",
      sep = "\n"
    ),
    x = "Patient",
    y = expression(R[section]^2)
  )

ggsave(
  "results/figures/section_effect_first10PCs.pdf",
  p_effect,
  width = 89 / 25.4,
  height = 72 / 25.4,
  device = cairo_pdf
)
ggsave(
  "results/figures/section_effect_first10PCs.png",
  p_effect,
  width = 89 / 25.4,
  height = 72 / 25.4,
  dpi = 400
)

# One spatial figure per patient. The colour cap is patient-specific because
# the inferential role is within-patient spatial heterogeneity, not a direct
# absolute comparison between patients.
for (patient_id in LETTERS[1:8]) {
  d <- erbb2_spots[erbb2_spots$patient == patient_id, ]
  cap <- as.numeric(quantile(d$ERBB2_normalized, 0.99, na.rm = TRUE))
  d$ERBB2_display <- pmin(d$ERBB2_normalized, cap)
  d$sample_id <- factor(d$sample_id, levels = sort(unique(d$sample_id)))

  p <- ggplot(d, aes(x = pixel_x, y = pixel_y, colour = ERBB2_display)) +
    geom_point(size = 0.65, alpha = 0.9) +
    scale_y_reverse() +
    scale_colour_viridis_c(
      option = "magma",
      limits = c(0, cap),
      name = "Normalized\nERBB2"
    ) +
    coord_fixed() +
    facet_wrap(~sample_id, ncol = 3) +
    theme_void(base_size = 7, base_family = "sans") +
    theme(
      plot.title = element_text(size = 8, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6.5, hjust = 0),
      strip.text = element_text(size = 7, face = "bold"),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      panel.spacing = grid::unit(1.2, "mm"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    ) +
    labs(
      title = paste0("Patient ", patient_id, ": spatial ERBB2 expression"),
      subtitle = "Log-normalized expression; colour capped at the patient-specific 99th percentile"
    )

  height_mm <- if (length(unique(d$sample_id)) > 3) 116 else 68
  base_name <- file.path("results", "figures", paste0("ERBB2_spatial_patient_", patient_id))
  ggsave(paste0(base_name, ".pdf"), p, width = 183 / 25.4, height = height_mm / 25.4, device = cairo_pdf)
  ggsave(paste0(base_name, ".png"), p, width = 183 / 25.4, height = height_mm / 25.4, dpi = 400)
}
