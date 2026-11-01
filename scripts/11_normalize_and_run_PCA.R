# Purpose: Normalize each patient and inspect section effects with PCA

library(Seurat)
library(ggplot2)
library(patchwork)

raw_files <- list.files(
  "data/processed",
  pattern = "HER2ST_patient_[A-H]_raw\\.rds$",
  full.names = TRUE
)

raw_files <- sort(raw_files)

stopifnot(length(raw_files) == 8)

pca_plot_list <- vector("list", length(raw_files))
normalization_manifest <- vector("list", length(raw_files))

for (i in seq_along(raw_files)) {
  raw_file <- raw_files[i]
  patient_id <- sub(
    ".*HER2ST_patient_([A-H])_raw\\.rds$",
    "\\1",
    raw_file
  )

  message("[", i, "/", length(raw_files), "] Normalizing patient ", patient_id)

  patient_object <- readRDS(raw_file)
  DefaultAssay(patient_object) <- "RNA"

  patient_object <- NormalizeData(
    patient_object,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  patient_object <- FindVariableFeatures(
    patient_object,
    selection.method = "vst",
    nfeatures = 3000,
    verbose = FALSE
  )

  patient_object <- ScaleData(
    patient_object,
    features = VariableFeatures(patient_object),
    verbose = FALSE
  )

  patient_object <- RunPCA(
    patient_object,
    features = VariableFeatures(patient_object),
    npcs = 30,
    verbose = FALSE
  )

  output_file <- file.path(
    "data", "processed",
    paste0("HER2ST_patient_", patient_id, "_lognorm_pca.rds")
  )
  saveRDS(patient_object, output_file)

  pca_coordinates <- as.data.frame(
    Embeddings(patient_object, reduction = "pca")[, 1:2, drop = FALSE]
  )
  pca_coordinates$sample_id <- patient_object$sample_id
  pca_coordinates$section <- factor(patient_object$section)

  pca_plot_list[[i]] <- ggplot(
    pca_coordinates,
    aes(x = PC_1, y = PC_2, color = section)
  ) +
    geom_point(size = 0.8, alpha = 0.7) +
    theme_classic() +
    labs(title = paste0("Patient ", patient_id), color = "Section")

  normalization_manifest[[i]] <- data.frame(
    patient = patient_id,
    number_of_spots = ncol(patient_object),
    number_of_genes = nrow(patient_object),
    variable_genes = length(VariableFeatures(patient_object)),
    PCA_dimensions = ncol(Embeddings(patient_object, reduction = "pca")),
    file = output_file
  )

  message("  Saved: ", output_file)
  rm(patient_object, pca_coordinates)
  gc()
}

normalization_manifest <- do.call(rbind, normalization_manifest)
print(normalization_manifest)

total_normalized_spots <- sum(normalization_manifest$number_of_spots)
cat("\nTotal normalized spots:", total_normalized_spots, "\n")
stopifnot(total_normalized_spots == 13302)

write.csv(
  normalization_manifest,
  "results/tables/normalization_manifest.csv",
  row.names = FALSE
)

pca_section_plot <- wrap_plots(pca_plot_list, ncol = 2)

ggsave(
  "results/figures/PCA_by_section_all_patients.png",
  pca_section_plot,
  width = 11,
  height = 14,
  dpi = 300
)
