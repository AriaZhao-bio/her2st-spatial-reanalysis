# 07_A1_quality_control.R
# Purpose: Calculate and visualize basic QC metrics for sample A1

library(ggplot2)
library(patchwork)

# Read expression counts
a1_counts <- read.delim(
  gzfile(
    "data/raw/count-matrices/count-matrices/A1.tsv.gz"
  ),
  check.names = FALSE
)

names(a1_counts)[1] <- "spot_id"

# Read matched spatial coordinates
a1_metadata <- read.csv(
  "data/interim/A1_spot_metadata.csv"
)

# Confirm that the spot order is still identical
stopifnot(
  identical(
    a1_counts$spot_id,
    a1_metadata$spot_id
  )
)

# Convert the gene-expression portion into a numeric matrix
count_matrix <- as.matrix(
  a1_counts[, -1]
)

storage.mode(count_matrix) <- "numeric"

# Calculate QC metrics
total_counts <- rowSums(count_matrix)

detected_genes <- rowSums(
  count_matrix > 0
)

# Check whether mitochondrial genes are available
mitochondrial_genes <- grepl(
  "^MT-",
  colnames(count_matrix)
)

number_of_mitochondrial_genes <- sum(
  mitochondrial_genes
)

cat(
  "Number of mitochondrial genes:",
  number_of_mitochondrial_genes,
  "\n"
)

if (number_of_mitochondrial_genes > 0) {
  
  percent_mitochondrial <- 100 *
    rowSums(
      count_matrix[, mitochondrial_genes, drop = FALSE]
    ) /
    total_counts
  
} else {
  
  percent_mitochondrial <- rep(
    NA_real_,
    nrow(count_matrix)
  )
}

# Extract raw ERBB2 counts
stopifnot(
  "ERBB2" %in% colnames(count_matrix)
)

erbb2_counts <- count_matrix[, "ERBB2"]

# Combine QC values with spatial coordinates
a1_qc <- data.frame(
  spot_id = a1_counts$spot_id,
  total_counts = total_counts,
  detected_genes = detected_genes,
  percent_mitochondrial = percent_mitochondrial,
  ERBB2_raw_counts = erbb2_counts,
  pixel_x = a1_metadata$pixel_x,
  pixel_y = a1_metadata$pixel_y
)

# Print numerical summaries
cat("\nTotal-count summary:\n")
print(summary(a1_qc$total_counts))

cat("\nDetected-gene summary:\n")
print(summary(a1_qc$detected_genes))

cat("\nERBB2 raw-count summary:\n")
print(summary(a1_qc$ERBB2_raw_counts))

cat(
  "\nSpots with detectable ERBB2:",
  sum(a1_qc$ERBB2_raw_counts > 0),
  "of",
  nrow(a1_qc),
  "\n"
)

cat("\nTen spots with the lowest total counts:\n")

print(
  head(
    a1_qc[
      order(a1_qc$total_counts),
      c(
        "spot_id",
        "total_counts",
        "detected_genes",
        "ERBB2_raw_counts"
      )
    ],
    10
  )
)

# Save the QC table
write.csv(
  a1_qc,
  "data/interim/A1_QC_metrics.csv",
  row.names = FALSE
)

# Distribution plots
p_total_counts <- ggplot(
  a1_qc,
  aes(x = total_counts)
) +
  geom_histogram(
    bins = 30,
    fill = "#377EB8",
    color = "white"
  ) +
  scale_x_log10() +
  theme_classic() +
  labs(
    title = "Total counts per spot",
    x = "Total counts (log10 scale)",
    y = "Number of spots"
  )

p_detected_genes <- ggplot(
  a1_qc,
  aes(x = detected_genes)
) +
  geom_histogram(
    bins = 30,
    fill = "#4DAF4A",
    color = "white"
  ) +
  theme_classic() +
  labs(
    title = "Detected genes per spot",
    x = "Number of detected genes",
    y = "Number of spots"
  )

p_counts_vs_genes <- ggplot(
  a1_qc,
  aes(
    x = total_counts,
    y = detected_genes
  )
) +
  geom_point(
    color = "#444444",
    alpha = 0.7
  ) +
  scale_x_log10() +
  theme_classic() +
  labs(
    title = "Counts versus detected genes",
    x = "Total counts (log10 scale)",
    y = "Detected genes"
  )

qc_distribution_plot <-
  (p_total_counts | p_detected_genes) /
  p_counts_vs_genes

ggsave(
  "results/figures/A1_QC_distributions.png",
  qc_distribution_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# Spatial QC maps without the H&E background
p_spatial_counts <- ggplot(
  a1_qc,
  aes(
    x = pixel_x,
    y = pixel_y,
    color = log10(total_counts + 1)
  )
) +
  geom_point(size = 2.6) +
  scale_y_reverse() +
  scale_color_viridis_c() +
  coord_fixed() +
  theme_void() +
  labs(
    title = "A1: total counts",
    color = "log10(counts + 1)"
  )

p_spatial_genes <- ggplot(
  a1_qc,
  aes(
    x = pixel_x,
    y = pixel_y,
    color = detected_genes
  )
) +
  geom_point(size = 2.6) +
  scale_y_reverse() +
  scale_color_viridis_c() +
  coord_fixed() +
  theme_void() +
  labs(
    title = "A1: detected genes",
    color = "Genes"
  )

qc_spatial_plot <- p_spatial_counts | p_spatial_genes

ggsave(
  "results/figures/A1_QC_spatial.png",
  qc_spatial_plot,
  width = 10,
  height = 5,
  dpi = 300
)