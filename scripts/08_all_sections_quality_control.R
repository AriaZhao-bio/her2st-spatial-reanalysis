# 08_all_sections_quality_control.R
# Purpose: Calculate comparable QC metrics for all 36 sections

library(ggplot2)
library(patchwork)

count_files <- list.files(
  path = "data/raw/count-matrices",
  pattern = "\\.tsv\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

count_files <- sort(count_files)

all_qc_list <- vector(
  mode = "list",
  length = length(count_files)
)

for (i in seq_along(count_files)) {
  
  count_file <- count_files[i]
  
  sample_id <- sub(
    "\\.tsv\\.gz$",
    "",
    basename(count_file)
  )
  
  message(
    "[",
    i,
    "/",
    length(count_files),
    "] Reading ",
    sample_id
  )
  
  count_table <- read.delim(
    gzfile(count_file),
    check.names = FALSE
  )
  
  spot_ids <- count_table[[1]]
  
  count_matrix <- as.matrix(
    count_table[, -1]
  )
  
  storage.mode(count_matrix) <- "numeric"
  
  total_counts <- rowSums(count_matrix)
  
  detected_genes <- rowSums(
    count_matrix > 0
  )
  
  all_qc_list[[i]] <- data.frame(
    sample_id = sample_id,
    patient = substr(sample_id, 1, 1),
    section = as.integer(
      substring(sample_id, 2)
    ),
    spot_id = spot_ids,
    total_counts = total_counts,
    detected_genes = detected_genes
  )
  
  rm(
    count_table,
    count_matrix
  )
}

all_spot_qc <- do.call(
  rbind,
  all_qc_list
)

row.names(all_spot_qc) <- NULL

# Provisional flag only; no spots are removed here
all_spot_qc$provisional_low_qc <-
  all_spot_qc$total_counts < 500 |
  all_spot_qc$detected_genes < 300

cat(
  "Number of sections:",
  length(unique(all_spot_qc$sample_id)),
  "\n"
)

cat(
  "Total number of spots:",
  nrow(all_spot_qc),
  "\n"
)

cat(
  "Provisionally flagged spots:",
  sum(all_spot_qc$provisional_low_qc),
  "\n"
)

# Produce one summary row per section
section_groups <- split(
  all_spot_qc,
  all_spot_qc$sample_id
)

section_summary <- do.call(
  rbind,
  lapply(
    section_groups,
    function(section_data) {
      
      data.frame(
        sample_id = section_data$sample_id[1],
        patient = section_data$patient[1],
        section = section_data$section[1],
        number_of_spots = nrow(section_data),
        median_total_counts = median(
          section_data$total_counts
        ),
        minimum_total_counts = min(
          section_data$total_counts
        ),
        median_detected_genes = median(
          section_data$detected_genes
        ),
        minimum_detected_genes = min(
          section_data$detected_genes
        ),
        provisional_low_qc_spots = sum(
          section_data$provisional_low_qc
        )
      )
    }
  )
)

row.names(section_summary) <- NULL

section_summary <- section_summary[
  order(
    section_summary$patient,
    section_summary$section
  ),
]

print(section_summary)

# Save both detailed and section-level tables
write.csv(
  all_spot_qc,
  "data/interim/all_spot_QC_metrics.csv",
  row.names = FALSE
)

write.csv(
  section_summary,
  "results/tables/section_QC_summary.csv",
  row.names = FALSE
)

# Preserve sample order in plots
sample_order <- section_summary$sample_id

all_spot_qc$sample_id <- factor(
  all_spot_qc$sample_id,
  levels = sample_order
)

# QC distribution across sections
p_counts <- ggplot(
  all_spot_qc,
  aes(
    x = sample_id,
    y = total_counts
  )
) +
  geom_boxplot(
    outlier.size = 0.4,
    fill = "#80B1D3"
  ) +
  scale_y_log10() +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(
    title = "Total counts across 36 sections",
    x = "Section",
    y = "Total counts per spot (log10 scale)"
  )

p_genes <- ggplot(
  all_spot_qc,
  aes(
    x = sample_id,
    y = detected_genes
  )
) +
  geom_boxplot(
    outlier.size = 0.4,
    fill = "#B3DE69"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(
    title = "Detected genes across 36 sections",
    x = "Section",
    y = "Detected genes per spot"
  )

all_sections_qc_plot <- p_counts / p_genes

ggsave(
  "results/figures/all_sections_QC_distributions.png",
  all_sections_qc_plot,
  width = 14,
  height = 9,
  dpi = 300
)