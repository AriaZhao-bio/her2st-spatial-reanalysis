# 09_finalize_QC_filter.R
# Purpose: Apply and document the final spot-level QC rule

library(ggplot2)

all_spot_qc <- read.csv(
  "data/interim/all_spot_QC_metrics.csv"
)

# Final rule based on the original HER2ST analysis
minimum_detected_genes <- 300

all_spot_qc$keep_after_qc <-
  all_spot_qc$detected_genes >= minimum_detected_genes

all_spot_qc$qc_decision <- ifelse(
  all_spot_qc$keep_after_qc,
  "retained",
  "removed: fewer than 300 detected genes"
)

total_spots <- nrow(all_spot_qc)

retained_spots <- sum(
  all_spot_qc$keep_after_qc
)

removed_spots <- sum(
  !all_spot_qc$keep_after_qc
)

removed_percent <- 100 *
  removed_spots /
  total_spots

cat("Total spots:", total_spots, "\n")
cat("Retained spots:", retained_spots, "\n")
cat("Removed spots:", removed_spots, "\n")
cat(
  "Removed percentage:",
  round(removed_percent, 2),
  "%\n"
)

# Section-level filtering summary
section_groups <- split(
  all_spot_qc,
  all_spot_qc$sample_id
)

filter_summary <- do.call(
  rbind,
  lapply(
    section_groups,
    function(section_data) {
      
      number_before <- nrow(section_data)
      
      number_retained <- sum(
        section_data$keep_after_qc
      )
      
      number_removed <- sum(
        !section_data$keep_after_qc
      )
      
      data.frame(
        sample_id = section_data$sample_id[1],
        patient = section_data$patient[1],
        section = section_data$section[1],
        spots_before_QC = number_before,
        spots_retained = number_retained,
        spots_removed = number_removed,
        percent_removed = round(
          100 * number_removed / number_before,
          2
        )
      )
    }
  )
)

row.names(filter_summary) <- NULL

filter_summary <- filter_summary[
  order(
    filter_summary$patient,
    filter_summary$section
  ),
]

print(filter_summary)

# Save without overwriting the original QC table
write.csv(
  all_spot_qc,
  "data/interim/all_spot_QC_filter_decisions.csv",
  row.names = FALSE
)

write.csv(
  filter_summary,
  "results/tables/QC_filter_summary.csv",
  row.names = FALSE
)

# Visualize how filtering differs among sections
filter_summary$sample_id <- factor(
  filter_summary$sample_id,
  levels = filter_summary$sample_id
)

qc_removal_plot <- ggplot(
  filter_summary,
  aes(
    x = sample_id,
    y = percent_removed,
    fill = patient
  )
) +
  geom_col() +
  geom_hline(
    yintercept = removed_percent,
    linetype = "dashed",
    color = "#333333"
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
    title = "Spot removal after QC",
    subtitle = paste0(
      "Final rule: at least ",
      minimum_detected_genes,
      " detected genes"
    ),
    x = "Section",
    y = "Spots removed (%)",
    fill = "Patient"
  )

ggsave(
  "results/figures/QC_spot_removal_by_section.png",
  qc_removal_plot,
  width = 11,
  height = 5,
  dpi = 300
)