# 04_match_A1_coordinates.R
# Purpose: Match A1 expression spots to spatial coordinates

library(ggplot2)

# Read the A1 expression matrix
a1_counts <- read.delim(
  gzfile(
    "data/raw/count-matrices/count-matrices/A1.tsv.gz"
  ),
  check.names = FALSE
)

# Give the unnamed first column a clear name
names(a1_counts)[1] <- "spot_id"

# Read the A1 spatial-coordinate table
a1_coordinates <- read.delim(
  gzfile(
    "data/raw/spot-selections/A1_selection.tsv.gz"
  ),
  check.names = FALSE
)

# Construct spot IDs such as 10x13
a1_coordinates$spot_id <- paste0(
  a1_coordinates$x,
  "x",
  a1_coordinates$y
)

cat("Expression spots:", nrow(a1_counts), "\n")
cat("Coordinate spots:", nrow(a1_coordinates), "\n")

# Check for duplicated IDs
cat(
  "Duplicated expression IDs:",
  anyDuplicated(a1_counts$spot_id),
  "\n"
)

cat(
  "Duplicated coordinate IDs:",
  anyDuplicated(a1_coordinates$spot_id),
  "\n"
)

# Compare the two ID sets
expression_without_coordinates <- setdiff(
  a1_counts$spot_id,
  a1_coordinates$spot_id
)

coordinates_without_expression <- setdiff(
  a1_coordinates$spot_id,
  a1_counts$spot_id
)

cat(
  "Expression spots without coordinates:",
  length(expression_without_coordinates),
  "\n"
)

cat(
  "Coordinate spots without expression:",
  length(coordinates_without_expression),
  "\n"
)

print(coordinates_without_expression)

# Reorder coordinates to exactly match the expression matrix
match_index <- match(
  a1_counts$spot_id,
  a1_coordinates$spot_id
)

a1_spatial_metadata <- a1_coordinates[match_index, ]

cat(
  "Spot order matches:",
  identical(
    a1_counts$spot_id,
    a1_spatial_metadata$spot_id
  ),
  "\n"
)

cat(
  "Missing coordinate rows after matching:",
  sum(is.na(a1_spatial_metadata$spot_id)),
  "\n"
)

# Save the matched coordinate table
write.csv(
  a1_spatial_metadata,
  "data/interim/A1_spot_metadata.csv",
  row.names = FALSE
)

# Draw the spot layout without the H&E image
a1_spot_plot <- ggplot(
  a1_spatial_metadata,
  aes(x = pixel_x, y = pixel_y)
) +
  geom_point(
    color = "#2C7FB8",
    size = 1.8
  ) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  labs(
    title = "A1 spatial spot layout",
    subtitle = "346 expression spots with matched image coordinates"
  )

print(a1_spot_plot)

ggsave(
  filename = "results/figures/A1_spot_layout.png",
  plot = a1_spot_plot,
  width = 6,
  height = 5,
  dpi = 300
)