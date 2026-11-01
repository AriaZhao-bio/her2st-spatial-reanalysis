# 06_overlay_A1_spots.R
# Purpose: Overlay A1 expression spots on the H&E tissue image

library(jpeg)

# File paths
a1_image_file <- paste0(
  "data/raw/images/images/HE/",
  "A1.jpg"
)

a1_metadata_file <- paste0(
  "data/interim/",
  "A1_spot_metadata.csv"
)

output_file <- paste0(
  "results/figures/",
  "A1_HE_spot_overlay.png"
)

# Read the large JPEG using a memory-efficient format
a1_image <- readJPEG(
  a1_image_file,
  native = TRUE
)

a1_spatial_metadata <- read.csv(
  a1_metadata_file
)

# JPEG dimensions
image_height <- dim(a1_image)[1]
image_width <- dim(a1_image)[2]

cat("Image width:", image_width, "\n")
cat("Image height:", image_height, "\n")
cat("Number of expression spots:", nrow(a1_spatial_metadata), "\n")

# Check whether all spot centers are inside the image
x_inside <- all(
  a1_spatial_metadata$pixel_x >= 0 &
    a1_spatial_metadata$pixel_x <= image_width
)

y_inside <- all(
  a1_spatial_metadata$pixel_y >= 0 &
    a1_spatial_metadata$pixel_y <= image_height
)

cat("All x coordinates inside image:", x_inside, "\n")
cat("All y coordinates inside image:", y_inside, "\n")

if (!x_inside || !y_inside) {
  stop("Some spot coordinates fall outside the image.")
}

# Preserve the original image aspect ratio
output_width <- 2400
output_height <- round(
  output_width * image_height / image_width
)

png(
  filename = output_file,
  width = output_width,
  height = output_height,
  res = 250,
  bg = "white"
)

par(
  mar = c(0, 0, 0, 0),
  xaxs = "i",
  yaxs = "i"
)

plot.new()

plot.window(
  xlim = c(0, image_width),
  ylim = c(image_height, 0),
  asp = 1
)

rasterImage(
  image = a1_image,
  xleft = 0,
  ybottom = image_height,
  xright = image_width,
  ytop = 0
)

points(
  x = a1_spatial_metadata$pixel_x,
  y = a1_spatial_metadata$pixel_y,
  pch = 21,
  bg = adjustcolor("#00FFFF", alpha.f = 0.45),
  col = "#111111",
  cex = 0.8,
  lwd = 0.4
)

dev.off()

cat("Saved overlay:", output_file, "\n")