# 05_download_extract_images.R
# Purpose: Download, verify, and extract the HER2 tissue images

options(timeout = 3600)

image_url <- paste0(
  "https://zenodo.org/records/4751624/files/",
  "images.zip?download=1"
)

image_zip <- "data/raw/images.zip"

# Download only when the ZIP file does not already exist
if (!file.exists(image_zip)) {
  
  message("Downloading images.zip ...")
  
  download.file(
    url = image_url,
    destfile = image_zip,
    mode = "wb",
    method = "libcurl"
  )
  
} else {
  
  message("images.zip already exists; skipping download.")
}

# Report downloaded size
image_size_mb <- file.info(image_zip)$size / 1024^2

cat(
  "Downloaded size:",
  round(image_size_mb, 2),
  "MB\n"
)

# Verify the digital fingerprint
observed_md5 <- unname(tools::md5sum(image_zip))
expected_md5 <- "383e23ac398790dbdb408ec7b8b5d3b2"

cat("Observed MD5:", observed_md5, "\n")
cat("Expected MD5:", expected_md5, "\n")
cat(
  "MD5 matches:",
  identical(observed_md5, expected_md5),
  "\n"
)

if (!identical(observed_md5, expected_md5)) {
  stop(
    "The image ZIP did not pass the MD5 check. ",
    "Do not extract it."
  )
}

# Locate the available 7-Zip program
seven_zip <- paste0(
  "C:/Program Files/NVIDIA Corporation/",
  "NVIDIA app/7z.exe"
)

if (!file.exists(seven_zip)) {
  stop("7-Zip executable was not found.")
}

seven_zip_command <- utils::shortPathName(seven_zip)

# Images use the same password as the count matrices
image_password <- "zNLXkYk3Q9znUseS"

image_output_dir <- "data/raw/images"

dir.create(
  image_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

exit_status <- system2(
  command = seven_zip_command,
  args = c(
    "x",
    image_zip,
    paste0("-o", image_output_dir),
    paste0("-p", image_password),
    "-y"
  )
)

if (exit_status != 0) {
  stop("Image extraction failed.")
}

# Inventory all extracted JPEG files
image_files <- list.files(
  path = image_output_dir,
  pattern = "\\.jpg$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

cat("Number of extracted JPG files:", length(image_files), "\n")

print(
  data.frame(
    file = basename(image_files),
    folder = dirname(image_files)
  )
)