# 02_extract_data.R
# Purpose: Extract the encrypted HER2 spatial data archives

seven_zip <- "C:/Program Files/NVIDIA Corporation/NVIDIA app/7z.exe"

if (!file.exists(seven_zip)) {
  stop("7-Zip executable was not found.")
}

# Convert to a Windows short path to avoid problems caused by spaces
seven_zip_command <- utils::shortPathName(seven_zip)

# Public password provided on the official Zenodo dataset page
zip_passwords <- c(
  "zNLXkYk3Q9znUseS",
  "yUx44SzG6NdB32gY",
  "yUx44SzG6NdB32gY"
)

extraction_plan <- data.frame(
  archive = file.path(
    "data",
    "raw",
    c(
      "count-matrices.zip",
      "meta.zip",
      "spot-selections.zip"
    )
  ),
  output_dir = file.path(
    "data",
    "raw",
    c(
      "count-matrices",
      "meta",
      "spot-selections"
    )
  )
)

for (i in seq_len(nrow(extraction_plan))) {
  
  archive_path <- extraction_plan$archive[i]
  output_path <- extraction_plan$output_dir[i]
  
  dir.create(
    output_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  message("Extracting: ", archive_path)
  
  exit_status <- system2(
    command = seven_zip_command,
    args = c(
      "x",
      archive_path,
      paste0("-o", output_path),
      paste0("-p", zip_passwords[i]),
      "-y"
    )
  )
  
  if (exit_status != 0) {
    stop("Extraction failed for: ", archive_path)
  }
}

extraction_report <- data.frame(
  folder = extraction_plan$output_dir,
  number_of_files = vapply(
    extraction_plan$output_dir,
    function(folder_path) {
      length(list.files(folder_path, recursive = TRUE))
    },
    integer(1)
  )
)

print(extraction_report)