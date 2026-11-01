# 01_download_data.R
# Purpose: Download the first set of public HER2 spatial data files

options(timeout = 600)

base_url <- "https://zenodo.org/records/4751624/files/"

data_files <- c(
  "count-matrices.zip",
  "meta.zip",
  "spot-selections.zip"
)

for (file_name in data_files) {
  
  destination <- file.path("data", "raw", file_name)
  download_url <- paste0(base_url, file_name, "?download=1")
  
  if (!file.exists(destination)) {
    message("Downloading: ", file_name)
    
    download.file(
      url = download_url,
      destfile = destination,
      mode = "wb",
      method = "libcurl"
    )
  } else {
    message("Already exists: ", file_name)
  }
}

download_report <- data.frame(
  file = data_files,
  size_MB = round(
    file.info(file.path("data", "raw", data_files))$size / 1024^2,
    2
  )
)

print(download_report)