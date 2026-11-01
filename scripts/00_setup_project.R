# 00_setup_project.R
# Purpose: Create the directory structure for the HER2 spatial project

project_folders <- c(
  "data/raw",
  "data/interim",
  "data/processed",
  "scripts",
  "results/tables",
  "results/figures",
  "docs"
)

for (folder in project_folders) {
  dir.create(
    path = folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

list.files(
  path = ".",
  recursive = TRUE,
  include.dirs = TRUE
)