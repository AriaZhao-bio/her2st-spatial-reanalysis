# 03_inspect_data_structure.R
# Purpose: Match count matrices, spatial coordinates, and pathology labels

count_files <- list.files(
  path = "data/raw/count-matrices",
  pattern = "\\.tsv\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

selection_files <- list.files(
  path = "data/raw/spot-selections",
  pattern = "_selection\\.tsv\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

meta_files <- list.files(
  path = "data/raw/meta",
  pattern = "_labeled_coordinates\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

# Extract section IDs from filenames
count_ids <- sub("\\.tsv\\.gz$", "", basename(count_files))

selection_ids <- sub(
  "_selection\\.tsv\\.gz$",
  "",
  basename(selection_files)
)

meta_ids <- sub(
  "_labeled_coordinates\\.tsv$",
  "",
  basename(meta_files)
)

# Print the number of files
cat("Count matrices:", length(count_files), "\n")
cat("Spot-selection files:", length(selection_files), "\n")
cat("Pathology metadata files:", length(meta_files), "\n")

# Check whether every expression matrix has a coordinate file
cat(
  "Count and coordinate IDs match:",
  setequal(count_ids, selection_ids),
  "\n"
)

# Build a section-level inventory
count_map <- setNames(count_files, count_ids)
selection_map <- setNames(selection_files, selection_ids)

section_inventory <- data.frame(
  sample_id = sort(count_ids),
  patient = substr(sort(count_ids), 1, 1),
  section = as.integer(substring(sort(count_ids), 2)),
  count_file = unname(count_map[sort(count_ids)]),
  selection_file = unname(selection_map[sort(count_ids)]),
  has_pathology_label = sort(count_ids) %in% meta_ids
)

print(section_inventory)

# Count sections per patient
patient_summary <- aggregate(
  sample_id ~ patient,
  data = section_inventory,
  FUN = length
)

names(patient_summary)[2] <- "number_of_sections"

print(patient_summary)

# Save the inventory
write.csv(
  section_inventory,
  file = "results/tables/file_inventory.csv",
  row.names = FALSE
)

# Read the first expression matrix: sample A1
a1_file <- count_map[["A1"]]

a1_counts <- read.delim(
  gzfile(a1_file),
  check.names = FALSE
)

cat("A1 matrix dimensions:", dim(a1_counts), "\n")

print(names(a1_counts)[1:10])
print(a1_counts[1:5, 1:6])