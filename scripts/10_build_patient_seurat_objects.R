# 10_build_patient_seurat_objects.R
# Purpose: Build one QC-filtered raw Seurat object per patient

library(Seurat)
library(Matrix)

qc_decisions <- read.csv(
  "data/interim/all_spot_QC_filter_decisions.csv"
)

qc_decisions$keep_after_qc <- as.logical(
  qc_decisions$keep_after_qc
)

# Map sample IDs to count files
count_files <- list.files(
  "data/raw/count-matrices",
  pattern = "\\.tsv\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

count_ids <- sub(
  "\\.tsv\\.gz$",
  "",
  basename(count_files)
)

count_map <- setNames(
  count_files,
  count_ids
)

# Map sample IDs to coordinate files
coordinate_files <- list.files(
  "data/raw/spot-selections",
  pattern = "_selection\\.tsv\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

coordinate_ids <- sub(
  "_selection\\.tsv\\.gz$",
  "",
  basename(coordinate_files)
)

coordinate_map <- setNames(
  coordinate_files,
  coordinate_ids
)

# Function: build one section-level Seurat object
build_section_object <- function(
    sample_id,
    patient_id
) {
  
  message("  Reading section ", sample_id)
  
  count_table <- read.delim(
    gzfile(count_map[[sample_id]]),
    check.names = FALSE
  )
  
  names(count_table)[1] <- "spot_id"
  
  section_decisions <- qc_decisions[
    qc_decisions$sample_id == sample_id &
      qc_decisions$keep_after_qc,
  ]
  
  retained_spot_ids <-
    section_decisions$spot_id
  
  count_index <- match(
    retained_spot_ids,
    count_table$spot_id
  )
  
  stopifnot(!anyNA(count_index))
  
  numeric_counts <- data.matrix(
    count_table[
      count_index,
      -1,
      drop = FALSE
    ]
  )
  
  unique_spot_ids <- paste(
    sample_id,
    retained_spot_ids,
    sep = "_"
  )
  
  section_counts <- Matrix(
    t(numeric_counts),
    sparse = TRUE
  )
  
  rownames(section_counts) <-
    names(count_table)[-1]
  
  colnames(section_counts) <-
    unique_spot_ids
  
  # Match spatial coordinates
  coordinate_table <- read.delim(
    gzfile(coordinate_map[[sample_id]]),
    check.names = FALSE
  )
  
  coordinate_table$spot_id <- paste0(
    coordinate_table$x,
    "x",
    coordinate_table$y
  )
  
  coordinate_index <- match(
    retained_spot_ids,
    coordinate_table$spot_id
  )
  
  stopifnot(!anyNA(coordinate_index))
  
  matched_coordinates <- coordinate_table[
    coordinate_index,
  ]
  
  section_metadata <- data.frame(
    original_spot_id = retained_spot_ids,
    patient = patient_id,
    sample_id = sample_id,
    section = as.integer(
      substring(sample_id, 2)
    ),
    array_x = matched_coordinates$x,
    array_y = matched_coordinates$y,
    pixel_x = matched_coordinates$pixel_x,
    pixel_y = matched_coordinates$pixel_y,
    total_counts_before_Seurat =
      section_decisions$total_counts,
    detected_genes_before_Seurat =
      section_decisions$detected_genes,
    row.names = unique_spot_ids
  )
  
  section_object <- CreateSeuratObject(
    counts = section_counts,
    project = paste0(
      "HER2ST_",
      sample_id
    ),
    assay = "RNA",
    min.cells = 0,
    min.features = 0,
    meta.data = section_metadata
  )
  
  return(section_object)
}

patients <- LETTERS[1:8]

manifest_list <- vector(
  "list",
  length(patients)
)

for (patient_index in seq_along(patients)) {
  
  patient_id <- patients[patient_index]
  
  message(
    "\nBuilding patient ",
    patient_id,
    " ..."
  )
  
  patient_samples <- sort(
    unique(
      qc_decisions$sample_id[
        qc_decisions$patient == patient_id
      ]
    )
  )
  
  section_objects <- lapply(
    patient_samples,
    build_section_object,
    patient_id = patient_id
  )
  
  # Seurat aligns differing gene sets by gene name
  if (length(section_objects) == 1) {
    
    patient_object <- section_objects[[1]]
    
  } else {
    
    patient_object <- merge(
      x = section_objects[[1]],
      y = section_objects[-1],
      merge.data = FALSE,
      project = paste0(
        "HER2ST_patient_",
        patient_id
      )
    )
    
    # Combine the section-specific count layers
    patient_object <- JoinLayers(
      patient_object
    )
  }
  
  stopifnot(
    ncol(patient_object) ==
      sum(
        qc_decisions$patient == patient_id &
          qc_decisions$keep_after_qc
      )
  )
  
  output_file <- file.path(
    "data",
    "processed",
    paste0(
      "HER2ST_patient_",
      patient_id,
      "_raw.rds"
    )
  )
  
  saveRDS(
    patient_object,
    output_file
  )
  
  manifest_list[[patient_index]] <- data.frame(
    patient = patient_id,
    number_of_sections =
      length(patient_samples),
    number_of_spots =
      ncol(patient_object),
    number_of_genes =
      nrow(patient_object),
    object_size_MB = round(
      as.numeric(object.size(patient_object)) /
        1024^2,
      1
    ),
    file = output_file
  )
  
  message(
    "  Saved patient ",
    patient_id,
    ": ",
    ncol(patient_object),
    " spots and ",
    nrow(patient_object),
    " genes"
  )
  
  rm(
    section_objects,
    patient_object
  )
  
  gc()
}

object_manifest <- do.call(
  rbind,
  manifest_list
)

print(object_manifest)

cat(
  "\nTotal spots in Seurat objects:",
  sum(object_manifest$number_of_spots),
  "\n"
)

write.csv(
  object_manifest,
  "results/tables/Seurat_object_manifest.csv",
  row.names = FALSE
)