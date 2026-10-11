############################################################
# Single Cell RNA-seq Analysis
# Dataset: GSE155698
# Cancer: Pancreatic Ductal Adenocarcinoma (PDAC)
#
# Goal:
# Load tumor and adjacent normal pancreatic tissue samples
# Create Seurat objects and merge into one dataset
############################################################


##############################
# 1. Load libraries
##############################

library(Seurat)
library(dplyr)
library(Matrix)


##############################
# 2. Set working directory
##############################

# Locate the repository root so the script runs from any working directory.
local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  root <- if (length(f)) dirname(dirname(normalizePath(f)))
          else if (file.exists("scripts")) normalizePath(".")
          else normalizePath("..")
  setwd(root)
})


##############################
# 3. Read sample metadata
##############################

metadata <- read.csv(
  "data/sample_metadata.csv"
)


print(metadata)



##############################
# 4. Create empty Seurat list
##############################

seurat_list <- list()



##############################
# 5. Load each sample
##############################

for(i in 1:nrow(metadata)){


  sample_name <- metadata$sample[i]

  sample_path <- metadata$path[i]


  print(
    paste(
      "Loading:",
      sample_name
    )
  )



  ###########################
  # Check 10X format
  ###########################


  # New 10X format
  # filtered_feature_bc_matrix

  if(
    dir.exists(sample_path)
  ){

    counts <- Read10X(
      data.dir = sample_path
    )


  } else {


    stop(
      paste(
        "Missing folder:",
        sample_path
      )
    )

  }



  ###########################
  # Create Seurat object
  ###########################

  obj <- CreateSeuratObject(

    counts = counts,

    project = sample_name,

    min.features = 200

  )



  ###########################
  # Add metadata
  ###########################

  obj$sample <- sample_name

  obj$condition <- metadata$condition[i]



  ###########################
  # Store object
  ###########################

  seurat_list[[sample_name]] <- obj


}




##############################
# 6. Check cell numbers
##############################

cell_numbers <- lapply(
  seurat_list,
  ncol
)


print(cell_numbers)



##############################
# 7. Merge all samples
##############################

pdac_combined <- merge(

  seurat_list[[1]],

  y = seurat_list[-1],

  add.cell.ids = names(seurat_list)

)



##############################
# 8. Check merged object
##############################

print(pdac_combined)



##############################
# 9. Save raw merged object
##############################

saveRDS(

  pdac_combined,

  file =
    "data/processed/PDAC_merged_raw.rds"

)



print(
  "Finished loading and merging PDAC scRNA-seq dataset"
)
