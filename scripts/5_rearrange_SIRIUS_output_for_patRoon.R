# 17 March 2026
# script for re-arranging SIRIUS output through patRoon so that it is compatible with MLinvitroTox
# R version 4.5.2


### ChatGPT code, slightly modified

# Set the root directory
root_dir <- "\\mlinvitrotox\\data\\sirius_hilic_compounds"

# List parent folders directly under root
parent_folders <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)

for (parent in parent_folders) {
  
  # List immediate subfolders of this parent
  subfolders <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
  
  for (sub in subfolders) {
    
    # List files in the subfolder (no recursion)
    files <- list.files(sub, full.names = TRUE, recursive = FALSE)
    
    if (length(files) == 0) next
    
    # Destination paths in the parent folder
    dest_files <- file.path(parent, basename(files))
    
    #  # Copy files (no overwrite by default)
    # file.copy(from = files, to = dest_files, overwrite = FALSE)
    
    # Move files (to reduce folder size)
    file.rename(from = files, to = dest_files)
  }
}
