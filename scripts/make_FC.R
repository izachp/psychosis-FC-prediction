library(tidyverse)

# -------------------------------- SETUP --------------------------------
repo = 'path/to/psychosis-FC-prediction'
prepro = 'dt_AROMA_8Phys-4GMR_bpf'      # Supps also used 'dt_AROMA_8Phys_bpf'
parcellation = '419_parc'               # Supps also used '328_parc'
# -----------------------------------------------------------------------

if (parcellation == '328_parc') {
  # Read in list of low-signal regions to be excluded
  excluded_regions <- read.table(paste0(repo, '/data/low_signal_regions.txt'), header = F) %>%
    select(2) %>%
    unlist()
}

# Loop through fMRI sessions ('1' is baseline, '2' is 3-month follow-up)
for (ses in 1:2) {
  ts_dir <- paste0(repo, '/data/timeseries/', prepro, '/', parcellation, '/ses-', ses)
  out_dir <- paste0(repo, '/data/conmats/', prepro, '/', parcellation, '/ses-', ses, '/')
  
  ts_list <- list.files(ts_dir, 
                        pattern = "*.txt", 
                        recursive = TRUE,
                        full.names = T)
  
  # Loop through subjects
  for (s in ts_list) {
    
    subtitle <- substr(s,nchar(ts_dir)+18,nchar(ts_dir)+30)  # e.g: 'sub-036_ses-1'
    ts <- as.matrix(read.table(s), header = F)
    
    if (parcellation == '328_parc') {
      ts <- ts[-excluded_regions, ]
   }
    
    # Calculate FC
    conmat <- cor(t(ts))
    
    # Save FC
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    write.table(round(conmat,3), 
                file = paste0(out_dir, 'conmat_', subtitle, '.txt'), 
                row.names = F, 
                col.names = F,
                quote = F)
  }
}

# --------------------------------- FC CHANGE MATRICES ---------------------------------
# Use subjects with ses-2 FC matrices
conmat_dir2 <- paste0(repo, '/data/conmats/', prepro, '/', parcellation, '/ses-2')

conmat_list2 <- list.files(conmat_dir2,
                           pattern = "*ses-2.txt",
                           recursive = TRUE,
                           full.names = T)

# Make output directory
out_dir <- paste0(repo, '/data/conmats/', prepro, '/', parcellation, '/change/')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Loop through subjects
for (s in conmat_list2) {
  
  subtitle = substr(s,nchar(conmat_dir2)+22,nchar(conmat_dir2)+28)  # e.g: 'sub-036'
  
  # Read in baseline and 3 month FC matrices
  fc_bl <- read.table(str_replace(str_replace(s, 'ses-2', 'ses-1'), 'ses-2', 'ses-1'))
  fc_3m <- read.table(s)
  
  # Calculate change matrix
  fc_change <- fc_3m - fc_bl
  
  # Save change matrix
  write.table(round(fc_change,3),
              file = paste0(out_dir, 'conmat_', subtitle, '_change.txt'),
              row.names = F,
              col.names = F,
              quote = F)
}