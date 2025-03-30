library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(stringr)

# Set repo directory
repo <- '/path/to/repo'

# Set preprocessing pipeline
prepro <- 'dt_AROMA_8Phys-4GMR_bpf'

breaksList <- seq(-1, 1, by = 0.1)  # Set upper and lower limits of plot
color <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(length(breaksList)) # Defines the vector of colors for the legend (it has to be of the same length of breaksList)

# Read in list of low-signal regions to be excluded from 332-region parcellation
excluded_regions <- read.table(paste0(repo, '/low_signal_regions.txt'), header = F) %>%
  select(2) %>%
  unlist()

# Loop through baseline and 3-month data for 332-region parcellation and baseline data for 419-region parcellation
for (ses in c('1', '2', '1_419-reg')) {
  
  # Find and loop through BOLD time series
  ts_list <- list.files(paste0(repo, '/timeseries/', prepro, '/ses-', ses), 
                        full.names = T)
  
  # Make folder to store FC matrices
  dir.create(paste0(repo, '/conmats/', prepro, '/ses-', ses, '/png'), recursive = T, showWarnings = T)
  
  for (s in ts_list) {
    
    # Read in timeseries and exclude low-signal regions
    ts <- as.matrix(read.table(s), header = F) %>%
    .[-excluded_regions,]
    
    # Check dimensions are correct for both parcellations and 234 timepoints
    if (ses == '1_419-reg') {
      if (dim(ts)[1] != 419 | dim(ts)[2] != 234) {
        print(paste0('Error: ', str_sub(s,-17,-5), ' has dimensions ', dim(ts)))
        next
      }
    } else {
    
      if (dim(ts)[1] != 332-length(excluded_regions) | dim(ts)[2] != 234) {
        print(paste0('Error: ', str_sub(s,-17,-5), ' has dimensions ', dim(ts)))
        next
      }
    }
    
    # Calculate FC matrix by Pearson's correlation
    conmat <- cor(t(ts))
    
    # Save FC matrix values and visualisation
    write.table(round(conmat,3), 
                file = paste0(repo, '/conmats/', prepro, '/ses-', ses, '/conmat_', str_sub(s,-17,-5), '.txt'), 
                row.names = F, 
                col.names = F,
                quote = F)
    
    png(file = paste0(repo, '/conmats/', prepro, '/ses-', ses, '/png/conmat_', str_sub(s,-17,-5), '.png'))
    print(pheatmap::pheatmap(conmat, treeheight_row = 0, treeheight_col = 0,
                             main = str_sub(s,-17,-5), 
                             cluster_rows=FALSE, 
                             cluster_cols=FALSE, 
                             breaks = breaksList, 
                             color = color))
    dev.off()
  }
}

## Now make matrices for 3-month change in FC (332-region parcellation) ##

# Make folder to store change in FC matrices
dir.create(paste0(repo, '/conmats/', prepro, '/change/png'), recursive = T, showWarnings = T)

# Find and loop through 3-month FC matrices
conmat_list_ses2 <- list.files(paste0(repo, '/conmats/', prepro, '/ses-2'),
                           pattern = ".txt",
                           full.names = T)

for (s in conmat_list_ses2) {
  
  # Read in 3 month FC matrix and corresponding baseline FC matrix
  fc_3m <- read.table(s)
  fc_bl <- read.table(str_replace(str_replace(s, 'ses-2', 'ses-1'), 'ses-2', 'ses-1'))
  
  # Create matrix for change in FC by element-wise subtraction
  fc_change <- fc_3m - fc_bl
  
  # Save change matrix
  write.table(round(fc_change,3),
              file = paste0(repo, '/conmats/', prepro, '/change/conmat_', str_sub(s,-17,-5), '_change.txt'),
              row.names = F,
              col.names = F,
              quote = F)
  
  # Plot change matrix
  png(file = paste0(repo, '/conmats/', prepro, '/change/png/conmat_', str_sub(s,-17,-5), '_change.png'))
  print(pheatmap::pheatmap(fc_change, treeheight_row = 0, treeheight_col = 0,
                           main = str_sub(s,-17,-5),
                           cluster_rows=FALSE,
                           cluster_cols=FALSE,
                           breaks = breaksList,
                           color = color))
  dev.off()
  
}