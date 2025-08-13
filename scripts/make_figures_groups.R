library(ggplot2)
library(ggpubr)
library(tidyverse)

# --------------------------------- SETUP ---------------------------------
repo <- 'path/to/psychosis-FC-prediction'

# Specify finer details of models
p_thresh <- '0.01'                   # supps also used p<0.05 and p<0.001
preproc <- 'dt_AROMA_8Phys-4GMR_bpf' # supps also used 'dt_AROMA_8Phys_bpf'
parc <- '419_parc'                   # supps also used '328_parc'
# -------------------------------------------------------------------------

# Initialise data frames to store stats for each algorithm
cpm_stats  <- data.frame('model' = character(), 'r_mean' = numeric(), 'p' = numeric(), 'p_FWE' = numeric())
krr_stats  <- data.frame('model' = character(), 'r_mean' = numeric(), 'p' = numeric(), 'p_FWE' = numeric())
meta_stats <- data.frame('model' = character(), 'r_mean' = numeric(), 'p' = numeric(), 'p_FWE' = numeric())

# Initialise data frames to store null r_mean values for each algorithm
cpm_nulls  <- data.frame(matrix(ncol = 0, nrow = 1000))
krr_nulls  <- data.frame(matrix(ncol = 0, nrow = 1000))
meta_nulls <- data.frame(matrix(ncol = 0, nrow = 1000))

# Calculate stats and make/save plots of prediction performance and significance
for (group in c('_placebo', '_medication')) {
  for (fc_type in c('fc_baseline', 'fc_change_bl3m')) {

    # ------------------------------- CONNECTOME-BASED PREDICTIVE MODELLING -------------------------------
    # Load results
    r_vals_sofas <- read.table(paste0(repo, '/results/cpm/sofas/4/', fc_type, '/r_vals_',
                                      preproc, '_', parc, '_', p_thresh, group, '.txt'))
    r_vals_bprs  <- read.table(paste0(repo, '/results/cpm/bprs/4/', fc_type, '/r_vals_',
                                      preproc, '_', parc, '_', p_thresh, group, '.txt'))
    
    # Plot prediction performancer
    r_vals <- cbind(_vals_sofas, r_vals_bprs) %>%
      setNames(c('pos_Functioning', 'neg_Functioning', 'pos_Symptoms', 'neg_Symptoms')) %>%
      gather(key = 'model', value = 'value') %>%
      mutate(predtype = ifelse(grepl('pos', model), 'pos', 'neg')) %>%
      mutate(model = gsub('pos_|neg_', '', model)) %>%
      mutate(model = factor(model, levels = c('Functioning', 'Symptoms')))
    
    fig <- ggplot(data = r_vals, aes(x = model, y = value,)) +
      geom_point(aes(color = predtype, group = predtype),
                 position = position_jitterdodge(
                   jitter.width = 0.2,
                   dodge.width  = 0.6),
                 size           = 1.5,
                 alpha          = 0.3) +
      geom_boxplot(aes(fill=predtype),
                   outlier.shape = NA,
                   width=0.6,
                   alpha= 0,
                   fatten=1) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.8, linetype = "dashed") +
      scale_color_brewer(palette = "Set1") +
      theme_classic2() +
      ylim(-1, 1) +
      scale_x_discrete(expand = c(0.1, 0)) +
      theme(axis.text=element_text(size=18),
            axis.title=element_blank(),
            legend.position = 'none')
    
    # Save plot
    dir.create(paste0(repo, '/results/figures/cpm/groups/', fc_type),
               recursive = TRUE,
               showWarnings = FALSE)    
    ggsave(paste0(repo, '/results/figures/cpm/groups/', fc_type, '/rvals_change-scores_',
                  preproc, '_', parc, '_', p_thresh, group, '.jpg'),
           plot = fig, width = 7, height = 7)
    
    # Calculate & plot significance if the null models have been run
    if (file.exists(paste0(repo, '/results/cpm/sofas/4/', fc_type,
                           '/null_r_means_', preproc, '_', parc, '_', p_thresh, group, '.txt')) &&
        file.exists(paste0(repo, '/results/cpm/bprs/4/', fc_type,
                           '/null_r_means_', preproc, '_', parc, '_', p_thresh, group, '.txt'))
    ) {
      
      # Grab mean r value for each of the 1000 permutations
      null_r_means_sofas  <- read.table(paste0(repo, '/results/cpm/sofas/4/', fc_type, '/null_r_means_',
                                               preproc, '_', parc, '_', p_thresh, group, '.txt')) %>%
        na.omit() %>%
        head(1000) %>%
        select(-1)
      
      null_r_means_bprs <- read.table(paste0(repo, '/results/cpm/bprs/4/', fc_type, '/null_r_means_',
                                             preproc, '_', parc, '_', p_thresh, group, '.txt')) %>%
        na.omit() %>%
        head(1000) %>%
        select(-1)
      
      # Calculate p values as proportion of null r values greater than observed r values
      p_sofas_pos <- sum(mean(r_vals_sofas$V1) < null_r_means_sofas$V2) / length(null_r_means_sofas$V2)
      p_sofas_neg <- sum(mean(r_vals_sofas$V2) < null_r_means_sofas$V3) / length(null_r_means_sofas$V3)
      p_bprs_pos  <- sum(mean(r_vals_bprs$V1)  < null_r_means_bprs$V2)  / length(null_r_means_bprs$V2)
      p_bprs_neg  <- sum(mean(r_vals_bprs$V2)  < null_r_means_bprs$V3)  / length(null_r_means_bprs$V3)
      
      cpm_nulls <- cbind(cpm_nulls, null_r_means_sofas, null_r_means_bprs)
      
    } else {
      p_sofas_pos <- p_sofas_neg <- p_bprs_pos <- p_bprs_neg <- NA
    }
    
    # Add to stats
    cpm_stats <- cpm_stats %>%
      add_row(model = paste0(fc_type, group, '_SOFAS_6mo_pos'),
              r_mean = mean(r_vals_sofas$V1),
              p = p_sofas_pos) %>%
      add_row(model = paste0(fc_type, group, '_SOFAS_6mo_neg'),
              r_mean = mean(r_vals_sofas$V2),
              p = p_sofas_neg) %>%
      add_row(model = paste0(fc_type, group, '_BPRS_6mo_pos'),
              r_mean = mean(r_vals_bprs$V1),
              p = p_bprs_pos) %>%
      add_row(model = paste0(fc_type, group, '_BPRS_6mo_neg'),
              r_mean = mean(r_vals_bprs$V2),
              p = p_bprs_neg)
    
    # ---------------------------------- KERNEL RIDGE REGRESSION ----------------------------------
    # Load results
    r_vals_sofas  <- read.table(paste0(repo, '/results/krr/sofas/4/', fc_type, '/r_vals_',
                                       preproc, '_', parc, group, '.txt'))
    r_vals_bprs <- read.table(paste0(repo, '/results/krr/bprs/4/', fc_type, '/r_vals_',
                                     preproc, '_', parc, group, '.txt'))
    
    # Plot prediction performance
    r_vals <- cbind(r_vals_sofas, r_vals_bprs) %>%
      setNames(c('Functioning', 'Symptoms')) %>%
      gather(key = 'timepoint', value = 'value') %>%
      mutate(timepoint = factor(timepoint, levels = c('Functioning', 'Symptoms')))
    
    fig <- ggplot(data = r_vals, aes(x = timepoint, y = value)) +
      geom_point(aes(fill = timepoint),
                 position = position_jitterdodge(
                   jitter.width = 0.7,
                   dodge.width  = 0.3),
                 size           = 1.5,
                 alpha          = 0.3,
                 color = 'darkgreen') +
      geom_boxplot(outlier.shape = NA,
                   width=0.45,
                   alpha= 0,
                   fatten=1) +
      theme_classic2() +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.8, linetype = "dashed") +
      ylim(-1, 1) +
      scale_x_discrete(expand = c(0.1, 0)) +
      theme(axis.text=element_text(size=18),
            axis.title=element_blank(),
            legend.position = 'none')
    
    # Save plot
    dir.create(paste0(repo, '/results/figures/krr/groups/', fc_type),
               recursive = TRUE,
               showWarnings = FALSE)
    ggsave(paste0(repo, '/results/figures/krr/groups/', fc_type, '/rvals_change-scores_',
                  preproc, '_', parc, group, '.jpg'),
           plot = fig, width = 4, height = 7)
    
    # Calculate & plot significance if the null models have been run
    if (file.exists(paste0(repo, '/results/krr/sofas/4/', fc_type,
                           '/null_r_means_', preproc, '_', parc, group, '.txt')) &&
        file.exists(paste0(repo, '/results/krr/bprs/4/', fc_type,
                           '/null_r_means_', preproc, '_', parc, group, '.txt'))) {
      
      # Grab mean r value for each of the 1000 permutations
      null_r_means_sofas  <- read.table(paste0(repo, '/results/krr/sofas/4/', fc_type,
                                             '/null_r_means_', preproc, '_', parc, group, '.txt'))
      null_r_means_bprs <- read.table(paste0(repo, '/results/krr/bprs/4/', fc_type,
                                             '/null_r_means_', preproc, '_', parc, group, '.txt'))
      
      # Calculate p values as proportion of null r values greater than observed r values
      p_sofas  <- sum(mean(r_vals_sofas$V1)  < null_r_means_sofas$V1)  / dim(null_r_means_sofas)[1]
      p_bprs <- sum(mean(r_vals_bprs$V1) < null_r_means_bprs$V1) / dim(null_r_means_bprs)[1]
      
      krr_nulls <- cbind(krr_nulls, null_r_means_sofas, null_r_means_bprs)
      
    } else {
      p_sofas  <- p_bprs <- NA
    }
    
    # Add to stats
    krr_stats <- krr_stats %>%
      add_row(model = paste0(fc_type, group, '_SOFAS_6mo'),
              r_mean = mean(r_vals_sofas$V1),
              p = p_sofas) %>%
      add_row(model = paste0(fc_type, group, '_BPRS_6mo'),
              r_mean = mean(r_vals_bprs$V1),
              p = p_bprs)
  }
  
  # --------------------------------------------- META-MATCHING ---------------------------------------------  
  # Load results
  r_vals_sofas  <- read.table(paste0(repo, '/results/meta-matching/sofas/4/r_vals_',
                                     preproc, group, '.txt'))
  r_vals_bprs <- read.table(paste0(repo, '/results/meta-matching/bprs/4/r_vals_',
                                   preproc, group, '.txt'))
  
  # Plot prediction performance
  r_vals <- cbind(r_vals_sofas, r_vals_bprs) %>%
    setNames(c('Functioning', 'Symptoms')) %>%
    gather(key = 'timepoint', value = 'value') %>%
    mutate(timepoint = factor(timepoint, levels = c('Functioning', 'Symptoms')))
  
  fig <- ggplot(data = r_vals, aes(x = timepoint, y = value)) +
    geom_point(aes(fill = timepoint),
               position = position_jitterdodge(
                 jitter.width = 0.7,
                 dodge.width  = 0.3),
               size           = 1.5,
               alpha          = 0.3,
               color = 'purple3') +
    geom_boxplot(outlier.shape = NA,
                 width=0.45,
                 alpha= 0,
                 fatten=1) +
    theme_classic2() +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8, linetype = "dashed") +
    ylim(-1, 1) +
    scale_x_discrete(expand = c(0.1, 0)) +
    theme(axis.text=element_text(size=18),
          axis.title=element_blank(),
          legend.position = 'none')
  
  # Save plot
  dir.create(paste0(repo, '/results/figures/meta-matching/groups'),
             recursive = TRUE,
             showWarnings = FALSE)
  ggsave(paste0(repo, '/results/figures/meta-matching/groups/rvals_change-scores_', preproc, group, '.jpg'),
         plot = fig, width = 4, height = 7)
  
  # Calculate & plot significance if the null models have been run
  if (file.exists(paste0(repo, '/results/meta-matching/sofas/4/null_r_means_', preproc, group, '.txt')) &&
      file.exists(paste0(repo, '/results/meta-matching/bprs/4/null_r_means_', preproc, group, '.txt'))) {
    
    # Grab mean r value for each of the 1000 permutations
    null_r_means_sofas  <- read.table(paste0(repo, '/results/meta-matching/sofas/4/null_r_means_',
                                             preproc, group, '.txt')) %>%
      select(-1)
    null_r_means_bprs <- read.table(paste0(repo, '/results/meta-matching/bprs/4/null_r_means_',
                                           preproc, group, '.txt')) %>%
      select(-1)
    
    # Calculate p values as proportion of null r values greater than observed r values
    p_sofas  <- sum(mean(r_vals_sofas$V1)  < null_r_means_sofas$V2)  / dim(null_r_means_sofas)[1]
    p_bprs <- sum(mean(r_vals_bprs$V1) < null_r_means_bprs$V2) / dim(null_r_means_bprs)[1]
    
    meta_nulls <- cbind(meta_nulls, null_r_means_sofas, null_r_means_bprs)
    
  } else {
    p_sofas <- p_bprs <- NA
  }
  
  # Add to stats
  meta_stats <- meta_stats %>%
    add_row(model = paste0('fc_baseline', group, '_SOFAS_6mo'),
            r_mean = mean(r_vals_sofas$V1),
            p = p_sofas) %>%
    add_row(model = paste0('fc_baseline', group, '_BPRS_6mo'),
            r_mean = mean(r_vals_bprs$V1),
            p = p_bprs)
}

# ------------------------------------------------------- FAMILY-WISE ERROR CORRECTION -------------------------------------------------------
# Apply the Westfall-Young method: retain the highest r_mean across the 16 CPM models / 8 KRR models / 4 meta-matching models
cpm_nulls  <- apply(cpm_nulls, 1, max)
krr_nulls  <- apply(krr_nulls, 1, max)
meta_nulls  <- apply(meta_nulls, 1, max)

# Initialise looping values for finding r_mean in stats dataframes
c <- 0
k <- 0
m <- 0

for (group in c('_placebo', '_medication')) {
  for (fc_type in c('fc_baseline', 'fc_change_bl3m')) {
    # ------------------------------------------------- CONNECTOME-BASED PREDICTIVE MODELLING -------------------------------------------------
    # Recalculate p-values
    cpm_stats[c+1,4] <- sum(cpm_stats[c+1,2] < cpm_nulls) / length(cpm_nulls)
    cpm_stats[c+2,4] <- sum(cpm_stats[c+2,2] < cpm_nulls) / length(cpm_nulls)
    cpm_stats[c+3,4] <- sum(cpm_stats[c+3,2] < cpm_nulls) / length(cpm_nulls)
    cpm_stats[c+4,4] <- sum(cpm_stats[c+4,2] < cpm_nulls) / length(cpm_nulls)
    
    # Plot model performances against the FWE-corrected null distribution
    vline_data <- data.frame(xintercept = c(cpm_stats[c+1,2], cpm_stats[c+2,2], cpm_stats[c+3,2], cpm_stats[c+4,2]),
                             model = c('6mo_pos', '6mo_neg', '12mo_pos', '12mo_neg'))
    
    fig <- ggplot(data = data.frame(cpm_nulls = cpm_nulls), aes(x = cpm_nulls)) +
      geom_histogram(binwidth = 0.025, fill = 'grey', color='black', alpha = 0.5) +
      geom_vline(data = vline_data, aes(xintercept = xintercept, color = model, linetype = model), linewidth = 1.5) +
      scale_color_manual(values = c('6mo_pos' = '#0000b3', '6mo_neg' = '#b30000', '12mo_pos' = '#0000b3', '12mo_neg' = '#b30000')) +
      scale_linetype_manual(values = c('6mo_pos' = 'solid', '6mo_neg' = 'solid', '12mo_pos' = 'dashed', '12mo_neg' = 'dashed')) +
      theme_minimal() +
      ylim(0, 90) +
      xlim(-0.7, 0.7) +
      theme(axis.text=element_text(size=16),
            axis.title=element_blank(),
            legend.position = 'none')
    
    # Save plot
    ggsave(paste0(repo, '/results/figures/cpm/groups/', fc_type, '/nulls_', preproc, '_', parc, '_', p_thresh, group, '.jpg'),
           plot = fig, width = 8, height = 6)
    
    # ----------------------------------------------- KERNEL RIDGE REGRESSION -----------------------------------------------
    # Recalculate p-values
    krr_stats[k+1,4] <- sum(krr_stats[k+1,2] < krr_nulls) / length(krr_nulls)
    krr_stats[k+2,4] <- sum(krr_stats[k+2,2] < krr_nulls) / length(krr_nulls)
    
    # Plot model performances against the FWE-corrected null distribution
    vline_data <- data.frame(xintercept = c(krr_stats[k+1,2], krr_stats[k+2,2]),
                             model = c('6mo', '12mo'))
    
    fig <- ggplot(data = data.frame(krr_nulls = krr_nulls), aes(x = krr_nulls)) +
      geom_histogram(binwidth = 0.025, fill = 'grey', color='black', alpha = 0.5) +
      geom_vline(data = vline_data, aes(xintercept = xintercept, linetype = model), linewidth = 1.5, color = 'darkgreen') +
      scale_linetype_manual(values = c('6mo' = 'solid', '12mo' = 'dashed')) +
      theme_minimal() +
      ylim(0, 90) +
      xlim(-0.6, 0.8) +
      theme(axis.text=element_text(size=16),
            axis.title=element_blank(),
            legend.position = 'none')
    
    # Save plot
    ggsave(paste0(repo, '/results/figures/krr/groups/', fc_type, '/nulls_', preproc, '_', parc, group, '.jpg'),
           plot = fig, width = 8, height = 6)
    
    # Bump up looping values
    c <- c + 4
    k <- k + 2
  }
  
  # --------------------------------------------------- META-MATCHING ---------------------------------------------------
  # Recalculate p-values
  meta_stats[m+1,4] <- sum(meta_stats[m+1,2] < meta_nulls) / length(meta_nulls)
  meta_stats[m+2,4] <- sum(meta_stats[m+2,2] < meta_nulls) / length(meta_nulls)
  
  # Plot model performances against the FWE-corrected null distribution
  vline_data <- data.frame(xintercept = c(meta_stats[m+1,2], meta_stats[m+2,2]),
                           model = c('6mo', '12mo'))
  
  fig <- ggplot(data = data.frame(meta_nulls = meta_nulls), aes(x = meta_nulls)) +
    geom_histogram(binwidth = 0.025, fill = 'grey', color='black', alpha = 0.5) +
    geom_vline(data = vline_data, aes(xintercept = xintercept, linetype = model), linewidth = 1.5, color = 'purple3') +
    scale_linetype_manual(values = c('6mo' = 'solid', '12mo' = 'dashed')) +
    theme_minimal() +
    ylim(0, 70) +
    xlim(-0.25, 0.85) +
    theme(axis.text=element_text(size=16),
          axis.title=element_blank(),
          legend.position = 'none')
  
  # Save plot
  ggsave(paste0(repo, '/results/figures/meta-matching/groups/nulls_', preproc, group, '.jpg'),
         plot = fig, width = 8, height = 6)
  
  # Bump up looping value
  m = m + 2
  }

# For any NA p-values, set p_FWE to NA
cpm_stats$p_FWE[is.na(cpm_stats$p)] <- NA
krr_stats$p_FWE[is.na(krr_stats$p)] <- NA
meta_stats$p_FWE[is.na(meta_stats$p)] <- NA

# Save stats
dir.create(paste0(repo, '/results/stats'), recursive = T, showWarnings = F)
write.csv(cpm_stats,  paste0(repo, '/results/stats/cpm_change-scores_', preproc, '_', parc, '_', p_thresh, '_groups.csv'), row.names=FALSE)
write.csv(krr_stats,  paste0(repo, '/results/stats/krr_change-scores_', preproc, '_', parc, '_groups.csv'), row.names=FALSE)
write.csv(meta_stats,  paste0(repo, '/results/stats/meta-matching_change-scores_', preproc, '_groups.csv'), row.names=FALSE)