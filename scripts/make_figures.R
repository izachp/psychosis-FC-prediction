library(ggplot2)
library(ggpubr)
library(tidyverse)

# --------------------------------- SETUP ---------------------------------
repo <- 'path/to/psychosis-FC-prediction'

# Specify finer details of models
p_thresh <- '0.01'                   # supps also used p<0.05 and p<0.001
preproc <- 'dt_AROMA_8Phys-4GMR_bpf' # supps also used 'dt_AROMA_8Phys_bpf'
parc <- '419_parc'                   # supps also used '328_parc'
scales <- c('sofas', 'bprs')         # supps also used 'bprs_pos'
# For the placebo/medication and slope outcome supplementary analyses,
#   use `make_figures_groups.R` and `make_figures_lmslopes.R`
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
for (scale in scales) {
  for (fc_type in c('fc_baseline', 'fc_change_bl3m')) {
    
    # ---------------------------------- CONNECTOME-BASED PREDICTIVE MODELLING ----------------------------------
    # Load results
    r_vals_6mo  <- read.table(paste0(repo, '/results/cpm/', scale, '/4/', fc_type, '/r_vals_',
                                     preproc, '_', parc, '_', p_thresh, '.txt'))
    
    r_vals_12mo <- read.table(paste0(repo, '/results/cpm/', scale, '/5/', fc_type, '/r_vals_',
                                     preproc, '_', parc, '_', p_thresh, '.txt'))
    
    # Plot prediction performance
    r_vals <- cbind(r_vals_6mo, r_vals_12mo) %>%
      setNames(c('pos_6 months', 'neg_6 months', 'pos_12 months', 'neg_12 months')) %>%
      gather(key = 'model', value = 'value') %>%
      mutate(predtype = ifelse(grepl('pos', model), 'pos', 'neg')) %>%
      mutate(model = gsub('pos_|neg_', '', model)) %>%
      mutate(model = factor(model, levels = c('6 months', '12 months')))
    
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
      ylim(-0.8, 0.8) +
      scale_x_discrete(expand = c(0.1, 0)) +
      theme(axis.text=element_text(size=18),
            axis.title=element_blank(),
            legend.position = 'none')
    
    # Save plot
    dir.create(paste0(repo, '/results/figures/cpm/', scale, '/', fc_type),
               recursive = T, showWarnings = F)
    ggsave(paste0(repo, '/results/figures/cpm/', scale, '/', fc_type,
                  '/rvals_', preproc, '_', parc, '_', p_thresh,  '.jpg'),
           plot = fig, width = 7, height = 7)
    
    # Calculate & plot significance if the null models have been run
    if (file.exists(paste0(repo, '/results/cpm/', scale, '/4/', fc_type,
                                      '/null_r_means_', preproc, '_', parc, '_', p_thresh, '.txt')) &&
        file.exists(paste0(repo, '/results/cpm/', scale, '/5/', fc_type,
                           '/null_r_means_', preproc, '_', parc, '_', p_thresh, '.txt'))
        ) {
      
      # Grab mean r value for each of the 1000 permutations
      null_r_means_6mo  <- read.table(paste0(repo, '/results/cpm/', scale, '/4/', fc_type,
                                             '/null_r_means_', preproc, '_', parc, '_', p_thresh, '.txt')) %>%
        na.omit() %>%
        head(1000) %>%
        select(-1)
      
      null_r_means_12mo <- read.table(paste0(repo, '/results/cpm/', scale, '/5/', fc_type,
                                             '/null_r_means_', preproc, '_', parc, '_', p_thresh, '.txt')) %>%
        na.omit() %>%
        head(1000) %>%
        select(-1)
  
      # Calculate p values as proportion of null r values greater than observed r values
      p_6mo_pos  <- sum(mean(r_vals_6mo$V1)  < null_r_means_6mo$V2)  / length(null_r_means_6mo$V2)
      p_6mo_neg  <- sum(mean(r_vals_6mo$V2)  < null_r_means_6mo$V3)  / length(null_r_means_6mo$V3)
      p_12mo_pos <- sum(mean(r_vals_12mo$V1) < null_r_means_12mo$V2) / length(null_r_means_12mo$V2)
      p_12mo_neg <- sum(mean(r_vals_12mo$V2) < null_r_means_12mo$V3) / length(null_r_means_12mo$V3)
      
      cpm_nulls <- cbind(cpm_nulls, null_r_means_6mo, null_r_means_12mo)
    
    } else {
      p_6mo_pos <- p_6mo_neg <- p_12mo_pos <- p_12mo_neg <- NA
    }
    
    # Add to stats
    cpm_stats <- cpm_stats %>%
      add_row(model = paste0(fc_type, '_', scale, '_6mo_pos'),
              r_mean = mean(r_vals_6mo$V1),
              p = p_6mo_pos) %>%
      add_row(model = paste0(fc_type, '_', scale, '_6mo_neg'),
              r_mean = mean(r_vals_6mo$V2),
              p = p_6mo_neg) %>%
      add_row(model = paste0(fc_type, '_', scale, '_12mo_pos'),
              r_mean = mean(r_vals_12mo$V1),
              p = p_12mo_pos) %>%
      add_row(model = paste0(fc_type, '_', scale, '_12mo_neg'),
              r_mean = mean(r_vals_12mo$V2),
              p = p_12mo_neg)
    
    # ---------------------------------- KERNEL RIDGE REGRESSION ----------------------------------    
    # Load results
    r_vals_6mo  <- read.table(paste0(repo, '/results/krr/', scale, '/4/', fc_type, '/r_vals_',
                                     preproc, '_', parc, '.txt'))
    r_vals_12mo <- read.table(paste0(repo, '/results/krr/', scale, '/5/', fc_type, '/r_vals_',
                                     preproc, '_', parc, '.txt'))
  
    # Plot prediction performance
    r_vals <- cbind(r_vals_6mo, r_vals_12mo) %>%
      setNames(c('6 months', '12 months')) %>%
      gather(key = 'timepoint', value = 'value') %>%
      mutate(timepoint = factor(timepoint, levels = c('6 months', '12 months')))
    
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
      ylim(-0.8, 0.8) +
      scale_x_discrete(expand = c(0.1, 0)) +
      theme(axis.text=element_text(size=18),
            axis.title=element_blank(),
            legend.position = 'none')
    
    # Save plot
    dir.create(paste0(repo, '/results/figures/krr/', scale, '/', fc_type),
               recursive = T, showWarnings = F)
    ggsave(paste0(repo, '/results/figures/krr/', scale, '/', fc_type,
                  '/rvals_', preproc, '_', parc, '.jpg'),
           plot = fig, width = 4, height = 7)
    
    # Calculate & plot significance if the null models have been run
    if (file.exists(paste0(repo, '/results/krr/', scale, '/4/', fc_type,
                           '/null_r_means_', preproc, '_', parc, '.txt')) &&
        file.exists(paste0(repo, '/results/krr/', scale, '/5/', fc_type,
                           '/null_r_means_', preproc, '_', parc, '.txt'))) {
    
      # Grab mean r value for each of the 1000 permutations
      null_r_means_6mo  <- read.table(paste0(repo, '/results/krr/', scale, '/4/', fc_type,
                                             '/null_r_means_', preproc, '_', parc, '.txt'))
      null_r_means_12mo <- read.table(paste0(repo, '/results/krr/', scale, '/5/', fc_type,
                                             '/null_r_means_', preproc, '_', parc, '.txt'))
      
      # Calculate p values as proportion of null r values greater than observed r values
      p_6mo  <- sum(mean(r_vals_6mo$V1)  < null_r_means_6mo$V1)  / dim(null_r_means_6mo)[1]
      p_12mo <- sum(mean(r_vals_12mo$V1) < null_r_means_12mo$V1) / dim(null_r_means_12mo)[1]
      
      krr_nulls <- cbind(krr_nulls, null_r_means_6mo, null_r_means_12mo)
    
    } else {
      p_6mo  <- p_12mo <- NA
    }
    
    # Add to stats
    krr_stats <- krr_stats %>%
      add_row(model = paste0(fc_type, '_', scale, '_6mo'),
              r_mean = mean(r_vals_6mo$V1),
              p = p_6mo) %>%
      add_row(model = paste0(fc_type, '_', scale, '_12mo'),
              r_mean = mean(r_vals_12mo$V1),
              p = p_12mo)
  }
  
  # --------------------------------------------- META-MATCHING ---------------------------------------------
  # Load results
  r_vals_6mo  <- read.table(paste0(repo, '/results/meta-matching/', scale, '/4/r_vals_',
                                   preproc, '.txt'))
  r_vals_12mo <- read.table(paste0(repo, '/results/meta-matching/', scale, '/5/r_vals_',
                                   preproc, '.txt'))
  
  # Plot prediction performance
  r_vals <- cbind(r_vals_6mo, r_vals_12mo) %>%
    setNames(c('6 months', '12 months')) %>%
    gather(key = 'timepoint', value = 'value') %>%
    mutate(timepoint = factor(timepoint, levels = c('6 months', '12 months')))
  
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
    ylim(-0.8, 0.8) +
    scale_x_discrete(expand = c(0.1, 0)) +
    theme(axis.text=element_text(size=18),
          axis.title=element_blank(),
          legend.position = 'none')
  
  # Save plot
  dir.create(paste0(repo, '/results/figures/meta-matching/', scale),
             recursive = T, showWarnings = F)
  ggsave(paste0(repo, '/results/figures/meta-matching/', scale, '/rvals_', preproc, '.jpg'),
         plot = fig, width = 4, height = 7)
  
  # Calculate & plot significance if the null models have been run
  if (file.exists(paste0(repo, '/results/meta-matching/', scale, '/4/null_r_means_', preproc, '.txt')) &&
      file.exists(paste0(repo, '/results/meta-matching/', scale, '/5/null_r_means_', preproc, '.txt'))) {
  
    # Grab mean r value for each of the 1000 permutations
    null_r_means_6mo  <- read.table(paste0(repo, '/results/meta-matching/', scale, '/4/null_r_means_',
                                           preproc, '.txt')) %>%
      select(-1)
    null_r_means_12mo <- read.table(paste0(repo, '/results/meta-matching/', scale, '/5/null_r_means_',
                                           preproc, '.txt')) %>%
      select(-1)

    # Calculate p values as proportion of null r values greater than observed r values
    p_6mo  <- sum(mean(r_vals_6mo$V1)  < null_r_means_6mo$V2)  / dim(null_r_means_6mo)[1]
    p_12mo <- sum(mean(r_vals_12mo$V1) < null_r_means_12mo$V2) / dim(null_r_means_12mo)[1]
    
    meta_nulls <- cbind(meta_nulls, null_r_means_6mo, null_r_means_12mo)
  
  } else {
    p_6mo <- p_12mo <- NA
  }
  
  # Add to stats
  meta_stats <- meta_stats %>%
    add_row(model = paste0('fc_baseline_', scale, '_6mo'),
            r_mean = mean(r_vals_6mo$V1),
            p = p_6mo) %>%
    add_row(model = paste0('fc_baseline_', scale, '_12mo'),
            r_mean = mean(r_vals_12mo$V1),
            p = p_12mo)
  
  # -------------------------- PLOT OBSERVED OUTCOMES --------------------------
  for (timepoint in c('4', '5')) {
    # Load data
    obs_fcbl <- read.csv(paste0(repo, '/results/cpm/', scale, '/', timepoint,
                                '/fc_baseline/outcomes_observed.csv'))

    obs_fcbl3m <- read.csv(paste0(repo, '/results/cpm/', scale, '/', timepoint,
                                  '/fc_change_bl3m/outcomes_observed.csv'))

    # Plot distributions
    fig <- ggplot() +
      geom_histogram(data = obs_fcbl, aes(x = outcome),
                     binwidth = 0.05,
                     fill = 'mediumpurple1',
                     color='black',
                     alpha = 0.5) +
      geom_histogram(data = obs_fcbl3m, aes(x = outcome),
                     binwidth = 0.05,
                     fill = 'yellow3',
                     color='black',
                     alpha = 0.5) +
      theme_minimal() +
      xlim(-1.2, 1.2) +
      ylim(0, 8) +
      theme(axis.text=element_text(size=16),
            axis.title=element_blank())

    # Save plot
    dir.create(paste0(repo, '/results/figures/observed', scale),
               recursive = T, showWarnings = F)
    ggsave(paste0(repo, '/results/figures/observed/', scale,
                  '/outcomes_observed_', timepoint, '.jpg'),
           plot = fig, width = 8, height = 6)
  }
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

for (scale in scales) {
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
      ylim(0, 120) +
      xlim(-0.5, 0.75) +
      theme(axis.text=element_text(size=16),
            axis.title=element_blank(),
            legend.position = 'none')
    
    # Save plot
    ggsave(paste0(repo, '/results/figures/cpm/', scale, '/', fc_type, '/nulls_', preproc, '_', parc, '_', p_thresh, '.jpg'),
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
      ylim(0, 120) +
      xlim(-0.5, 0.75) +
      theme(axis.text=element_text(size=16),
            axis.title=element_blank(),
            legend.position = 'none')
    
    # Save plot
    ggsave(paste0(repo, '/results/figures/krr/', scale, '/', fc_type, '/nulls_', preproc, '_', parc, '.jpg'),
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
    ylim(0, 120) +
    xlim(-0.5, 0.75) +
    theme(axis.text=element_text(size=16),
          axis.title=element_blank(),
          legend.position = 'none')
  
  # Save plot
  ggsave(paste0(repo, '/results/figures/meta-matching/', scale, '/nulls_', preproc, '.jpg'),
         plot = fig, width = 8, height = 6)
  
  # Bump up looping value
  m <- m + 2
  
  # Test for group differences between dropouts and 3-month completers
  for (timepoint in c('4', '5')) {
    obs_fcbl3m <- read.csv(paste0(repo, '/results/cpm/', scale, '/', timepoint,
                                  '/fc_change_bl3m/outcomes_observed.csv'))
    
    completers <- obs_fcbl3m %>%
      select(1) %>%
      unlist()
    
    obs_fcbl3m <- obs_fcbl3m %>%
      select(2) %>%
      unlist()
    
    obs_dropouts <- read.csv(paste0(repo, '/results/cpm/', scale, '/', timepoint,
                                    '/fc_baseline/outcomes_observed.csv')) %>%
      filter(!BIDS_ID %in% completers) %>%  # remove completers
      select(2) %>%
      unlist()
    
    # Perform t-test
    ttest <- t.test(obs_fcbl3m, obs_dropouts)
    print(paste0('t-test comparing ', scale, ' change scores at timepoint ', timepoint,
                 ' between patients with usable 3-month FC data and those without: p = ',
                 ttest[3]))
  }
}

# For any NA p-values, set p_FWE to NA
cpm_stats$p_FWE[is.na(cpm_stats$p)] <- NA
krr_stats$p_FWE[is.na(krr_stats$p)] <- NA
meta_stats$p_FWE[is.na(meta_stats$p)] <- NA

# Save stats
dir.create(paste0(repo, '/results/stats'), recursive = T, showWarnings = F)
write.csv(cpm_stats,  paste0(repo, '/results/stats/cpm_change-scores_', preproc, '_', parc, '_', p_thresh, '.csv'), row.names=FALSE)
write.csv(krr_stats,  paste0(repo, '/results/stats/krr_change-scores_', preproc, '_', parc, '.csv'), row.names=FALSE)
write.csv(meta_stats,  paste0(repo, '/results/stats/meta-matching_change-scores_', preproc, '.csv'), row.names=FALSE)
