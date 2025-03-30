library(ggplot2)
library(readxl)
library(dplyr)
library(ggpubr)
library(writexl)
library(gtools)
library(tidyverse)

# Specify path to this repo
repo <- 'path/to/repo'

# Initialise data frames with column names
cpm_stats  <- data.frame('model' = character(), 'r_mean' = numeric(), 'p' = numeric(), 'p_FWE' = numeric())
krr_stats  <- data.frame('model' = character(), 'r_mean' = numeric(), 'p' = numeric(), 'p_FWE' = numeric())
meta_stats <- data.frame('model' = character(), 'r_mean' = numeric(), 'p' = numeric(), 'p_FWE' = numeric())

# Initialise null data frames with 1000 rows
cpm_nulls  <- data.frame(matrix(ncol = 0, nrow = 1000))
krr_nulls  <- data.frame(matrix(ncol = 0, nrow = 1000))
meta_nulls <- data.frame(matrix(ncol = 0, nrow = 1000))

for (scale in c('sofas', 'bprs')) {

  # Plot prediction performance & calculate uncorrected p-values
  for (fc_type in c('fc_baseline', 'fc_change_bl3m')) {
    
    dir.create(paste0(repo, '/figures/', scale, '/', fc_type), recursive = T, showWarnings = T)
    
    # ============================================= CPM =============================================
    
    # Load data
    r_vals_6mo  <- read_excel(paste0(repo, '/cpm/', scale, '/4/', fc_type, '/r_vals.xlsx'), col_names = TRUE)
    r_vals_12mo <- read_excel(paste0(repo, '/cpm/', scale, '/5/', fc_type, '/r_vals.xlsx'), col_names = TRUE)
    
    # Plot r values for 4 models
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
      ylim(-0.6, 0.6) +
      scale_x_discrete(expand = c(0.1, 0)) +
      theme(axis.text=element_text(size=16),
            axis.title=element_blank(),
            legend.position = 'none')
    
    ggsave(paste0(repo, '/figures/', scale, '/', fc_type, '/cpm_rvals.jpg'), plot = fig, width = 7, height = 7)

    # Grab the null model's mean r values
    null_r_means_6mo  <- read_excel(paste0(repo, '/cpm/', scale, '/4/', fc_type, '/null_r_means.xlsx'), col_names = TRUE) %>%
      select(-1) %>%
      na.omit() %>%
      head(1000)
    null_r_means_12mo <- read_excel(paste0(repo, '/cpm/', scale, '/5/', fc_type, '/null_r_means.xlsx'), col_names = TRUE) %>%
      select(-1) %>%
      na.omit() %>%
      head(1000)

    # Calculate p values as proportion of null r values greater than observed r values
    p_6mo_pos  <- sum(mean(r_vals_6mo$pos)  < null_r_means_6mo$pos)  / length(null_r_means_6mo$pos)
    p_6mo_neg  <- sum(mean(r_vals_6mo$neg)  < null_r_means_6mo$neg)  / length(null_r_means_6mo$neg)
    p_12mo_pos <- sum(mean(r_vals_12mo$pos) < null_r_means_12mo$pos) / length(null_r_means_12mo$pos)
    p_12mo_neg <- sum(mean(r_vals_12mo$neg) < null_r_means_12mo$neg) / length(null_r_means_12mo$neg)
    
    cpm_stats <- cpm_stats %>%
      add_row(model = paste0(fc_type, '_', scale, '_6mo_pos'),
              r_mean = mean(r_vals_6mo$pos),
              p = p_6mo_pos) %>%
      add_row(model = paste0(fc_type, '_', scale, '_6mo_neg'),
              r_mean = mean(r_vals_6mo$neg),
              p = p_6mo_neg) %>%
      add_row(model = paste0(fc_type, '_', scale, '_12mo_pos'),
              r_mean = mean(r_vals_12mo$pos),
              p = p_12mo_pos) %>%
      add_row(model = paste0(fc_type, '_', scale, '_12mo_neg'),
              r_mean = mean(r_vals_12mo$neg),
              p = p_12mo_neg)

    # Retain only the highest mean r value across the 16 CPM models for family-wise error correction
    cpm_nulls <- cbind(cpm_nulls, null_r_means_6mo, null_r_means_12mo)
    
    # ============================================= KRR =============================================
    
    # Load the data
    r_vals_6mo  <- read_excel(paste0(repo, '/krr/', scale, '/4/', fc_type, '/r_vals.xlsx'), col_names = FALSE)
    r_vals_12mo <- read_excel(paste0(repo, '/krr/', scale, '/5/', fc_type, '/r_vals.xlsx'), col_names = FALSE)
  
    r_vals <- cbind(r_vals_6mo, r_vals_12mo) %>%
      setNames(c('6 months', '12 months')) %>%
      gather(key = 'timepoint', value = 'value') %>%
      mutate(timepoint = factor(timepoint, levels = c('6 months', '12 months')))
    
    # Plot KRR r values
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
      ylim(-0.6, 0.6) +
      scale_x_discrete(expand = c(0.1, 0)) +
      theme(axis.text=element_text(size=16),
            axis.title=element_blank(),
            legend.position = 'none')
    
    ggsave(paste0(repo, '/figures/', scale, '/', fc_type, '/krr_rvals.jpg'), plot = fig, width = 4, height = 7)
    
    null_r_means_6mo  <- read_excel(paste0(repo, '/krr/', scale, '/4/', fc_type, '/null_r_means.xlsx'), col_names = FALSE) %>%
      na.omit()
    null_r_means_12mo <- read_excel(paste0(repo, '/krr/', scale, '/5/', fc_type, '/null_r_means.xlsx'), col_names = FALSE) %>%
      na.omit()
    
    # Calculate p values as proportion of null r values greater than observed r values
    p_6mo  <- sum(mean(r_vals_6mo$...1)  < null_r_means_6mo[,1])  / dim(null_r_means_6mo)[1]
    p_12mo <- sum(mean(r_vals_12mo$...1) < null_r_means_12mo[,1]) / dim(null_r_means_12mo)[1]
    
    krr_stats <- krr_stats %>%
      add_row(model = paste0(fc_type, '_', scale, '_6mo'),
              r_mean = mean(r_vals_6mo$...1),
              p = p_6mo) %>%
      add_row(model = paste0(fc_type, '_', scale, '_12mo'),
              r_mean = mean(r_vals_12mo$...1),
              p = p_12mo)
    
    # Keep highest r value across models (6mo/12mo) for family-wise error correction
    krr_nulls <- cbind(krr_nulls, null_r_means_6mo, null_r_means_12mo)
  }
  
  # ======================================== META-MATCHING ========================================
  
  # Load the data
  r_vals_6mo  <- read_excel(paste0(repo, '/meta-matching/', scale, '/4/r_vals.xlsx'), col_names = TRUE)
  r_vals_12mo <- read_excel(paste0(repo, '/meta-matching/', scale, '/5/r_vals.xlsx'), col_names = TRUE)
  
  r_vals <- cbind(r_vals_6mo, r_vals_12mo) %>%
    setNames(c('6 months', '12 months')) %>%
    gather(key = 'timepoint', value = 'value') %>%
    mutate(timepoint = factor(timepoint, levels = c('6 months', '12 months')))
  
  # Plot meta-matching r values
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
    theme(axis.text=element_text(size=16),
          axis.title=element_blank(),
          legend.position = 'none')
  
  ggsave(paste0(repo, '/figures/', scale, '/fc_baseline/meta-matching_rvals.jpg'), plot = fig, width = 4, height = 7)
  
  null_r_means_6mo  <- read_excel(paste0(repo, '/meta-matching/', scale, '/4/null_r_means.xlsx'), col_names = TRUE) %>%
    select(-1) %>%
    na.omit()
  null_r_means_12mo <- read_excel(paste0(repo, '/meta-matching/', scale, '/5/null_r_means.xlsx'), col_names = TRUE) %>%
    select(-1) %>%
    na.omit()

  # Calculate p values as proportion of null r values greater than observed r values
  p_6mo  <- sum(mean(r_vals_6mo$`r value`)  < null_r_means_6mo[,1])  / dim(null_r_means_6mo)[1]
  p_12mo <- sum(mean(r_vals_12mo$`r value`) < null_r_means_12mo[,1]) / dim(null_r_means_12mo)[1]
  
  meta_stats <- meta_stats %>%
    add_row(model = paste0('fc_baseline_', scale, '_6mo'),
            r_mean = mean(r_vals_6mo$`r value`),
            p = p_6mo) %>%
    add_row(model = paste0('fc_baseline_', scale, '_12mo'),
            r_mean = mean(r_vals_12mo$`r value`),
            p = p_12mo)
  
  # Keep highest r value across models (6mo/12mo) for family-wise error correction
  meta_nulls <- cbind(meta_nulls, null_r_means_6mo, null_r_means_12mo)
  
  # ====================================== OBSERVED OUTCOMES ======================================
  for (timepoint in c('4', '5')) {
    
    obs_fcbl <- read_excel(paste0(repo, '/cpm/', scale, '/', timepoint, '/fc_baseline/change_observed.xlsx'), col_names = TRUE) %>%
      select(2) %>%
      as.data.frame()
    
    obs_fcbl3m <- read_excel(paste0(repo, '/cpm/', scale, '/', timepoint, '/fc_change_bl3m/change_observed.xlsx'), col_names = TRUE) %>%
      select(2) %>%
      as.data.frame()
    
    # Make plots of observed change scores
    fig <- ggplot() +
      geom_histogram(data = obs_fcbl, aes(x = change_score),
                     binwidth = 0.05,
                     fill = 'mediumpurple1',
                     color='black',
                     alpha = 0.5) +
      geom_histogram(data = obs_fcbl3m, aes(x = change_score),
                     binwidth = 0.05,
                     fill = 'yellow3',
                     color='black',
                     alpha = 0.5) +
      theme_minimal() +
      xlim(-1, 1) +
      ylim(0, 8) +
      theme(axis.text=element_text(size=16),
            axis.title=element_blank())
    
    ggsave(paste0(repo, '/figures/', scale, '/observed_change_scores_', timepoint, '.jpg'), plot = fig, width = 8, height = 6)
  }
}

# Family-wise error correction: keep the highest r_mean across the 16 CPM models / 8 KRR models / 4 meta-matching models
cpm_nulls  <- apply(cpm_nulls, 1, max)
krr_nulls  <- apply(krr_nulls, 1, max)
meta_nulls  <- apply(meta_nulls, 1, max)

# Initialise looping values for finding r_means in cpm_stats, krr_stats, meta_stats
c <- 0
k <- 0
m <- 0

for (scale in c('sofas', 'bprs')) {
  
  for (fc_type in c('fc_baseline', 'fc_change_bl3m')) {
    
    # ============================================= CPM =============================================
    
    # Recalculate p-values
    cpm_stats[c+1,4] <- sum(cpm_stats[c+1,2] < cpm_nulls) / length(cpm_nulls)
    cpm_stats[c+2,4] <- sum(cpm_stats[c+2,2] < cpm_nulls) / length(cpm_nulls)
    cpm_stats[c+3,4] <- sum(cpm_stats[c+3,2] < cpm_nulls) / length(cpm_nulls)
    cpm_stats[c+4,4] <- sum(cpm_stats[c+4,2] < cpm_nulls) / length(cpm_nulls)
    
    vline_data <- data.frame(xintercept = c(cpm_stats[c+1,2], cpm_stats[c+2,2], cpm_stats[c+3,2], cpm_stats[c+4,2]),
                             model = c('6mo_pos', '6mo_neg', '12mo_pos', '12mo_neg'))
    
    # Plot r_means against the FWE-corrected null distribution
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
    
    ggsave(paste0(repo, '/figures/', scale, '/', fc_type, '/cpm_nulls.jpg'), plot = fig, width = 8, height = 6)
    
    # ============================================= KRR =============================================
    
    # Recalculate p-values
    krr_stats[k+1,4] <- sum(krr_stats[k+1,2] < krr_nulls) / length(krr_nulls)
    krr_stats[k+2,4] <- sum(krr_stats[k+2,2] < krr_nulls) / length(krr_nulls)
    
    vline_data <- data.frame(xintercept = c(krr_stats[k+1,2], krr_stats[k+2,2]),
                             model = c('6mo', '12mo'))
    
    # Plot r_means against the FWE-corrected null distribution
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
    
    ggsave(paste0(repo, '/figures/', scale, '/', fc_type, '/krr_nulls.jpg'), plot = fig, width = 8, height = 6)
    
    c <- c + 4
    k <- k + 2
  }
  
  # ============================================= META-MATCHING =============================================
  
  # Recalculate p-values
  meta_stats[m+1,4] <- sum(meta_stats[m+1,2] < meta_nulls) / length(meta_nulls)
  meta_stats[m+2,4] <- sum(meta_stats[m+2,2] < meta_nulls) / length(meta_nulls)
  
  vline_data <- data.frame(xintercept = c(meta_stats[m+1,2], meta_stats[m+2,2]),
                           model = c('6mo', '12mo'))
  
  # Plot r_means against the FWE-corrected null distribution
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
  
  ggsave(paste0(repo, '/figures/', scale, '/fc_baseline/meta-matching_nulls.jpg'), plot = fig, width = 8, height = 6)
  
  m <- m + 2
  
  for (timepoint in c('4', '5')) {
    # Ad hoc test for group differences between dropouts and 3-month completers
    obs_fcbl3m <- read_excel(paste0(repo, '/cpm/', scale, '/', timepoint, '/fc_change_bl3m/change_observed.xlsx'), col_names = TRUE) %>%
      as.data.frame()
    
    completers <- obs_fcbl3m %>%
      select(1) %>%
      unlist()
    
    obs_fcbl3m <- obs_fcbl3m %>%
      select(2) %>%
      unlist()
    
    obs_dropouts <- read_excel(paste0(repo, '/cpm/', scale, '/', timepoint, '/fc_baseline/change_observed.xlsx'), col_names = TRUE) %>%
      as.data.frame() %>%
      # remove completers
      filter(!BIDS_ID %in% completers) %>%
      select(2) %>%
      unlist()
    
    # Perform t-test
    ttest <- t.test(obs_fcbl3m, obs_dropouts)
    
    print(paste0('t-test comparing ', scale, ' change scores at timepoint ', timepoint, ' between patients with usable 3-month FC data and those without: p = ', ttest[3]))
  }
}

# Save stats
write_xlsx(cpm_stats,  paste0(repo, '/cpm/cpm_stats.xlsx'))
write_xlsx(krr_stats,  paste0(repo, '/krr/krr_stats.xlsx'))
write_xlsx(meta_stats, paste0(repo, '/meta-matching/meta-matching_stats.xlsx'))