########################################
### Custom functions for DA analysis ###
###                 -Christian Ayala ###
########################################

# -------------------------------------------------------------------------
global.norm <- function(matrix, transform_data = TRUE){
  # This function will perform normalization based in the global AUC of each sample
  # and the median of such intensities across samples
  
  colsum <- colSums(matrix, na.rm = TRUE)
  colsum.median <- median(colsum)
  norm.matrix <- data.frame(matrix(NA, nrow = nrow(matrix), ncol = ncol(matrix)))
  for(col in 1:ncol(matrix)){
    norm.matrix[,col] <- (matrix[,col] / colsum[col]) * colsum.median
  }
  colnames(norm.matrix) <- colnames(matrix)
  rownames(norm.matrix) <- rownames(matrix)
  if(transform_data == TRUE){
    norm.matrix <- log2(norm.matrix + 1)
  }
  
  return(norm.matrix)
}


# -------------------------------------------------------------------------
median.norm <- function(matrix, transform_data = TRUE){
  # This function will perform data normalization based in the  median AUC of each sample
  
  colmedian <- apply(matrix, 2, FUN = median, na.rm = TRUE)
  colmedian.mean <- mean(colmedian)
  norm.matrix <- data.frame(matrix(NA, nrow = nrow(matrix), ncol = ncol(matrix)))
  for(col in 1:ncol(matrix)){
    norm.matrix[,col] <- (matrix[,col] / colmedian[col]) * colmedian.mean
  }
  colnames(norm.matrix) <- colnames(matrix)
  rownames(norm.matrix) <- rownames(matrix)
  if(transform_data == TRUE){
    norm.matrix <- log2(norm.matrix + 1)
  }
  return(norm.matrix)
} 

# -------------------------------------------------------------------------
mean.norm <- function(matrix, transform_data = TRUE){
  # This function will perform data normalization based in the  mean AUC of each sample
  
  colmean <- colMeans(matrix, na.rm = TRUE)
  colmean.mean <- mean(colmean)
  norm.matrix <- data.frame(matrix(NA, nrow = nrow(matrix), ncol = ncol(matrix)))
  for(col in 1:ncol(matrix)){
    norm.matrix[,col] <- (matrix[,col] / colmean[col]) * colmean.mean
  }
  colnames(norm.matrix) <- colnames(matrix)
  rownames(norm.matrix) <- rownames(matrix)
  if(transform_data == TRUE){
    norm.matrix <- log2(norm.matrix + 1)
  }
  return(norm.matrix)
} 

# -------------------------------------------------------------------------
vsn.norm <- function(matrix){
  # This functions tries to adjust the data to the vsn normalization
  
  norm.matrix <- suppressMessages(vsn::justvsn(as.matrix(matrix)))
  norm.matrix <- as.data.frame(norm.matrix)
  return(norm.matrix)
}
# -------------------------------------------------------------------------
cycloess.norm <- function(matrix){
  # This functions tries to adjust the data to the vsn normalization
  
  norm.matrix <- log2(matrix)
  norm.matrix <- limma::normalizeCyclicLoess(norm.matrix, method = 'fast')
  norm.matrix <- as.data.frame(norm.matrix)
  rownames(norm.matrix) <- rownames(matrix)
  return(norm.matrix)
}

# -------------------------------------------------------------------------
plot_boxplot <- function(df, my_x, my_y){
  ggplot(df,
         aes(x = {{my_x}},
             y = {{my_y}})) +
           geom_boxplot() +
           theme_bw() +
           theme(plot.title = element_text(face = 'bold', hjust = 0.5))
}

# -------------------------------------------------------------------------
normalize_by_all <- function(df){
  
  no_norm.df <- gather(df, key = 'SampleID', value = 'AUC')
  no_norm.plot <- plot_boxplot(no_norm.df, SampleID, AUC) +
    labs(title = 'No normalization') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  
  gi_norm.df <- global.norm(df)
  gi_norm.df <- gather(gi_norm.df, key = 'SampleID', value = 'AUC')
  gi_norm.plot <- plot_boxplot(gi_norm.df, SampleID, AUC) +
    labs(title = 'Normalization by global AUC',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  
  
  mean_norm.df <- mean.norm(df)
  mean_norm.df <- gather(mean_norm.df, key = 'SampleID', value = 'AUC')
  mean_norm.plot <- plot_boxplot(mean_norm.df, SampleID, AUC) +
    labs(title = 'Normalization by mean',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  
  
  median_norm.df <- median.norm(df)
  median_norm.df <- gather(median_norm.df, key = 'SampleID', value = 'AUC')
  median_norm.plot <- plot_boxplot(median_norm.df, SampleID, AUC) +
    labs(title = 'Normalization by median',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  
  
  vsn_norm.df <- vsn.norm(df)
  vsn_norm.df <- gather(vsn_norm.df, key = 'SampleID', value = 'AUC')
  vsn_norm.plot <- plot_boxplot(vsn_norm.df, SampleID, AUC) +
    labs(title = 'VSN',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  
  
  cycloess_norm.df <- cycloess.norm(df)
  cycloess_norm.df <- gather(cycloess_norm.df, key = 'SampleID', value = 'AUC')
  cycloess_norm.plot <- plot_boxplot(cycloess_norm.df, SampleID, AUC) +
    labs(title = 'LOESS normalization',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  
  
  all_norm.plot <- ggarrange(no_norm.plot, gi_norm.plot, mean_norm.plot,
                                median_norm.plot, vsn_norm.plot, cycloess_norm.plot,
                                nrow = 2,
                                ncol = 3)
  
  return(all_norm.plot)
  
}

# -------------------------------------------------------------------------
get_diff_table <- function(auc_matrix, control.sample_list, treatment.sample_list, log2_transformed = FALSE){
  
  # Get the AUC values per each sample and calculate the means per feature
  temp.df_control <- auc_matrix %>% 
    select(all_of(control.sample_list))
  control_means <- rowMeans(temp.df_control, na.rm = TRUE)
  
  temp.df_treatment <- auc_matrix %>% 
    select(all_of(treatment.sample_list))
  treatment_means <- rowMeans(temp.df_treatment, na.rm = TRUE)
  
  diff_table <- as.data.frame(cbind(control_means, treatment_means))
  
  if(log2_transformed == TRUE){
    diff_table <- diff_table %>% 
      mutate(log2FC = treatment_means - control_means)
  } else {
    diff_table <- diff_table %>% 
      mutate(ratio = treatment_means/control_means) %>% # get the control/treatment ratio
      mutate(log2FC = log2(ratio)) # calculate log2FC
  }
  
  rownames(diff_table) <- rownames(auc_matrix)
  
  # Initialize pvalues matrix
  pvalues <- data.frame(row.names = rownames(auc_matrix), pval = rep(0, length(rownames(auc_matrix))))
  
  #Calculate pvalue per each of the features
  for(i in 1:nrow(pvalues)){
    t.test <- t.test(as.numeric(temp.df_control[i,]), as.numeric(temp.df_treatment[i,]), paired = FALSE)
    pvalues$pval[i] <- t.test$p.value
  }
  pvalues <- pvalues %>% 
    rownames_to_column(var = 'FeatureID')
  
  diff_table <- diff_table %>% 
    rownames_to_column(var = 'FeatureID') %>% 
    left_join(pvalues, by = 'FeatureID')
  
  diff_table$pval.adj <- p.adjust(diff_table$pval, method = 'fdr')
  
  return(diff_table)
  
}

# -------------------------------------------------------------------------
plot_volcano <- function(df, column_with_log2FC, column_with_pval, log2FC.threshold, pval.threshold){
  
  #Generate label for the plot
  
  significant_points <- df %>% 
    select(FeatureID, {{column_with_log2FC}}, {{column_with_pval}}) %>% 
    filter(abs({{column_with_log2FC}}) > log2FC.threshold,
           -log10({{column_with_pval}}) > -log10(0.05)) %>% 
    pull(FeatureID)
  
  plot <- df %>%
    mutate(color4plot = case_when(FeatureID %in% significant_points & {{column_with_log2FC}} > 0 ~ 'upregulated',
                                  FeatureID %in% significant_points & {{column_with_log2FC}} < 0 ~ 'downregulated',
                                  TRUE ~ 'non-significant')) %>% 
    ggplot(aes(x = {{column_with_log2FC}},
               y = -log10({{column_with_pval}}))) +
    geom_point(aes(color = color4plot)) +
    scale_color_manual(values = c("upregulated" = 'red', 'non-significant' = 'gray', 'downregulated' = 'blue')) +
    geom_vline(xintercept = c(-{{log2FC.threshold}}, {{log2FC.threshold}}),
               linetype = 'dotted',
               size = 1,
               color = 'blue') +
    geom_hline(yintercept = -log10({{pval.threshold}}),
               linetype = 'dotted',
               size = 1,
               color = 'blue') +
    theme_bw() +
    labs(title = 'Volcano plot',
         x = expression("Log"[2]*" Fold Change"),
         y = expression("-Log"[10]*" pvalue")) +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = 'bold'),
          plot.subtitle = element_text(hjust = 0.5,
                                       face = 'bold'),
          legend.title = element_blank())
  
  return(plot)
}