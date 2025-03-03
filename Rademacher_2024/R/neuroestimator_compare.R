# scripts to annotate and compare Neuroestimator predicted activity results
library(dplyr)
library(ggplot2)
library(ggpubr)
library(plyr)
library(tidyr)
library(stringr)
library(Seurat)
library(semla)

#------------ Inputs -----------#
setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/neuroestimator/analysis")
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/neuroestimator/analysis/neuroestimator_res_291124.Rdata")


#------------ annnotate neuroestimator results -----------#
setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/neuroestimator/241105_CountsForZac")

count_files <- list.files(pattern=".rds")
res2 <- list()

run_index <- 1:length(count_files)
run_index <- run_index[!run_index %in% c(9,14,21)]

for (i in run_index) {
  print(i)
  counts <- readRDS(count_files[i])

  tmp <- t(as.data.frame(unlist(strsplit(gsub("transcripts", "", count_files[i]), "_")[[1]])))
  
  
  row.names(tmp) <- NULL
  tmp2 <- cbind(res[[i]],tmp[rep(seq_len(nrow(tmp)), each = nrow(res[[i]])), 1:3])
  colnames(tmp2) <- c("predicted_activity","sample","replicate","region")
  
  res2[[i]] <- tmp2
}


#------------ evaluate neuroestimator results - targeted spots -----------#
setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/neuroestimator/analysis")

# format
dat <- do.call(rbind,res2)
dat$sample <- factor(dat$sample, levels = c("CNOOnly","GqVEH","GqCNO"))

# Perform KS tests and generate plots per region
generate_region_plot <- function(data, region_name) {
  # Subset the data for the current region
  region_data <- data %>% filter(region == region_name)
  
  # Perform pairwise KS tests
  ks_results <- combn(levels(region_data$sample), 2, function(samples) {
    group1 <- region_data %>% filter(sample == samples[1]) %>% pull(predicted_activity)
    group2 <- region_data %>% filter(sample == samples[2]) %>% pull(predicted_activity)
    p_value <- ks.test(group1, group2)$p.value
    data.frame(
      sample1 = samples[1],
      sample2 = samples[2],
      p_value = p_value,
      region = region_name
    )
  }, simplify = FALSE)
  
  ks_results <- do.call(rbind, ks_results)
  
  # Calculate positions for annotation lines
  ks_results <- ks_results %>%
    mutate(
      x1 = as.numeric(factor(sample1, levels = levels(region_data$sample))),
      x2 = as.numeric(factor(sample2, levels = levels(region_data$sample))),
      y = max(region_data$predicted_activity, na.rm = TRUE) * 1.05 + 0.1 * (seq_len(nrow(ks_results)) - 1)
    )
  
  # Create boxplot
  plot <- ggplot(region_data, aes(x = sample, y = predicted_activity, fill = sample)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = paste("Region:", region_name),
         x = "Sample",
         y = "NEUROeSTIMator Output (Predicted Activity)") +
    theme_minimal() +
    theme(legend.position = "none") + ylim(0,1.2) + xlab("") +
    scale_x_discrete(guide = guide_axis(angle = 45)) 
  
  # Add p-values and lines
  plot <- plot +
    geom_segment(data = ks_results,
                 aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, color = "black") +
    geom_text(data = ks_results,
              aes(x = (x1 + x2) / 2, y = y + 0.05, label = paste0("p = ", signif(p_value, 3))),
              inherit.aes = FALSE, size = 3)
  
  return(plot)
}

# Generate plots for all regions
regions <- unique(dat$region)
plots <- lapply(regions, function(region) generate_region_plot(dat, region))

# Combine plots using ggarrange
combined_plot <- ggarrange(plotlist = plots, ncol = 3, nrow = 2)
print(combined_plot)

# save
ggsave("boxplots_neuroestimator_KS_mm_031224.png", combined_plot, width = 10, height = 10)


#------------ evaluate neuroestimator results - targeted spots -----------#
setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/neuroestimator/analysis")

# Step 1: Calculate the mean of predicted_activity in CP for each sample group
dat$sample <- factor(dat$sample, levels = c("CNOOnly","GqVEH","GqCNO"))

cp_means <- dat %>%
  filter(region == "CP") %>%
  dplyr::group_by(sample) %>%
  dplyr::summarize(cp_mean_activity = mean(predicted_activity, na.rm = TRUE), .groups = 'drop')

# Normalize SN and VTA by CP mean
dat_normalized <- dat %>%
  left_join(cp_means, by = "sample") %>%
  mutate(
    normalized_activity = ifelse(region %in% c("SN", "VTA"),
                                 predicted_activity / cp_mean_activity,
                                 predicted_activity)
  )
# Initialize an empty list to store results
ks_results <- list()

# Perform pairwise comparisons explicitly
pairs <- list(
  c("CNOOnly", "GqCNO"),
  c("CNOOnly", "GqVEH"),
  c("GqCNO", "GqVEH")
)

for (pair in pairs) {
  group1 <- region_data %>% filter(sample == pair[1]) %>% pull(normalized_activity) %>% na.omit()
  group2 <- region_data %>% filter(sample == pair[2]) %>% pull(normalized_activity) %>% na.omit()
  print(group1)
  print(group2)
  
  if (length(group1) > 0 && length(group2) > 0) {
    p_value <- ks.test(group1, group2)$p.value
  } else {
    p_value <- NA
  }
  print(pair[1])
  print(pair[2])
  
  ks_results[[paste(pair[1], pair[2], sep = "_vs_")]] <- data.frame(
    sample1 = pair[1],
    sample2 = pair[2],
    p_value = p_value,
    region = region_name
  )
}

# Combine all results into a single data frame
ks_results_df <- do.call(rbind, ks_results)

# Print results for verification
print(ks_results_df)

# Step 3: Update the function to use normalized_activity for SN and VTA
generate_region_plot <- function(data, region_name) {
  # Subset the data for the current region
  region_data <- data %>% filter(region == region_name)
  
  # Perform pairwise KS tests
  ks_results <- combn(levels(region_data$sample), 2, function(samples) {
    group1 <- region_data %>% filter(sample == samples[1]) %>% pull(normalized_activity)
    group2 <- region_data %>% filter(sample == samples[2]) %>% pull(normalized_activity)
    p_value <- ks.test(group1, group2)$p.value
    data.frame(
      sample1 = samples[1],
      sample2 = samples[2],
      p_value = p_value,
      region = region_name
    )
  }, simplify = FALSE)
  
  ks_results <- do.call(rbind, ks_results)
  
  # Calculate positions for annotation lines
  ks_results <- ks_results %>%
    mutate(
      x1 = as.numeric(factor(sample1, levels = levels(region_data$sample))),
      x2 = as.numeric(factor(sample2, levels = levels(region_data$sample))),
      y = max(region_data$normalized_activity, na.rm = TRUE) * 1.05 + 0.1 * (seq_len(nrow(ks_results)) - 1)
    )
  
  # Create boxplot
  plot <- ggplot(region_data, aes(x = sample, y = normalized_activity, fill = sample)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = paste("Region:", region_name),
         x = "Sample",
         y = "CP Normalized Activity") +
    theme_minimal() +
    theme(legend.position = "none") + ylim(0,2) + xlab("") +
    scale_x_discrete(guide = guide_axis(angle = 45)) 
  
  # Add p-values and lines
  plot <- plot +
    geom_segment(data = ks_results,
                 aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, color = "black") +
    geom_text(data = ks_results,
              aes(x = (x1 + x2) / 2, y = y + 0.05, label = paste0("p = ", signif(p_value, 3))),
              inherit.aes = FALSE, size = 3)
  
  return(plot)
  dev.off()
}

# Generate plots for all regions
regions <- unique(dat$region)
plots <- lapply(regions, function(region) generate_region_plot(dat_normalized, region))

# Combine plots using ggarrange
combined_plot <- ggarrange(plotlist = list(plots[[2]],plots[[3]]), ncol = 3, nrow = 2)
print(combined_plot)

# Save the plot
ggsave("boxplots_neuroestimator_KS_normalized_101224.png", combined_plot, width = 10, height = 10)



#------------ evaluate neuroestimator results - total Visium DREAD model -----------#
# load previos files that define the spots wthin the SN, VTA and CP
setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/neuroestimator/241105_CountsForZac")
count_files <- list.files(pattern=".rds")

# Function to extract the prefix (everything before the last "_")
extract_prefix <- function(filename) {
  str_extract(filename, "^[^_]+_[^_]+")
}

# Extract prefixes
file_groups <- split(count_files, sapply(count_files, extract_prefix))

# # Read and bind files row-wise for each prefix
# merged_counts <- lapply(file_groups, function(files) {
#   message("Merging: ", paste(files, collapse = ", "))
#   counts_list <- lapply(files, readRDS)  # Read each file
#   do.call(cbind, counts_list)            # Bind row-wise
# })

# Example: Access merged object for a specific prefix
merged_counts[["transcriptsCNOOnly_B1"]]












setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/neuroestimator/241105_CountsForZac")

count_files <- list.files(pattern=".rds")
res2 <- list()

run_index <- 1:length(count_files)
run_index <- run_index[!run_index %in% c(9,14,21)]

for (i in run_index) {
  print(i)
  counts <- readRDS(count_files[i])
  
  spots <- as.data.frame(colnames(counts))
  tmp <- t(as.data.frame(unlist(strsplit(gsub("transcripts", "", count_files[i]), "_")[[1]])))
  spots$sample <- tmp[1]
  spots$replicate <- tmp[2]
  spots$region<- tmp[3]
  
  spots <- spots[!spots[,1] == "FeatureName", ]

  res2[[i]] <- spots
}

meta1 <- do.call(rbind,res2)
colnames(meta1) <- c("spot","sample","replicate","region") 


# load new neuroestimator results
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/neuroestimator/250110_ForZac/neuroestimator_res_130125.Rdata")

# read in seurat files 
setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/neuroestimator/250114_ForZac")

list_of_paths <- list.files(pattern=".rds")
seurat_files <- lapply(list_of_paths, readRDS)

# update meta for each seurat object with neuroestimator predictions and spot region annotation
for (i in 1:length(seurat_files)) {
  print(i)
  
  # get sample information from file names
  tmp <- t(as.data.frame(unlist(strsplit(gsub("transcripts", "", list_of_paths[i]), "_")[[1]])))
  row.names(tmp) <- NULL
  tmp2 <- cbind(res[[i]],tmp[rep(seq_len(nrow(tmp)), each = nrow(res[[i]])), 1:3])
  colnames(tmp2) <- c("predicted_activity","sample","replicate","region")
  
  # Add previous positions of SN/ VTA/ CP annotations
  tmp2$region <- NA
  tmp2$region <- mapvalues(row.names(tmp2),
                           from =  meta1$spot,
                           to  = meta1$region, warn_missing = F)
  tmp2$region[!row.names(tmp2) %in% meta1$spot ] <- NA
  
  
  # Add metadata to Seurat object
  seurat_files[[i]] <- AddMetaData(seurat_files[[i]],tmp2)

}


## Get spot neighbours using semla
semla_files <- list()

for (i in seq_along(seurat_files)) {
  print(i)
  
  # Convert to Staffli
  se <- UpdateSeuratForSemla(seurat_files[[i]])
  spatial_data <- GetStaffli(se)
  se <- LoadImages(se)
  
  # Check available regions
  available_regions <- unique(se@meta.data$region)
  
  # Define target labels
  target_labels <- c("CP", "SN", "VTA")
  
  # Run RegionNeighbors only if label exists in 'region'
  for (label in target_labels) {
    if (label %in% available_regions) {
      se <- RegionNeighbors(se, column_name = "region",
                            column_labels = label, mode = "all_inner_outer")
    } else {
      warning(paste("Skipping", label, "as it is not found in 'region'"))
    }
  }
  
  semla_files[[i]] <- se
}


# loop through to plot enrichment of slice
for (i in 1:length(seurat_files)){
  plot2 <- SpatialFeaturePlot(seurat_files[[i]], features = "predicted_activity", 
                              images = "slice1", pt.size.factor = 2.5, image.alpha = 1) + 
    theme(
      legend.title = element_text(color = "white", size = 12),
      legend.key.size = unit(0.8, 'cm'),
      #ggtitle(meta_to_plot[i]),
      plot.background = element_rect(fill = "black"),
      panel.background = element_rect(fill = "black"),
      panel.grid = element_blank(),  # Remove all grid lines
      # Change legend
      legend.position = "right",
      legend.background = element_rect(fill = "black", color = NA),
      legend.text = element_text(color = "white", size = 12)
    ) +
    geom_point(shape = 21, stroke = NA, fill = NA, size = 0.5)
  
  ggsave(plot2,file=paste0("spatplot_neuroest_activity_",
                           unique(paste0(seurat_files[[i]]$sample,"_",seurat_files[[i]]$replicate))
                           ,"_black.png")) 
}


# # merge all semla objects
merged_seurat <- Reduce(function(x, y) merge(x, y), semla_files)
merged_seurat$sample_rep <- paste0(merged_seurat$sample,"_",merged_seurat$replicate)


## comparisons
dat <- merged_seurat@meta.data

# format
dat$nb_to_SN[!is.na(dat$nb_to_SN)] <- "SN" 
dat$nb_to_VTA[!is.na(dat$nb_to_VTA)] <- "VTA" 
dat$nb_to_CP[!is.na(dat$nb_to_CP)] <- "CP" 
dat$sample <- factor(dat$sample, levels = c("CNOOnly","GqVEH","GqCNO"))
dat$region_single <- dat$region

# Perform KS tests and generate plots per region
generate_region_plot <- function(data, region_name) {
  # Subset the data for the current region
  region_data <- data %>% filter(region == region_name)
  
  # Perform pairwise KS tests
  ks_results <- combn(levels(region_data$sample), 2, function(samples) {
    group1 <- region_data %>% filter(sample == samples[1]) %>% pull(predicted_activity)
    group2 <- region_data %>% filter(sample == samples[2]) %>% pull(predicted_activity)
    p_value <- ks.test(group1, group2)$p.value
    data.frame(
      sample1 = samples[1],
      sample2 = samples[2],
      p_value = p_value,
      region = region_name
    )
  }, simplify = FALSE)
  
  ks_results <- do.call(rbind, ks_results)
  
  # Calculate positions for annotation lines
  ks_results <- ks_results %>%
    mutate(
      x1 = as.numeric(factor(sample1, levels = levels(region_data$sample))),
      x2 = as.numeric(factor(sample2, levels = levels(region_data$sample))),
      y = max(region_data$predicted_activity, na.rm = TRUE) * 1.05 + 0.1 * (seq_len(nrow(ks_results)) - 1)
    )
  
  # Create boxplot
  plot <- ggplot(region_data, aes(x = sample, y = predicted_activity, fill = sample)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = paste("Region:", region_name),
         x = "Sample",
         y = "Neuronal Activity") +
    theme_minimal() +
    theme(legend.position = "none") + ylim(0,1.2) + xlab("") +
    scale_x_discrete(guide = guide_axis(angle = 45)) 
  
  # Add p-values and lines
  plot <- plot +
    geom_segment(data = ks_results,
                 aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, color = "black") +
    geom_text(data = ks_results,
              aes(x = (x1 + x2) / 2, y = y + 0.05, label = paste0("p = ", signif(p_value, 3))),
              inherit.aes = FALSE, size = 3)
  
  return(plot)
}

# gererate plots
dat$region <- dat$nb_to_SN
plot_nb_SN <- generate_region_plot(dat,"SN")
dat$region <- dat$nb_to_VTA
plot_nb_VTA <- generate_region_plot(dat, "VTA")
dat$region <- dat$nb_to_CP
plot_nb_CP <- generate_region_plot(dat, "CP")

# Combine plots using ggarrange
combined_plot <- ggarrange(plotlist = list(plot_nb_SN,plot_nb_VTA,plot_nb_CP), ncol = 3, nrow = 2)
print(combined_plot)

# save
ggsave("boxplots_neuroestimator_KS_neighbour1.png", combined_plot, width = 10, height = 10)


## evaluate activity of total slices
# Create the boxplot with statistical comparisons

# Define all pairwise comparisons
comparisons <- combn(unique(dat$sample_rep), 2, simplify = FALSE)

# Create the boxplot with t-test comparisons
plot1 <- ggplot(dat, aes(x = sample_rep, y = predicted_activity, fill = sample_rep)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Total",
       x = "Sample",
       y = "NEUROeSTIMator Output (Predicted Activity)") +
  theme_minimal() +
  theme(legend.position = "none") +
  ylim(0, 1.2) +
  xlab("") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif")  # Adds p-value brackets

ggsave("boxplots_neuroestimator_per_slice.png", plot1, width = 10, height = 10)


## Attempt to normalize to the activity of entire slice
# Step 1: Calculate the mean of predicted_activity inof each slice
dat$sample <- factor(dat$sample, levels = c("CNOOnly","GqVEH","GqCNO"))
dat$region <- dat$region_single

slice_means <- dat %>%
  dplyr::group_by(sample_rep) %>%
  dplyr::summarize(slice_mean_activity = mean(predicted_activity, na.rm = TRUE), .groups = 'drop')

# Normalize SN and VTA by CP mean
dat_normalized <- dat %>%
  left_join(slice_means , by = "sample_rep") %>%
  mutate(
    normalized_activity = ifelse(region %in% c("SN", "VTA","CP"),
                                 predicted_activity / slice_mean_activity,
                                 predicted_activity)
  )


# Perform KS tests and generate plots per region
generate_region_plot_norm  <- function(data, region_name) {
  # Define custom colors
  sample_colors <- c("CNOOnly" = "dodgerblue", "GqVEH" = "lightblue", "GqCNO" = "red")
  
  # Subset the data for the current region
  region_data <- data %>% filter(region == region_name)
  
  # Perform pairwise KS tests
  ks_results <- combn(levels(region_data$sample), 2, function(samples) {
    group1 <- region_data %>% filter(sample == samples[1]) %>% pull(normalized_activity)
    group2 <- region_data %>% filter(sample == samples[2]) %>% pull(normalized_activity)
    p_value <- ks.test(group1, group2)$p.value
    
    # Format p-values
    if (p_value < 0.001) {
      p_value_formatted <- bquote(.(formatC(p_value, format = "e", digits = 1)))
    } else {
      p_value_formatted <- format(round(p_value, 3), nsmall = 3)
    }
    
    data.frame(
      sample1 = samples[1],
      sample2 = samples[2],
      p_value = p_value,
      p_value_formatted = p_value_formatted,
      region = region_name
    )
  }, simplify = FALSE)
  
  ks_results <- do.call(rbind, ks_results)
  
  # Ensure proper formatting of scientific notation for ggplot2 parsing
  ks_results$p_value_formatted <- ifelse(
    ks_results$p_value < 0.001, 
    paste0("italic(p) == ", gsub("e", " %*% 10^", formatC(ks_results$p_value, format = "e", digits = 1))),
    paste0("italic(p) == ", ks_results$p_value_formatted)
  )
  
  # Calculate positions for annotation lines
  ks_results <- ks_results %>%
    mutate(
      x1 = as.numeric(factor(sample1, levels = levels(region_data$sample))),
      x2 = as.numeric(factor(sample2, levels = levels(region_data$sample))),
      y = max(region_data$normalized_activity, na.rm = TRUE) * 1.05 + 0.1 * (seq_len(nrow(ks_results)) - 1)
    )
  
  # Create boxplot with custom fill colors
  plot <- ggplot(region_data, aes(x = sample, y = normalized_activity, fill = sample)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = region_name,
         x = "Sample",
         y = "Predicted Activity (Normalized)") +
    theme_minimal() +
    theme(legend.position = "none") + ylim(0,1.75) + xlab("") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_manual(values = sample_colors)  # Apply custom colors
  
  # Add p-values and lines
  plot <- plot +
    geom_segment(data = ks_results,
                 aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, color = "black") +
    geom_text(data = ks_results,
              aes(x = (x1 + x2) / 2, y = y + 0.05, label = p_value_formatted),
              inherit.aes = FALSE, size = 3, parse = TRUE)  # Correctly parse scientific notation
  
  return(plot)
}

# gererate plots
dat_normalized$region <- dat_normalized$nb_to_SN
plot_norm_SN <- generate_region_plot_norm(dat_normalized,"SN")
dat_normalized$region <- dat_normalized$nb_to_VTA
plot_norm_VTA <- generate_region_plot_norm(dat_normalized,"VTA")
# dat_normalized$region <- dat_normalized$nb_to_CP
# plot_norm_CP <- generate_region_plot_norm(dat_normalized,"CP")

# Combine plots using ggarrange
combined_plot <- ggarrange(plotlist = list(plot_norm_SN,plot_norm_VTA), ncol = 3, nrow = 2)

# Save the plot
ggsave("boxplots_neuroestimator_KS_normalized_neighbours.pdf", combined_plot, width = 10, height = 10)



















