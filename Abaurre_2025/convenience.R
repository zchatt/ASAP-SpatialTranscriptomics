# convenience scripts for ST validation of Abaurre references

library(patchwork)
library(spacexr)
library(RColorBrewer)

## color palletes

Cell_col = c("Astro_CYP4F12" = brewer.pal(n = 9, name = "Reds")[7],
  "Astro_GBP2_SPOCD1" = brewer.pal(n = 9, name = "Reds")[7],
  "Astro_GJB6_OXTR" = brewer.pal(n = 9, name = "Reds")[7],
  "Astro_GLYATL2" = brewer.pal(n = 9, name = "Reds")[7],
  "Astro_GUCY1A2" = brewer.pal(n = 9, name = "Reds")[7],
  "Astro_SERPINA3" = brewer.pal(n = 9, name = "Reds")[7],
  "Astro_SIDT1" = brewer.pal(n = 9, name = "Reds")[7],
  "Astro_VIM_TNFSRF12A" = brewer.pal(n = 9, name = "Reds")[7],
  "CALB1_CBLN4_PAX5" = brewer.pal(n = 9, name = "Greens")[7],
  "CALB1_CRYM" = brewer.pal(n = 9, name = "Greens")[7],
  "CALB1_CRYM_CALCR" = brewer.pal(n = 9, name = "Greens")[7],
  "CALB1_NEUROD6_PPP1R17" = brewer.pal(n = 9, name = "Greens")[7],
  "CALB1_PRLHR_TRHR" = brewer.pal(n = 9, name = "Greens")[7],
  "CALB1_SEMA3D_RSPO3" = brewer.pal(n = 9, name = "Greens")[7],
  "CALB1_VIP_NPPC" = brewer.pal(n = 9, name = "Greens")[7],
  "Endo_COL6A3" = brewer.pal(n = 9, name = "Greys")[7],
  "Endo_DCN_ABCC9" = brewer.pal(n = 9, name = "Greys")[7],
  "Endo_IL27RA" = brewer.pal(n = 9, name = "Greys")[7],
  "Endo_NOTCH3_PLK2" = brewer.pal(n = 9, name = "Greys")[7],
  "Endo_SLIT3" = brewer.pal(n = 9, name = "Greys")[7],
  "Endo_SNTG2" = brewer.pal(n = 9, name = "Greys")[7],
  "Ex_CYP2J2" = brewer.pal(n = 9, name = "Blues")[7],
  "Ex_EBF2_CTC-552D5.1" = brewer.pal(n = 9, name = "Blues")[7],
  "Ex_LAMP5_BAIAP3" = brewer.pal(n = 9, name = "Blues")[7],
  "Ex_LAMP5_NTNG2" = brewer.pal(n = 9, name = "Blues")[7],
  "Ex_MYO5B" = brewer.pal(n = 9, name = "Blues")[7],
  "Ex_OPRD1" = brewer.pal(n = 9, name = "Blues")[7],
  "Ex_POSTN" = brewer.pal(n = 9, name = "Blues")[7],
  "Ex_PPP1R1C" = brewer.pal(n = 9, name = "Blues")[7],
  "Ex_SATB2" = brewer.pal(n = 9, name = "Blues")[7],
  "Ex_VWA5B1_CALB1" = brewer.pal(n = 9, name = "Blues")[7],
  "GAD2_CALCRL_KCNK13" = brewer.pal(n = 9, name = "Oranges")[7],
  "GAD2_EBF2_NPSR1" = brewer.pal(n = 9, name = "Oranges")[7],
  "Inh_IGFBP5" = brewer.pal(n = 9, name = "Blues")[7],
  "Inh_INHBA" = brewer.pal(n = 9, name = "Blues")[7],
  "Inh_OTX2_CASR" = brewer.pal(n = 9, name = "Blues")[7],
  "Inh_PAX5_CCBE1" = brewer.pal(n = 9, name = "Blues")[7],
  "Inh_PAX5_VCAN" = brewer.pal(n = 9, name = "Blues")[7],
  "Inh_PRLR_RP11-384J4.2" = brewer.pal(n = 9, name = "Blues")[7],
  "Inh_SIX3" = brewer.pal(n = 9, name = "Purples")[7],
  "Macro_CD200R1" = brewer.pal(n = 9, name = "Reds")[7],
  "MG_CCL3" = brewer.pal(n = 9, name = "Reds")[7],
  "MG_CECR2_FGL1" = brewer.pal(n = 9, name = "Reds")[7],
  "MG_FOSL2" = brewer.pal(n = 9, name = "Reds")[7],
  "MG_GPNMB_LPL" = brewer.pal(n = 9, name = "Reds")[7],
  "MG_GPNMB_SULT1C2" = brewer.pal(n = 9, name = "Reds")[7],
  "MG_GPNMB_SUSD1" = brewer.pal(n = 9, name = "Reds")[7],
  "MG_MGAM" = brewer.pal(n = 9, name = "Reds")[7],
  "MG_MKI67" = brewer.pal(n = 9, name = "Reds")[7],
  "MG_OPRM1" = brewer.pal(n = 9, name = "Reds")[7],
  "MG_SPON1" = brewer.pal(n = 9, name = "Reds")[7],
  "MG_TSPO_VIM" = brewer.pal(n = 9, name = "Reds")[7],
  "Olig_ENPP6_ACTN2" = brewer.pal(n = 9, name = "RdPu")[7],
  "Olig_ENPP6_EMILIN2" = brewer.pal(n = 9, name = "RdPu")[7],
  "Olig_ENPP6_LUCAT1" = brewer.pal(n = 9, name = "RdPu")[7],
  "Olig_PLXDC2" = brewer.pal(n = 9, name = "RdPu")[7],
  "Olig_PLXDC2_KCNAB1" = brewer.pal(n = 9, name = "RdPu")[7],
  "Olig_PLXDC2_KCNK10" = brewer.pal(n = 9, name = "RdPu")[7],
  "Olig_PLXDC2_SFRP1" = brewer.pal(n = 9, name = "RdPu")[7],
  "OPC_CACNG4" = brewer.pal(n = 9, name = "RdPu")[7],
  "OPC_HOXD3" = brewer.pal(n = 9, name = "RdPu")[7],
  "OPC_KIAA0040" = brewer.pal(n = 9, name = "RdPu")[7],
  "SOX6_AGTR1_NOX4" = brewer.pal(n = 9, name = "Purples")[7],
  "SOX6_GFRA2_TBC1D8B" = brewer.pal(n = 9, name = "Purples")[7],
  "SOX6_SMOC1_LPL" = brewer.pal(n = 9, name = "Purples")[7]
)

DA_col = c("CALB1_CBLN4_PAX5" = brewer.pal(n = 9, name = "Greens")[3],
             "CALB1_CRYM" = brewer.pal(n = 9, name = "Greens")[4],
             "CALB1_CRYM_CALCR" = brewer.pal(n = 9, name = "Greens")[5],
             "CALB1_NEUROD6_PPP1R17" = brewer.pal(n = 9, name = "Greens")[6],
             "CALB1_PRLHR_TRHR" = brewer.pal(n = 9, name = "Greens")[7],
             "CALB1_SEMA3D_RSPO3" = brewer.pal(n = 9, name = "Greens")[8],
             "CALB1_VIP_NPPC" = brewer.pal(n = 9, name = "Greens")[9],
             "GAD2_CALCRL_KCNK13" = brewer.pal(n = 9, name = "Oranges")[5],
             "GAD2_EBF2_NPSR1" = brewer.pal(n = 9, name = "Oranges")[7],
             "SOX6_AGTR1_NOX4" = brewer.pal(n = 9, name = "Purples")[8],
             "SOX6_GFRA2_TBC1D8B" = brewer.pal(n = 9, name = "Purples")[7],
             "SOX6_SMOC1_LPL" = brewer.pal(n = 9, name = "Purples")[6]
)

DA_col2 = c("CALB1_CBLN4" = brewer.pal(n = 9, name = "Greens")[3],
             "CALB1_CRYM" = brewer.pal(n = 9, name = "Greens")[4],
             "CALB1_CRYM_CALCR" = brewer.pal(n = 9, name = "Greens")[5],
             "CALB1_NEUROD6_PPP1R17" = brewer.pal(n = 9, name = "Greens")[6],
             "CALB1_TRHR" = brewer.pal(n = 9, name = "Greens")[7],
             "CALB1_SEMA3D" = brewer.pal(n = 9, name = "Greens")[8],
             "CALB1_VIP" = brewer.pal(n = 9, name = "Greens")[9],
             "GAD2_CALCRL" = brewer.pal(n = 9, name = "Oranges")[5],
             "GAD2_EBF2" = brewer.pal(n = 9, name = "Oranges")[7],
             "SOX6_SMOC1" = brewer.pal(n = 9, name = "Purples")[4],
             "SOX6_AGTR1" = brewer.pal(n = 9, name = "Purples")[8],
             "SOX6_PART1" = brewer.pal(n = 9, name = "Purples")[7],
             "SOX6_LPL" = brewer.pal(n = 9, name = "Purples")[6]
)

##  Version 1
DA_col_cytograph = c('#efa2fe', '#0074db', '#983e00', '#4b005b', '#005b30', '#2acd47', '#fecb98', '#7f7f7f', '#93feb4', '#8e7b00', '#9ccb00', '#c10087')
names(DA_col_cytograph) <- c("CALB1_CBLN4_PAX5","CALB1_CRYM","CALB1_CRYM_CALCR","CALB1_NEUROD6_PPP1R17","CALB1_PRLHR_TRHR","CALB1_SEMA3D_RSPO3","CALB1_VIP_NPPC","GAD2_CALCRL_KCNK13","GAD2_EBF2_NPSR1","SOX6_AGTR1_NOX4","SOX6_GFRA2_TBC1D8B","SOX6_SMOC1_LPL")


##  Version 2
# DA_col_cytograph = c('#efa2fe', '#0074db', '#983e00', '#4b005b', '#005b30', '#2acd47', '#fecb98', '#7f7f7f', '#93feb4', '#8e7b00', '#9ccb00', '#c10087')
# names(DA_col_cytograph) <- c("CALB1_CBLN4","CALB1_CRYM","CALB1_CRYM_CALCR","CALB1_NEUROD6_PPP1R17",
#                              "CALB1_PRLHR_TRHR","CALB1_SEMA3D_RSPO3","CALB1_VIP_NPPC","GAD2_CALCRL_KCNK13","GAD2_EBF2_NPSR1","SOX6_AGTR1_NOX4","SOX6_GFRA2_TBC1D8B","SOX6_SMOC1_LPL")
# 

## Version 3
# ## get new cytograph - Version 3
# # load sc
# sc_seurat_obj = "/Users/zacc/USyd/spatial_transcriptomics/reports/Abaurre_DA_manuscript/20250122-abaurre_ref_seurat.rds"
# sc_obj <- LoadSeuratRds(file=sc_seurat_obj)
# 
# # read in colors
# colors_df <- read.csv("/Users/zacc/USyd/spatial_transcriptomics/reports/Abaurre_DA_manuscript/20250122-abaurre_cluster_colors.csv")
# 
# # da neurons
# da_cells <- names(table(sc_obj@meta.data$cell_type_merge[is.na(sc_obj@meta.data$dataset_merge)]))
# 
# # Extracting clusters
# sc_clusters <- as.numeric(unique(sc_obj$Clusters[sc_obj$cell_type_merge %in% da_cells])) + 1
# names(sc_clusters) <- unique(sc_obj$cell_type_merge[sc_obj$cell_type_merge %in% da_cells])
# 
# # Match clusters to colors using the provided color map
# DA_col_cytograph <- colors_df$Color[match(sc_clusters, colors_df$Cluster)]
# 
# # Assign names for better visualization
# names(DA_col_cytograph) <- names(sc_clusters)

# DA_col_cytograph = c( 
#   "SOX6_LPL"                 = "#efa2fe", 
#   "SOX6_SMOC1"               = "#0075db", 
#   "SOX6_AGTR1"               = "#983f00", 
#   "SOX6_PART1_GFRA2"         = "#4c005c", 
#   "CALB1_CRYM"               = "#005c31", 
#   "CALB1_CRYM_CALCR"         = "#2bcd48", 
#   "GAD2_EBF2"                = "#fecb98", 
#   "GAD2_CALCRL"              = "#808080", 
#   "CALB1_NEUROD6_PPP1R17"     = "#93feb4", 
#   "CALB1_PAX5"               = "#8e7c00", 
#   "CALB1_TRHR"               = "#9ccb00", 
#   "CALB1_SEMA3D"             = "#c10087", 
#   "CALB1_VIP"                = "#003380", 
#   "CALB1_NPW"                = "#fea305"
# )



## functions
create.RCTD_logFC<- function (spatialRNA, reference, max_cores = 4, test_mode = FALSE, 
          gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 2e-04, 
          fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 2e+07, counts_MIN = 10, 
          UMI_min_sigma = 300, class_df = NULL, CELL_MIN_INSTANCE = 25, 
          cell_type_names = NULL, MAX_MULTI_TYPES = 4, keep_reference = F, 
          cell_type_profiles = NULL, CONFIDENCE_THRESHOLD = 5, DOUBLET_THRESHOLD = 20) 
{
  config <- list(gene_cutoff = gene_cutoff, fc_cutoff = fc_cutoff, 
                 gene_cutoff_reg = gene_cutoff_reg, fc_cutoff_reg = fc_cutoff_reg, 
                 UMI_min = UMI_min, UMI_min_sigma = UMI_min_sigma, max_cores = max_cores, 
                 N_epoch = 8, N_X = 50000, K_val = 100, N_fit = 1000, 
                 N_epoch_bulk = 30, MIN_CHANGE_BULK = 1e-04, MIN_CHANGE_REG = 0.001, 
                 UMI_max = UMI_max, counts_MIN = counts_MIN, MIN_OBS = 3, 
                 MAX_MULTI_TYPES = MAX_MULTI_TYPES, CONFIDENCE_THRESHOLD = CONFIDENCE_THRESHOLD, 
                 DOUBLET_THRESHOLD = DOUBLET_THRESHOLD)
  if (test_mode) 
    config <- list(gene_cutoff = 0.00125, fc_cutoff = 0.5, 
                   gene_cutoff_reg = 0.002, fc_cutoff_reg = 0.75, UMI_min = 1000, 
                   N_epoch = 1, N_X = 50000, K_val = 100, N_fit = 50, 
                   N_epoch_bulk = 4, MIN_CHANGE_BULK = 1, MIN_CHANGE_REG = 0.001, 
                   UMI_max = 2e+05, MIN_OBS = 3, max_cores = 1, counts_MIN = 5, 
                   UMI_min_sigma = 300, MAX_MULTI_TYPES = MAX_MULTI_TYPES, 
                   CONFIDENCE_THRESHOLD = CONFIDENCE_THRESHOLD, DOUBLET_THRESHOLD = DOUBLET_THRESHOLD)
  if (is.null(cell_type_profiles)) {
    if (is.null(cell_type_names)) 
      cell_type_names <- levels(reference@cell_types)
    cell_type_info <- list(info = spacexr:::process_cell_type_info(reference, 
                                                         cell_type_names = cell_type_names, CELL_MIN = CELL_MIN_INSTANCE), 
                           renorm = NULL)
  }
  else {
    cell_type_names <- colnames(cell_type_profiles)
    cell_type_info <- list(info = list(cell_type_profiles, 
                                       cell_type_names, length(cell_type_names)), renorm = NULL)
  }
  if (!keep_reference) 
    reference <- spacexr:::create_downsampled_data(reference, n_samples = 5)
  puck.original = restrict_counts(spatialRNA, rownames(spatialRNA@counts), 
                                  UMI_thresh = config$UMI_min, UMI_max = config$UMI_max, 
                                  counts_thresh = config$counts_MIN)
  message("create.RCTD: getting regression differentially expressed genes: ")
  gene_list_reg = get_de_genes_logFC(cell_type_info$info, puck.original, 
                               fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, 
                               MIN_OBS = config$MIN_OBS)
  return(gene_list_reg)
  # if (length(gene_list_reg) < 10) 
  #   stop("create.RCTD: Error: fewer than 10 regression differentially expressed genes found")
  # message("create.RCTD: getting platform effect normalization differentially expressed genes: ")
  # gene_list_bulk = get_de_genes_logFC(cell_type_info$info, puck.original, 
  #                               fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff, 
  #                               MIN_OBS = config$MIN_OBS)
  # if (length(gene_list_bulk) < 10) 
  #   stop("create.RCTD: Error: fewer than 10 bulk differentially expressed genes found")
  # puck = restrict_counts(puck.original, gene_list_bulk, UMI_thresh = config$UMI_min, 
  #                        UMI_max = config$UMI_max, counts_thresh = config$counts_MIN)
  # puck = restrict_puck(puck, colnames(puck@counts))
  # if (is.null(class_df)) 
  #   class_df <- data.frame(cell_type_info$info[[2]], row.names = cell_type_info$info[[2]])
  # colnames(class_df)[1] = "class"
  # internal_vars <- list(gene_list_reg = gene_list_reg, gene_list_bulk = gene_list_bulk, 
  #                       proportions = NULL, class_df = class_df, cell_types_assigned = F)
  # new("RCTD", spatialRNA = puck, originalSpatialRNA = puck.original, 
  #     reference = reference, config = config, cell_type_info = cell_type_info, 
  #     internal_vars = internal_vars)
}


get_de_genes_logFC <- function (cell_type_info, puck, fc_thresh = 1.25, expr_thresh = 0.00015, 
          MIN_OBS = 3) 
{
  total_gene_list = c()
  epsilon = 1e-09
  bulk_vec = rowSums(puck@counts)
  gene_list = rownames(cell_type_info[[1]])
  prev_num_genes <- min(length(gene_list), length(names(bulk_vec)))
  if (length(grep("mt-", gene_list)) > 0) 
    gene_list = gene_list[-grep("mt-", gene_list)]
  gene_list = intersect(gene_list, names(bulk_vec))
  if (length(gene_list) == 0) 
    stop("get_de_genes: Error: 0 common genes between SpatialRNA and Reference objects. Please check for gene list nonempty intersection.")
  gene_list = gene_list[bulk_vec[gene_list] >= MIN_OBS]
  if (length(gene_list) < 0.1 * prev_num_genes) 
    stop("get_de_genes: At least 90% of genes do not match between the SpatialRNA and Reference objects. Please examine this. If this is intended, please remove the missing genes from the Reference object.")
  
  deg_info_list <- list() # Z.C modified to extract the logFC and exp for DEG used in RCTD analysis
  for (cell_type in cell_type_info[[2]]) {
    if (cell_type_info[[3]] > 2) 
      other_mean = rowMeans(cell_type_info[[1]][gene_list, 
                                                cell_type_info[[2]] != cell_type])
    else {
      other_mean <- cell_type_info[[1]][gene_list, cell_type_info[[2]] != 
                                          cell_type]
      names(other_mean) <- gene_list
    }
    logFC = log(cell_type_info[[1]][gene_list, cell_type] + 
                  epsilon) - log(other_mean + epsilon)
    type_gene_list = which((logFC > fc_thresh) & (cell_type_info[[1]][gene_list, 
                                                                      cell_type] > expr_thresh))
    message(paste0("get_de_genes: ", cell_type, " found DE genes: ", 
                   length(type_gene_list)))
    total_gene_list = union(total_gene_list, type_gene_list)
    
  # Z.C modified to extract the logFC and exp for DEG used in RCTD analysis
    deg_info_list[[which(cell_type == cell_type_info[[2]])]] <-  cbind(Genes = gene_list[type_gene_list], 
                logFC = logFC[type_gene_list], 
                exp = cell_type_info[[1]][gene_list,cell_type][type_gene_list],
                cell_type = rep(cell_type,length(type_gene_list))
                )
  }
  total_gene_list = gene_list[total_gene_list]
  message(paste0("get_de_genes: total DE genes: ", length(total_gene_list)))
  #return(total_gene_list)
  return(deg_info_list)  # Z.C modified to extract the logFC and exp for DEG used in RCTD analysis
}

SpatialFeaturePlotBlend_2 <- function(cells_obj, column_1, column_2, slice,ratio, combine = TRUE)  {
  
  # Convert decimal number to hexadecimal. Pad with 0s if only a single
  # character following conversion.
  as_hex <- function(num) {
    hex_str <- as.character(as.hexmode(num))
    if (nchar(hex_str) == 1) {
      hex_str <- paste0("0", hex_str)
    }
    
    return(hex_str)
  }
  
  metadata_to_hexadecimal <- function(in_dat) {
    apply(in_dat, 2,
          function(x) {
            # Make minimum 0
            x - min(x)
          }) %>%
      apply(2,
            function(x) {
              # Constrain to range [0, 255]
              round(255 * (x / max(x)))
            }) %>%
      apply(1,
            function(x) {
              # Convert to hexadecimal codes
              toupper(paste0("#", as_hex(x[1]), as_hex(x[2]), "00"))
            })
  }
  
  blend_plot_theme <- theme(legend.position = "none",
                            plot.title = element_text(hjust = 0.5))
  
  plot_list <- lapply(c(column_1, column_2),
                      function(column) {
                        max_color <- ifelse(column == column_1,
                                            "#FF0000", "#00FF00")
                        SpatialFeaturePlot(cells_obj, column) +
                          scale_fill_gradient(low = "#000000",
                                              high = max_color) +
                          ggtitle(column) +
                          blend_plot_theme
                      })
  
  dat <- FetchData(cells_obj, c(column_1, column_2))
  #dat <- FetchData(cells_obj, c(column_1, column_2), layer = "scale.data")
  dat <- scale(dat)
  colors <- as.matrix(dat) %>% metadata_to_hexadecimal()
  #print(head(colors))
  
  new_md_column <- paste0(column_1, "_vs_", column_2)
  cells_obj[[new_md_column]] <- colors
  names(colors) <- as.character(colors)
  
  plot_list[[1]] <- SpatialDimPlot(cells_obj, new_md_column, cols = colors, images = slice, pt.size.factor = 2.5, image.alpha = 0) +
    ggtitle(paste0(column_1, "/", column_2)) +
    blend_plot_theme 
  
  # Custom legend
  custom_legend <- data.frame(
    x = c(0, 1, 1,30),
    y = c(0, 5, 6,8),
    color = c("white","#FF0000", "#00FF00","white"),
    label = c(NA,column_1, column_2,NA)
  )
  
  plot_list[[2]] <- ggplot(custom_legend, aes(x = x, y = y, color = color)) +
    geom_point(shape = 19 , size = 5) +
    scale_color_identity() +
    geom_text(aes(label = label), vjust = -1) +
    theme_void() +
    theme(legend.position = "none", aspect.ratio = 1)
  
  if (combine == FALSE) {
    return(plot_list)
  } else {
    p <- wrap_plots(plot_list, nrow = 1,
                    widths = c(0.9, 0.9))
    return(p)
  }

}

SpatialFeaturePlotBlendMask_2 <- function(cells_obj, column_1, column_2,column_mask, slice,ratio, combine = TRUE)  {
  
  # Convert decimal number to hexadecimal. Pad with 0s if only a single
  # character following conversion.
  as_hex <- function(num) {
    hex_str <- as.character(as.hexmode(num))
    if (nchar(hex_str) == 1) {
      hex_str <- paste0("0", hex_str)
    }
    
    return(hex_str)
  }
  
  metadata_to_hexadecimal <- function(in_dat) {
    apply(in_dat, 2,
          function(x) {
            # Make minimum 0
            x - min(x)
          }) %>%
      apply(2,
            function(x) {
              # Constrain to range [0, 255]
              round(255 * (x / max(x)))
            }) %>%
      apply(1,
            function(x) {
              # Convert to hexadecimal codes
              toupper(paste0("#", as_hex(x[1]), as_hex(x[2]), "00"))
            })
  }
  
  blend_plot_theme <- theme(legend.position = "none",
                            plot.title = element_text(hjust = 0.5))
  
  plot_list <- lapply(c(column_1, column_2),
                      function(column) {
                        max_color <- ifelse(column == column_1,
                                            "#FF0000", "#00FF00")
                        SpatialFeaturePlot(cells_obj, column) +
                          scale_fill_gradient(low = "#000000",
                                              high = max_color) +
                          ggtitle(column) +
                          blend_plot_theme
                      })
  
  dat <- FetchData(cells_obj, c(column_1, column_2))
  #dat <- FetchData(cells_obj, c(column_1, column_2), layer = "scale.data")
  dat <- scale(dat)
  colors <- as.matrix(dat) %>% metadata_to_hexadecimal()
  #print(head(colors))
  
  # Add mask "grey"
  dat2 <- FetchData(cells_obj, c(column_mask))
  colors[dat2[,1] == 0] <- "grey"
  
  new_md_column <- paste0(column_1, "_vs_", column_2)
  cells_obj[[new_md_column]] <- colors
  names(colors) <- as.character(colors)
  
  plot_list[[1]] <- SpatialDimPlot(cells_obj, new_md_column, cols = colors, images = slice, pt.size.factor = 2.5, image.alpha = 0) +
    ggtitle(paste0(column_1, "/", column_2)) +
    blend_plot_theme 
  
  # Custom legend
  custom_legend <- data.frame(
    x = c(0, 1, 1,30),
    y = c(0, 5, 6,8),
    color = c("white","#FF0000", "#00FF00","white"),
    label = c(NA,column_1, column_2,NA)
  )
  
  plot_list[[2]] <- ggplot(custom_legend, aes(x = x, y = y, color = color)) +
    geom_point(shape = 19 , size = 5) +
    scale_color_identity() +
    geom_text(aes(label = label), vjust = -1) +
    theme_void() +
    theme(legend.position = "none", aspect.ratio = 1)
  
  if (combine == FALSE) {
    return(plot_list)
  } else {
    p <- wrap_plots(plot_list, nrow = 1,
                    widths = c(0.9, 0.9))
    return(p)
  }
  
}

SpatialFeaturePlotBlend_3 <- function(cells_obj, column_1, column_2, column_3, slice,ratio, combine = TRUE)  {
  
  # Convert decimal number to hexadecimal. Pad with 0s if only a single
  # character following conversion.
  as_hex <- function(num) {
    hex_str <- as.character(as.hexmode(num))
    if (nchar(hex_str) == 1) {
      hex_str <- paste0("0", hex_str)
    }
    
    return(hex_str)
  }
  
  metadata_to_hexadecimal <- function(in_dat) {
    apply(in_dat, 2,
          function(x) {
            # Make minimum 0
            x - min(x)
          }) %>%
      apply(2,
            function(x) {
              # Constrain to range [0, 255]
              round(255 * (x / max(x)))
            }) %>%
      apply(1,
            function(x) {
              # Convert to hexadecimal codes
              toupper(paste0("#", as_hex(x[1]), as_hex(x[2]), as_hex(x[3])))
            })
  }
  
  blend_plot_theme <- theme(legend.position = "none",
                            plot.title = element_text(hjust = 0.5),
                            aspect.ratio = ratio)
  
  plot_list <- list()
  dat <- FetchData(cells_obj, c(column_1, column_2, column_3))
  #dat <- FetchData(cells_obj, c(column_1, column_2, column_3), layer = "scale.data")
  dat <- scale(dat)
  colors <- as.matrix(dat) %>% metadata_to_hexadecimal()
  
  new_md_column <- paste0(column_1, "_vs_", column_2, "_vs_", column_3)
  cells_obj[[new_md_column]] <- colors
  names(colors) <- as.character(colors)
  
  plot_list[[1]] <- SpatialDimPlot(cells_obj, new_md_column, cols = colors, images = slice, pt.size.factor = 2.5, image.alpha = 0) +
    ggtitle(paste0(column_1, "/", column_2, "/", column_3)) +
    blend_plot_theme 
  
  # Custom legend
  custom_legend <- data.frame(
    x = c(0,1, 1, 1,30),
    y = c(0,4, 5, 6,8),
    color = c("white","#FF0000", "#00FF00", "#0000FF","white"),
    label = c(NA,column_1, column_2, column_3,NA)
  )
  
  plot_list[[2]] <- ggplot(custom_legend, aes(x = x, y = y, color = color)) +
    geom_point(shape = 19 , size = 5) +
    scale_color_identity() +
    geom_text(aes(label = label), vjust = -1) +
    theme_void() +
    theme(legend.position = "none", aspect.ratio = 1)

  if (combine == FALSE) {
    return(plot_list)
  } else {
    p <- wrap_plots(plot_list, nrow = 1,
                    widths = c(0.9, 0.9))
    return(p)
  }
  
}

SpatialFeaturePlotBlendMask_3 <- function(cells_obj, column_1, column_2, column_3, column_mask, slice,ratio, combine = TRUE)  {
  
  # Convert decimal number to hexadecimal. Pad with 0s if only a single
  # character following conversion.
  as_hex <- function(num) {
    hex_str <- as.character(as.hexmode(num))
    if (nchar(hex_str) == 1) {
      hex_str <- paste0("0", hex_str)
    }
    
    return(hex_str)
  }
  
  metadata_to_hexadecimal <- function(in_dat) {
    apply(in_dat, 2,
          function(x) {
            # Make minimum 0
            x - min(x)
          }) %>%
      apply(2,
            function(x) {
              # Constrain to range [0, 255]
              round(255 * (x / max(x)))
            }) %>%
      apply(1,
            function(x) {
              # Convert to hexadecimal codes
              toupper(paste0("#", as_hex(x[1]), as_hex(x[2]), as_hex(x[3])))
            })
  }
  
  blend_plot_theme <- theme(legend.position = "none",
                            plot.title = element_text(hjust = 0.5),
                            aspect.ratio = ratio)
  
  plot_list <- list()
  dat <- FetchData(cells_obj, c(column_1, column_2, column_3))
  #dat <- FetchData(cells_obj, c(column_1, column_2, column_3), layer = "scale.data")
  dat <- scale(dat)
  colors <- as.matrix(dat) %>% metadata_to_hexadecimal()
  
  # Add mask "grey"
  dat2 <- FetchData(cells_obj, c(column_mask))
  colors[dat2[,1] == 0] <- "grey"
  
  new_md_column <- paste0(column_1, "_vs_", column_2, "_vs_", column_3)
  cells_obj[[new_md_column]] <- colors
  names(colors) <- as.character(colors)
  
  plot_list[[1]] <- SpatialDimPlot(cells_obj, new_md_column, cols = colors, images = slice, pt.size.factor = 2.5, image.alpha = 0) +
    ggtitle(paste0(column_1, "/", column_2, "/", column_3)) +
    blend_plot_theme 
  
  # Custom legend
  custom_legend <- data.frame(
    x = c(0,1, 1, 1,30),
    y = c(0,4, 5, 6,8),
    color = c("white","#FF0000", "#00FF00", "#0000FF","white"),
    label = c(NA,column_1, column_2, column_3,NA)
  )
  
  plot_list[[2]] <- ggplot(custom_legend, aes(x = x, y = y, color = color)) +
    geom_point(shape = 19 , size = 5) +
    scale_color_identity() +
    geom_text(aes(label = label), vjust = -1) +
    theme_void() +
    theme(legend.position = "none", aspect.ratio = 1)
  
  if (combine == FALSE) {
    return(plot_list)
  } else {
    p <- wrap_plots(plot_list, nrow = 1,
                    widths = c(0.9, 0.9))
    return(p)
  }
  
}
