# linear regression of GeoMx gene and cell-types x Neuromelanin quantification.

library(readxl)
library(edgeR)
library(dplyr)
library(plyr)
library(EnhancedVolcano)
library(ggplot2)
library(ggridges)
library(MASS)
library(ggplot2)
library(viridis)
library(ggpubr)
library(EnvStats)
library(rstatix)

theme_set(theme_minimal())


## inputs
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis"
meta_path <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/Annotations_NM.xlsx"
geomx_annot_path = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/Annotations_remapped.xlsx"
NM_data_dir = "/Users/zacc/USyd/spatial_transcriptomics/analysis/NM_data_new"
#NM_data_dir = "/Users/zacc/USyd/spatial_transcriptomics/analysis/NM_data"

##############################################################
### Part 1: Extract NM quantification's and add to ST data ###
##############################################################
# list of NM quant files
file.names <- list.files(path=NM_data_dir,full.names=T)
ag_list <- list()
for (i in 1:length(file.names)){
  print(length(file.names)-i)
  print(file.names[i])
  df <- read_xlsx(file.names[i], 1)
  # remove rows without ROI eg. summary rows.
  df <- df[!is.na(df$ROI),]
  # add Brainbank_ID
  df$Brainbank_ID <- rep(gsub("_.*","",basename(file.names[i])),nrow(df))
  # select cols
  #df <- df[,c("Area [µm²]","ROI","Perimeter [µm]","Mean (Radius) [µm]","Mean (Gray Intensity Value)","Brainbank_ID")]
  df <- df[,c("Area [µm²]","ROI","Mean (Gray Intensity Value)","Brainbank_ID")]

  # process totals 
  df_roi <- read_xlsx(file.names[i], 3)
  df_roi <-  df_roi[!is.na(df_roi$ROI),]
  # select cols
  df_roi <- df_roi[,c("ROI","Sum (Area) [µm²]","Object Area Fraction ROI [%]","Mean (Mean (Gray Intensity Value))","Mean (Shape Factor)","ROI Area [µm²]")]
  # aggregate over ROI
  tmp <- aggregate(cbind(`Sum (Area) [µm²]`,`Mean (Mean (Gray Intensity Value))`,`ROI Area [µm²]`) ~ ROI, data = df_roi, FUN = sum, na.rm = TRUE)
  # mapvalues
  df$`Sum (Area) [µm²]` <- mapvalues(df$ROI,tmp$ROI,tmp$`Sum (Area) [µm²]`)
  df$`Mean (Mean (Gray Intensity Value))` <- mapvalues(df$ROI,tmp$ROI,tmp$`Mean (Mean (Gray Intensity Value))`)
  df$`ROI Area [µm²]` <- mapvalues(df$ROI,tmp$ROI,tmp$`ROI Area [µm²]`)
  
  # add to list
  ag_list[[i]] <- df
}

# collapse to dataframe
df_agg <- as.data.frame(do.call("rbind", ag_list))

# # read in meta data
# setwd(analysis_dir)
# load(rdata)
# meta <- gxdat_s@meta.data
meta <- read_excel(meta_path)

# check if any metadata is missing
unique(df_agg$Brainbank_ID)[!unique(df_agg$Brainbank_ID) %in% meta$Brainbank_ID]

# # create brainbank id and region identifiers
# meta$brainbank_region <- paste0(meta$Brainbank_ID,sep="_",meta$Brainregion_2)
# df_agg$brainbank_region <- paste0(df_agg$Brainbank_ID,sep="_",df_agg$ROI)

# # merge and format
# meta_simple <- meta[!duplicated(meta$brainbank_region),] # get just simple meta for NM assessments first
df_agg_merged <- merge(df_agg,meta,by="Brainbank_ID",all.x=TRUE)
df_agg_merged$Diagnosis <- factor(df_agg_merged$Diagnosis,levels=c("CTR","ILBD", "ePD","lPD"))
df_agg_merged <- df_agg_merged[df_agg_merged$ROI != "ROI 1",]
df_agg_merged$log10_Area <- log10(df_agg_merged$`Area [µm²]`)
  
# format age group
df_agg_merged$Age <- as.numeric(df_agg_merged$Age)
dplot <- df_agg_merged
dplot <- dplot %>% 
  mutate(
    # Create categories
    age_group = dplyr::case_when(
      Age <= 70            ~ "<70",
      Age > 70 & Age <= 80 ~ "70-80",
      Age > 80 & Age <= 90 ~ "80-90",
      Age > 90             ~ "> 90"
    ),
    # Convert to factor
    age_group = factor(
      age_group,
      level = c("<70", "70-80","80-90", "> 90")
    )
  )
df_agg_merged <- dplot

# histograms of cohort
dplot <- df_agg_merged
#dplot$Age <- factor(dplot$Age,levels=as.character(unique(sort(as.numeric(dplot$Age)))))
pa <- ggplot(dplot, aes(x = log10(`Area [µm²]`))) +
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins =100)+
  geom_density(alpha=.2, fill="grey") + xlim(0,4)

dplot <- df_agg_merged
#dplot$Age <- factor(dplot$Age,levels=as.character(unique(sort(as.numeric(dplot$Age)))))
p1 <- ggplot(dplot, aes(x = log10(`Area [µm²]`), y = age_group)) +
  geom_density() + theme_minimal() +  theme(legend.position = "none") +
  geom_density_ridges(aes(fill = age_group))  + ylab("Age Group")

dplot <- df_agg_merged
p2 <- ggplot(dplot, aes(x = log10(`Area [µm²]`), y=Diagnosis)) +
  geom_density() + theme_minimal() + theme(legend.position = "none") +
  geom_density_ridges(aes(fill = Diagnosis)) 

dplot <- df_agg_merged
p3 <- ggplot(dplot, aes(x = log10(`Area [µm²]`), y=ROI)) +
  geom_density() + theme_minimal() + ylab("Brain Region") + theme(legend.position = "none") +
  geom_density_ridges(aes(fill = ROI)) 

# colour palettes
Diagnosis_col = c("grey","#00AFBB", "#E7B800","red")
age_group_col = viridis(4)
ROI_col = viridis(6)

# Violin plot function
violin_plot_function <- function(data_table, x_variable, y_variable, x_lab, y_lab, colour_palette) {
  # make violin plot
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none")
  
  # perform t-test
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  stat.test <- data_table %>% t_test(y ~ x) %>% add_significance("p.adj")
  stat.test <- stat.test[stat.test$p.adj < 0.05,]
  print(stat.test)
  
  # add p-values to plot and save the result
  bxp <- bxp + stat_pvalue_manual(stat.test,
                                  y.position = max(data_table$y) * 1.4, step.increase = 0.1,
                                  label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
  
  # return object
  return(bxp)
}

# plots
v1 <- violin_plot_function(dplot,"Diagnosis","log10_Area","Diagnosis","Log10(Area [µm²])", Diagnosis_col)
v2 <- violin_plot_function(dplot,"age_group","log10_Area","age_group","Log10(Area [µm²])", age_group_col)
v3 <- violin_plot_function(dplot,"ROI","log10_Area","Brain Region","Log10(Area [µm²])", ROI_col)


arrange <- ggarrange(plotlist=list(pa,p1,p2,p3), nrow=2, ncol=2, widths = c(2,2))
ggsave("distributions_NMarea.png", arrange)
arrange <- ggarrange(plotlist=list(pa,v1,v2,v3), nrow=2, ncol=2, widths = c(2,2))
ggsave("violin_NMarea_DxAgeRegion.png", arrange)


# run anova
fit1 <- lm(log10(`Area [µm²]`) ~ Age + Diagnosis + ROI + Sex, data = df_agg_merged)
a <- anova(fit1)
nfac <- length(a[,1])-1
maxval = 100

names_2 = c("Age","Diagnosis","Brainregion","Sex")

pdf(file="anova_con.var_area.pdf")
op <- par(cex.main = 1.5, mar = c(18, 6, 4, 3) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1 , font.lab = 1, cex.axis = 1, bty = "n", las = 1)
barplot(100*a$"Sum Sq"[1:nfac]/sum(a$"Sum Sq"[1:nfac]),names.arg=names_2,ylim=c(0,maxval),las=3, 
        ylab="Contribution to Variance (%)", border = FALSE,cex.names=0.9,cex.axis=0.9)
abline(h=2,lty=2, col="red",lwd=1.5)
dev.off()


#######################################################
### Evalutate the number of modes in a distribution ###
#######################################################
#var_interest <- unique(df_agg_merged$ROI)
#var_interest <- unique(df_agg_merged$Diagnosis)
var_interest = "total"

for (i in 1:length(var_interest)) {
  print(i)
  x <- log10(df_agg_merged$`Area [µm²]`)
  #x <- log10(df_agg_merged$`Area [µm²]`[df_agg_merged$ROI == var_interest[i]])
  #x <- log10(df_agg_merged$`Area [µm²]`[df_agg_merged$Diagnosis == var_interest[i]])
  
  # Find the modes of a KDE: it assumes no mode spans more than one x value.)
  findmodes <- function(kde) {
    kde$x[which(c(kde$y[-1],NA) < kde$y & kde$y > c(NA,kde$y[-length(kde$y)]))]
  }

  # Compute the mode trace by varying the bandwidth within a factor of 10 of
  # the default bandwidth.  Track the modes as the bandwidth is decreased from
  # its largest to its smallest value.
  # This calculation is fast, so we can afford a detailed search.

  m <- mean(x)
  id <- 1
  bw <- density(x)$bw * 10^seq(1,-1, length.out=101) 
  modes.lst <- lapply(bw, function(h) {
    m.new <- sort(findmodes(density(x, bw=h)))
    # -- Associate each previous mode with a nearest new mode.
    if (length(m.new)==1) delta <- Inf else delta <- min(diff(m.new))/2
    d <- outer(m.new, m, function(x,y) abs(x-y))
    i <- apply(d, 2, which.min)
    g <- rep(NA_integer_, length(m.new))
    g[i] <- id[1:ncol(d)]
    #-- Create new ids for new modes that appear.
    k <- is.na(g)
    g[k] <- (sum(!k)+1):length(g)
    id <<- g
    m <<- m.new
    data.frame(bw=h, Mode=m.new, id=g)
  })
  X <- do.call(rbind, args=modes.lst)
  X$id <- factor(X$id)

  # Locate the modes at the most vertical portions of their traces.
  minslope <- function(x, y) {
    f <- splinefun(x, y)
    e <- diff(range(x)) * 1e-4
    df2 <- function(x) ((f(x+e)-f(x-e)) / (2*e))^2 # Numerical derivative, squared
    v <- optimize(df2, c(min(x),max(x)))
    c(bw=v$minimum, slope=v$objective, Mode=f(v$minimum))
  }
  # Retain the desired modes.
  n.modes <- 2 # USER SELECTED: Not automatic
  bw.max <- max(subset(X, id==n.modes)$bw)
  modes <- sapply(1:n.modes, function(i) {
    Y <- subset(X, id==i & bw <= bw.max)
    minslope(Y$bw, Y$Mode)
  })
  #
  print(modes)
  # Plot the results.
  g1 <- ggplot(X, aes(bw, Mode)) +
    geom_line(aes(col=id), size=1.2, show.legend=FALSE) +
    geom_point(aes(bw, Mode), data=as.data.frame(t(modes)), size=3, col="Black", alpha=1/2) +
    scale_x_log10() + ylab("Mode") + xlab("Bandwidth") +
    coord_flip() + ylim(0.5,3.5) +
    ggtitle(var_interest[i])
  
  
  g2 <- ggplot(data.frame(x), aes(x, ..density..)) +
    geom_histogram(aes(y=..density..), colour="black", fill="white", bins =100)+
    geom_density(alpha=.2, fill="grey") + xlim(0.5,3.5) +
    geom_vline(data=as.data.frame(t(modes)),
               mapping=aes(xintercept=Mode), col="#D18A4e", size=1) +
    xlab("log10(`Area [µm²]`)") + ylab("Density") 
  
  arrange <- ggarrange(plotlist=list(g1,g1,g2,g2), nrow=2, ncol=2, widths = c(2,2))
  ggsave(paste0("Mode_Trace_",var_interest[i],".png"), arrange)
}

# select 2 modes/ clusters
#######################################################

#############################################
### Assigning intra and extra-cellular NM ###
#############################################
x <- log10(df_agg_merged$`Area [µm²]`)

# k-means clustering on log10 reads
set.seed(20)
clusters <- kmeans(x, 2)
df_agg_merged$cluster <- clusters$cluster

# calculate 95% confidence intervals for normal distribution of highest K-cluster by log10(unique reads)
# fit normal distribution
extra_NM_cluster = 2

fit <- fitdistr(x[df_agg_merged$cluster == extra_NM_cluster], "normal")
para <- fit$estimate
a <- para[[1]]
s <- para[[2]]

# select 95% CI threshold
left_extra <- a-(1.96 * s)
right_extra <- a+(1.96 * s)

# fit normal distribution
intra_NM_cluster = 1

fit <- fitdistr(x[df_agg_merged$cluster == intra_NM_cluster], "normal")
para <- fit$estimate
a <- para[[1]]
s <- para[[2]]

# select 95% CI threshold
left_intra <- a-(1.96 * s)
right_intra <- a+(1.96 * s)

# assign clusters
df_agg_merged$intra.extra[df_agg_merged$cluster == 2 & x < left_intra] <- "eNM"
df_agg_merged$intra.extra[df_agg_merged$cluster == 1 & x > right_extra] <- "iNM"

# plot data

p3 <- df_agg_merged %>%
  ggplot(aes(log10(`Area [µm²]`), ..scaled.., fill=factor(intra.extra))) +
  geom_histogram(bins=20,aes(y=..count../sum(..count..))) +
  scale_fill_manual(values=c("grey","red","black")) +
  theme_minimal() + 
  xlim(0.5,3.5) + 
  labs(fill="") + geom_vline(xintercept = left_intra, color = "red", lty=2) +
  geom_vline(xintercept = right_extra, color = "grey", lty=2) +
  geom_density(aes(group=1)) 

ggsave("dist_iNMeNM.png", p3)


# outliers - Rosner test
test <- rosnerTest(log10(df_agg_merged$`Area [µm²]`),
                   k = 3
)
test
test <- rosnerTest(log10(df_agg_merged$`Area [µm²]`[df_agg_merged$intra.extra == "iNM"]),
                   k = 3
)
test
test <- rosnerTest(log10(df_agg_merged$`Area [µm²]`[df_agg_merged$intra.extra == "eNM"]),
                   k = 3
)
test


############################
#### COUNTS of NM
############################
dplot <- df_agg_merged %>% dplyr::count(intra.extra,Diagnosis,Brainbank_ID, ROI,age_group,`ROI Area [µm²]`, sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$n_per.µm2 <- dplot$n / as.numeric(dplot$`ROI Area [µm²]`)

# colour palettes
Diagnosis_col = c("grey","#00AFBB", "#E7B800","red")
age_group_col = viridis(4)
ROI_col = viridis(6)

# plots
d2 <- dplot[dplot$intra.extra == "iNM",]
v1 <- violin_plot_function(d2,"Diagnosis","n","Diagnosis","iNM (n)", Diagnosis_col)
v2 <- violin_plot_function(d2,"ROI","n","Brain Region","iNM (n)", ROI_col)
v3 <- violin_plot_function(d2,"age_group","n","Age Group","iNM (n)", age_group_col)

d2 <- dplot[dplot$intra.extra == "eNM",]
v4 <- violin_plot_function(d2,"Diagnosis","n","Diagnosis","eNM (n)", Diagnosis_col)
v5 <- violin_plot_function(d2,"ROI","n","Brain Region","eNM (n)", ROI_col)
v6 <- violin_plot_function(d2,"age_group","n","Age Group","eNM (n)", age_group_col)

arrange <- ggarrange(plotlist=list(v1,v2,v3), nrow=2, ncol=2, widths = c(2,2))
ggsave("violin_iNM_counts.png", arrange)

arrange <- ggarrange(plotlist=list(v4,v5,v6), nrow=2, ncol=2, widths = c(2,2))
ggsave("violin_eNM_counts.png", arrange)

d2 <- dplot[dplot$intra.extra == "iNM" & dplot$Diagnosis == "CTR",]
v2 <- violin_plot_function(d2,"ROI","n","Brain Region","iNM (n)", ROI_col)
v2b <- violin_plot_function(d2,"ROI","n_per.µm2","Brain Region","iNM (n/µm²)", ROI_col)

d2 <- dplot[dplot$intra.extra == "eNM" & dplot$Diagnosis == "CTR",]
v5 <- violin_plot_function(d2,"ROI","n","Brain Region","eNM (n)", ROI_col)
v5b <- violin_plot_function(d2,"ROI","n_per.µm2","Brain Region","eNM (n/µm²)", ROI_col)

arrange <- ggarrange(plotlist=list(v2,v5,v2b,v5b), nrow=2, ncol=2, widths = c(2,2))
ggsave("violin_iNMeNM_counts_controls.png", arrange)


# plots - normalised to Area
d2 <- dplot[dplot$intra.extra == "iNM",]
v1b <- violin_plot_function(d2,"Diagnosis","n_per.µm2","Diagnosis","iNM (n/µm²)", Diagnosis_col)
v2b <- violin_plot_function(d2,"ROI","n_per.µm2","Brain Region","iNM (n/µm²)", ROI_col)
v3b <- violin_plot_function(d2,"age_group","n_per.µm2","Age Group","iNM (n/µm²)", age_group_col)

d2 <- dplot[dplot$intra.extra == "eNM",]
v4b <- violin_plot_function(d2,"Diagnosis","n_per.µm2","Diagnosis","eNM (n/µm²)", Diagnosis_col)
v5b <- violin_plot_function(d2,"ROI","n_per.µm2","Brain Region","eNM (n/µm²)", ROI_col)
v6b <- violin_plot_function(d2,"age_group","n_per.µm2","Age Group","eNM (n/µm²)", age_group_col)


# ratio intra/extra
tmp <- df_agg_merged[df_agg_merged$intra.extra == "iNM",] %>% dplyr::count(Diagnosis,Brainbank_ID, ROI,age_group, sort = TRUE)
tmp$brainbank_region <- paste0(tmp$Brainbank_ID,"_",tmp$ROI)
tmp1 <- tmp

tmp <- df_agg_merged[df_agg_merged$intra.extra == "eNM",] %>% dplyr::count(Diagnosis,Brainbank_ID, ROI,age_group, sort = TRUE)
tmp$brainbank_region <- paste0(tmp$Brainbank_ID,"_",tmp$ROI)
tmp2 <- tmp

tmp3 <- merge(tmp1,tmp2,by="brainbank_region")
tmp3$intra.extra.ratio.µm2 <- tmp3$n.x/tmp3$n.y
tmp3 <- tmp3[!is.na(tmp3$Diagnosis.x),]

v7 <- violin_plot_function(tmp3,"Diagnosis.x","intra.extra.ratio.µm2","Diagnosis","iNM (n/µm²) / eNM (n/µm²)", Diagnosis_col)
v8 <- violin_plot_function(tmp3,"ROI.x","intra.extra.ratio.µm2","Brain Region","iNM (n/µm²) / eNM (n/µm²)", ROI_col)
v9 <- violin_plot_function(tmp3,"age_group.x","intra.extra.ratio.µm2","Age Group","iNM (n/µm²) / eNM (n/µm²)", age_group_col)


arrange <- ggarrange(plotlist=list(v1b,v2b,v3b), nrow=2, ncol=2, widths = c(2,2))
ggsave("violin_iNM_counts.µm2.png", arrange)

arrange <- ggarrange(plotlist=list(v4b,v5b,v6b), nrow=2, ncol=2, widths = c(2,2))
ggsave("violin_eNM_counts.µm2.png", arrange)

arrange <- ggarrange(plotlist=list(v7,v8,v9), nrow=2, ncol=2, widths = c(2,2))
ggsave("violin_iNM.eNMratio_counts.png", arrange)


############################################################################################
###### Part 2: Area and Optical density of iNM and eNM
############################################################################################

dplot <- df_agg_merged[df_agg_merged$intra.extra == "iNM",]
dplot <- dplot[!is.na(dplot$intra.extra),]

v1 <- violin_plot_function(dplot,"Diagnosis","log10_Area","Diagnosis","iNM Log10(Area [µm²])", Diagnosis_col)
v2 <- violin_plot_function(dplot,"age_group","log10_Area","age_group","iNM Log10(Area [µm²])", age_group_col)
v3 <- violin_plot_function(dplot,"ROI","log10_Area","Brain Region","iNM Log10(Area [µm²])", ROI_col)

dplot <- df_agg_merged[df_agg_merged$intra.extra == "eNM",]
dplot <- dplot[!is.na(dplot$intra.extra),]

v4 <- violin_plot_function(dplot,"Diagnosis","log10_Area","Diagnosis","eNM Log10(Area [µm²])", Diagnosis_col)
v5 <- violin_plot_function(dplot,"age_group","log10_Area","age_group","eNM Log10(Area [µm²])", age_group_col)
v6 <- violin_plot_function(dplot,"ROI","log10_Area","Brain Region","eNM Log10(Area [µm²])", ROI_col)

arrange <- ggarrange(plotlist=list(v1,v2,v3), nrow=2, ncol=2, widths = c(2,2))
ggsave("violin_iNM_area.png", arrange)

arrange <- ggarrange(plotlist=list(v4,v5,v6), nrow=2, ncol=2, widths = c(2,2))
ggsave("violin_eNM_area.png", arrange)

##

dplot <- df_agg_merged[df_agg_merged$intra.extra == "iNM",]
dplot <- dplot[!is.na(dplot$intra.extra),]
colnames(dplot) <- make.names(colnames(dplot))

v1 <- violin_plot_function(dplot,"Diagnosis","Mean..Gray.Intensity.Value.","Diagnosis","iNM Optical Density", Diagnosis_col)
v2 <- violin_plot_function(dplot,"age_group","Mean..Gray.Intensity.Value.","age_group","iNM Optical Density", age_group_col)
v3 <- violin_plot_function(dplot,"ROI","Mean..Gray.Intensity.Value.","Brain Region","iNM Optical Density", ROI_col)

dplot <- df_agg_merged[df_agg_merged$intra.extra == "eNM",]
dplot <- dplot[!is.na(dplot$intra.extra),]
colnames(dplot) <- make.names(colnames(dplot))

v4 <- violin_plot_function(dplot,"Diagnosis","Mean..Gray.Intensity.Value.","Diagnosis","eNM Optical Density", Diagnosis_col)
v5 <- violin_plot_function(dplot,"age_group","Mean..Gray.Intensity.Value.","age_group","eNM Optical Density", age_group_col)
v6 <- violin_plot_function(dplot,"ROI","Mean..Gray.Intensity.Value.","Brain Region","eNM Optical Density", ROI_col)

arrange <- ggarrange(plotlist=list(v1,v2,v3), nrow=2, ncol=2, widths = c(2,2))
ggsave("violin_iNM_opticaldensity.png", arrange)

arrange <- ggarrange(plotlist=list(v4,v5,v6), nrow=2, ncol=2, widths = c(2,2))
ggsave("violin_eNM_opticaldensity.png", arrange)



############################################################################################
###### Add NM quantifications to GeoMx Annotations
############################################################################################
annot_geomx <- read_excel(geomx_annot_path)

# counts per µm²
dplot <- df_agg_merged %>% dplyr::count(intra.extra,Diagnosis,Brainbank_ID,ROI,age_group,`ROI Area [µm²]`, sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$n_per.µm2 <- dplot$n / as.numeric(dplot$`ROI Area [µm²]`)
dplot$brainbank_region <- paste0(dplot$Brainbank_ID,"_",dplot$ROI)

iNM_n_per.µm2 <- dplot[dplot$intra.extra == "iNM",]
eNM_n_per.µm2 <- dplot[dplot$intra.extra == "eNM",]

# mean area
dplot <- df_agg_merged

tmp <- as.data.frame(dplot %>% dplyr::group_by(Brainbank_ID,ROI, intra.extra) %>%
  dplyr::summarize(NM_mean.area = mean(`Area [µm²]`, na.rm=TRUE)))

tmp <- tmp[!is.na(tmp$intra.extra),]

iNM_NM_mean.area <- tmp[tmp$intra.extra == "iNM",]
eNM_NM_mean.area <- tmp[tmp$intra.extra == "eNM",]

# mean optical density
tmp <- as.data.frame(dplot %>% dplyr::group_by(Brainbank_ID,ROI, intra.extra) %>%
                       dplyr::summarize(NM_mean.optical.density = mean(`Mean (Gray Intensity Value)`, na.rm=TRUE)))

tmp <- tmp[!is.na(tmp$intra.extra),]

iNM_NM_mean.optical.density <- tmp[tmp$intra.extra == "iNM",]
eNM_NM_mean.optical.density <- tmp[tmp$intra.extra == "eNM",]

# combine
iNM_n_per.µm2$Brainbank_ID == eNM_n_per.µm2$Brainbank_ID
tmp <- merge()


#put all data frames into list
df_list <- list(iNM_n_per.µm2,eNM_n_per.µm2,
                iNM_NM_mean.area,eNM_NM_mean.area,
                iNM_NM_mean.optical.density,eNM_NM_mean.optical.density)

#merge all data frames in list
tmp <- Reduce(function(x, y) merge(x, y, by=c("Brainbank_ID","ROI"),all=TRUE), df_list)
tmp <- tmp[,c("Brainbank_ID","ROI",
              "ROI Area [µm²].x","n.x","n_per.µm2.x",
              "ROI Area [µm²].y","n.y","n_per.µm2.y",
              "NM_mean.area.x","NM_mean.area.y",
              "NM_mean.optical.density.x","NM_mean.optical.density.y")]
colnames(tmp) <- gsub("\\.x",".iNM",colnames(tmp))
colnames(tmp) <- gsub("\\.y",".eNM",colnames(tmp))

# add to annotatations and write to file
tmp_m <- merge(annot_geomx, tmp, by = c("Brainbank_ID","ROI"), all.x = TRUE)
tmp_m <- tmp_m[order(tmp_m$Sample_ID),]
writexl::write_xlsx(tmp_m,path = paste0(analysis_dir,"/Annotations_nm.xlsx"))


