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

theme_set(theme_minimal())


## inputs
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis"
rdata = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis/geomx_oct2023_seurat.Rdata"
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

# read in geomx data
setwd(analysis_dir)
load(rdata)
meta <- gxdat_s@meta.data

# check if any brainbank ID data missing
unique(df_agg$Brainbank_ID)[!unique(df_agg$Brainbank_ID) %in% meta$Brainbank_ID]
table(meta$Brainbank_ID[!meta$Brainbank_ID %in% df_agg$Brainbank_ID])
table(meta$Diagnosis[!meta$Brainbank_ID %in% df_agg$Brainbank_ID],meta$Brainbank_ID[!meta$Brainbank_ID %in% df_agg$Brainbank_ID])


# create brainbank id and region identifiers
meta$brainbank_region <- paste0(meta$Brainbank_ID,sep="_",meta$Brainregion_2)
df_agg$brainbank_region <- paste0(df_agg$Brainbank_ID,sep="_",df_agg$ROI)

# merge
meta_simple <- meta[!duplicated(meta$brainbank_region),] # get just simple meta for NM assessments first
df_agg_merged <- merge(df_agg,meta_simple,by="brainbank_region")
df_agg_merged$Diagnosis <- factor(df_agg_merged$Diagnosis,levels=c("CTR","ILBD", "ePD","lPD"))
 
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
p3 <- ggplot(dplot, aes(x = log10(`Area [µm²]`), y=Brainregion_2)) +
  geom_density() + theme_minimal() + ylab("Brain Region") + theme(legend.position = "none") +
  geom_density_ridges(aes(fill = Brainregion_2)) 

dplot <- df_agg_merged
p4 <- ggplot(dplot, aes(x = Diagnosis, y=log10(`Area [µm²]`), fill=Diagnosis)) +
  geom_boxplot(position=position_dodge(1)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  scale_fill_brewer(palette='viridis') + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) + 
  stat_compare_means( label = "p.signif") 


dplot <- df_agg_merged[!duplicated(df_agg_merged$Brainbank_ID.x),]
p5 <- ggplot(dplot, aes(x = Diagnosis, y=Age, fill=Diagnosis)) +
  geom_boxplot(position=position_dodge(1)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  scale_fill_brewer(palette='viridis') + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) + 
  stat_compare_means( label = "p.signif") 


arrange <- ggarrange(plotlist=list(pa,p1,p2,p3), nrow=2, ncol=2, widths = c(2,2))
ggsave("distributions_NMarea.png", arrange)
arrange <- ggarrange(plotlist=list(p4,p5), nrow=2, ncol=2, widths = c(2,2))
ggsave("boxplot_NMarea_DxAgeRegion.png", arrange)


# run anova
fit1 <- lm(log10(`Area [µm²]`) ~ Age + Diagnosis + Brainregion_2 + Sex, data = df_agg_merged)
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

#var_interest <- unique(df_agg_merged$Brainregion_2)
var_interest <- unique(df_agg_merged$Diagnosis)
var_interest = "total"

for (i in 1:length(var_interest)) {
  print(i)
  #x <- log10(df_agg_merged$`Area [µm²]`)
  #x <- log10(df_agg_merged$`Area [µm²]`[df_agg_merged$Brainregion_2 == var_interest[i]])
  x <- log10(df_agg_merged$`Area [µm²]`[df_agg_merged$Diagnosis == var_interest[i]])
  
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

# aggregate counts
library(ggpubr)
dplot <- df_agg_merged %>% dplyr::count(intracellular, Diagnosis, brainbank_region, Brainregion_2,sort = TRUE)
dplot <- dplot[dplot$intracellular == TRUE,]
p4 <- ggplot(dplot, aes(x = Brainregion_2, y=n, fill=Diagnosis)) +
  geom_boxplot(position=position_dodge(1)) + 
  #facet_grid( ~ Brainregion_2) +
  scale_fill_brewer(palette='viridis') + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) + 
  stat_compare_means( label = "p.signif") + ggtitle("intracellular")

dplot <- df_agg_merged %>% dplyr::count(extracellular, Diagnosis, brainbank_region, Brainregion_2,sort = TRUE)
dplot <- dplot[dplot$extracellular == TRUE,]
p5 <- ggplot(dplot, aes(x = Brainregion_2, y=n, fill=Diagnosis)) +
  geom_boxplot(position=position_dodge(1)) + 
  #facet_grid( ~ Brainregion_2) +
  scale_fill_brewer(palette='viridis') + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) + 
  stat_compare_means( label = "p.signif") + ggtitle("extracellular")

# ratio intra/extra
tmp <- df_agg_merged %>% dplyr::count(intracellular, Diagnosis, brainbank_region, Brainregion_2,sort = TRUE)
tmp <- tmp[tmp$intracellular == "TRUE",]
tmp2 <- df_agg_merged %>% dplyr::count(extracellular, Diagnosis, brainbank_region, Brainregion_2,sort = TRUE)
tmp2 <- tmp2[tmp2$extracellular == "TRUE",]
tmp3 <- merge(tmp,tmp2,by="brainbank_region")
tmp3$intra.extra.ratio <- tmp3$n.x/tmp3$n.y
dplot <- tmp3

p6 <- ggplot(dplot, aes(x = Brainregion_2.x, y=intra.extra.ratio, fill=Diagnosis.x)) +
  geom_boxplot(position=position_dodge(1)) + 
  #facet_grid( ~ Brainregion_2) +
  scale_fill_brewer(palette='viridis') + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) + 
  stat_compare_means( label = "p.signif") + ggtitle("intra:extracellular ratio")

ggsave("NM_dist_age.png",p1,device = "png")
ggsave("NM_dist_dx.png",p2,device = "png")
ggsave("NM_size_dx.png",p2b,device = "png")
ggsave("NM_define_intra.extra.png",p3,device = "png")
ggsave("NM_quant_compare_intra.png",p4,device = "png")
ggsave("NM_quant_compare_extra.png",p5,device = "png")
ggsave("NM_quant_compare_intra.extra.png",p6,device = "png")



############################################################################################
###### Part 2: Regression analysis of C-SIDE results; brain-region x disease stage
############################################################################################
# get C-SIDE and PAGE results and format
setwd(analysis_dir)
load(rdata)
barcodes <- colnames(myRCTD@spatialRNA@counts)
norm_weights <- normalize_weights(myRCTD@results$weights)
meta <- gxdat_s@meta.data[gxdat_s@meta.data$segment.y != "TH" & gxdat_s$Diagnosis != "NTC",]
meta <- cbind(meta,norm_weights[row.names(meta),])

# setup design matrix and perform contrast
brain_regions <- c("A10","A9","A6")
cell_types <- colnames(norm_weights)
res_list <- list()
for (i in 1:length(brain_regions)){
  run_index <- meta$Brainregion %in% brain_regions[i]
  targ <- meta[run_index,]
  colnames(targ) <- make.names(colnames(targ))
  y <- t(meta[run_index,cell_types]) # C-SIDE cell-types
  #y <- y + 30 # adjust so not negative
  
  targ$DV200 <- as.numeric(targ$DV200)
  targ$Age <- as.numeric(targ$Age)
  targ$area <- as.numeric(targ$area)
  targ$IHC.score <- as.numeric(targ$IHC.score)
  targ$Dx_cont <- recode(targ$Diagnosis,"CTR" = 1,"ePD" = 2, "ILBD" = 3,"lPD" = 4 )
  
  ## Create design matrix
  design <- model.matrix(~ Dx_cont + Age + Sex + DV200 + area ,data=targ)
  colnames(design) <- make.names(colnames(design))
  v <- voom(y,design)
  vfit <- lmFit(v)
  
  # Perform LIMMA contrasts
  cont.matrix <- makeContrasts(A=Dx_cont,levels=design)
  fit2 <- contrasts.fit(vfit, cont.matrix)
  vfit2 <- eBayes(fit2)
  options(digits=3)
  
  # Select significant k-mers for each contrast and bind together
  topA <- topTable(vfit2,coef="A",number=Inf,sort.by="P")
  
  res_list[[i]] <- topA
}

# plot logFC vs p-value
pdf("volcano_limma_csidecell_full_res1.pdf")
EnhancedVolcano(toptable = res_list[[1]], x = "logFC",y = "adj.P.Val",
                lab = rownames(res_list[[1]]),pCutoff = 0.05,
                title = brain_regions[1],FCcutoff = 0,xlim = c(-0.05,0.05), ylim = c(0,10),
                subtitle = "")
dev.off()
pdf("volcano_limma_csidecell_full_res2.pdf")
EnhancedVolcano(toptable = res_list[[2]], x = "logFC",y = "adj.P.Val",
                lab = rownames(res_list[[2]]),pCutoff = 0.05,
                title = brain_regions[2],FCcutoff = 0,xlim = c(-0.05,0.05), ylim = c(0,10),
                subtitle = "")
dev.off()
pdf("volcano_limma_csidecell_full_res3.pdf")
EnhancedVolcano(toptable = res_list[[3]], x = "logFC",y = "adj.P.Val",
                lab = rownames(res_list[[3]]),pCutoff = 0.05,
                title = brain_regions[3],FCcutoff = 0,xlim = c(-0.05,0.05), ylim = c(0,10),
                subtitle = "")

dev.off()