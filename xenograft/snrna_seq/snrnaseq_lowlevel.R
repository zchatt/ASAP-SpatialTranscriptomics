---
  title: 'cellranger and kallisto snRNA-seq integration and normalization'
output:
  html_document:
  theme: united
pdf_document: default
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
---
  
  ```{r setup, include=FALSE}
all_times <- list()  # store the time for each chunk
knitr::knit_hooks$set(time_it = local({
  now <- NULL
  function(before, options) {
    if (before) {
      now <<- Sys.time()
    } else {
      res <- difftime(Sys.time(), now, units = "secs")
      all_times[[options$label]] <<- res
    }
  }
}))
knitr::opts_chunk$set(
  tidy = TRUE,
  tidy.opts = list(width.cutoff = 95),
  fig.width = 10,
  message = FALSE,
  warning = FALSE,
  time_it = TRUE
)
```

```{r data}
library(Seurat)
#library(SeuratData)
library(patchwork)
```

```{r define.paths}
# path of cellranger output
cellranger_out="/Users/zacc/USyd/ASAP/snrna_seq/cellranger_out"

# path of kallisto output
kallisto_out="/Users/zacc/USyd/ASAP/snrna_seq/kallistor_out"
```



```{r init, results='hide', message=FALSE, fig.keep='none'}
# load dataset

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- SCTransform(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)
```