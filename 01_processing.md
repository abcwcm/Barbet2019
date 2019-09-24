scRNA-seq data was generated on the 10X Chromium platform (see the manuscript for details) and processed in `R`.
Both gene expression and VDJ sequencing was performed for two sample types (KO and WT).
After initial processing with the CellRanger pipeline, the data was processed following the [recommendations of the Satija Lab](https://satijalab.org/seurat/v3.1/integration.html).
For adding the information of the VDJ sequencing to the Seurat object, we used the code by Jared Andrews and Atakan Ekiz as shared on [biostars](ttps://www.biostars.org/p/384640/).

```
library(Seurat)
library(magrittr)
library(data.table)
```

## Processing of the scRNA-seq expression data

```
## READ-IN & PROCESSING =====================================================
KO.data <- Read10X(data.dir = "data/KO_Seurat/")
colnames(KO.data) <- paste0("KO_", colnames(KO.data))

WT.data <- Read10X(data.dir = "data/WT_Seurat/")
colnames(WT.data) <- paste0("WT_", colnames(WT.data))

# make a list of Seurat objects
srt_list <- list()
srt_list$WT <- CreateSeuratObject(counts = WT.data,
    project = "KO_explore",
    names.field = 1, 
    min.cells = 3, min.features = 500)
srt_list$WT$genotype <- "WT"

srt_list$KO <- CreateSeuratObject(counts = KO.data,
    project = "KO_explore",
    names.field = 1, 
    min.cells = 3, min.features = 500)
srt_list$KO$genotype <- "KO"

## NORMALIZATION ==========================================
for (i in seq_len(length(srt_list)) ) {
    print(names(srt_list[i]))
    srt_list[[i]] <- SCTransform(srt_list[[i]], verbose = FALSE)
}

## INTEGRATION ============================================
# select features for downstream integration
# run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.
int_features <- SelectIntegrationFeatures(object.list = srt_list, nfeatures = 3000)
# for PrepSCTIntegration, we need to adjust the memory allocation
options(future.globals.maxSize = 4000 * 1024^2)
srt_list <- PrepSCTIntegration(object.list = srt_list,
    anchor.features = int_features, 
    verbose = FALSE)

# now identify anchors and integrate the datasets.
# Commands are identical to the standard workflow,
# but make sure to set normalization.method = 'SCT':
int_anchors <- FindIntegrationAnchors(object.list = srt_list,
    normalization.method = "SCT", 
    anchor.features = int_features, verbose = FALSE)

# Combine the two data sets using the information from the anchors
combined <- IntegrateData(anchorset = int_anchors,
    normalization.method = "SCT", 
    verbose = FALSE)

rm(WT.data); rm(KO.data); rm(srt_list)

## DIMENSIONALITY REDUCTION =========================================
combined <- RunPCA(combined, verbose = FALSE)

png("pca_elbowPlot.png", width = 500)
ElbowPlot(combined)
dev.off()

combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

save(int_anchors, int_features, file = "integration_anchors_and_features.rda")
saveRDS(combined, file = "seurat_integrated_WT-KO.rds")

## MARKER GENE DETECTION =============================================
DefaultAssay(combined) <- "SCT"
markers_clustersRes05 <- FindAllMarkers(combined) 

## cell cycle score -------------------------------------------------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
combined <- CellCycleScoring(combined, s.features = s.genes,
    g2m.features = g2m.genes)

## counting cells per group etc. ------------------------------------
ci <- as.data.table(combined@meta.data)

count_cells <-  ci[, .N, by = c("genotype","cluster.group")]
count_cells[, all_per_sample := sum(N), by = "genotype"]
count_cells[, pct := N/all_per_sample * 100]

count_cells.clsts <-  ci[, .N, by = c("genotype","seurat_clusters")]
count_cells.clsts[, all_per_sample := sum(N), by = "genotype"]
count_cells.clsts[, pct := N/all_per_sample * 100]
```

## Adding the VDJ sequencing info

To add the clonotype information obtained from V(D)J-sequencing to the Seurat object to enable their visualization, we followed recommendations from [Jared Andrews and Atakan Ekiz](https://www.biostars.org/p/384640/).

```
### WT VDJ info
wt_bar <- read.csv("WT-TCR/outs/filtered_contig_annotations.csv")
wt_clone <- read.csv("WT-TCR/outs/clonotypes.csv")
barcodes_clonal_type_col_number <<- which(colnames(wt_bar) == "raw_clonotype_id")
barcodes_barcode_col_number <<- which(colnames(wt_bar) == "barcode")
wt_bar <- subset(unique(wt_bar[,c(barcodes_clonal_type_col_number,barcodes_barcode_col_number)]), raw_clonotype_id != "None")
wt_clone <- merge(wt_clone, wt_bar, by.x = "clonotype_id", by.y="raw_clonotype_id")

### KO VDJ info
tap_bar <- read.csv("KO-TCR/outs/filtered_contig_annotations.csv")
tap_clone <- read.csv("KO-TCR/outs/clonotypes.csv")
barcodes_clonal_type_col_number <<- which(colnames(tap_bar) == "raw_clonotype_id")
barcodes_barcode_col_number <<- which(colnames(tap_bar) == "barcode")
tap_bar <- subset(unique(tap_bar[,c(barcodes_clonal_type_col_number,barcodes_barcode_col_number)]),
    raw_clonotype_id != "None")
tap_clone <- merge(tap_clone, tap_bar, by.x = "clonotype_id", by.y="raw_clonotype_id")


### Create a function to trim unwanted "-1" and append sample information before barcodes
barcoder <- function(df, prefix,  trim="\\-1"){
    
    df$barcode <- gsub(trim, "", df$barcode)
    df$barcode <- paste0(prefix, df$barcode)
    
    df
}

## Override barcode data with the suitable replacement
wt_clone <- barcoder(wt_clone, prefix = "WT_")
tap_clone <- barcoder(tap_clone, prefix = "KO_")

## Combine VDJ data from the separate experiments
all_cln <- rbind(wt_clone,tap_clone)

rownames(all_cln) <- all_cln$barcode
all_cln <- all_cln[,-c(which(colnames(all_cln) == "barcode"))]

## Add back to Seurat object
seurat.obj.with.cln <- AddMetaData(seurat.obj, metadata = all_cln)
saveRDS(seurat.obj.with.cln, "seurat.obj.with.cln.RDS", compress=F)
```


## TRM signature score

Following the description of [Chihara et al.](Ref: https://www.nature.com/articles/s41586-018-0206-z).

```
## defining signature genes -------------------------------------------------
genes_up_in_trm <- c("Prdm1","Cd244","Cdh1","Chn2","Ctla4","Hpgds","Hspa1a",
    "Icos","Inpp4b","Itga1","Itgae","Litaf","Zfp683","Nr4a1","Nr4a2","Qpct",
    "Rgs1","Rgs2","Sik1","Skil","Tmem123","Vps37b","Xcl1","Cd44","Cd69","Il2rb",
    "Il7r","Il2rg","Csf1","Gzmb","Tnfsf10","Notch1","Klf6","Ets1","Zfp36l1",
    "Gpbp1","Nfe2l2","Rela","Fos","Nfat5","Zc3h12a","Cited2","Nr3c1","Rel",
    "Ahr","Nfkb1","Crem","Nr4a3","Nfkb2","Batf","Irf4","Rbpj","Atf3","Fosb",
    "Egr2")

genes_down_in_trm <- c("Cmah","Elovl7","Eomes","Ripor2","Fgf13","Klre1",
    "Ly6c2","Rasgrp2","S1pr1","S1pr5","Sidt1","Slamf6","Tlr1","Usp33")
# these genes that belong to downregulated set were not in the dataset 
genes_down_in_trm[!c(genes_down_in_trm %in% row.names(scale_data))]

# remove the genes that are not in the dataset
genes_down_in_trm <- genes_down_in_trm[c(genes_down_in_trm %in% row.names(scale_data))]

## calculate the signature score ---------------------------------------
signature_res <- lapply(colnames(scale_data), function(cell){
    cell = cell
    df <- data.frame(cell=scale_data[order(scale_data[,cell], 
        decreasing=F),cell])
    df$rank <- seq(1, nrow(df))
    sum_up = sum(df[genes_up_in_trm,]$rank)
    sum_down = sum(df[genes_down_in_trm,]$rank)
    return_df <- data.frame(
        sum_up = sum_up,
        sum_down=sum_down,
        score=sum_up-sum_down)
    return_df$cell <- cell
    return_df
}) %>% rbindlist()

signature_res <- as.data.frame(signature_res)
row.names(signature_res) <- signature_res$cell
signature_res <- signature_res[,-c(which(colnames(signature_res) == "cell"))]
```
