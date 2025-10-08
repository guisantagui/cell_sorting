################################################################################
# Age Sorter: Preprocess the alignments performed with Cell Ranger on DLPFC.   #
# experiment 1 data, and include in the metadata info of cell type annotation, #
# age, neurodegeneration status etc.                                           #
################################################################################
if (!require("BiocManager", quietly = T)){
        install.packages("BiocManager",
                         repos='http://cran.us.r-project.org')
}
if (!require("devtools", quietly = T)){
        install.packages("devtools",
                         repos='http://cran.us.r-project.org')
}
if (!require("ggplot2", quietly = T)){
        install.packages("ggplot2",
                         repos='http://cran.us.r-project.org')
}
if (!require("patchwork", quietly = T)){
        install.packages("patchwork",
                         repos='http://cran.us.r-project.org')
}

if(!require("Seurat", quietly = T)) remotes::install_github("satijalab/seurat",
                                                            "seurat5",
                                                            quiet = TRUE,
                                                            upgrade = "never")
if(!require("DoubletFinder", quietly = T)){
        remotes::install_github("chris-mcginnis-ucsf/DoubletFinder",
                                quiet = T,
                                upgrade = "never")
}
if(!require("plotUtils", quietly = T)){
        remotes::install_github("guisantagui/plotUtils",
                                quiet = T,
                                upgrade = "never")
}
if (!require("ggExtra", quietly = T)){
        install.packages("ggExtra",
                         repos='http://cran.us.r-project.org')
}
if (!require("biomaRt", quietly = T)){
        BiocManager::install("biomaRt", update = F)
}

library(Seurat)
library(SeuratWrappers)
library(plotUtils)
library(patchwork)
library(ggplot2)
library(ggExtra)
library(dplyr)
library(biomaRt)
library(DoubletFinder)
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Parse Cell Ranger outputs of ROSMAP data.")

parser <- add_argument(parser = parser,
                       arg = c("--counts_dir",
                               "--clin_metdat",
                               "--expr_metdat",
                               "--cell_annot",
                               "--cell_cyc_genes",
                               "--outDir"),
                       help = c("Directory where outputs of cell ranger are",
                                "Clinical metadata.",
                                "Experiment metadata.",
                                "Cell annotation",
                                "Cell cycle genes",
                                "Output directory where Seurat objects will be stored"),
                       flag = c(F, F, F, F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
counts_dir <- parsed$counts_dir
clin_mdat_f <- parsed$clin_metdat
expr_mdat_f <- parsed$expr_metdat
cell_annot_f <- parsed$cell_annot
cellCycGenesFile <- parsed$cell_cyc_genes

#seur_out_dir <- sprintf("%s/seur_objs/", dirname(counts_dir))
seur_out_dir <- parsed$outDir



#clin_mdat_f <- "/work/projects/age_sorter/data/counts/ROSMAP_sc/metadata/ROSMAP_clinical.csv"
#cell_annot_f <- "/work/projects/age_sorter/data/counts/ROSMAP_sc/metadata/cell-annotation.full-atlas.csv"
#expr_mdat_f <- "/work/projects/age_sorter/data/counts/ROSMAP_sc/metadata/ROSMAP_assay_scrnaSeq_metadata.csv"
cellCycGenesFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/utility_files/CellCycleGenes_Human.csv"

expr_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/cell_ranger_outs/DLPFC_exp1/ROSMAP_assay_scrnaSeq_metadata.csv"
clin_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/cell_ranger_outs/ROSMAP_clinical.csv"
samp_metdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/cell_ranger_outs/DLPFC_exp1/outs/SYNAPSE_METADATA_MANIFEST.tsv"
cell_metdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/cell_ranger_outs/DLPFC_exp1/ROSMAP_Brain.snRNAseq_metadata_cells_20230420.csv"
barcodes_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/cell_ranger_outs/DLPFC_exp1/outs/barcodes/" #Cellranger was ran with all samples at once, so no info of samp ID in metadata. These barcodes were obtained from the FASTQ files.
counts_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/cell_ranger_outs/DLPFC_exp1/outs/filtered_feature_bc_matrix/"
seur_out_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/seur_objs/DLPFCexp1/"


create_dir_if_not(seur_out_dir)


plot_dir <- sprintf("%s/plots_preproc/", seur_out_dir)
create_dir_if_not(plot_dir)

# Functions
################################################################################
filterCells <- function(seur, mad.coeff = 3,pass = -1, org = "HUMAN", plot.path = "./"){
        seur@meta.data$Barcodes <- rownames(seur@meta.data)
        # Remove cells with zero counts
        seur <- seur[, WhichCells(seur, expression = nCount_RNA != 0)]
        
        #Calculate percent.mito and percent.ribo metadata columns if they are not there
        if(!any(colnames(seur@meta.data) == "percent.mito")){
                if(org == "HUMAN"){
                        seur[["percent.mito"]] <- PercentageFeatureSet(seur,
                                                                       pattern = "^MT-")
                } else if(org == "MOUSE"){
                        seur[["percent.mito"]] <- PercentageFeatureSet(seur,
                                                                       pattern = "^mt-")
                } else {
                        stop("The specified organism is not supported")
                }
        }
        if(!any(colnames(seur@meta.data) == "percent.ribo")){
                if(org == "HUMAN"){
                        seur[["percent.ribo"]] <- PercentageFeatureSet(seur,
                                                                       pattern = "(^RPL|^RPS|^MRP)")
                } else if(org == "MOUSE"){
                        seur[["percent.ribo"]] <- PercentageFeatureSet(seur,
                                                                       pattern = "(^Rpl|^Rps|^Mrp)")
                } else {
                        stop("The specified organism is not supported")
                }
        }
        # Filtering cells based on percentage of mitochondrial transcripts
        cell.QC.stat <- seur@meta.data
        
        max.mito.thr <- median(cell.QC.stat$percent.mito) + mad.coeff*mad(cell.QC.stat$percent.mito)
        min.mito.thr <- median(cell.QC.stat$percent.mito) - mad.coeff*mad(cell.QC.stat$percent.mito)
        
        p1 <- ggplot(cell.QC.stat, aes(x=nFeature_RNA, y=percent.mito)) +
                geom_point() +
                geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
                geom_hline(aes(yintercept = min.mito.thr), colour = "red", linetype = 2) +
                annotate(geom = "text",
                         label = paste0(as.numeric(table(cell.QC.stat$percent.mito > max.mito.thr | cell.QC.stat$percent.mito < min.mito.thr)[2]),
                                        " cells removed\n",
                                        as.numeric(table(cell.QC.stat$percent.mito > max.mito.thr | cell.QC.stat$percent.mito < min.mito.thr)[1]),
                                        " cells remain"),
                         x = 6000,
                         y = -10)
        
        p <- ggMarginal(p1, type = "histogram", fill="lightgrey", bins=100) 
        ggsave(paste0(plot.path,"Mitofilter_Marginal_Pass",pass,".png"),plot = p)
        
        cell.QC.stat <- cell.QC.stat %>%
                dplyr::filter(percent.mito <= max.mito.thr) %>% 
                dplyr::filter(percent.mito >= min.mito.thr)
        
        # Filtering cells based on number of genes and transcripts detected
        # Set low and hight thresholds on the number of detected genes
        min.features.thr <- median(log10(cell.QC.stat$nFeature_RNA)) - mad.coeff*mad(log10(cell.QC.stat$nFeature_RNA))
        max.features.thr <- median(log10(cell.QC.stat$nFeature_RNA)) + mad.coeff*mad(log10(cell.QC.stat$nFeature_RNA))
        
        # Set hight threshold on the number of transcripts
        max.count.thr <- median(log10(cell.QC.stat$nCount_RNA)) + mad.coeff*mad(log10(cell.QC.stat$nCount_RNA))
        
        p1 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
                geom_point() +
                geom_smooth(method="lm") +
                geom_hline(aes(yintercept = min.features.thr), colour = "green", linetype = 2) +
                geom_hline(aes(yintercept = max.features.thr), colour = "green", linetype = 2) +
                geom_vline(aes(xintercept = max.count.thr), colour = "red", linetype = 2)
        
        p <- ggMarginal(p1, type = "histogram", fill="lightgrey")
        ggsave(paste0(plot.path,"FeatureAndCountFilter_1_Pass",pass,".png"),plot = p)
        
        # Filter cells base on both metrics
        cell.QC.stat <- cell.QC.stat %>% 
                dplyr::filter(log10(nFeature_RNA) > min.features.thr) %>%
                dplyr::filter(log10(nFeature_RNA) < max.features.thr) %>%
                dplyr::filter(log10(nCount_RNA) < max.count.thr)
        
        lm.model <- lm(data = cell.QC.stat, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
        
        p2 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
                geom_point() +
                geom_smooth(method="lm") +
                geom_hline(aes(yintercept = min.features.thr), colour = "green", linetype = 2) +
                geom_hline(aes(yintercept = max.features.thr), colour = "green", linetype = 2) +
                geom_vline(aes(xintercept = max.count.thr), colour = "red", linetype = 2) +
                geom_abline(intercept = lm.model$coefficients[1] - 0.09 ,
                            slope = lm.model$coefficients[2],
                            color="orange") +
                annotate(geom = "text", label = paste0(dim(cell.QC.stat)[1],
                                                       " QC passed cells"), x = 4, y = 3.8)
        
        p <- ggMarginal(p2, type = "histogram", fill="lightgrey")
        ggsave(paste0(plot.path,
                      "FeatureAndCountOutlier_2_Pass",
                      pass,
                      ".png"),
               plot = p)
        
        # Cells to exclude lie below an intercept offset of -0.09
        cell.QC.stat$validCells <- log10(cell.QC.stat$nFeature_RNA) > (log10(cell.QC.stat$nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] - 0.09))
        
        p3 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
                geom_point(aes(colour = validCells), fill = "black",pch=21) +
                scale_color_manual(breaks = c("TRUE", "FALSE"), 
                                   values = c("black","firebrick1")) +
                geom_smooth(method="lm") +
                geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") + 
                theme(legend.position="none") +
                annotate(geom = "text", label = paste0(as.numeric(table(cell.QC.stat$validCells)[2]), " QC passed cells\n",
                                                       as.numeric(table(cell.QC.stat$validCells)[1]), " QC filtered"), x = 4, y = 3.8)
        
        p <- ggMarginal(p3, type = "histogram", fill="lightgrey")
        ggsave(paste0(plot.path,"FeatureAndCountOutlier_3_Pass",pass,".png"),plot = p)
        
        # Remove invalid cells
        cell.QC.stat <- cell.QC.stat %>% dplyr::filter(validCells)
        
        seur <- subset(seur, subset = Barcodes %in% cell.QC.stat$Barcodes)
        return(seur)
}

# Find the optimal number of PCs
findNumPCs <- function(seur){
        pct <- seur[["pca"]]@stdev / sum(seur[["pca"]]@stdev) * 100
        cumu <- cumsum(pct)
        
        co1 <- which(cumu > 90 & pct < 5)[1]
        
        co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05),
                    decreasing = T)[1] + 1
        
        return(min(co1,co2))
}

# Load data 
################################################################################
samp_metdat <- read_table_fast(samp_metdat_f)
cell_metdat <- read.csv(cell_metdat_f)

clin_mdat <- read.csv(clin_mdat_f, row.names = 1)
expr_mdat <- read.csv(expr_mdat_f)

# Remove samples with ambiguous ages from metadata
clin_mdat <- clin_mdat[clin_mdat$age_death != "90+", ]
clin_mdat <- clin_mdat[clin_mdat$age_death != "", ]
clin_mdat$age_death <- as.numeric(clin_mdat$age_death)
clin_mdat <- clin_mdat[!is.na(clin_mdat$cogdx), ]

cellcyclegenes <- read.csv(cellCycGenesFile)

mart <- useDataset("hsapiens_gene_ensembl",
                   useMart("ensembl"))
ensToSymb <- getBM(filters= "ensembl_gene_id",
                   attributes= c("ensembl_gene_id","hgnc_symbol"),
                   values=cellcyclegenes$geneID,mart= mart)
cellcyclegenes <- merge(cellcyclegenes,ensToSymb,by.x = 2, by.y = 1, all.x = T)

# Load the sc counts into a seurat object
################################################################################
count_files <- list.files(counts_dir)

# Cell ranger was ran with all the samples at the same time, so because of this
# it doesnt store sample info into the barcodes file it generates. So we can
# retrieve the barcodes from the fastqs per fastq file (R1) and store them in
# txt files. We then parse the files into a dataframe of barcodes per sample,
# so we're able to assign a sample to each cell in the seurat object.
barcode_files <- list.files(barcodes_dir, full.names = T)

barcodes_per_samp <- data.frame(matrix(nrow = 0, ncol = 2,
                                       dimnames = list(NULL,
                                                       c("barcode", "sample"))))
for (i in seq_along(barcode_files)){
        samp <- basename(barcode_files[i])
        samp <- gsub("\\..*", "", samp)
        samp <- gsub("\\_.*", "", samp)
        bc <- read_table_fast(barcode_files[i])
        bc <- data.frame(barcode = paste0(bc$V2, "-1"),
                         sample = rep(samp, length(bc$V2)))
        barcodes_per_samp <- rbind.data.frame(barcodes_per_samp, bc)
}

# Create a seurat object
################################################################################
expression_matrix <- Read10X(data.dir = counts_dir)
# Create the seurat object
seurat_object = CreateSeuratObject(counts = expression_matrix)
seurat_object@meta.data$specimenID = barcodes_per_samp$sample[match(rownames(seurat_object@meta.data),
                                                                    barcodes_per_samp$barcode)]
table(seurat_object@meta.data$specimenID)
# Remove samples that have very low cell count
low_cell_samps <- names(table(seurat_object@meta.data$specimenID))[table(seurat_object@meta.data$specimenID) < 100]

seurat_object <- seurat_object[, !seurat_object@meta.data$specimenID %in% low_cell_samps]

seurat_object@meta.data$individualID = samp_metdat$individualID[match(seurat_object@meta.data$specimenID,
                                                                      samp_metdat$specimenID)]
table(seurat_object@meta.data$specimenID)
# Remove the cells that come from individuals that are not in clinical metadata,
# since we removed from there those without valid ages.
seurat_object <- seurat_object[, seurat_object@meta.data$individualID %in% clin_mdat$individualID]
seurat_object@meta.data$age_death <- clin_mdat$age_death[match(seurat_object@meta.data$individualID,
                                                               clin_mdat$individualID)]
seurat_object@meta.data$braaksc <- clin_mdat$braaksc[match(seurat_object@meta.data$individualID,
                                                           clin_mdat$individualID)]
seurat_object@meta.data$ceradsc <- clin_mdat$ceradsc[match(seurat_object@meta.data$individualID,
                                                           clin_mdat$individualID)]
seurat_object@meta.data$cogdx <- clin_mdat$cogdx[match(seurat_object@meta.data$individualID,
                                                       clin_mdat$individualID)]
seurat_object@meta.data$libraryBatch <- expr_mdat$libraryBatch[match(seurat_object$specimenID,
                                                                     expr_mdat$specimenID)]
seurat_object@meta.data$sequencingBatch <- expr_mdat$sequencingBatch[match(seurat_object$specimenID,
                                                                           expr_mdat$specimenID)]
head(seurat_object@meta.data)
# Preprocess each sample (specimenID) separately
uniq_samps <- unique(seurat_object@meta.data$specimenID)
for (i in seq_along(uniq_samps)){
        i <- 1
        samp <- uniq_samps[i]
        seur_samp <- seurat_object[, seurat_object@meta.data$specimenID == samp]
        minCells <- 3
        minFeats <- 700
        nfeatures <- Matrix::colSums(x = seur_samp@assays$RNA$counts > 0)
        nCellsPassed <- length(which(x = nfeatures >= minFeats))
        cellsPassed <- names(nfeatures)[nfeatures >= minFeats]
        seur_samp <- seur_samp[, cellsPassed]
        plot_samp_dir <- sprintf("%s%s/", plot_dir, samp)
        create_dir_if_not(plot_samp_dir)
        seur_samp <- filterCells(seur_samp,
                                 mad.coeff = 3,
                                 pass = -1,
                                 org = "HUMAN",
                                 plot.path = plot_samp_dir)
        # Log-normalize data
        seur_samp <- NormalizeData(seur_samp)
        # Find variable features --> Identify features that are outliers on a mean
        # variability plot
        seur_samp <- FindVariableFeatures(seur_samp,
                                          selection.method = "vst",
                                          nfeatures = 3000)
        # Scale data --> Scales and centers features in the dataset
        seur_samp <- ScaleData(seur_samp,
                               features = VariableFeatures(object = seur_samp))
        # Run PCA
        if(ncol(seur_samp) < 100){
                nPCs <- ncol(seur_samp) - 1
        }else{
                nPCs <- 100
        }
        seur_samp <- RunPCA(seur_samp,
                            features = VariableFeatures(object = seur_samp),
                            npcs = nPCs)
        
        g2m_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "G2/M"]
        g2m_genes <- intersect(g2m_genes, rownames(seur_samp))
        s_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "S"]
        s_genes <- intersect(s_genes, rownames(seur_samp))
        
        seur_samp <- CellCycleScoring(seur_samp,
                                      s.features = s_genes,
                                      g2m.features = g2m_genes,
                                      set.ident = F)
        
        ggsave(paste0(plot_samp_dir,
                      "PCA_CellCycleGenes.png"),
               DimPlot(seur_samp,
                       group.by = "Phase"))
        
        seur_samp$CC.Difference <- seur_samp$S.Score - seur_samp$G2M.Score
        
        seur_samp <- SCTransform(seur_samp,
                                 vars.to.regress = "CC.Difference",
                                 vst.flavor = "v2",
                                 verbose = T)
        
        seur_samp <- RunPCA(seur_samp,
                            assay = "SCT",
                            npcs = nPCs,
                            verbose = T)
        
        numPCs <- findNumPCs(seur_samp) #+
        #round(findNumPCs(tissSampFilt)/2)
        
        elbPlot <- ElbowPlot(seur_samp,
                             ndims = nPCs) +
                geom_vline(aes(xintercept = numPCs),
                           colour = "red",
                           linetype = "dashed")
        
        ggsave(paste0(plot_samp_dir,
                      "ElbowPlot.png"),
               elbPlot)
        
        # paramSweep_v3 finds pN and pK values for removing the doublets
        # later on
        sweep.list <- paramSweep(seur_samp,
                                 PCs = 1:numPCs,
                                 sct = T)
        sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        
        ## Estimate expected percentage of doublets from 10X Genomics 
        # estimates from 3' Gene Expression v3.1 assay##
        estDoublets <- c(0.4,0.8,1.6,2.4,3.2,4,4.8,5.6,6.4,7.2,8)
        numCellsRec <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
        
        scatter.smooth(numCellsRec, estDoublets) #Looks linear
        lm_doublets <- lm(estDoublets ~ numCellsRec)
        summary(lm_doublets) #Perfect linear relationship r2 = 1
        
        nExp <- round(ncol(seur_samp) * (unname(predict(lm_doublets,
                                                        data.frame(numCellsRec = ncol(seur_samp))))/100))
        
        pK <- as.numeric(levels(bcmvn$pK)[bcmvn$BCmetric == max(bcmvn$BCmetric)])
        seur_samp <- doubletFinder(seu = seur_samp,
                                   PCs = 1:numPCs,
                                   pN = 0.25, #default
                                   pK = pK,
                                   nExp = nExp,
                                   reuse.pANN = NULL,
                                   sct = TRUE,
                                   annotations = NULL)
        
        seur_samp$Doublets <- "Singlet"
        
        seur_samp$Doublets[seurat_object[[paste0("pANN_0.25_",
                                                 pK,
                                                 "_",
                                                 nExp)]] >= 0.5] <- "Doublet"
        
        seur_samp$Doublets <- factor(seur_samp$Doublets,
                                     levels = c("Doublet",
                                                "Singlet"))
        
        p <- ggplot(seur_samp@meta.data,
                    aes(x=log10(nCount_RNA),
                        y=log10(nFeature_RNA))) +
                geom_point(aes(colour = Doublets), fill = "black",pch=21) +
                scale_color_manual(breaks = c("Singlet", "Doublet"), 
                                   values = c("black","firebrick1")) +
                geom_smooth(method="lm") +
                theme(legend.position="none") +
                annotate(geom = "text",
                         label = paste0(as.numeric(table(seur_samp@meta.data$Doublets)[2]),
                                        " Singlets\n",
                                        as.numeric(table(seur_samp@meta.data$Doublets)[1]),
                                        " Doublets"),
                         x = 4,
                         y = 3.8)
        
        ggsave(paste0(plot_samp_dir,"Doublets.png"),plot = p)
        
        seur_samp <- subset(seur_samp,
                            subset = Doublets == "Singlet")
        seur_samp <- RunUMAP(seur_samp,
                             dims = 1:numPCs,
                             n.neighbors = 20)
        
        cds <- as.cell_data_set(seur_samp)
        cds <- cluster_cells(cds,
                             k = 10,
                             reduction_method = "UMAP",
                             cluster_method = "leiden",
                             #resolution = 0.0001,
                             num_iter = 5
        )
        
        seur_samp$seurat_clusters <- clusters(cds)
        Idents(seur_samp) <- clusters(cds)
        
        umap_seur_clusts <- DimPlot(seur_samp,
                                    reduction = "umap",
                                    group.by = "seurat_clusters")
        ggsave(filename = sprintf("%sumap_seur_clusts.pdf",
                                  plot_samp_dir),
               plot = umap_seur_clusts)
}