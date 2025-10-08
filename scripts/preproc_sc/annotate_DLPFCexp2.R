################################################################################
# Age Sorter: Preprocess the alignments performed with Cell Ranger on DLPFC.   #
# experiment 2 data, and include in the metadata info of cell type annotation, #
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
                               "--cell_cyc_genes"),
                       help = c("Directory where outputs of cell ranger are",
                                "Clinical metadata.",
                                "Experiment metadata.",
                                "Cell annotation",
                                "Cell cycle genes"),
                       flag = c(F, F, F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
clin_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/ROSMAP_sc/metadata/ROSMAP_clinical.csv"
cell_annot_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/ROSMAP_sc/metadata/cell-annotation.full-atlas.csv"
#cell_annot_n424_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/ROSMAP_sc/metadata/cell-annotation.n424.csv"
#cell_annot_n427_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/ROSMAP_sc/metadata/cell-annotation.n437.csv"
expr_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/ROSMAP_sc/metadata/ROSMAP_assay_scrnaSeq_metadata.csv"
cellCycGenesFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/utility_files/CellCycleGenes_Human.csv"

counts_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/ROSMAP_sc/counts/"

counts_dir <- parsed$counts_dir
clin_mdat_f <- parsed$clin_metdat
expr_mdat_f <- parsed$expr_metdat
cell_annot_f <- parsed$cell_annot
cellCycGenesFile <- parsed$cell_cyc_genes

seur_out_dir <- sprintf("%s/seur_objs/", dirname(counts_dir))
create_dir_if_not(seur_out_dir)


plot_dir <- sprintf("%s/plots_preproc/", dirname(counts_dir))
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
                geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") +
                annotate(geom = "text", label = paste0(dim(cell.QC.stat)[1], " QC passed cells"), x = 4, y = 3.8)
        
        p <- ggMarginal(p2, type = "histogram", fill="lightgrey")
        ggsave(paste0(plot.path,"FeatureAndCountOutlier_2_Pass",pass,".png"),plot = p)
        
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
clin_mdat <- read.csv(clin_mdat_f, row.names = 1)
expr_mdat <- read.csv(expr_mdat_f)

# Remove samples with ambiguous ages from metadata
clin_mdat <- clin_mdat[clin_mdat$age_death != "90+", ]
clin_mdat <- clin_mdat[clin_mdat$age_death != "", ]
clin_mdat$age_death <- as.numeric(clin_mdat$age_death)
clin_mdat <- clin_mdat[!is.na(clin_mdat$cogdx), ]

cell_annot <- read_table_fast(cell_annot_f)
cell_annot <- cell_annot[, 2:ncol(cell_annot)]
colnames(cell_annot) <- cell_annot[1, ]
cell_annot <- cell_annot[2:nrow(cell_annot), ]
cell_annot

cell_annot_n424 <- read_table_fast(cell_annot_n424_f)
cell_annot_n427 <- read_table_fast(cell_annot_n427_f)

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
samples <- unique(gsub("\\..*", "", count_files))

# Create a temp dir for loading the matrixes into seurat, as names need to be
# changed
temp_dir <- sprintf("%s/temp/", dirname(counts_dir))
create_dir_if_not(temp_dir)

temp_dir_files <- list.files(temp_dir, full.names = T)

barcodes <- read_table_fast(temp_dir_files[1])

# Create a seurat object
for (s in samples){
        s <- samples[1]
        # Copy sample files into temp dir
        system(sprintf("cp %s%s* %s", counts_dir, s, temp_dir))
        # Change the names to remove reference to sample
        system(sprintf("mv %s%s.barcodes.tsv.gz %sbarcodes.tsv.gz", temp_dir, s, temp_dir))
        system(sprintf("mv %s%s.features.tsv.gz %sfeatures.tsv.gz", temp_dir, s, temp_dir))
        system(sprintf("mv %s%s.matrix.mtx.gz %smatrix.mtx.gz", temp_dir, s, temp_dir))
        # Load the expression matrix
        expression_matrix <- Read10X(data.dir = temp_dir)
        # Create the seurat object
        seurat_object = CreateSeuratObject(counts = expression_matrix)
        cell_annot_samp <- cell_annot[grep(s, cell_annot$cell), ]
        cell_annot_samp$cell <- gsub(paste0(s, "_"), "", cell_annot_samp$cell)
        
        cell_annot_n424_samp <- cell_annot_n424[grep(s, cell_annot_n424$barcode), ]
        cell_annot_n427_samp <- cell_annot_n427[grep(s, cell_annot_n427$cell), ]
        
        colnames(cell_annot_n424)
        
        seurat_object@meta.data$cell_type <- cell_annot_samp$cell.type[match(rownames(seurat_object@meta.data),
                                                                             cell_annot_samp$cell)]
        seurat_object@meta.data$individualID <- cell_annot_samp$individualID[match(rownames(seurat_object@meta.data),
                                                                                   cell_annot_samp$cell)]
        seurat_object@meta.data$age_death <- clin_mdat$age_death[match(seurat_object@meta.data$individualID,
                                                                       clin_mdat$individualID)]
        seurat_object@meta.data$braaksc <- clin_mdat$braaksc[match(seurat_object@meta.data$individualID,
                                                                   clin_mdat$individualID)]
        seurat_object@meta.data$ceradsc <- clin_mdat$ceradsc[match(seurat_object@meta.data$individualID,
                                                                   clin_mdat$individualID)]
        seurat_object@meta.data$cogdx <- clin_mdat$cogdx[match(seurat_object@meta.data$individualID,
                                                               clin_mdat$individualID)]
        cells_wAllInfo <- rownames(seurat_object@meta.data[!is.na(seurat_object@meta.data$age_death), ])
        seurat_object <- seurat_object[, cells_wAllInfo]
        
        # Remove cells that have less than 700 features
        minCells <- 3
        minFeats <- 700
        colnames(seurat_object@assays$RNA)
        nfeatures <- Matrix::colSums(x = seurat_object@assays$RNA$counts > 0)
        nCellsPassed <- length(which(x = nfeatures >= minFeats))
        cellsPassed <- names(nfeatures)[nfeatures >= minFeats]
        seurat_object <- seurat_object[, cellsPassed]
        
        plot_samp_dir <- sprintf("%s%s/", plot_dir, s)
        create_dir_if_not(plot_samp_dir)
        seurat_object <- filterCells(seurat_object,
                                     mad.coeff = 3,
                                     pass = -1,
                                     org = "HUMAN",
                                     plot.path = plot_samp_dir)
        # Log-normalize data
        seurat_object <- NormalizeData(seurat_object)
        # Find variable features --> Identify features that are outliers on a mean
        # variability plot
        seurat_object <- FindVariableFeatures(seurat_object,
                                              selection.method = "vst",
                                              nfeatures = 3000)
        # Scale data --> Scales and centers features in the dataset
        seurat_object <- ScaleData(seurat_object,
                                   features = VariableFeatures(object = seurat_object))
        
        # Run PCA
        if(ncol(seurat_object) < 100){
                nPCs <- ncol(seurat_object) - 1
        }else{
                nPCs <- 100
        }
        seurat_object <- RunPCA(seurat_object,
                                features = VariableFeatures(object = seurat_object),
                                npcs = nPCs)
        g2m_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "G2/M"]
        g2m_genes <- intersect(g2m_genes,rownames(seurat_object))
        s_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "S"]
        s_genes <- intersect(s_genes,rownames(seurat_object))
        
        seurat_object <- CellCycleScoring(seurat_object,
                                          s.features = s_genes,
                                          g2m.features = g2m_genes,
                                          set.ident = F)
        
        ggsave(paste0(plot_samp_dir,
                      "PCA_CellCycleGenes.png"),
               DimPlot(seurat_object,
                       group.by = "Phase"))
        
        
        seurat_object$CC.Difference <- seurat_object$S.Score - seurat_object$G2M.Score
        
        seurat_object <- SCTransform(seurat_object,
                                     vars.to.regress = "CC.Difference",
                                     vst.flavor = "v2",
                                     verbose = T)
        
        seurat_object <- RunPCA(seurat_object,
                                assay = "SCT",
                                npcs = nPCs,
                                verbose = T)
        
        numPCs <- findNumPCs(seurat_object) #+
        #round(findNumPCs(tissSampFilt)/2)
        
        elbPlot <- ElbowPlot(seurat_object,
                             ndims = nPCs) +
                geom_vline(aes(xintercept = numPCs),
                           colour = "red",
                           linetype = "dashed")
        
        ggsave(paste0(plot_samp_dir,
                      "ElbowPlot.png"),
               elbPlot)
        
        # paramSweep_v3 finds pN and pK values for removing the doublets
        # later on
        sweep.list <- paramSweep_v3(seurat_object,
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
        
        nExp <- round(ncol(seurat_object) * (unname(predict(lm_doublets,
                                                            data.frame(numCellsRec = ncol(seurat_object))))/100))
        
        pK <- as.numeric(levels(bcmvn$pK)[bcmvn$BCmetric == max(bcmvn$BCmetric)])
        seurat_object <- doubletFinder_v3(seu = seurat_object,
                                          PCs = 1:numPCs,
                                          pN = 0.25, #default
                                          pK = pK,
                                          nExp = nExp,
                                          reuse.pANN = FALSE,
                                          sct = TRUE,
                                          annotations = NULL)
        
        seurat_object$Doublets <- "Singlet"
        
        seurat_object$Doublets[seurat_object[[paste0("pANN_0.25_",
                                                     pK,
                                                     "_",
                                                     nExp)]] >= 0.5] <- "Doublet"
        
        seurat_object$Doublets <- factor(seurat_object$Doublets,
                                         levels = c("Doublet",
                                                    "Singlet"))
        
        p <- ggplot(seurat_object@meta.data,
                    aes(x=log10(nCount_RNA),
                        y=log10(nFeature_RNA))) +
                geom_point(aes(colour = Doublets), fill = "black",pch=21) +
                scale_color_manual(breaks = c("Singlet", "Doublet"), 
                                   values = c("black","firebrick1")) +
                geom_smooth(method="lm") +
                theme(legend.position="none") +
                annotate(geom = "text",
                         label = paste0(as.numeric(table(seurat_object@meta.data$Doublets)[2]),
                                        " Singlets\n",
                                        as.numeric(table(seurat_object@meta.data$Doublets)[1]),
                                        " Doublets"),
                         x = 4,
                         y = 3.8)
        
        ggsave(paste0(plot_samp_dir,"Doublets.png"),plot = p)
        
        seurat_object <- subset(seurat_object,
                                subset = Doublets == "Singlet")
        seurat_object <- RunUMAP(seurat_object,
                                 dims = 1:numPCs,
                                 n.neighbors = 20)
        #Let's see how it looks like
        umap_cell_types <- DimPlot(seurat_object,
                                   reduction = "umap",
                                   group.by = "cell_type")
        ggsave(filename = sprintf("%sumap_celltype.pdf", plot_samp_dir),
               plot = umap_cell_types)
        saveRDS(seurat_object, file = sprintf("%s%s.rds", seur_out_dir, s))
        print(sprintf("%s.rds saved at %s.", s, seur_out_dir))
}
