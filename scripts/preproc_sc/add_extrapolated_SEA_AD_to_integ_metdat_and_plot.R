library(plotUtils)

SEA_AD_mDat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/SEA_AD/metadata/SEA-AD_individual_metadata.xlsx"
integ_mDat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/results/integ_metadata/metdat_DLPFCexp2_ageAnno_GSE25469.csv"
STAB2_mDat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/STAB2/metadata/human_cell_meta.csv"
STAB2_srcs_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/STAB2/metadata/human_source.xlsx"
STAB2_autism_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/STAB2/metadata/autism/meta.tsv"
STAB2_morabito_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/STAB2/metadata/morabito/snRNA_metadta.csv"
STAB2_franjic_cellmeta_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/STAB2/metadata/franjic/GSE186538_Human_cell_meta.txt"
STAB2_franjic_indivmeta_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/STAB2/metadata/franjic/NIHMS1758771-supplement-Table_S1.xlsx"
STAB2_lake_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/STAB2/metadata/lake/41587_2018_BFnbt4038_MOESM11_ESM.xlsx"
STAB2_jakel_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/STAB2/metadata/jakel/GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt"
STAB2_jakel_indiv_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/STAB2/metadata/jakel/41586_2019_903_MOESM3_ESM.xlsx"
STAB2_absinta_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/STAB2/metadata/absinta/GSE180759_annotation.txt"
STAB2_absinta_indiv_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/STAB2/metadata/absinta/41586_2021_3892_MOESM4_ESM.xlsx"
STAB2_herring_indiv_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/STAB2/metadata/herring/ScienceDirect_files_31Mar2025_13-55-36.176/1-s2.0-S0092867422012582-mmc1.xlsx"
STAB2_zhong_cells_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/STAB2/metadata/zhong/GSE119212_hippocampus_aggr_8_barcodes.tsv.gz"

GSM6657986_cell_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/GSM6657986/metadata/GSM6657986_cell_annotation.csv.gz"
GSM6657986_indiv_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/GSM6657986/metadata/41588_2023_1572_MOESM3_ESM.xlsx"

GSE193688_indiv_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/GSE193688/metadata/43587_2024_583_MOESM3_ESM.xlsx"

out_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/results/integ_metadata/"

integ_mDat <- read_table_fast(integ_mDat_f, row.names = 1)
SEA_AD_mDat <- as.data.frame(readxl::read_xlsx(SEA_AD_mDat_f))
STAB2_mDat <- read_table_fast(STAB2_mDat_f, row.names = 1)
STAB2_srcs <- as.data.frame(readxl::read_xlsx(STAB2_srcs_f))
colnames(STAB2_srcs) <- STAB2_srcs[2, ]
STAB2_srcs <- STAB2_srcs[3:nrow(STAB2_srcs), ]

STAB2_autism <- read_table_fast(STAB2_autism_f, row.names = 1)
# Keep only controls
STAB2_autism <- STAB2_autism[STAB2_autism$diagnosis == "Control", ]
# Keep only adults (we dont want to pick up development)
STAB2_autism <- STAB2_autism[STAB2_autism$age >= 18, ]

# Keep only 10X
STAB2_mDat <- STAB2_mDat[STAB2_mDat$platform == "10X", ]


STAB2_morabito <- read_table_fast(STAB2_morabito_f, row.names = 1)
table(STAB2_morabito$Diagnosis)

STAB2_franjic_cellmeta <- read_table_fast(STAB2_franjic_cellmeta_f, row.names = 1)
STAB2_franjic_indivmeta <- data.frame(readxl::read_xlsx(STAB2_franjic_indivmeta_f))
colnames(STAB2_franjic_indivmeta) <- STAB2_franjic_indivmeta[1, ]
STAB2_franjic_indivmeta <- STAB2_franjic_indivmeta[2:nrow(STAB2_franjic_indivmeta), ]
# Keep only human
STAB2_franjic_indivmeta <- STAB2_franjic_indivmeta[STAB2_franjic_indivmeta$Species == "Human", ]
STAB2_franjic_cellmeta <- STAB2_franjic_cellmeta[STAB2_franjic_cellmeta$samplename %in% STAB2_franjic_indivmeta$Case, ]
STAB2_franjic_cellmeta$age_death <- STAB2_franjic_indivmeta$`Age (Y)`[match(STAB2_franjic_cellmeta$samplename,
                                                                            STAB2_franjic_indivmeta$Case)]
STAB2_franjic_cellmeta$age_death <- as.numeric(STAB2_franjic_cellmeta$age_death)


STAB2_lake_indiv <- data.frame(readxl::read_xlsx(STAB2_lake_f, sheet = 1))
colnames(STAB2_lake_indiv) <- STAB2_lake_indiv[2, ]
STAB2_lake_indiv <- STAB2_lake_indiv[3:nrow(STAB2_lake_indiv), 2:ncol(STAB2_lake_indiv)]
STAB2_lake_indiv <- STAB2_lake_indiv[1:(which(is.na(STAB2_lake_indiv$AGE)) - 1), ]
STAB2_lake_indiv$`Patient UMB#`
STAB2_lake_indiv$AGE <- as.numeric(gsub(" years", "", STAB2_lake_indiv$AGE))

STAB2_lake_cell <- data.frame(readxl::read_xlsx(STAB2_lake_f, sheet = 2))
View(STAB2_lake_cell)

colnames(STAB2_lake_cell) <- STAB2_lake_cell[4, ]
STAB2_lake_cell <- STAB2_lake_cell[5:nrow(STAB2_lake_cell), 1:6]
colnames(STAB2_lake_cell) <- gsub("Patient UMB#", "individualID", colnames(STAB2_lake_cell))
STAB2_lake_cell$age_death <- STAB2_lake_indiv$AGE[match(STAB2_lake_cell$individualID,
                                                        STAB2_lake_indiv$`Patient UMB#`)]
STAB2_jakel <- read_table_fast(STAB2_jakel_f, row.names = 1)
STAB2_jakel <- STAB2_jakel[STAB2_jakel$Condition == "Ctrl", ]
STAB2_jakel_indiv <- data.frame(readxl::read_xlsx(STAB2_jakel_indiv_f))
colnames(STAB2_jakel_indiv) <- STAB2_jakel_indiv[3, ]
STAB2_jakel_indiv <- STAB2_jakel_indiv[5:13, ]
STAB2_jakel$age_death <- STAB2_jakel_indiv$`Age at death`[match(STAB2_jakel$Sample,
                                                                STAB2_jakel_indiv$ID)]
STAB2_jakel$age_death <- as.numeric(STAB2_jakel$age_death)

STAB2_absinta <- read_table_fast(STAB2_absinta_f, row.names = 1)
STAB2_absinta <- STAB2_absinta[STAB2_absinta$pathology == "control_white_matter", ]
STAB2_absinta_indiv <- data.frame(readxl::read_xlsx(STAB2_absinta_indiv_f, sheet = 2))
View(STAB2_absinta_indiv)
colnames(STAB2_absinta_indiv) <- STAB2_absinta_indiv[2, ]
STAB2_absinta_indiv <- STAB2_absinta_indiv[21:26, ]

STAB2_absinta_indiv_parsed <- data.frame(individualID = STAB2_absinta_indiv[1:6 %% 2 == 0, ]$`Case #`,
                                         age_death = STAB2_absinta_indiv[1:6 %% 2 == 1, ]$Age)
STAB2_absinta_indiv_parsed$age_death <- as.numeric(STAB2_absinta_indiv_parsed$age_death)
STAB2_absinta$NBB_case <- gsub("11_69", "11_069", STAB2_absinta$NBB_case)
STAB2_absinta$age_death <- STAB2_absinta_indiv_parsed$age_death[match(STAB2_absinta$NBB_case,
                                                                      STAB2_absinta_indiv_parsed$individualID)]

STAB2_herring_indiv <- data.frame(readxl::read_xlsx(STAB2_herring_indiv_f, sheet = 2))
View(STAB2_herring_indiv)
STAB2_herring_indiv <- STAB2_herring_indiv[!is.na(STAB2_herring_indiv$BrainBank), ]
STAB2_herring_indiv <- STAB2_herring_indiv[STAB2_herring_indiv$Stage == "Adult", ]

GSM6657986_cell <- read_table_fast(GSM6657986_cell_f, row.names = 1)
GSM6657986_cell <- GSM6657986_cell[!is.na(GSM6657986_cell$Individual_ID), ]
GSM6657986_indiv <- data.frame(readxl::read_xlsx(GSM6657986_indiv_f, sheet = 14))
View(GSM6657986_indiv)
colnames(GSM6657986_indiv) <- GSM6657986_indiv[1, ]
GSM6657986_indiv <- GSM6657986_indiv[2:nrow(GSM6657986_indiv), ]
GSM6657986_indiv$AGE <- as.numeric(GSM6657986_indiv$AGE)

GSM6657986_cell$age_death <- GSM6657986_indiv$AGE[match(GSM6657986_cell$Individual_ID, GSM6657986_indiv$Human_ID)]

GSE193688_indiv <- data.frame(readxl::read_xlsx(GSE193688_indiv_f, sheet = 1))
colnames(GSE193688_indiv) <- GSE193688_indiv[1, ]
GSE193688_indiv <- GSE193688_indiv[2:nrow(GSE193688_indiv), ]
GSE193688_indiv <- GSE193688_indiv[GSE193688_indiv$Diagnosis == "Unaffected Control", ]

head(SEA_AD_mDat)
tail(SEA_AD_mDat)
SEA_AD_mDat_toBind <- SEA_AD_mDat[, c("individualID", "ageDeath", "diagnosis")]
colnames(SEA_AD_mDat_toBind) <- gsub("ageDeath", "age_death",
                                     colnames(SEA_AD_mDat_toBind))
colnames(SEA_AD_mDat_toBind) <- gsub("diagnosis", "diagn",
                                     colnames(SEA_AD_mDat_toBind))
SEA_AD_mDat_toBind <- SEA_AD_mDat_toBind[SEA_AD_mDat_toBind$age_death != "90+", ]
SEA_AD_mDat_toBind$age_death <- as.numeric(SEA_AD_mDat_toBind$age_death)
SEA_AD_mDat_toBind$diagn[grepl("control", tolower(SEA_AD_mDat_toBind$diagn))] <- "Control"
SEA_AD_mDat_toBind$diagn[grepl("not applicable", tolower(SEA_AD_mDat_toBind$diagn))] <- "Control"
# Remove cancer sample
SEA_AD_mDat_toBind <- SEA_AD_mDat_toBind[!grepl("Cancer", SEA_AD_mDat_toBind$diagn), ]
SEA_AD_mDat_toBind$diagn[SEA_AD_mDat_toBind$diagn != "Control"] <- "Neurodegeneration"
table(SEA_AD_mDat_toBind$diagn)
# We dont know how many cells are in SEA_AD samples. So, let's count the
# average of cells that are on a single sample, and then infer basedd on that.

uniq_indivs <- unique(integ_mDat$individualID)
uniq_cell_types <- unique(integ_mDat$cell_type)
ctype_ind_df <- data.frame(matrix(nrow = 0, ncol = length(uniq_cell_types),
                                  dimnames = list(NULL, uniq_cell_types)))

for (i in seq_along(uniq_indivs)){
        indiv <- uniq_indivs[i]
        ind_cells <- table(integ_mDat$cell_type[integ_mDat$individualID == indiv])
        names(ind_cells)[names(ind_cells) == ""] <- "V17"
        toBind <- data.frame(matrix(rep(0, ncol(ctype_ind_df)),
                                    nrow = 1, ncol = ncol(ctype_ind_df),
                                    dimnames = list(indiv,
                                                    colnames(ctype_ind_df))))
        toBind[1, match(make.names(names(ind_cells)), colnames(toBind))] <- as.vector(ind_cells)
        ctype_ind_df <- rbind.data.frame(ctype_ind_df, toBind)
}

ctype_medians <- round(apply(ctype_ind_df, 2, median))

library(dplyr)
library(tidyr)
cell_counts_df <- tibble(cell_type = gsub(".", " ", names(ctype_medians), fixed = T),
                         count = ctype_medians)
cell_counts_df$cell_type <- gsub("CD8  T Cells", "CD8+ T Cells", cell_counts_df$cell_type)

SEA_AD_mDat_toBind <- as_tibble(SEA_AD_mDat_toBind)

expanded_df <- SEA_AD_mDat_toBind %>%
        crossing(cell_counts_df) %>%  # Cartesian join: repeat each sample for each cell type
        uncount(count)
expanded_df <- as.data.frame(expanded_df)
expanded_df <- expanded_df[, c("individualID", "age_death", "cell_type", "diagn")]
expanded_df$study <- "SEA_AD"

# Expand as well Herring and GSE193688, as we dont have annotation data
# available

STAB2_herring_indiv_toBind <- STAB2_herring_indiv
STAB2_herring_indiv_toBind <- STAB2_herring_indiv_toBind[, c("snRNA.seq.library.ID", "Age")]
STAB2_herring_indiv_toBind$diagn <- "Control"
colnames(STAB2_herring_indiv_toBind) <- gsub("Age", "age_death", colnames(STAB2_herring_indiv_toBind))
STAB2_herring_indiv_toBind$age_death <- as.numeric(gsub("yr", "", STAB2_herring_indiv_toBind$age_death))
colnames(STAB2_herring_indiv_toBind) <- gsub("snRNA.seq.library.ID", "individualID", colnames(STAB2_herring_indiv_toBind))
STAB2_herring_indiv_toBind <- as_tibble(STAB2_herring_indiv_toBind)
herring_expanded <- STAB2_herring_indiv_toBind %>%
        crossing(cell_counts_df) %>%
        uncount(count)
herring_expanded <- as.data.frame(herring_expanded)
herring_expanded <- herring_expanded[, c("individualID", "age_death", "cell_type", "diagn")]
herring_expanded$study <- "STAB2_herring"

GSE193688_indiv_toBind <- GSE193688_indiv
head(GSE193688_indiv_toBind)
GSE193688_indiv_toBind <- GSE193688_indiv_toBind[, c("Sample ID", "Age", "Diagnosis")]
colnames(GSE193688_indiv_toBind) <- gsub("Sample ID", "individualID", colnames(GSE193688_indiv_toBind))
colnames(GSE193688_indiv_toBind) <- gsub("Age", "age_death", colnames(GSE193688_indiv_toBind))
colnames(GSE193688_indiv_toBind) <- gsub("Diagnosis", "diagn", colnames(GSE193688_indiv_toBind))
GSE193688_indiv_toBind$diagn[grep("Control", GSE193688_indiv_toBind$diagn)] <- "Control"
GSE193688_indiv_toBind$age_death <- as.numeric(GSE193688_indiv_toBind$age_death)
GSE193688_indiv_toBind <- as_tibble(GSE193688_indiv_toBind)
GSE193688_expanded <- GSE193688_indiv_toBind %>%
        crossing(cell_counts_df) %>%
        uncount(count)
GSE193688_expanded <- as.data.frame(GSE193688_expanded)
GSE193688_expanded <- GSE193688_expanded[, c("individualID", "age_death", "cell_type", "diagn")]
GSE193688_expanded$study <- "GSE193688"

head(STAB2_autism)
STAB2_autism <- STAB2_autism[, c("individual", "age", "cluster", "diagnosis")]
colnames(STAB2_autism) <- gsub("individual", "individualID", colnames(STAB2_autism))
colnames(STAB2_autism) <- gsub("age", "age_death", colnames(STAB2_autism))
colnames(STAB2_autism) <- gsub("cluster", "cell_type", colnames(STAB2_autism))
colnames(STAB2_autism) <- gsub("diagnosis", "diagn", colnames(STAB2_autism))
table(expanded_df$cell_type)
table(STAB2_autism$cell_type)
STAB2_autism$cell_type <- gsub("OPC", "OPCs", STAB2_autism$cell_type)
STAB2_autism$cell_type <- gsub("AST-FB|AST-PP", "Astrocytes", STAB2_autism$cell_type)
STAB2_autism$cell_type <- gsub("IN-PV|IN-SST|IN-SV2C|IN-VIP", "Inhibitory Neurons", STAB2_autism$cell_type)
STAB2_autism$cell_type <- gsub("L2/3|L4|L5/6|L5/6-CC|Neu-NRGN-I| Neu-NRGN-II", "Excitatory Neurons", STAB2_autism$cell_type)
STAB2_autism$cell_type <- gsub("Excitatory NeuronsI", "Excitatory Neurons", STAB2_autism$cell_type)
STAB2_autism$cell_type <- gsub("Neu-mat", "Maturing Neurons", STAB2_autism$cell_type)
STAB2_autism$study <- "STAB2_velmeshev"

head(STAB2_morabito)
STAB2_morabito <- STAB2_morabito[, c("Sample.ID", "Age", "celltype", "Diagnosis")]
colnames(STAB2_morabito) <- gsub("Sample.ID", "individualID", colnames(STAB2_morabito))
colnames(STAB2_morabito) <- gsub("Age", "age_death", colnames(STAB2_morabito))
colnames(STAB2_morabito) <- gsub("celltype", "cell_type", colnames(STAB2_morabito))
colnames(STAB2_morabito) <- gsub("Diagnosis", "diagn", colnames(STAB2_morabito))
table(STAB2_morabito$cell_type)
STAB2_morabito$cell_type <- gsub("ASC", "Astrocytes", STAB2_morabito$cell_type)
STAB2_morabito$cell_type <- gsub("EX", "Excitatory Neurons", STAB2_morabito$cell_type)
STAB2_morabito$cell_type <- gsub("INH", "Inhibitory Neurons", STAB2_morabito$cell_type)
STAB2_morabito$cell_type <- gsub("MG", "Microglia", STAB2_morabito$cell_type)
STAB2_morabito$cell_type <- gsub("ODC", "Oligodendrocytes", STAB2_morabito$cell_type)
STAB2_morabito$cell_type <- gsub("OPC", "OPCs", STAB2_morabito$cell_type)
STAB2_morabito$cell_type <- gsub("PER.END", "Endothelial", STAB2_morabito$cell_type)
STAB2_morabito$diagn[STAB2_morabito$diagn != "Control"] <- "Neurodegeneration"
STAB2_morabito$study <- "STAB2_morabito"

STAB2_franjic <- STAB2_franjic_cellmeta
STAB2_franjic$diagn <- "Control"
STAB2_franjic <- STAB2_franjic[, c("samplename", "age_death", "cluster", "diagn")]
colnames(STAB2_franjic) <- gsub("samplename", "individualID", colnames(STAB2_franjic))
colnames(STAB2_franjic) <- gsub("cluster", "cell_type", colnames(STAB2_franjic))
unique(STAB2_franjic$cell_type)
STAB2_franjic$cell_type <- gsub("CA3 CFAP299 SYN3|CA1 dorsal GRIK1 GRM3|CA2 CFAP299 HGF|CA1 ventral ACVR1C SYT13|DG GC PROX1 SGCZ|DG GC PROX1 PDLIM5|DG MC ARHGAP24 DLC1|SUB proximal ROBO1 SEMA3E|SUB distal FN1 NTNG1|SUB proximal ROBO1 COL5A2",
                                "Excitatory Neurons",
                                STAB2_franjic$cell_type)
STAB2_franjic$cell_type[grep("EC L", STAB2_franjic$cell_type)] <- "Excitatory Neurons"
STAB2_franjic$cell_type[grep("InN", STAB2_franjic$cell_type)] <- "Inhibitory Neurons"
STAB2_franjic$cell_type[grep("Astro", STAB2_franjic$cell_type)] <- "Astrocytes"
STAB2_franjic$cell_type[grep("Oligo", STAB2_franjic$cell_type)] <- "Oligodendrocytes"
STAB2_franjic$cell_type[grep("OPC", STAB2_franjic$cell_type)] <- "OPCs"
STAB2_franjic$cell_type[grep("Endo", STAB2_franjic$cell_type)] <- "Endothelial"
STAB2_franjic$cell_type[grep("Micro", STAB2_franjic$cell_type)] <- "Microglia"
STAB2_franjic$cell_type[grep("Myeloid", STAB2_franjic$cell_type)] <- "Myeloid Cells"
STAB2_franjic$cell_type[grep("Macro", STAB2_franjic$cell_type)] <- "Macrophages"
STAB2_franjic$cell_type[grep("SMC", STAB2_franjic$cell_type)] <- "SMC"
STAB2_franjic$cell_type[grep("COP", STAB2_franjic$cell_type)] <- "OPCs"
STAB2_franjic$cell_type[grep("CR", STAB2_franjic$cell_type)] <- "Maturing Neurons"
STAB2_franjic$cell_type[grep("T SKAP", STAB2_franjic$cell_type)] <- "CD8+ T Cells"
STAB2_franjic$cell_type[grep("PC", STAB2_franjic$cell_type)] <- "Pericytes"
STAB2_franjic$cell_type[grep("VLMC", STAB2_franjic$cell_type)] <- "Fibroblast"
table(STAB2_franjic$cell_type)
STAB2_franjic$study <- "STAB2_franjic"

STAB2_lake <- STAB2_lake_cell
STAB2_lake <- STAB2_lake[, c("individualID", "age_death", "Identity")]
STAB2_lake$diagn <- "Control"
STAB2_lake$study <- "STAB2_lake"
colnames(STAB2_lake) <- gsub("Identity", "cell_type", colnames(STAB2_lake))
head(STAB2_lake)
table(STAB2_lake$cell_type)
STAB2_lake$cell_type[grep("Ast", STAB2_lake$cell_type)] <- "Astrocytes"
STAB2_lake$cell_type[grep("End", STAB2_lake$cell_type)] <- "Endothelial"
STAB2_lake$cell_type[grep("Ex", STAB2_lake$cell_type)] <- "Excitatory Neurons"
STAB2_lake$cell_type[grep("In", STAB2_lake$cell_type)] <- "Inhibitory Neurons"
STAB2_lake$cell_type[grep("Mic", STAB2_lake$cell_type)] <- "Microglia"
STAB2_lake$cell_type[grep("Oli", STAB2_lake$cell_type)] <- "Oligodendrocytes"
STAB2_lake$cell_type[grep("OPC", STAB2_lake$cell_type)] <- "OPCs"
STAB2_lake$cell_type[grep("Per", STAB2_lake$cell_type)] <- "Pericytes"
STAB2_lake$cell_type[grep("Purk", STAB2_lake$cell_type)] <- "Inhibitory Neurons"
STAB2_lake <- STAB2_lake[STAB2_lake$cell_type != "NoN", ]
STAB2_lake <- STAB2_lake[STAB2_lake$cell_type != "Gran", ]

head(STAB2_jakel)
STAB2_jakel <- STAB2_jakel[, c("Sample", "age_death", "Celltypes", "Condition")]
colnames(STAB2_jakel) <- gsub("Sample", "individualID", colnames(STAB2_jakel))
colnames(STAB2_jakel) <- gsub("Celltypes", "cell_type", colnames(STAB2_jakel))
colnames(STAB2_jakel) <- gsub("Condition", "diagn", colnames(STAB2_jakel))
STAB2_jakel$diagn[STAB2_jakel$diagn == "Ctrl"] <- "Control"
table(STAB2_jakel$cell_type)
STAB2_jakel$cell_type[grep("Astrocytes", STAB2_jakel$cell_type)] <- "Astrocytes"
STAB2_jakel$cell_type[grep("Endothelial", STAB2_jakel$cell_type)] <- "Endothelial"
STAB2_jakel$cell_type[grep("COPs", STAB2_jakel$cell_type)] <- "OPCs"
STAB2_jakel$cell_type[grep("Oligo", STAB2_jakel$cell_type)] <- "Oligodendrocytes"
STAB2_jakel$cell_type[grep("ImOlGs", STAB2_jakel$cell_type)] <- "Oligodendrocytes"
STAB2_jakel$cell_type[grep("Vasc", STAB2_jakel$cell_type)] <- "SMC"
STAB2_jakel$cell_type[grep("Neuron1|Neuron2|Neuron3", STAB2_jakel$cell_type)] <- "Excitatory Neurons"
STAB2_jakel$cell_type[grep("Neuron4|Neuron5", STAB2_jakel$cell_type)] <- "Inhibitory Neurons"
STAB2_jakel$cell_type[grep("Microglia", STAB2_jakel$cell_type)] <- "Microglia"
STAB2_jakel$study <- "STAB2_jakel"

head(STAB2_absinta)
STAB2_absinta <- STAB2_absinta[, c("NBB_case", "age_death", "cell_type", "pathology")]
colnames(STAB2_absinta) <- gsub("NBB_case", "individualID", colnames(STAB2_absinta))
colnames(STAB2_absinta) <- gsub("pathology", "diagn", colnames(STAB2_absinta))
STAB2_absinta$diagn <- "Control"
table(STAB2_absinta$cell_type)
STAB2_absinta$cell_type[grep("astrocytes", STAB2_absinta$cell_type)] <- "Astrocytes"
STAB2_absinta$cell_type[grep("oligodendrocytes", STAB2_absinta$cell_type)] <- "Oligodendrocytes"
STAB2_absinta$cell_type[grep("opc", STAB2_absinta$cell_type)] <- "OPCs"
STAB2_absinta$cell_type[grep("lymphocytes", STAB2_absinta$cell_type)] <- "CD8+ T Cells"
STAB2_absinta$cell_type[grep("immune", STAB2_absinta$cell_type)] <- "Microglia"
STAB2_absinta$cell_type[grep("vascular", STAB2_absinta$cell_type)] <- "Endothelial"
STAB2_absinta$cell_type[grep("neurons", STAB2_absinta$cell_type)] <- "Excitatory Neurons"
STAB2_absinta$study <- "STAB2_absinta"

head(GSM6657986_cell)
GSM6657986_cell <- GSM6657986_cell[, c("Individual_ID", "age_death", "Cell_type", "Condition")]
colnames(GSM6657986_cell) <- gsub("Individual_ID", "individualID", colnames(GSM6657986_cell))
colnames(GSM6657986_cell) <- gsub("Cell_type", "cell_type", colnames(GSM6657986_cell))
colnames(GSM6657986_cell) <- gsub("Condition", "diagn", colnames(GSM6657986_cell))
GSM6657986_cell$diagn[GSM6657986_cell$diagn == "WT"] <- "Control"
GSM6657986_cell$diagn[GSM6657986_cell$diagn == "AD"] <- "Neurodegeneration"
table(GSM6657986_cell$cell_type)
GSM6657986_cell$cell_type[grep("Oligodendrocyte progenitor cells", GSM6657986_cell$cell_type)] <- "OPCs"
GSM6657986_cell$cell_type[grep("Cortical projection neurons", GSM6657986_cell$cell_type)] <- "Excitatory Neurons"
GSM6657986_cell$cell_type[grep("Interneurons", GSM6657986_cell$cell_type)] <- "Inhibitory Neurons"
GSM6657986_cell$cell_type[grep("Vascular leptomeningeal cells", GSM6657986_cell$cell_type)] <- "Fibroblast"
GSM6657986_cell <- GSM6657986_cell[GSM6657986_cell$cell_type != "Choroid plexus epithelial cells", ]
GSM6657986_cell$cell_type[grep("Endothelial", GSM6657986_cell$cell_type)] <- "Endothelial"
GSM6657986_cell <- GSM6657986_cell[GSM6657986_cell$cell_type != "Ependymal cells", ]
GSM6657986_cell$study <- "GSM6657986"

integ_mDat_all <- rbind.data.frame(integ_mDat, expanded_df,
                                   STAB2_autism,
                                   STAB2_morabito,
                                   STAB2_franjic,
                                   STAB2_lake,
                                   STAB2_jakel,
                                   STAB2_absinta,
                                   GSM6657986_cell,
                                   herring_expanded,
                                   GSE193688_expanded)

get_cell_bplot <- function(met_dat, samps = "healthy"){
        if (samps == "healthy"){
                bplot_df <- data.frame(cell_type = names(table(met_dat$cell_type[met_dat$diagn == "Control"])),
                                       counts = as.vector(table(met_dat$cell_type[met_dat$diagn == "Control"])))
        }else{
                bplot_df <- data.frame(cell_type = names(table(met_dat$cell_type[met_dat$diagn != "Control"])),
                                       counts = as.vector(table(met_dat$cell_type[met_dat$diagn != "Control"])))
        }
        bplot_df$cell_type <- factor(bplot_df$cell_type, levels = bplot_df$cell_type[order(bplot_df$counts)])
        cells_bplot <- ggplot(data = bplot_df,
                              mapping = aes(y = cell_type,
                                            x = counts)) +
                geom_bar(stat = "identity") +
                labs(x = "Number of cells", y = "Cell type") +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA,
                                                  linewidth = 1),
                      #panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        return(cells_bplot)
}

get_age_hist <- function(met_dat, samps){
        if (samps == "healthy"){
                aplot_df <- data.frame(individualID = unique(met_dat$individualID[met_dat$diagn == "Control"]))
        }else{
                aplot_df <- data.frame(individualID = unique(met_dat$individualID[met_dat$diagn != "Control"]))
        }
        aplot_df$age_death <- met_dat$age_death[match(aplot_df$individualID,
                                                      met_dat$individualID)]
        ages_hist <- ggplot(aplot_df, aes(round(age_death))) +
                geom_histogram(binwidth = 1) +
                labs(x = "Age (years)", y = "Subject count") +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA,
                                                  linewidth = 1),
                      #panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        return(ages_hist)
}

integ_mDat_all <- integ_mDat_all[integ_mDat_all$cell_type != "", ]
cell_bplot_healthy <- get_cell_bplot(met_dat = integ_mDat_all, samps = "healthy") +
        scale_x_log10()
cell_bplot_neurodg <- get_cell_bplot(met_dat = integ_mDat_all, samps = "neurodeg") +
        scale_x_log10()

age_hist_healthy <- get_age_hist(met_dat = integ_mDat_all, samps = "healthy")
age_hist_neurodg <- get_age_hist(met_dat = integ_mDat_all, samps = "neurodeg")

plt <- ggarrange(cell_bplot_healthy, age_hist_healthy,
                 cell_bplot_neurodg, age_hist_neurodg)

plt_healthy <- ggarrange(age_hist_healthy, cell_bplot_healthy)

plt_neurodg <- ggarrange(age_hist_neurodg, cell_bplot_neurodg)

ggsave(plot = plt, filename = sprintf("%shealth_and_dis_cells_and_ages.pdf", out_dir),
       height = 5, width = 7)

ggsave(plot = plt_healthy, filename = sprintf("%shealth_ages_and_cells.pdf", out_dir),
       height = 2.5, width = 6)

ggsave(plot = plt_neurodg, filename = sprintf("%sneurodg_ages_and_cells.pdf", out_dir),
       height = 2.5, width = 6)
