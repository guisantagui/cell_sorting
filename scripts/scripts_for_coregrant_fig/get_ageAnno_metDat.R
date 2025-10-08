ageAnno <- readRDS("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/ageAnno_brain_annot_badam/Annotated_sce.rds")

colnames(ageAnno@assays@data@listData$counts)
ageAnno_metDat <- as.data.frame(ageAnno@colData)

metDatSCFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/aging/data/gkac847_supplemental_files/Supplementary Tables.xlsx"
metDatSC <- as.data.frame(readxl::read_xlsx(metDatSCFile, sheet = 1))
metDatSC <- metDatSC[(which(metDatSC$`Table S1. Age information of the sample in AgeAnno` == "scRNA") + 1):(which(metDatSC$`Table S1. Age information of the sample in AgeAnno` == "scATAC") - 1), ]
colnames(metDatSC) <- unlist(as.vector(metDatSC[1, ]))
metDatSC <- metDatSC[2:(nrow(metDatSC) - 1), ]
metDatSC <- metDatSC[grepl("brain",
                           metDatSC$sample), colnames(metDatSC) != "organs"]
rownames(metDatSC) <- 1:nrow(metDatSC)


# The numbered mid1, mid2, mid3 etc of the counts_sc DF correspond to the
# metdat samples in the same order, so let's match them
metDatSC$samp <- tolower(metDatSC$`age group`)
for(grp in unique(metDatSC$samp)){
        metDatSC$samp[metDatSC$samp == grp] <- paste0(grp,
                                                      1:sum(metDatSC$samp == grp))
}
metDatSC$samp == unique(metDatSC$samp)[1]
head(ageAnno_metDat)

ageAnno_metDat$age_death <- as.numeric(metDatSC$`age (year)`[match(ageAnno_metDat$orig.ident, metDatSC$samp)])
ageAnno_metDat <- ageAnno_metDat[!is.na(ageAnno_metDat$age_death), ]
write.csv(ageAnno_metDat, file = "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/ageAnno_brain_annot_badam/met_dat.csv")


