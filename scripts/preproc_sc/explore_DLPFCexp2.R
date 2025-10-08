################################################################################
# Plot number of cells and age distributions of the seurat objects from SC     #  
# rosmap exp2.                                                                 #
################################################################################

if (!require("ggplot2", quietly = T)){
        install.packages("ggplot2", repos='http://cran.us.r-project.org')
}
if (!require("ggpubr", quietly = T)){
        install.packages("ggpubr", repos='http://cran.us.r-project.org')
}
if(!require("plotUtils", quietly = T)){
        remotes::install_github("guisantagui/plotUtils",
                                quiet = T,
                                upgrade = "never")
}

library(ggplot2)
library(ggpubr)
library(plotUtils)
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Plot counts of each cell type and age distribution of the Seurat objects.")

parser <- add_argument(parser = parser,
                       arg = c("--seur_dir",
                               "--out_dir"),
                       help = c("Directory with the Seurat objects.",
                                "Output directory."),
                       flag = c(F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
seur_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/ROSMAP_sc/seur_objs/"
out_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/results/counts_data/ROSMAP_sc/exploration_plots/"

seur_dir <- add_slash_if_not(parsed$seur_dir)
out_dir <- add_slash_if_not(parsed$out_dir)

create_dir_if_not(out_dir)


# Get metadata of the Seurat objects in the seur_dir
################################################################################
seur_files <- list.files(seur_dir)

metDat <- data.frame()
for(s in seur_files){
        seur <- readRDS(sprintf("%s%s", seur_dir, s))
        to_bind <- seur@meta.data
        colnames(to_bind)[grepl("pANN|DF.classifications",
                                colnames(to_bind))] <- gsub("\\_.*",
                                                            "",
                                                            colnames(to_bind)[grepl("pANN|DF.classifications",
                                                                                    colnames(to_bind))])
        if (nrow(metDat) == 0){
                metDat <- to_bind
        }else{
                metDat <- rbind.data.frame(metDat, to_bind)
        }
}

get_cell_bplot <- function(met_dat, samps = "healthy"){
        if (samps == "healthy"){
                bplot_df <- data.frame(cell_type = names(table(met_dat$cell_type[met_dat$cogdx == 1])),
                                       counts = as.vector(table(met_dat$cell_type[met_dat$cogdx == 1])))
        }else{
                bplot_df <- data.frame(cell_type = names(table(met_dat$cell_type[met_dat$cogdx != 1])),
                                       counts = as.vector(table(met_dat$cell_type[met_dat$cogdx != 1])))
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
                aplot_df <- data.frame(individualID = unique(met_dat$individualID[met_dat$cogdx == 1]))
        }else{
                aplot_df <- data.frame(individualID = unique(met_dat$individualID[met_dat$cogdx != 1]))
        }
        aplot_df$age_death <- met_dat$age_death[match(aplot_df$individualID,
                                                      met_dat$individualID)]
        ages_hist <- ggplot(aplot_df, aes(round(age_death))) +
                geom_histogram(binwidth = 1) +
                labs(x = "Counts", y = "Age (years)") +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA,
                                                  linewidth = 1),
                      #panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        return(ages_hist)
}

cell_bplot_healthy <- get_cell_bplot(met_dat = metDat, samps = "healthy")
cell_bplot_neurodg <- get_cell_bplot(met_dat = metDat, samps = "neurodeg")

age_hist_healthy <- get_age_hist(met_dat = metDat, samps = "healthy")
age_hist_neurodg <- get_age_hist(met_dat = metDat, samps = "neurodeg")

plt <- ggarrange(cell_bplot_healthy, age_hist_healthy,
                 cell_bplot_neurodg, age_hist_neurodg)

ggsave(filename = sprintf("%shealth_and_dis_cells_and_ages.pdf", out_dir),
       height = 5, width = 6)