if (!require("synapser", quietly = T)){
        install.packages("synapser",
                         repos=c("http://ran.synapse.org",
                                 "https://cloud.r-project.org"))
}
if (!require("plotUtils", quietly = T)){
        devtools::install_github("guisantagui/plotUtils", upgrade = "never")
}
library(synapser)
library(plotUtils)

parent_synID <- "syn51792520"
filt_fastqs <- T
out_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/SEA_AD/synIDs/"
create_dir_if_not(out_dir)

files <- synGetChildren(parent_synID)
files <- as.list(files)
files_df <- data.frame(file_name = unlist(lapply(files, function(x) x$name)),
                       synID = unlist(lapply(files, function(x) x$id)))
if (filt_fastqs){
        files_df <- files_df[grepl("fastq", files_df$file_name), ]
}

write_table_fast(files_df,
                 f = sprintf("%s%s_children.csv",
                             out_dir,
                             parent_synID))
