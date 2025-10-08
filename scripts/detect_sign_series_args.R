################################################################################
# Distinguish signal from noise in time series.                                #
################################################################################

if (!require(plotUtils)){
        devtools::install_github("guisantagui/plotUtils",
                                 upgrade = "never")
}
if(!require("Kendall", quietly = T)){
        install.packages("Kendall")
}
if(!require("ggplot2", quietly = T)){
        install.packages("ggplot2")
}
if(!require("boot", quietly = T)){
        install.packages("boot")
}
if(!require("boot", quietly = T)){
        install.packages("boot")
}
if(!require("ggpubr", quietly = T)){
        install.packages("ggpubr")
}
if(!require("tseries", quietly = T)){
        install.packages("tseries")
}

library(Kendall)
library(plotUtils)
library(argparser)
library(boot)
library(ggpubr)
library(tseries)

# Argument parser
################################################################################
parser <- arg_parser("Distinguish true signal from noise in time series.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--do_plots",
                               "--smooth_loess",
                               "--nboot",
                               "--outDir"),
                       help = c("CSV file, with series in columns and timepoints in rows.",
                                "If plots should be generated.",
                                "If series should be LOESS smoother before computing the p-values",
                                "Number of block bootstrap replicates",
                                "Output directory"),
                       flag = c(F, T, T, F, F))

parsed <- parse_args(parser)

# Inputs and directory
################################################################################
series_df_f <- parsed$input
do_plots <- parsed$do_plots
smooth_loess <- parsed$smooth_loess
nboot <- as.numeric(parsed$nboot)
outDir <- parsed$outDir

series_df_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/test_data/sample_series.csv"
do_plots <- T
smooth_loess <- T
nboot <- 5000
outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/results/series_stats/"

outDir <- add_slash_if_not(outDir)
create_dir_if_not(outDir)

# Functions
################################################################################
plot_prog <- function(x, y, get_pval = F){
        df <- data.frame(x = x,
                         y = y)
        plt <- ggplot(data = df, mapping = aes(x = x, y = y)) +
                geom_line()
        if (get_pval){
                box_pval <- Box.test(y,
                                     lag = floor(sqrt(length(y))),
                                     type = "Ljung-Box")$p.value
                mk_pval <- as.numeric(MannKendall(y)$sl)
                plt <- plt +
                        labs(caption = sprintf("Ljung-Box p-value = %s; Mann-Kendall p-value = %s",
                                               box_pval,
                                               mk_pval))
                out <- list(plot = plt, ljung_box_pval = box_pval,
                            mann_kendall_pval = mk_pval)
        }else{
                out <- list(plot = plt)
        }
        return(out)
}

mk_test <- function(series){
        res <- MannKendall(series)$tau
        return(res)
}

loess_smooth <- function(series, span = 0.2){
        #series <- series_df[, i] 
        time <- 1:length(series)
        loess_fit <- loess(series ~ time, span = span)
        smoothed_vals <- predict(loess_fit)
        return(smoothed_vals)
}

#plot_prog(x = 1:length(series), y = loess_smooth(series, span = .1))

#MannKendall(loess_smooth(series, span = .1))
# Load data
################################################################################
series_df <- read.csv(series_df_f, row.names = 1)

# Do tests
################################################################################
stats_df <- data.frame(matrix(nrow = 0, ncol = 6,
                              dimnames = list(NULL, c("gene",
                                                      "bj_pval",
                                                      "mk_pval",
                                                      "mk_boot_pval",
                                                      "adf_pval_stat",
                                                      "adf_pval_expl"))))
plot_list <- list()
for(i in 1:ncol(series_df)){
        ser <- series_df[, i]
        ser_presmooth <- ser
        if (smooth_loess){
                ser <- loess_smooth(ser, span = 0.1)
        }
        gene <- colnames(series_df)[i]
        bj_pval <- Box.test(ser,
                            lag = floor(sqrt(length(ser))),
                            type = "Ljung-Box")$p.value
        mk_pval <- as.numeric(MannKendall(ser)$sl)
        mk_boot <- tsboot(ser,
                          mk_test,
                          R = nboot,
                          l = round(length(ser)^(1/3)),
                          sim = "fixed")
        mk_boot_pval <- mean(abs(mk_boot$t) > abs(mk_test(ser)))
        adf_pval_stat <- adf.test(ser, alternative = "stationary")$p.value
        adf_pval_expl <- adf.test(ser, alternative = "explosive")$p.value
        
        to_bind <- data.frame(gene = gene,
                              bj_pval = bj_pval,
                              mk_pval = mk_pval,
                              mk_boot_pval = mk_boot_pval,
                              adf_pval_stat = adf_pval_stat,
                              adf_pval_expl = adf_pval_stat)
        stats_df <- rbind.data.frame(stats_df, to_bind)
        if (do_plots){
                plot_dir <- sprintf("%sseries_plots/",
                                    outDir)
                create_dir_if_not(plot_dir)
                if (smooth_loess){
                        plt_smooth <- plot_prog(x = 1:length(ser),
                                                y = ser,
                                                get_pval = T)$plot
                        plt_presmooth <- plot_prog(x = 1:length(ser_presmooth),
                                                   y = ser_presmooth,
                                                   get_pval = T)$plot
                        plt <- ggarrange(plt_presmooth, plt_smooth, common.legend = T)
                        plt_width <- 10
                }else{
                        plt <- plot_prog(x = 1:length(ser),
                                         y = ser,
                                         get_pval = T)$plot
                        plt_width <- 5
                }
                ggsave(filename = sprintf("%s%s.pdf",
                                          plot_dir,
                                          gene),
                       plot = plt,
                       height = 5,
                       width = plt_width)
                plot_list[[gene]] <- plt
        }
}

stats_df$bj_padj <- p.adjust(stats_df$bj_pval, method = "BH")
stats_df$mk_padj <- p.adjust(stats_df$mk_pval, method = "BH")
stats_df$mk_boot_padj <- p.adjust(stats_df$mk_boot_pval, method = "BH")
stats_df$adf_stat_padj <- p.adjust(stats_df$adf_pval_stat, method = "BH")
stats_df$adf_expl_padj <- p.adjust(stats_df$adf_pval_expl, method = "BH")

# Save results dataframe
################################################################################
write.csv(stats_df, file = sprintf("%sseries_stats.csv", outDir))
