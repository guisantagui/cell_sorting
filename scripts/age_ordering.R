library(Seurat)
if (!require(transport, quietly = T)) install.packages("transport")
library(transport)
if (!require(emdist, quietly = T)) install.packages("emdist")
library(emdist)
if (!require(philentropy, quietly = T)) install.packages("philentropy")
library(philentropy)
if (!require(irlba, quietly = T)) install.packages("irlba")
library(irlba)
if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes")
}
if (!require(SeuratDisk, quietly = T)) remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
if (!require(transformGamPoi, quietly = T)) BiocManager::install("transformGamPoi", update = F)
library(transformGamPoi)
if (!require(biomaRt, quietly = T)) BiocManager::install("biomaRt", update = F)
require(biomaRt)
if (!require(igraph, quietly = T)) install.packages("igraph")
library(igraph)

if(!require(energy, quietly = T)) install.packages("energy")
library(energy)

library(diptest)

# Datadir contains the matrixes extraacted from the h5ad anndata file downloaded from GEO
data_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/GSE254569/counts/extracted_data/"
mdat_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/GSE254569/metadata/metadata/metadata.csv"
target_cell <- "Exc_L2-3"
clock_file <- "/Users/guillem.santamaria/Documents/postdoc/manuscripts/brain_clock/submission_advanced_science/round_3/supplementary_material/supplementary_tables.xlsx"


clock_coefs <- as.data.frame(readxl::read_xlsx(clock_file))

test <- Seurat::Read10X(data.dir = data_dir)
test <- CreateSeuratObject(counts = test)
metdat <- plotUtils::read_table_fast(mdat_file)
feats <- read.table(list.files(data_dir, full.names = T)[2])

rownames(feats) <- feats[, 2]

test@meta.data$celltype_final <- metdat$celltypes_final[match(rownames(test@meta.data),
                                                               metdat$index)]
test@meta.data$celltype_major <- metdat$major_celltypes[match(rownames(test@meta.data),
                                                              metdat$index)]
test@meta.data$donor <- metdat$Donor[match(rownames(test@meta.data),
                                           metdat$index)]
test@meta.data$age <- metdat$Age[match(rownames(test@meta.data),
                                       metdat$index)]
test[["RNA"]]@meta.features <- feats

#test <- readRDS("~/Downloads/b752cc64-2ea2-4dc1-ade2-2c8e8fea949d.rds") #Read Seurat object
test <- subset(test, subset = celltype_final == target_cell) #Subset Seurat object to one cell type


toy_samps <- data.frame(donor = paste0("S", 1:length(rep(seq(1, 100, by = 10), each = 3))),
                        age = rep(seq(1, 100, by = 10), each = 3))
toy_centroids <- matrix(data = NA, nrow = 5, ncol = nrow(toy_samps),
                        dimnames = list(paste0("G", 1:5),
                                        toy_samps$donor))

gene_func1 <- function(age, sd = 2, intercept = 5, slope = 2){
        expr <- intercept + age * slope + rnorm(length(age), sd = sd)
        return(expr)
}

gene_func2 <- function(age, sd = 2, intercept = 5, slope = 2, thrshld = 60, slopemult = 2){
        func <- function(age_one, sd, intercept, slope, thrshld, slopemult){
                if (age_one <= thrshld){
                        expr <- intercept + age_one * slope + rnorm(length(age_one), sd = sd)
                }else{
                        y_thres_1 <- intercept + thrshld * slope
                        y_thres_2 <- intercept + thrshld * slope * slopemult
                        expr <- intercept - (y_thres_2 - y_thres_1) + age_one * (slope * slopemult) + rnorm(length(age_one), sd = sd)
                }
                return(expr)
        }
        out <- sapply(age, function(x) func(x, sd = sd, intercept = intercept,
                                            slope = slope, thrshld = thrshld,
                                            slopemult = slopemult))
        return(out)
}

gene_func3 <- function(age, sd = 2, slope = 2, age_start = 50){
        #age_start_lev <- age_start * slope
        expr <- -age_start * slope + age * slope + rnorm(length(age), sd = sd)
        expr[expr < 0] <- 0 + rnorm(length(expr[expr < 0]), sd = sd)
        return(expr)
}

toy_centroids[1, ] <- gene_func1(toy_samps$age, slope = 2, intercept = 5, sd = 10)
toy_centroids[2, ] <- gene_func1(toy_samps$age, slope = 3, intercept = 10, sd = 10)
toy_centroids[3, ] <- gene_func1(toy_samps$age, slope = 1, intercept = 10, sd = 10)
toy_centroids[4, ] <- gene_func1(toy_samps$age, slope = 2, intercept = 10, sd = 10)
toy_centroids[5, ] <- gene_func1(toy_samps$age, slope = 4.5, intercept = 10, sd = 10)

toy_centroids[1, 1:ncol(toy_centroids) %% 3 == 0] <- gene_func2(toy_samps$age[1:ncol(toy_centroids) %% 3 == 0], slope = 2, intercept = 5, sd = 10,
                                                                thrshld = 40, slopemult = 5)
toy_centroids[2, 1:ncol(toy_centroids) %% 3 == 0] <- gene_func2(toy_samps$age[1:ncol(toy_centroids) %% 3 == 0], slope = 3, intercept = 10, sd = 10,
                                                                thrshld = 40, slopemult = 5)
toy_centroids[3, 1:ncol(toy_centroids) %% 3 == 0] <- gene_func2(toy_samps$age[1:ncol(toy_centroids) %% 3 == 0], slope = 1, intercept = 10, sd = 10,
                                                                thrshld = 40, slopemult = 5)
toy_centroids[4, 1:ncol(toy_centroids) %% 3 == 0] <- gene_func2(toy_samps$age[1:ncol(toy_centroids) %% 3 == 0], slope = 2, intercept = 10, sd = 10,
                                                                thrshld = 40, slopemult = 5)
toy_centroids[5, 1:ncol(toy_centroids) %% 3 == 0] <- gene_func2(toy_samps$age[1:ncol(toy_centroids) %% 3 == 0], slope = 4.5, intercept = 10, sd = 10,
                                                                thrshld = 40, slopemult = 5)


# Simulate genes that get activated at 50 (genes 4 and 5), in one sample out
# of each 3

toy_centroids[4, ] <- rnorm(length(toy_centroids[4, ]), sd = 10)
toy_centroids[5, ] <- rnorm(length(toy_centroids[5, ]), sd = 10)

toy_centroids[4, 1:ncol(toy_centroids) %% 3 == 0] <- gene_func3(toy_samps$age[1:ncol(toy_centroids) %% 3 == 0], slope = 2, age_start = 50, sd = 10)
toy_centroids[5, 1:ncol(toy_centroids) %% 3 == 0] <- gene_func3(toy_samps$age[1:ncol(toy_centroids) %% 3 == 0], slope = 4.5, age_start = 50, sd = 10)

toy_centroids

# Filter out low quality cells
################################################################################

#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#ensToSymb <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=rownames(test),mart= mart)
#mt.features <- ensToSymb$ensembl_gene_id[which(grepl("^mt-",ensToSymb$mgi_symbol))]
mt.features <- rownames(test)[which(grepl("^MT-", rownames(test)))]
#ribo.features <- ensToSymb$ensembl_gene_id[which(grepl("(^Mrpl|^Mrps|^mrp)",ensToSymb$mgi_symbol))]
ribo.features <- rownames(test)[which(grepl("(^MRPL|^MRPS|^MRP)",
                                            rownames(test)))]

# Percent of mitochondrial counts
test[["percent.mito"]] <- PercentageFeatureSet(test, features = mt.features)

# Percent of mitochondrial ribosomal
test[["percent.ribo"]] <- PercentageFeatureSet(test, features = ribo.features)
###############################

# Perform shift-log normalization
################################################################################
DefaultAssay(test) <- "RNA"
mat <- LayerData(test, layer = "counts")
mat <- as.matrix(mat)
mat <- shifted_log_transform(mat)

meta <- test@meta.data[which(rownames(test@meta.data) %in% colnames(mat)),] #Extract metadata from Seurat object
all(rownames(meta) == colnames(mat)) #Check if metadata is concordant with matrix

donors <- as.character(unique(meta$donor)) #Extract donors. Each donor will represent one node.

centroids <- lapply(donors,function(x){ #Compute centroids as mean expression for each donor
  rowMeans(mat[,rownames(meta)[which(meta$donor == x)],drop = F])
})

centroids_mat <- do.call("cbind",centroids)
colnames(centroids_mat) <- donors

centroids_mat <- centroids_mat[rownames(centroids_mat) %in% clock_coefs$symbol, ]

centroids_mat <- toy_centroids

# Create age dataframe and age vector
################################################################################

##### Make ages a named vector ######
ages_df <- meta[, c("donor", "age")]
ages_df <- ages_df[!duplicated(ages_df),]
ages_df$age <- as.numeric(gsub("m","",ages_df$age))

ages_df <- toy_samps

rownames(ages_df) <- ages_df$donor



ages_vec <- ages_df$age
names(ages_vec) <- ages_df$donor

# Add a bit of noise to the ages to break ties between equal ages

set.seed(1)
ages_vec <- ages_vec + rnorm(length(ages_vec), mean = 0, sd = 0.000001)


# Do variable selection
################################################################################

# Obtain distance covariance with age
dcors <- apply(centroids_mat, 1, dcor, y = ages_df$age)

dcor(centroids_mat[1, ], ages_df$age)

mat_no_don <- mat


# Get variance of the genes
vars <- apply(mat[, colnames(mat) %in% rownames(meta)[meta$celltype_final == target_cell]], 1, var)

keep <- vars[vars >= quantile(vars, probs = 0.9)]

keep <- dcors[dcors >= quantile(dcors, probs = 0.5)]


win_frac <- .2
step_frac <- .1

age_range <- range(ages_df$age)

win_width <- diff(age_range) * win_frac
step <- diff(age_range) * step_frac

starts <- seq(age_range[1], age_range[2]-win_width, by=step)
window_pvals <- numeric(length = length(starts))

gene_vec <- mat[1, ]

for(i in seq_along(starts))
{
        #i <- 1
        idx <- which(meta$age >= starts[i] & meta$age <= starts[i] + win_width)
        if (length(idx) >= 5)
        {
                window_pvals[i] <- dip.test(gene_vec[idx])$p.value
        }
}
library(plotly)

get_win_stats <- function(gene_vec,
                          metdat,
                          win_var = "age",
                          starts,
                          step,
                          win_width,
                          plotdens = F)
{
        window_pvals <- numeric(length = length(starts))
        window_ds <- numeric(length = length(starts))
        dens_list <- list()
        for(i in seq_along(starts))
        {
                idx <- which(metdat[, win_var] >= starts[i] & metdat[, win_var] <= starts[i] + win_width)
                if (length(idx) >= 5)
                {
                        if (plotdens){
                                d <- density(gene_vec[idx], n = 50)
                                d_df <- data.frame(age_mid = mean(c(starts[i],
                                                                     starts[i] + win_width)),
                                                   expr = d$x,
                                                   density = d$y)
                                dens_list[[i]] <- d_df
                        }
                        res <- dip.test(gene_vec[idx])
                        window_pvals[i] <- res$p.value
                        window_ds[i] <- res$statistic
                }
        }
        if (plotdens){
                dens_list <- do.call(rbind, dens_list)
                plt <- plot_ly(dens_list, x = ~age_mid, y = ~expr, z = ~density,
                        type = "mesh3d",
                        intensity = ~density,
                        colors = colorRamp(c("white","blue")))
        }
        win_names <- make.names(paste(starts,
                                      starts + win_width, sep = "_"))
        names(window_pvals) <- win_names
        names(window_ds) <- win_names
        out <- data.frame(pval = window_pvals,
                          d_stat = window_ds)
        return(out)
}

test_modality <- function(gene_vec,
                          metdat,
                          reg_var = NULL,
                          win_var = "age",
                          win_frac = 0.2,
                          step_frac = 0.1,
                          nperms = 1)
{
        gene_vec <- mat[1, ]
        metdat <- meta
        reg_var = "donor"
        win_var = "age"
        win_frac = 0.2
        step_frac = 0.1
        df <- metdat
        df$expression <- gene_vec
        
        if (!is.null(reg_var)){
                if (reg_var %in% colnames(metdat)){
                        message(sprintf("Regressing out %s from the expression vector...",
                                        reg_var))
                        gene_vec <- residuals(lm(expression ~ donor, data = df))
                        message("Done.")
                }else if(!reg_var %in% colnames(metdat)){
                        stop(sprintf("%s is not present in the metadata, so it cannot be regressed out.",
                                     reg_var),
                             call. = F)
                }
        }
        
        
        age_range <- range(metdat[, win_var])
        win_width <- diff(age_range) * win_frac
        step <- diff(age_range) * step_frac
        
        age_range <- range(metdat[, win_var])
        starts <- seq(age_range[1],
                      age_range[2] - win_width,
                      by=step)
        
        real_stats <- get_win_stats(gene_vec = gene_vec,
                                    metdat = metdat,
                                    win_var = win_var,
                                    starts = starts,
                                    win_width = win_width,
                                    step = step,
                                    plotdens = T)
        
        seed <- 123
        perm_list <- list()
        if (nperms > 1){
                message("Computing p-values through a permutation test...")
                pb <- txtProgressBar(min = 0, max = nperms, style = 3)
                for(i in 1:nperms)
                {
                        #set.seed(seed + i)
                        setTxtProgressBar(pb, i)
                        perm_stats <- get_win_stats(gene_vec = sample(gene_vec,
                                                                      size = length(gene_vec),
                                                                      replace = F),
                                                    metdat = metdat,
                                                    win_var = win_var,
                                                    starts = starts,
                                                    win_width = win_width,
                                                    step = step)
                        perm_d <- perm_stats$d_stat
                        perm_list[[i]] <- perm_d
                }
                close(pb)
                perm_list <- do.call(cbind, perm_list)
                perm_pvals <- sapply(1:nrow(perm_list),
                                     function(x) sum(perm_list[x, ] > real_stats$d_stat[x])/length(perm_list[x, ]))
                names(perm_pvals) <- rownames(perm_list)
                real_stats$pval_perm <- perm_pvals
        }
        return(real_stats)
}

test_modality(gene_vec = mat[1, ],
              metdat = meta,
              reg_var = "donor",
              nperms = 0)

test_modality(gene_vec = centroids_mat[1, ],
              metdat = ages_df,
              reg_var = NULL,
              nperms = 1000)

centroids_mat <- centroids_mat[rownames(centroids_mat) %in% names(keep), ]

dists <- as.matrix(dist(t(centroids_mat))) #Compute distance between centroids

##### This is a test of the wasserstein metric #####
# dists <- matrix(0,nrow = length(donors),ncol = length(donors))
# for(i in 1:nrow(dists)){
#   for(j in 1:ncol(dists)){
#     if(i == j){
#       dists[i,j] <- 0
#     }else{
#       dists[i,j] <- transport::wasserstein1d(centroids_mat[,i],centroids_mat[,j])
#     }
#   }
# }

colnames(dists) <- rownames(dists) <- colnames(centroids_mat)


####################################

#ages_vec <- ages_vec[ages_vec != 59]

#donors <- donors[donors %in% names(ages_vec)]

# Build directed graph with allowed edges (monotonically increasing attribute, no strict monotonicity required)
donors <- ages_df$donor
g <- make_empty_graph(n = length(donors), directed = TRUE)



for (i in 1:length(donors)) {
        for (j in 1:length(donors)) {
                if (ages_vec[donors[j]] >= ages_vec[donors[i]] & i != j) {
                        g <- add_edges(g, c(i, j), attr = list(weight = dists[donors[i], donors[j]]))
                }
        }
}

V(g)$name <- donors

is_dag(g)

edge_mat <- as.data.frame(ends(g, E(g)))
edge_mat$w <- E(g)$weight

edge_mat[edge_mat$V1 == "s54" & edge_mat$V2 == "s65", ]
edge_mat[edge_mat$V1 == "s65" & edge_mat$V2 == "s54", ]

## More efficient version
FindCycles = function(g) {
        Cycles = NULL
        pb <- txtProgressBar(min = 0, max = length(V(g)), style = 3)
        i <- 0
        for(v1 in V(g)) {
                setTxtProgressBar(pb, i)
                if(degree(g, v1, mode="in") == 0) { next }
                GoodNeighbors = neighbors(g, v1, mode="out")
                GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
                for(v2 in GoodNeighbors) {
                        TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"),
                                         function(p) c(v1,p))
                        TempCyc = TempCyc[which(sapply(TempCyc, length) > 3)]
                        TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
                        Cycles  = c(Cycles, TempCyc)
                }
                i <- i + 1
        }
        return(Cycles)
}

find_all_cycles <- function(g)
{
        cycles <- NULL
        pb <- txtProgressBar(min = 0, max = length(V(g)), style = 3)
        i <- 0
        for(v1 in V(g))
        {
                setTxtProgressBar(pb, i)
                for(v2 in neighbors(g, v1, mode = "out"))
                {
                        cycles <- c(cycles,
                                    lapply(all_simple_paths(g, v2, v1, mode = "out"),
                                           function(p) c(v1, p)))
                }
                i <- i + 1
        }
        return(cycles)
}

find_all_cycles(g)

constructTree <- function(g,attribute,root){
        tree <- make_empty_graph(n =0)
        tree <- add_vertices(tree,1,attr = list(name = root))

        notAdded <- setdiff(V(g)$name,
                            root)
        notAdded <- notAdded[order(attribute[notAdded],decreasing = F)]
        for(v in notAdded){
                leafs <- V(tree)$name
    
                avg_path_weight <- rep(-1,length(leafs))
                names(avg_path_weight) <- leafs
                for(l in leafs){
                        tmp <- tree
                        tmp <- add_vertices(tmp, 1, attr = list(name = v))
                        tmp <- add_edges(tmp,c(l,v),attr=list(weight = E(g)[.from(l) & .to(v)]$weight))
                        all_paths <- all_simple_paths(tmp,root,to = c(setdiff(leafs,l),v))
                        path_weights <- sapply(all_paths, function(p) {
                                edges <- E(tmp, path = p)
                                total <- sum(E(tmp)[edges]$weight)
                                avg <- total / (length(p) - 1)
                                return(avg)
                        })
                        avg_path_weight[l] <- sum(path_weights)/length(all_paths)
                }
                selectedLeaf <- names(which.min(avg_path_weight))
                tree <- add_vertices(tree,1,attr = list(name = v))
                tree <- add_edges(tree,c(selectedLeaf,v),attr=list(weight = E(g)[.from(selectedLeaf) & .to(v)]$weight))
        }
        return(tree)
}

min_avg_weight_path <- function(g, from, to, max_depth = 20) {
        all_paths <- all_simple_paths(g, from = from, to = to, mode = "out", cutoff = max_depth)
  
        if (length(all_paths) == 0) return(NULL)
  
        path_weights <- sapply(all_paths, function(p) {
                edges <- E(g, path = p)
                total <- sum(E(g)[edges]$weight)
                avg <- total / (length(p) - 1)
                return(avg)
        })
  
        best_idx <- which.min(path_weights)
        return(all_paths[[best_idx]])
}

min_avg_weight_path_dp <- function(g, source, target) {
        Vn <- vcount(g)
        En <- ecount(g)
        w <- E(g)$weight
        edges <- ends(g, E(g), names = FALSE)
        sources <- edges[,1]
        targets <- edges[,2]
        
        # DP table: rows = path length (0..Vn), cols = vertices. Stores weight
        # (cost) between source node and each other node in the graph at a
        # particular number of steps
        dp <- matrix(Inf, nrow = Vn + 1, ncol = Vn,
                     dimnames = list(0:Vn, names(V(g))))
        # DP table: rows = path length (0..Vn), cols = vertices. Stores previous
        # node to reach a particular node. Used to reconstruct the path
        prev <- matrix(NA_integer_, nrow = Vn + 1, ncol = Vn,
                       dimnames = list(0:Vn, names(V(g))))
        
        # start at source
        dp[1, source] <- 0
        
        # Fill DP: O(V * E)
        for (k in 2:Vn) {
                # start fresh: paths with exactly k edges are not known yet
                dp[k, ] <- Inf  
                
                for (e in seq_len(En)) {
                        u <- sources[e]
                        v <- targets[e]
                        
                        if (is.finite(dp[k-1, u])) {
                                new_cost <- dp[k-1, u] + w[e]   # total cost for exactly k edges
                                if (new_cost < dp[k, v]) {
                                        dp[k, v] <- new_cost
                                        prev[k, v] <- u
                                }
                        }
                }
        }
        # Find the best average cost path to target
        best_avg <- Inf
        best_len <- NA
        for (k in 2:Vn) {
                if (is.finite(dp[k, target])) {
                        avg <- dp[k, target] / (k-1)
                        if (avg < best_avg) {
                                best_avg <- avg
                                best_len <- k
                        }
                }
        }
        
        if (is.infinite(best_avg)) return(NULL)
        
        # Reconstruct path
        path <- integer()
        v <- target
        k <- best_len
        while (!is.na(v) && k > 1) {
                path <- c(v, path)
                v <- names(V(g))[prev[k, v]]
                k <- k - 1
        }
        path <- c(source, path)
        
        return(V(g)[path])
}

min_avg_weight_path_dp(g, "s55", "s12")


pairs <- t(combn(V(g)$name,2))
pairs <- rbind(pairs,pairs[,c(2,1)])
pairs <- data.frame(pairs)

# Keep only pairs where age of the first element is less or equal than age of
# second.
pairs$age_1 <- ages_vec[match(pairs[, 1], names(ages_vec))]
pairs$age_2 <- ages_vec[match(pairs[, 2], names(ages_vec))]
pairs <- as.matrix(pairs[pairs$age_1 <= pairs$age_2, c(1, 2)])


pairs <- split(pairs,1:nrow(pairs))

library(parallel)

detectCores()

all_paths <- mclapply(pairs, function(pair)
{
        min_avg_weight_path_dp(g, pair[1], pair[2])
}, mc.cores = detectCores())

# Collect all unique paths
pb <- txtProgressBar(min = 0, max = length(pairs), style = 3)
all_paths <- list()
for (i in seq_along(pairs)) {
        pair <- pairs[[i]]
        setTxtProgressBar(pb = pb, value = i)
        path <- min_avg_weight_path_dp(g, pair[1], pair[2])
        if (!is.null(path)) {
                all_paths <- append(all_paths, list(path))
        }
}
close(pb)

edge_list <- do.call(rbind, lapply(all_paths, function(p) {if (length(p) >= 2) {
    from <- p[-length(p)]$name
    to <- p[-1]$name
    return(cbind(from, to))
  } else {
    return(NULL)
  }
}))

# Remove duplicate edges
edge_list <- unique(edge_list)

# Build a subgraph with those edges
merged_subgraph <- subgraph.edges(g, E(g, P = t(edge_list)))

# Plot the merged subgraph
plot(merged_subgraph, vertex.label=V(merged_subgraph)$name,
     vertex.size = 5)

#df <- as_data_frame(merged_subgraph)

#ages <- ages_df$age

#breaks1 <- seq(min(ages), max(ages), by = 5)
#breaks2 <- seq(min(ages) + 2.5, max(ages), by = 5)

#bins1 <- cut(ages, breaks = breaks1, include.lowest = T, right = F)
#bins2 <- cut(ages, breaks = breaks2, include.lowest = T, right = F)

#overlap_bins <- lapply(levels(bins1), function(b) which(bins1 == b))
#overlap_bins2 <- lapply(levels(bins2), function(b) which(bins2 == b))

#all_bins <- c(overlap_bins, overlap_bins2)
#names(all_bins) <- make.names(c(levels(bins1), levels(bins2)))

#bracket_mat <- matrix(data = 0, ncol = length(all_bins), nrow = length(ages_df$donor),
#                      dimnames = list(ages_df$donor, names(all_bins)))
#for(i in 1:nrow(bracket_mat)){
#        samp <- rownames(bracket_mat)[i]
#        bracket_mat[i, ] <- as.numeric(unlist(lapply(all_bins,
#                                                     function(x) i %in% x)))
#}

jump_thrshld <- 15
k <- 2
library(dplyr)


sub_g <- merged_subgraph
all_edges <- as.data.frame(ends(sub_g, E(sub_g)))
colnames(all_edges) <- c("from", "to")
all_edges$w <- E(sub_g)$weight
all_edges$age_diff <- ages_df$age[match(all_edges$to,
                                        ages_df$donor)] -
        ages_df$age[match(all_edges$from,
                          ages_df$donor)]

all_edges <- all_edges[all_edges$age_diff <= jump_thrshld, ]

best_edges <- all_edges %>%
        group_by(from) %>%
        mutate(min_w = min(w),
               z = (w - min_w)/sd(w)) %>%
        filter(z <= 1)

best_edges <- as.data.frame(best_edges)

best_edges <- as.matrix(best_edges[, 1:2])

subgraph_pruned <- subgraph.edges(sub_g, E(sub_g, P = t(best_edges)))

rownames(ages_df) <- ages_df$donor

plot_graph <- function(g, node_info, color){
        #g <- subgraph_pruned
        #node_info <- ages_df
        #color <- "age"
        colFun <- colorRampPalette(c("blue", "red"))
        round_up <- function(x, digits = 0) {
                multiplier <- 10^digits
                rounded <- ceiling(x * multiplier) / multiplier
                return(rounded)
        }
        if (min(node_info[, color]) < 0 & max(node_info[, color]) > 0){
                minMaxVec <- c(abs(min(node_info[, color])),
                               abs(max(node_info[, color])))
                lims <- round_up(max(minMaxVec), digits = 2)
                minMaxScale <- seq(-lims, lims, .01)
        }else{
                minMaxVec <- c(min(node_info[, color]),
                               max(node_info[, color]))
                lims <- round_up(minMaxVec, digits = 2)
                minMaxScale <- seq(lims[1], lims[2], .1)                
        }
        colBins <- colFun(length(minMaxScale))
        colScaleDF <- data.frame(value = minMaxScale,
                                 color = colBins)
        
        colVecGrads <- sapply(node_info[, color], 
                              function(x) colScaleDF$color[max(which(x >= colScaleDF$value))])
        names(colVecGrads) <- rownames(node_info)
        V(g)$color <- colVecGrads[match(make.names(names(V(g))),
                                        make.names(names(colVecGrads)))]
        plot(g,
             layout=layout.fruchterman.reingold,
             vertex.color=V(g)$color,
             edge.arrow.size=.5,
             vertex.size = 8)
}

plot(subgraph_pruned)

for(i in seq_along(V(g))){
        i <- 5
        print(i)
        samp <- V(g)[i]$name
        samp_edges <- as.data.frame(ends(g, E(g)))
        colnames(samp_edges) <- c("from", "to")
        samp_edges$w <- E(g)$weight
        samp_edges <- samp_edges[samp_edges[, 1] == samp, ]
        samp_edges$age_diff <- ages_df$age[match(samp_edges$to,
                                                 ages_df$donor)] -
                ages_df$age[match(samp_edges$from,
                                  ages_df$donor)]
        samp_edges <- samp_edges[samp_edges$age_diff <= jump_thrshld, ]
        best_edges <- samp_edges %>%
                group_by(from) %>%
                mutate(min_w = min(w),
                       z = (w - min_w)/sd(w)) %>%
                filter(z <= 1)
        best_edges <- as.data.frame(best_edges)
        pruned_edges <- rbind(pruned_edges,
                              as.matrix(best_edges[, c("from", "to")]))
}

ages_df

write.table(df,file = "~/Downloads/Edges_test.txt",quote= F, sep = "\t",row.names = F, col.names = T)

# Optional: Visualize
plot(g, vertex.label=donors)

# Find a minimum spanning tree / arborescence
# Since we now have a DAG with weights, use shortest paths or minimal arborescence
# Find root candidate (lowest attribute)
root <- which.min(ages_vec)[1] #Need to set a root node. Not sure how to select if multiple donors have the same age

# Use Dijkstra's algorithm from root (handles weights)
paths <- shortest_paths(merged_subgraph, from = root, mode = "out", output = "both")

# Extract tree edges
tree_edges <- unlist(paths$epath)

# Tree subgraph
tree <- subgraph.edges(merged_subgraph, tree_edges)

# Plot the tree
plot(tree, vertex.label=donors)

plot_graph(tree, node_info = ages_df, color = "age")

get_mst <- function(g, node_info, variable, mode = "out", output = "both")
{
        #mode <- "out"
        #output <- "both"
        root <- rownames(node_info)[which.min(node_info[, variable])]
        paths <- shortest_paths(g, from = root, mode = mode, output = output)
        tree_edges <- unlist(paths$epath)
        tree <- subgraph.edges(g, tree_edges)
        return(tree)
}

library(RBGL)

get_opt_branching <- function(g)
{
        gNEL <- igraph::as_graphnel(g)
        arbo <- edmondsOptimumBranching(gNEL)
        edge_list <- t(arbo$edgeList)  # returns from-to as list
        tree <- graph_from_edgelist(edge_list, directed = TRUE)
        E(tree)$weight <- E(g)[match(apply(edge_list, 1, paste, collapse="-"),
                                     apply(as_edgelist(g),
                                           1,
                                           paste,
                                           collapse="-"))]$weight
        return(tree)
}

plot_graph(merged_subgraph, node_info = ages_df, color = "age")

plot_graph(get_opt_branching(merged_subgraph), node_info = ages_df, color = "age")
plot_graph(get_opt_branching(subgraph_pruned), node_info = ages_df, color = "age")


gNEL <- igraph::as_graphnel(merged_subgraph)


root <- rownames(node_info)[which.min(node_info[, age_variable])]

# Chu–Liu/Edmonds algorithm
arbo <- edmondsOptimumBranching(gNEL)


# get edge list from arborescence
edge_list <- t(arbo$edgeList)  # returns from-to as list

# convert to matrix for igraph
el <- do.call(rbind, lapply(edge_list, function(x) c(names(x[1]), names(x[2]))))

# create igraph object
tree <- graph_from_edgelist(edge_list, directed = TRUE)

# optionally add weights from original g
E(tree)$weight <- E(g)[match(apply(edge_list, 1, paste, collapse="-"),
                             apply(as_edgelist(g), 1, paste, collapse="-"))]$weight

plot_graph(tree, node_info = ages_df, color = "age")

plot_graph(get_mst(subgraph_pruned, node_info = ages_df, variable = "age"),
           node_info = ages_df, color = "age")

plot_graph(get_mst(merged_subgraph, node_info = ages_df, variable = "age"),
           node_info = ages_df, color = "age")

plot_graph(mst(g), node_info = ages_df, color = "age")
# Bayesian graphs (networks)
