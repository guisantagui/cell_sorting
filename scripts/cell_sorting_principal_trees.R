
if(!require(Iso, quietly = T)) install.packages("Iso")
if(!require(DDRTree, quietly = T)){
        devtools::install_github("cole-trapnell-lab/DDRTree",
                                 ref="simple-ppt-like")
}
library(DDRTree)
library(ggplot2)
library(igraph)
library(Iso)
library(Seurat)
library(transformGamPoi)

outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/results/cellsort/"
seur_file <- sprintf("%sGSE254569_ctrls_seur.rds", outDir)

cell_ids_in <- "celltype_final"
target_cell <- "Exc_L2-3"

# cellsorter functions
################################################################################

# Filter Seurat object to keep only cells of a given cell type.
subset_celltype <- function(seur, cell_ids_in, target_cell){
        if (cell_ids_in %in% colnames(seur@meta.data)){
                if (sum(seur@meta.data[[cell_ids_in]] == target_cell) == 0){
                        warning(sprintf("No cells labeled as %s were found",
                                        target_cell),
                                call. = F)
                }
                seur <- seur[, seur@meta.data[[cell_ids_in]] == target_cell]
        }else{
                stop(sprintf("%s is not a column in seur meta.data slot.",
                             cell_ids_in),
                     call. = F)
        }
        return(seur)
}

create_cellsort_obj <- function(seur,
                                cell_ids_in = NULL,
                                target_cell = NULL,
                                assay = "RNA",
                                layer = "counts",
                                #dim_red = NULL,
                                n_var_features = NULL,
                                shiftlog = T){
        #assay <- "RNA"
        #layer <- "counts"
        if(class(seur)[1] != "Seurat"){
                stop("'seur' is not a Seurat object")
        }
        if (!is.null(cell_ids_in) & !is.null(target_cell)){
                seur <- subset_celltype(seur, cell_ids_in, target_cell)
        }
        if (!is.null(n_var_features)){
                if (is.numeric(n_var_features)){
                        message(sprintf("Identifying %s most variable genes...",
                                        n_var_features))
                        seur <- FindVariableFeatures(seur,
                                                     nfeatures = n_var_features)
                        genes <- VariableFeatures(seur)
                }else{
                        stop("Invalid 'n_var_features' value", call. = F)
                }
        }else{
                genes <- rownames(seur)
        }
        
        DefaultAssay(seur) <- assay
        mat <- LayerData(seur, layer = "counts")
        #mat <- as.matrix(seur)
        if (shiftlog){
                mat <- shifted_log_transform(mat)
        }
        
        mat <- mat[genes, ]
        meta <- seur@meta.data
        obj <- list(cell_expr = mat,
                    cell_expr_stand = NULL,
                    #cent_expr = NULL,
                    #cent_expr_stand = NULL,
                    means = NULL,
                    sds = NULL,
                    cell_metadata = meta,
                    #centroid_metadata = NULL,
                    centroid_graph = NULL,
                    cell_graph = NULL)
        class(obj) <- "cellsort"
        return(obj)
}

get_cell_to_edges_dist <- function(cell, centroids_A, centroids_B) {
        v <- centroids_B - centroids_A
        u <- cell - centroids_A
        v2 <- colSums(v^2)
        t <- colSums(u * v) / v2
        t <- pmin(pmax(t, 0), 1)
        distances <- sqrt(colSums((u - v * matrix(t, nrow = nrow(v), ncol = length(t), byrow = TRUE))^2))
        return(data.frame(t = t, distance = distances))
}

# Given the normalized cell expression matrix, the normalized centroid
# expression matrix and the edge matrix of the skeleton graph, uses
# get_cell_to_edge_dist to comput the distances of each cell to each edge
# and assigns the cell to the edge with the minimum distance.
assign_cell_to_edge <- function(cell_mat, centroid_mat, edges,
                                parallelize = F,
                                cores = NULL){
        #cores <- detectCores() - 1
        #if (cores < 1){
        #        cores <- 1
        #}
        message("Assigning cells to best edge in backbone graph...")
        start_time <- Sys.time()
        
        
        
        centroids_A <- centroid_mat[, edges[, 1]]
        centroids_B <- centroid_mat[, edges[, 2]]
        
        cell_fun <- function(i) {
                cell_name <- colnames(cell_mat)[i]
                cell_vec <- cell_mat[, i]
                dists <- get_cell_to_edges_dist(cell_vec,
                                                centroids_A,
                                                centroids_B)
                j <- which.min(dists$distance)
                return(data.frame(from = edges[j, 1],
                                  to = edges[j, 2],
                                  dist = dists$distance[j],
                                  cell = cell_name))
        }
        if (!parallelize){
                cell_edge <- data.frame(matrix(nrow = ncol(cell_mat), ncol = 4,
                                               dimnames = list(NULL,
                                                               c("from",
                                                                 "to",
                                                                 "dist",
                                                                 "cell"))))
                pb <- txtProgressBar(min = 0, max = ncol(cell_mat), style = 3)
                for (i in 1:ncol(cell_mat)){
                        setTxtProgressBar(pb, i)
                        cell_name <- colnames(cell_mat)[i]
                        cell_vec <- cell_mat[, i]
                        dists <- get_cell_to_edges_dist(cell_vec,
                                                        centroids_A,
                                                        centroids_B)
                        tobind <- matrix(edges[which.min(dists$distance), ],
                                         nrow = 1,
                                         dimnames = list(cell_name,
                                                         c("from", "to")))
                        tobind <- data.frame(tobind)
                        tobind$dist <- min(dists$distance)
                        tobind$cell <- cell_name
                        cell_edge$from[i] <- tobind$from
                        cell_edge$to[i] <- tobind$to
                        cell_edge$dist[i] <- tobind$dist
                        cell_edge$cell[i] <- tobind$cell
                }
                close(pb)
        }else{
                if (is.null(cores)){
                        cores <- max(1, detectCores() - 1)
                }
                cell_edge <- mclapply(seq_len(ncol(cell_mat)),
                                      cell_fun,
                                      mc.cores = cores)
                cell_edge <- do.call(rbind, cell_edge)
        }
        end_time <- Sys.time()
        elapsed <- round(as.numeric(end_time - start_time, units = "mins"),
                         digits = 3)
        message(sprintf("Cells have been assigned to edges. %s mins elapsed.",
                        elapsed))
        
        return(cell_edge)
}

# Functions for creating the toy dataset

# Functions for creating toy dataset
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

# Plot the expression of two genes in the centroid matrix
plot_centroids <- function(centroid_mat, x, y){
        #x <- 1
        #y <- 2
        #centroid_mat <- t(scale(t(X)))
        centroid_mat <- as.matrix(centroid_mat)
        attributes(centroid_mat) <- attributes(centroid_mat)[c("dim", "dimnames")]
        plot_df <- data.frame(x = as.vector(centroid_mat[x, ]),
                              y = as.vector(centroid_mat[y, ]),
                              labels = colnames(centroid_mat))
        plt <- ggplot(plot_df,
                      mapping = aes(x = x, y = y, label = labels)) +
                geom_text() +
                xlab(rownames(centroid_mat)[x]) +
                ylab(rownames(centroid_mat)[y]) +
                theme_minimal()
        return(plt)
}

plot_cells <- function(cell_mat, x, y){
        #x <- 1
        #y <- 2
        #centroid_mat <- t(scale(t(X)))
        cell_mat <- as.matrix(cell_mat)
        attributes(cell_mat) <- attributes(cell_mat)[c("dim", "dimnames")]
        plot_df <- data.frame(x = as.vector(cell_mat[x, ]),
                              y = as.vector(cell_mat[y, ]),
                              labels = colnames(cell_mat),
                              donor = gsub("\\_.*", "", colnames(cell_mat)))
        plt <- ggplot(plot_df,
                      mapping = aes(x = x, y = y, col = donor)) +
                xlab(rownames(cell_mat)[x]) +
                ylab(rownames(cell_mat)[y]) +
                geom_point() +
                #geom_text() +
                theme_minimal()
        return(plt)
}

get_toy_cells <- function(centroid_mat, sd = 1, ncells = 100){
        cell_mat <- matrix(nrow = nrow(centroid_mat), ncol = 0,
                           dimnames = list(rownames(centroid_mat),
                                           NULL))
        for(j in 1:ncol(centroid_mat)){
                don_mat <- matrix(nrow = 0, ncol = ncells,
                                  dimnames = list(NULL,
                                                  paste(colnames(centroid_mat)[j],
                                                        paste0("C", 1:ncells),
                                                        sep = "_")))
                for(i in 1:nrow(centroid_mat)){
                        cell_vec <- rnorm(ncells, mean = centroid_mat[i, j], sd = sd)
                        to_bind <- matrix(cell_vec, nrow = 1,
                                          ncol = length(cell_vec),
                                          dimnames = list(rownames(centroid_mat)[i],
                                                          colnames(don_mat)))
                        don_mat <- rbind(don_mat, to_bind)
                }
                cell_mat <- cbind(cell_mat, don_mat)
        }
        return(cell_mat)
}

get_centroids <- function(cell_mat){
        donors <- unique(gsub("\\_.*", "", colnames(X_cells)))
        centroid_mat <- sapply(donors,
                               function(x) rowMeans(X_cells[, gsub("\\_.*",
                                                                   "",
                                                                   colnames(X_cells)) == x]))
        return(centroid_mat)
}

# Create toy dataset
################################################################################

toy_samps <- data.frame(donor = paste0("S", 1:length(rep(seq(1, 100, by = 10), each = 4))),
                        age = rep(seq(1, 100, by = 10), each = 4))
rownames(toy_samps) <- toy_samps$donor
toy_centroids <- matrix(data = NA, nrow = 2, ncol = nrow(toy_samps),
                        dimnames = list(paste0("G", 1:2),
                                        toy_samps$donor))

toy_centroids[1, ] <- gene_func1(toy_samps$age, slope = 2, intercept = 5, sd = 5)
toy_centroids[2, ] <- gene_func1(toy_samps$age, slope = 3, intercept = 10, sd = 5)

all_constant <- toy_centroids
activtline_3 <- toy_centroids
accelerate_3 <- toy_centroids

accelerate_3[1, 1:ncol(accelerate_3) %% 2 == 0] <- gene_func2(toy_samps$age[1:ncol(accelerate_3) %% 2 == 0],
                                                              slope = 2, intercept = 5, sd = 5,
                                                              thrshld = 20, slopemult = 8)
accelerate_3[2, 1:ncol(accelerate_3) %% 2 == 0] <- gene_func2(toy_samps$age[1:ncol(accelerate_3) %% 2 == 0],
                                                              slope = 3, intercept = 10, sd = 5,
                                                              thrshld = 20, slopemult = 8)

activtline_3[2, ] <- rnorm(length(activtline_3[2, ]), sd = 5)

activtline_3[2, 1:ncol(activtline_3) %% 2 == 0] <- gene_func3(toy_samps$age[1:ncol(activtline_3) %% 2 == 0],
                                                              slope = 2, age_start = 40, sd = 5)

X <- activtline_3

X_cells <- get_toy_cells(X, sd = 15, ncells = 10)

cell_metdat <- data.frame(cell_name = colnames(X_cells),
                          donor = gsub("\\_.*", "", colnames(X_cells)))

cell_metdat$age <- toy_samps$age[match(cell_metdat$donor,
                                       toy_samps$donor)]

#X <- plotUtils::stand(t(X))

plot_cells(X_cells, x = 1, y = 2)

X_back <- get_centroids(X_cells)

plot_centroids(X_back, 1, 2)

X <- X_back

samp_info <- toy_samps


tree <- DDRTree(
        X_cells,
        ncenter = length(unique(cell_metdat$donor))
)
library(DDRTree)

# Load real data
################################################################################

seur <- readRDS(seur_file)


cellsort <- create_cellsort_obj(seur = seur,
                                cell_ids_in = cell_ids_in,
                                target_cell = target_cell,
                                n_var_features = 5000)
start_time <- Sys.time()
pca <- irlba::prcomp_irlba(t(cellsort$cell_expr), n = 100, center = T, scale. = T)
end_time <- Sys.time()
elapsed <- round(as.numeric(end_time - start_time, unit = "mins"), digits = 3)
print(sprintf("PCA finished. %s mins elapsed.", elapsed))

sorted_age_PCs <- sort(abs(apply(pca$x, 2, function(x) cor(x, cellsort$cell_metadata$age))),
                       decreasing = T)

start_time <- Sys.time()
#tree <- DDRTree(as.matrix(cellsort$cell_expr), dimensions = 2, maxIter = 5, 
#                          ncenter = length(unique(cellsort$cell_metadata$donor)),
#                          verbose = T)

tree <- DDRTree(as.matrix(t(pca$x[, names(sorted_age_PCs)[1:2]])), dimensions = 2, maxIter = 5, 
                ncenter = length(unique(cellsort$cell_metadata$donor)),
                verbose = T,
                no_reduction = T)
end_time <- Sys.time()
elapsed <- round(as.numeric(end_time - start_time, unit = "mins"), digits = 3)
print(sprintf("DDRTree finished. %s mins elapsed.", elapsed))

plot(tree$Y[1, ], tree$Y[2, ])

Z <- t(tree$Z)  # cells × dims
Y <- t(tree$Y)  # centers × dims

# Compute squared distances
Z_sq <- rowSums(Z^2)         # 80606
Y_sq <- rowSums(Y^2)         # 33
cross <- Z %*% t(Y)         # 80606 x 33

dists_sq <- matrix(Z_sq, nrow=length(Z_sq), ncol=length(Y_sq)) +
        matrix(Y_sq, nrow=length(Z_sq), ncol=length(Y_sq), byrow=TRUE) -
        2 * cross

cell2center <- max.col(-dists_sq)

centage_tab <- table(cellsort$cell_metadata$age, cell2center)


cent_age <- apply(centage_tab, 2,
                  function(x) sum(as.numeric(rownames(centage_tab)) * x)/sum(x))

g <- graph_from_adjacency_matrix(as.matrix(dist(Y)),
                                 mode = "undirected",
                                 weighted = T,
                                 diag = T)

g <- mst(g)

leavs <- which(igraph::degree(g) == 1)

min_cent <- which.min(cent_age[names(leavs)])

get_pseudotime <- function(g, root_node){
        dists <- igraph::distances(g, 
                                   v = root_node, 
                                   to = V(g), 
                                   mode = "out")
        pseudo <- dists[1, ]
        pseudo <- (pseudo - min(pseudo, na.rm = TRUE)) / 
                (max(pseudo, na.rm = TRUE) - min(pseudo, na.rm = TRUE))
        return(pseudo)
}

pseudo <- get_pseudotime(g, names(min_cent))
ord <- order(pseudo)

plot(x = cent_age, y = pseudo)
isofit <- ufit(cent_age[ord], x = pseudo[ord], w = table(cell2center)[names(pseudo[ord])]/sum(table(cell2center)))

isofit <- ufit(cent_age[ord], x = pseudo[ord], w = table(cell2center)[names(pseudo[ord])]/sum(table(cell2center)))

cor(pseudo[names(pseudo) %in% names(cent_age)], cent_age, method = "spearman")


isofit$y

plot(isofit,type="l",
     xlab="pseudotime",ylab="centroid_age")
points(pseudo,cent_age,pch="+",col="red")
dev.off()

plot(y = isofit$y, x = pseudo)

isoresids <- cent_age[ord] - isofit$y



# compute Euclidean distance from each cell to each center
dmat <- as.matrix(dist(rbind(Z, Y)))
dmat <- dmat[1:nrow(Z), (nrow(Z)+1):(nrow(Z)+nrow(Y))]

# assign each cell to its closest center
cell2center <- max.col(-dmat)  # index of nearest Y for each cell


g <- graph_from_adjacency_matrix(
        as.matrix(tree_real_data$stree),
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE,
)

plot(tree$Z[1, ], tree$Z[2, ])
plot(tree_real_data$Y[1, ], tree_real_data$Y[2, ])

g <- graph_from_adjacency_matrix(
        as.matrix(tree$stree),
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE,
)
V(g)$name <- colnames(X_cells)
V(g)$age <- cell_metdat$age

edges <- as_edgelist(g)

endpoints <- data.frame(node = V(g)$name[degree(g) == 1],
                        age = V(g)$age[degree(g) == 1])



edges[edges[, 2] == endpoints[1], ]

root <- endpoints$node[which.min(endpoints$age)]
pseudotime <- distances(g, v = root)

df_reg <- data.frame(cell = names(pseudotime[1, ]),
                     pseudotime = pseudotime[1, ],
                     age = cell_metdat$age[match(names(pseudotime[1, ]),
                                                 cell_metdat$cell_name)])

ufit(df_reg$age, x = df_reg$pseudotime)

data("vigour")
par(mfrow=c(3,2),mar=c(4,4,3,1))
for(i in 2:6) {
        plot(ufit(vigour[,i],x=vigour[,1]),type="l",ylim=c(0,0.3),
             xlab="year",ylab="vigour",main=paste("stand",i-1),cex.main=1.5)
        points(vigour[,1],vigour[,i],pch="+",col="red")

}


plot(df_reg$pseudotime, df_reg$age)

as_edgelist(g)
plot(g)



data('iris')
subset_iris_mat <- as.matrix(t(iris[c(1, 2, 52, 103), 1:4])) #subset the data
#run DDRTree with ncenters equal to species number
DDRTree_res <- DDRTree(subset_iris_mat, dimensions = 2, maxIter = 5, sigma = 1e-2,
                       lambda = 1, ncenter = 3, param.gamma = 10, tol = 1e-2, verbose = FALSE)
Z <- DDRTree_res$Z #obatain matrix
Y <- DDRTree_res$Y
stree <- DDRTree_res$stree
plot(Z[1, ], Z[2, ], col = iris[c(1, 2, 52, 103), 'Species']) #reduced dimension
legend("center", legend = unique(iris[c(1, 2, 52, 103), 'Species']), cex=0.8,
       col=unique(iris[c(1, 2, 52, 103), 'Species']), pch = 1) #legend
title(main="DDRTree reduced dimension", col.main="red", font.main=4)
dev.off()
plot(Y[1, ], Y[2, ], col = 'blue', pch = 17) #center of the Z
title(main="DDRTree smooth principal curves", col.main="red", font.main=4)


#run DDRTree with ncenters equal to species number
DDRTree_res <- DDRTree(subset_iris_mat, dimensions = 2, maxIter = 5, sigma = 1e-3,
                       lambda = 1, ncenter = NULL,param.gamma = 10, tol = 1e-2, verbose = FALSE)
Z <- DDRTree_res$Z #obatain matrix
Y <- DDRTree_res$Y
stree <- DDRTree_res$stree
plot(Z[1, ], Z[2, ], col = iris[c(1, 2, 52, 103), 'Species']) #reduced dimension
legend("center", legend = unique(iris[c(1, 2, 52, 103), 'Species']), cex=0.8,
       col=unique(iris[c(1, 2, 52, 103), 'Species']), pch = 1) #legend
title(main="DDRTree reduced dimension", col.main="red", font.main=4)
dev.off()
plot(Y[1, ], Y[2, ], col = 'blue', pch = 2) #center of the Z
title(main="DDRTree smooth principal graphs", col.main="red", font.main=4)
