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

# Datadir contains the matrixes extraacted from the h5ad anndata file downloaded from GEO
data_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/GSE254569/counts/extracted_data/"
mdat_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/GSE254569/metadata/metadata/metadata.csv"


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
test <- subset(test, subset = celltype_final == "Exc_L2-3") #Subset Seurat object to one cell type

###### This is optional. You may want to exclude bad quality cells #######
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
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

#Perform shift-log normalization
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

colnames(dists) <- rownames(dists) <- donors

##### Make ages a named vector ######
ages_df <- meta[, c("donor", "age")]
ages_df <- ages_df[!duplicated(ages_df),]
ages_df$age <- as.numeric(gsub("m","",ages_df$age))

ages_vec <- ages_df$age
names(ages_vec) <- ages_df$donor

# Add a bit of noise to the ages to break ties between equal ages

set.seed(1)
ages_vec <- ages_vec + rnorm(length(ages_vec), mean = 0, sd = 0.000001)
####################################

#ages_vec <- ages_vec[ages_vec != 59]

#donors <- donors[donors %in% names(ages_vec)]

# Build directed graph with allowed edges (monotonically increasing attribute, no strict monotonicity required)
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

min_avg_weight_path_dp <- function(g, source, target, max_depth = 20) {
        #source <- "s55"
        #target <- "s12"
        nodes <- V(g)$name
        n <- length(nodes)
        
        dp <- matrix(Inf, nrow = n, ncol = max_depth + 1,
                     dimnames = list(nodes, 0:max_depth))
        pred <- matrix(NA, nrow = n, ncol = max_depth + 1,
                       dimnames = list(nodes, 0:max_depth))
        
        dp[source, "0"] <- 0
        for (k in 1:max_depth) {
                for (v in nodes) {
                        for (e in incident(g, v, mode = "in")) {
                                u <- ends(g, e)[1]
                                w <- E(g)[e]$weight
                                if (dp[u, as.character(k-1)] + w < dp[v, as.character(k)]) {
                                        dp[v, as.character(k)] <- dp[u, as.character(k-1)] + w
                                        pred[v, as.character(k)] <- u
                                }
                        }
                }
        }
        # Find minimum average
        k_vals <- 1:max_depth
        avg_vals <- dp[target, as.character(k_vals)] / k_vals
        best_k <- k_vals[which.min(avg_vals)]
        
        # Reconstruct path
        path <- target
        curr <- target
        for (k in best_k:1) {
                curr <- pred[curr, as.character(k)]
                path <- c(curr, path)
        }
        
        return(path)
}

min_avg_weight_path_dp <- function(g, source, target, max_depth = 20) {
        source <- "s55"
        target <- "s12"
        max_depth <- 20
        nodes <- V(g)$name
        n <- length(nodes)
        
        dp <- matrix(Inf, nrow = n, ncol = max_depth + 1,
                     dimnames = list(nodes, 0:max_depth))
        dp[source, "0"] <- 0
        
        prev <- matrix(NA, nrow = n, ncol = max_depth + 1,
                       dimnames = list(nodes, 0:max_depth))

        edges_mat <- ends(g, E(g))
        wghts <- E(g)$weight
        u_nodes <- edges_mat[, 1]
        v_nodes <- edges_mat[, 2]
        for (k in 1:max_depth) {
                #k <- 1
                candidate <- dp[u_nodes, as.character(k-1)] + wghts
                
                improved <- candidate < dp[v_nodes, k]
                dp[v_nodes[improved], as.character(k)] <- candidate[improved]
                
                # Record predecessor for path reconstruction
                prev[v_nodes[improved], as.character(k)] <- u_nodes[improved]
        }
        avg_weights <- dp[target, as.character(1:max_depth)] / (1:max_depth)
        # Find best path to end
        which.min(avg_weights)
        best_k <- which.min(dp[target, ])
        if (is.infinite(dp[target, best_k])) return(NULL)
        
        path <- character(best_k)
        path[best_k] <- target
        for (k in best_k:2) {
                path[k-1] <- prev[path[k], k]
        }
        return(path)
}


min_avg_weight_path_dp_vec <- function(g, start, end, max_depth = 20) {
        start <- "s55"
        end <- "s12"
        nodes <- V(g)$name
        n <- length(nodes)
        
        # Initialize DP: minimum total weight to reach node using k edges
        dp <- matrix(Inf, nrow = n, ncol = max_depth)
        rownames(dp) <- nodes
        dp[start, 1] <- 0
        
        # Predecessor matrix
        prev <- matrix(NA, nrow = n, ncol = max_depth)
        rownames(prev) <- nodes
        
        # Precompute edges
        edges_mat <- ends(g, E(g))
        weights <- E(g)$weight
        u_nodes <- edges_mat[,1]
        v_nodes <- edges_mat[,2]
        
        # DP over path lengths
        for (k in 2:max_depth) {
                candidate <- dp[u_nodes, k-1] + weights
                improved <- candidate < dp[v_nodes, k]
                dp[v_nodes[improved], k] <- candidate[improved]
                prev[v_nodes[improved], k] <- u_nodes[improved]
        }
        
        # Compute minimum average weight to 'end'
        avg_weights <- dp[end, ] / (1:max_depth)
        best_k <- which.min(avg_weights)
        if (is.infinite(dp[end, best_k])) return(NULL)
        
        # Reconstruct path
        path <- character(best_k)
        path[best_k] <- end
        for (k in best_k:2) {
                path[k-1] <- prev[path[k], k]
        }
        
        return(path)
}

g_full <- g
adjmat <- matrix(c(0, 2, 5, 0, 0, 1, 3, 0, 0),
                 ncol = 3, nrow = 3,
                 dimnames = list(c("A", "B", "C"),
                                 c("A", "B", "C")))

adjmat_2 <- matrix(c(0, 0, 0, 0, 0,
                     1, 0, 0, 0, 0,
                     0, 3, 0, 0, 0,
                     0, 8, 0, 0, 2,
                     0, 0, 1, 0, 0),
                   ncol = 5, nrow = 5,
                   dimnames = list(c("A", "B", "C", "D", "E"),
                                   c("A", "B", "C", "D", "E")))

g_test <- graph_from_adjacency_matrix(adjmatrix = adjmat_2,
                                      mode = "directed",
                                      weighted = T)

plot(g_test)

min_avg_weight_path_dp <- function(g, source, target) {
        #g <- g_sub
        #source <- "s55"
        #target <- "s12"
        
        Vn <- vcount(g)
        En <- ecount(g)
        w <- E(g)$weight
        edges <- ends(g, E(g), names = FALSE)
        sources <- edges[,1]
        targets <- edges[,2]
        
        # DP table: rows = path length (0..Vn), cols = vertices
        dp <- matrix(Inf, nrow = Vn + 1, ncol = Vn,
                     dimnames = list(0:Vn, names(V(g))))
        prev <- matrix(NA_integer_, nrow = Vn + 1, ncol = Vn,
                       dimnames = list(0:Vn, names(V(g))))
        
        # start at source
        dp[1, source] <- 0
        
        # Fill DP: O(V * E)
        i <- 1
        #for (k in 2:(Vn)) {
        #        dp[k, ] <- dp[k-1, ]
        #        for (e in seq_len(En)) {
        #                u <- sources[e]
        #                v <- targets[e]
        #                new_cost <- (dp[k-1, u] * i + w[e]) / (i + 1)
        #                
        #                #new_cost <- (dp[k-1, u] + w[e]) / i
        #                if (new_cost < dp[k, v]) {
        #                        dp[k, v] <- new_cost
        #                        prev[k, v] <- u
        #                }
        #        }
        #        i <- i + 1
        #}
        
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

min_avg_weight_path_dp <- function(g, source, target) {
        #g <- g_sub
        #g <- g_full
        #source <- "s55"
        #target <- "s12"
        Vn <- vcount(g)
        En <- ecount(g)
        w <- E(g)$weight
        edges <- ends(g, E(g), names = FALSE)
        sources <- edges[,1]
        targets <- edges[,2]
        
        # DP table: rows = path length (0..Vn), cols = vertices
        dp <- matrix(Inf, nrow = Vn + 1, ncol = Vn,
                     dimnames = list(0:Vn, names(V(g))))
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

v <- c(4, 5, 2, 6, 9)
avrg <- 5
for(i in 1:(length(v))){
        print(i)
        avrg <- (avrg * (i) + v[i]) / (i + 1)
        
}


min_avg_weight_path_dp(g, "s55", "s12")

scc <- components(g_full, mode = "strong")
scc$csize    # sizes of each SCC
scc$membership

table(scc$membership)

ages_df[match(names(scc$membership[scc$membership == 8]), ages_df$donor), ]
ages_df[match(names(scc$membership[scc$membership == 14]), ages_df$donor), ]
ages_df[match(names(scc$membership[scc$membership == 15]), ages_df$donor), ]
ages_df[match(names(scc$membership[scc$membership == 16]), ages_df$donor), ]
ages_df[match(names(scc$membership[scc$membership == 17]), ages_df$donor), ]

min_avg_weight_path_dp(g, "s55", "s12")

min_avg_weight_path_dp(g, "s35", "s9")
# Example: from node 1 to 5
result <- min_avg_weight_path(g, "s35", "s25")
print(result)



edge_mat <- ends(g, E(g))
edge_mat <- as.data.frame(edge_mat)
edge_mat$w <- E(g)$weight

#lapply(2:)

#pairs <- list(V(g)$name[c(1, 2)],
#              V(g)$name[c(1, 3)],
#              V(g)$name[c(1, 4)],
#              V(g)$name[c(1, 5)],
#              V(g)$name[c(1, 6)],
#              V(g)$name[c(1, 7)],
#              V(g)$name[c(1,8)],
#              V(g)$name[c(1, 9)],
#             V(g)$name[c(1, 10)],
#             V(g)$name[c(1,11)],
#             V(g)$name[c(1,12)],
#             V(g)$name[c(1,13)],
#             V(g)$name[c(1,14)],
#              V(g)$name[c(1,15)],
#              V(g)$name[c(1,16)],
#              V(g)$name[c(1,17)],
#              V(g)$name[c(1,18)],
#              V(g)$name[c(1,19)],
#              V(g)$name[c(1,20)],
#              V(g)$name[c(1,21)],
#              V(g)$name[c(1,22)],
#              V(g)$name[c(1,23)],
#             V(g)$name[c(1,24)],
#              V(g)$name[c(1,25)],
#              V(g)$name[c(1,26)],
#              V(g)$name[c(1,27)],
#              V(g)$name[c(1,28)],
#              V(g)$name[c(1,29)],
#              V(g)$name[c(1,30)],
#              V(g)$name[c(1,31)],
#              V(g)$name[c(1,32)],
#              V(g)$name[c(1,33)])

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

edge_list <- do.call(rbind, lapply(all_paths, function(p) {
  if (length(p) >= 2) {
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
plot(merged_subgraph, vertex.label=V(merged_subgraph)$name)

df <- as_data_frame(merged_subgraph)

write.table(df,file = "~/Downloads/Edges_test.txt",quote= F, sep = "\t",row.names = F, col.names = T)

# Optional: Visualize
plot(g, vertex.label=donors)

# Find a minimum spanning tree / arborescence
# Since we now have a DAG with weights, use shortest paths or minimal arborescence
# Find root candidate (lowest attribute)
root <- which.min(ages_vec)[1] #Need to set a root node. Not sure how to select if multiple donors have the same age

# Use Dijkstra's algorithm from root (handles weights)
paths <- shortest_paths(g, from = root, mode = "out", output = "both")

# Extract tree edges
tree_edges <- unlist(paths$epath)

# Tree subgraph
tree <- subgraph.edges(g, tree_edges)

# Plot the tree
plot(tree, vertex.label=donors)




# Bayesian graphs (networks)
