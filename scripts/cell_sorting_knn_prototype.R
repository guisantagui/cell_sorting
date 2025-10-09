################################################################################
# Age Ordering KNN.                                                            #
################################################################################
library(ggplot2)
library(RBGL)
library(graph)
library(igraph)
library(parallel)
library(Matrix)
# Functions
################################################################################
out_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/results/graphs/"
plotUtils::create_dir_if_not(out_dir)


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


# Find paths with maximum average weight (minimize distances between paths)

max_avg_weight_path_dp <- function(g, source, target) {
        Vn <- vcount(g)
        En <- ecount(g)
        w <- E(g)$weight
        edges <- ends(g, E(g), names = FALSE)
        sources <- edges[,1]
        targets <- edges[,2]
        
        # DP table: rows = path length (0..Vn), cols = vertices. Stores weight
        # (cost) between source node and each other node in the graph at a
        # particular number of steps
        dp <- matrix(-Inf, nrow = Vn + 1, ncol = Vn,
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
                dp[k, ] <- -Inf  
                
                for (e in seq_len(En)) {
                        u <- sources[e]
                        v <- targets[e]
                        #if (donor_info$age[v] < donor_info$age[u]) next
                        if (is.finite(dp[k-1, u])) {
                                new_cost <- dp[k-1, u] + w[e]   # total cost for exactly k edges
                                if (new_cost > dp[k, v]) {
                                        dp[k, v] <- new_cost
                                        prev[k, v] <- u
                                }
                        }
                }
        }
        # Find the best average cost path to target
        best_avg <- -Inf
        best_len <- NA
        for (k in 2:Vn) {
                if (is.finite(dp[k, target])) {
                        avg <- dp[k, target] / (k-1)
                        if (avg > best_avg) {
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
        names(path) <- path
        #return(V(g)[path])
        return(path)
}

get_all_max_avg_weight_g <- function(g, parallelize = T)
{
        pairs <- t(combn(V(g)$name,2))
        pairs <- rbind(pairs,pairs[,c(2,1)])
        pairs <- data.frame(pairs)
        
        ages <- data.frame(node = names(V(g)),
                           age = V(g)$age)
        
        pairs$age_1 <- ages$age[match(pairs[, 1], ages$node)]
        pairs$age_2 <- ages$age[match(pairs[, 2], ages$node)]
        pairs <- as.matrix(pairs[pairs$age_1 <= pairs$age_2, c(1, 2)])
        pairs <- split(pairs,
                       1:nrow(pairs))
        if (!parallelize){
                pb <- txtProgressBar(min = 0, max = length(pairs), style = 3)
                all_paths <- list()
                for (i in seq_along(pairs)) {
                        pair <- pairs[[i]]
                        setTxtProgressBar(pb = pb, value = i)
                        path <- max_avg_weight_path_dp(g, pair[1], pair[2])
                        if (!is.null(path)) {
                                all_paths <- append(all_paths, list(path))
                        }
                }
                close(pb)
        }else{
                all_paths <- mclapply(pairs, function(pair)
                {
                        max_avg_weight_path_dp(g, pair[1], pair[2])
                }, mc.cores = detectCores())
                all_paths <- all_paths[unlist(lapply((all_paths),
                                                     function(x) !is.null(x)))]
        }
        
        edge_list <- do.call(rbind,
                             lapply(all_paths,
                                    function(p) {
                                            if (length(p) >= 2) {
                                                    from <- p[-length(p)]#$name
                                                    to <- p[-1]#$name
                                                    return(cbind(from, to))
                                            } else {
                                                    return(NULL)
                                            }
                                    }))
        
        # Remove duplicate edges
        edge_list <- unique(edge_list)
        g_pruned <- subgraph.edges(g, E(g, P = t(edge_list)))
        return(g_pruned)
}

# Find optimal branching, based in Edmonds algorithm
get_opt_branching <- function(g)
{
        gNEL <- igraph::as_graphnel(g)
        arbo <- edmondsOptimumBranching(gNEL)
        edge_list <- t(arbo$edgeList)  # returns from-to as list
        tree <- subgraph.edges(g, E(g, P = t(edge_list)))
        #E(tree)$weight <- E(g)[match(apply(edge_list, 1, paste, collapse="-"),
        #                             apply(as_edgelist(g),
        #                                   1,
        #                                   paste,
        #                                   collapse="-"))]$weight
        V(tree)$age <- V(g)$age[match(names(V(tree)),
                                      names(V(g)))]
        return(tree)
}


# Plot the graph
plot_graph <- function(g, node_info = NULL, color = NULL,
                       plot_labels = T,
                       node_size = 8,
                       edge_arrow_size = 0.5,
                       graph_layout = "fr",
                       root = NULL,
                       coords){
        
        #g <- cell_graph
        #node_info <- NULL
        #color <- "age"
        #node_size = 8
        #edge_arrow_size = 0.5
        #plot_labels <- F
        
        if (is.null(node_info)){
                attribs <- vertex_attr_names(g)
                node_info <- do.call(cbind.data.frame,
                                     lapply(attribs,
                                            function(x) vertex_attr(g,
                                                                    name = x)))
                colnames(node_info) <- attribs
                rownames(node_info) <- node_info$name
        }
        
        colFun <- colorRampPalette(c("blue", "red"))
        round_up <- function(x, digits = 0) {
                multiplier <- 10^digits
                rounded <- ceiling(x * multiplier) / multiplier
                return(rounded)
        }
        if (min(node_info[, color]) < 0 & max(node_info[, color]) > 0){
                minMaxVec <- c(abs(min(node_info[, color])),
                               abs(max(node_info[, color])))
                #lims <- round_up(max(minMaxVec), digits = 2)
                lims <- max(minMaxVec)
                minMaxScale <- seq(-lims, lims, .01)
        }else{
                minMaxVec <- c(min(node_info[, color]),
                               max(node_info[, color]))
                #lims <- round_up(minMaxVec, digits = 2)
                lims <- minMaxVec
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
        if (graph_layout == "fr"){
                lyout <- layout.fruchterman.reingold
        }else if (graph_layout == "tree"){
                if (is.null(root)){
                        stop("If graph_layout is set to 'tree', a root must be indicated", call. = F)
                }
                lyout <- layout_as_tree(g, root = root)
        }else if (graph_layout == "coords"){
                if (is.null(coords)){
                        stop("If graph_layout is set to 'coords', the 2D coords must be indicated", call. = F)
                }
                lyout <- coords[match(V(g)$name, rownames(coords)), ]
        }else{
                stop("Invalid graph_layout.", call. = F)
        }
        if (plot_labels){
                plot(g,
                     layout = lyout,
                     vertex.color = V(g)$color,
                     edge.arrow.size = edge_arrow_size,
                     vertex.size = node_size)
        }else{
                plot(g,
                     layout = lyout,
                     vertex.color = V(g)$color,
                     edge.arrow.size = edge_arrow_size,
                     vertex.size = node_size,
                     vertex.label = NA)
        }
        
}

# Computes KNN and builds a directed graph based on it, with the inverted 
# distances (euclidean) between the points as weights.
build_knn_graph <- function(X,
                            k = 6,
                            samp_info,
                            standardize = T){
        if (standardize){
                X <- plotUtils::stand(t(X))
        }else{
                X <- t(X)
        }
        
        nn <- FNN::get.knn(X, k = k)
        n <- nrow(nn$nn.dist)
        edges <- cbind(samp_info$donor[rep(1:n, k)],
                       samp_info$donor[as.vector(nn$nn.index)])
        # Prune edges that don't satisfy age constraint
        edges <- as.data.frame(edges)
        edges$age_1 <- samp_info$age[match(edges[, 1], samp_info$donor)]
        edges$age_2 <- samp_info$age[match(edges[, 2], samp_info$donor)]
        edges$w <- 1 / (1 + as.vector(nn$nn.dist))
        edges <- edges[edges$age_2 >= edges$age_1, ]
        w <- edges$w
        edges <- as.matrix(edges[, 1:2])
        # Build the graph
        g <- graph_from_edgelist(edges, directed = T)
        #V(g)
        E(g)$weight <- w
        V(g)$age <- samp_info$age[match(names(V(g)),
                                        samp_info$donor)]
        return(g)
}

# For samples that are the same age, groups them based on Louvain or Leiden
# community detection algorithms. Intended to be used along
# build_collapsed_sample_info and get_mean_close_samps to obtain merged
# centroids and avoid the cycles.
get_close_sameage_samps <- function(g, resolution = 1, alg = "louvain"){
        
        
        dup_ages <- unique(V(g)$age[duplicated(V(g)$age)])
        grpvec <- c()
        for (age in dup_ages){
                #age <- 51
                sameage_g <- as.undirected(subgraph(g,
                                                    names(V(g))[V(g)$age == age]),
                                           mode = "collapse",
                                           edge.attr.comb = "mean")
                subg_names <- names(V(sameage_g))
                if (alg == "louvain"){
                        clusts <- cluster_louvain(sameage_g, resolution = resolution)
                }else if (alg == "leiden"){
                        clusts <- cluster_leiden(sameage_g,
                                                 resolution_parameter = resolution,
                                                 objective_function = "modularity")#$membership
                }
                sameage_clust <- data.frame(samps = subg_names,
                                            membership = clusts$membership)
                grp <- sapply(unique(sameage_clust$membership),
                              function(x) paste(sameage_clust$samps[sameage_clust$membership == x],
                                                collapse = "-"))
                grpvec <- c(grpvec, grp)
        }
        return(grpvec)
}

# Given the sample_info object and a set of close_samps, generates a new
# sample_info object collapsing the samples that are similar into a single row.
build_collapsed_sample_info <- function(samp_info, close_samps){
        
        close_samps_split <- unique(unlist(sapply(close_samps,
                                                  function(x) strsplit(x,
                                                                       split = "-")[[1]])))
        collapsed_sampinfo <- data.frame(donor = c(close_samps,
                                                   samp_info$donor[!samp_info$donor %in% close_samps_split]))
        collapsed_sampinfo$age <- samp_info$age[match(collapsed_sampinfo$donor, samp_info$donor)]
        collapsed_sampinfo$age[is.na(collapsed_sampinfo$age)] <- as.numeric(sapply(close_samps,
                                                                                   function(x) paste(unique(samp_info$age[match(strsplit(x, "-")[[1]], samp_info$donor)]),
                                                                                                     collapse = "-")))
        collapsed_sampinfo$donor <- make.names(collapsed_sampinfo$donor)
        rownames(collapsed_sampinfo) <- collapsed_sampinfo$donor
        return(collapsed_sampinfo)
}

# Given the close_samps object generated by get_close_sameage_samps, generates
# a new centroid matrix aveaging the samples that are same age and that are
# close in the base graph.
get_mean_close_samps <- function(X, close_samps){
        #close_samps <- close_sameage_samps
        collapsed_df <- X
        for(i in seq_along(close_samps)){
                #i <- 1
                samp_group <- close_samps[i]
                samp_group_vec <- strsplit(samp_group, split = "-")[[1]]
                collapsed_df <- collapsed_df[, !colnames(collapsed_df) %in% samp_group_vec]
                to_bind <- data.frame(matrix(apply(X[, samp_group_vec],
                                                   1,
                                                   mean),
                                             nrow = nrow(X), ncol = 1,
                                             dimnames = list(rownames(X),
                                                             samp_group)))
                collapsed_df <- cbind.data.frame(collapsed_df, to_bind)
        }
        return(collapsed_df)
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

# KNN workflow
################################################################################
plot_centroids(all_constant, 1, 2)
plot_centroids(activtline_3, 1, 2)
plot_centroids(accelerate_3, 1, 2)


k <- 10

#X <- all_constant
X <- activtline_3
#X <- accelerate_3

X_cells <- get_toy_cells(X, sd = 15, ncells = 10)

cell_metdat <- data.frame(cell_name = colnames(X_cells),
                          donor = gsub("\\_.*", "", colnames(X_cells)))

cell_metdat$age <- toy_samps$age[match(cell_metdat$donor,
                                       toy_samps$donor)]

#X <- plotUtils::stand(t(X))

plot_cells(X_cells, x = 1, y = 2)
ggsave(sprintf("%scell_gene_exp.pdf", out_dir), width = 7, height = 8)


# Double check that cell inferrence is done correctly
get_centroids <- function(cell_mat){
        donors <- unique(gsub("\\_.*", "", colnames(X_cells)))
        centroid_mat <- sapply(donors,
                               function(x) rowMeans(X_cells[, gsub("\\_.*",
                                                                   "",
                                                                   colnames(X_cells)) == x]))
        return(centroid_mat)
}



X_back <- get_centroids(X_cells)

plot_centroids(X_back, 1, 2)
ggsave(sprintf("%scentroid_gene_exp.pdf", out_dir), width = 6, height = 8)


X <- X_back

samp_info <- toy_samps

# Build base graph
g_base <- build_knn_graph(X = X,
                          samp_info = samp_info,
                          k = k,
                          standardize = T)

pdf(file = sprintf("%scentroid_KNN_graph.pdf", out_dir))
plot_graph(g_base, node_info = samp_info, color = "age")
dev.off()

if(!is_dag(g_base)){
        close_sameage_samps <- get_close_sameage_samps(g_base,
                                                       resolution = 1,
                                                       alg = "louvain")
        collapsed_sampinfo <- build_collapsed_sample_info(samp_info,
                                                          close_samps = close_sameage_samps)
        
        X <- get_mean_close_samps(X, close_sameage_samps)
        
        g <- build_knn_graph(X = X,
                             samp_info = collapsed_sampinfo, k = k,
                             standardize = T)
        plot_graph(g, node_info = collapsed_sampinfo, color = "age")
}else{
        g <- g_base
}

g_pruned <- get_all_max_avg_weight_g(g)

pdf(file = sprintf("%scentroid_KNN_graph_collapsed_max_avg_paths.pdf", out_dir))
plot_graph(g_pruned, node_info = collapsed_sampinfo, color = "age")
dev.off()

g_pruned <- get_opt_branching(g_pruned)

pdf(file = sprintf("%scentroid_KNN_graph_collapsed_max_avg_paths_edmonds.pdf", out_dir))
plot_graph(g_pruned, node_info = collapsed_sampinfo, color = "age")
dev.off()

# Identify common nodes in a path. Use: find trunk of the tree.
intersect_recursive <- function(paths){
        if (length(paths) == 1) return(paths[[1]])
        if (length(paths) == 2) return(intersect(paths[[1]], paths[[2]]))
        first <- paths[[1]]
        rest <- intersect_recursive(paths[-1])
        intersect(first, rest)
}

# Identify last element in a linear branch.
last_in_branch <- function(branch, g, ori){
        next_node <- names(neighbors(g, ori, mode = "out"))
        if(length(next_node) == 0) return(ori)
        if (!all(next_node %in% branch)){
                return(ori)
        }else{
                ori <- next_node
                return(last_in_branch(branch, g, ori))
        }
}

# Given a cell vector and two centroid vectors, identifies the projection of 
# the cell over the vector delimited by centroids A and B, and returns distance
# between the cell and the projection --> Basically finds the length of the 
# perpendicular line between the cell point and the edge.
get_cell_to_edge_dist <- function(cell, centroid_A, centroid_B){
        v <- centroid_B - centroid_A
        u <- cell - centroid_A
        t <- sum(u * v) / sum(v^2)  # scalar projection
        t_clipped <- pmin(pmax(t, 0), 1) # restrict to edge segment
        distance <- sqrt(sum((u - t_clipped * v)^2))
        return(data.frame(t = t_clipped, distance = distance))
}

# Given the normalized cell expression matrix, the normalized centroid
# expression matrix and the edge matrix of the skeleton graph, uses
# get_cell_to_edge_dist to comput the distances of each cell to each edge
# and assigns the cell to the edge with the minimum distance.
assign_cell_to_edge <- function(cell_mat, centroid_mat, edges){
        cell_edge <- data.frame(matrix(nrow = 0, ncol = 4,
                                       dimnames = list(NULL,
                                                       c("from",
                                                         "to",
                                                         "dist",
                                                         "cell"))))
        message("Assigning cells to best edge in backbone graph...")
        pb <- txtProgressBar(min = 0, max = ncol(cell_mat), style = 3)
        for (i in 1:ncol(cell_mat)){
                setTxtProgressBar(pb, i)
                cell_name <- colnames(cell_mat)[i]
                cell_vec <- cell_mat[, i]
                dists <- apply(edges,
                               1,
                               function(x) get_cell_to_edge_dist(cell_vec,
                                                                 centroid_mat[, x[1]],
                                                                 centroid_mat[, x[2]]))
                dists <- do.call(rbind,
                                 dists)
                tobind <- matrix(edges[which.min(dists$distance), ], nrow = 1,
                                 dimnames = list(cell_name,
                                                 c("from", "to")))
                tobind <- data.frame(tobind)
                tobind$dist <- min(dists$distance)
                tobind$cell <- cell_name
                cell_edge <- rbind.data.frame(cell_edge, tobind)
        }
        close(pb)
        return(cell_edge)
}

# Multi-core version
get_cell_to_edges_dist <- function(cell, centroids_A, centroids_B) {
        v <- centroids_B - centroids_A
        u <- cell - centroids_A
        v2 <- colSums(v^2)
        t <- colSums(u * v) / v2
        t <- pmin(pmax(t, 0), 1)
        distances <- sqrt(colSums((u - v * matrix(t, nrow = nrow(v), ncol = length(t), byrow = TRUE))^2))
        return(data.frame(t = t, distance = distances))
}
assign_cell_to_edge <- function(cell_mat, centroid_mat, edges,
                                parallelize = F){
        cores <- max(1, detectCores() - 1)
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
                dists <- get_cell_to_edges_dist(cell_vec, centroids_A, centroids_B)
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
                cell_edge <- mclapply(seq_len(ncol(cell_mat)),
                                      cell_fun,
                                      mc.cores = cores)
                cell_edge <- do.call(rbind, cell_edge)
        }
        end_time <- Sys.time()
        message(sprintf("Cells have been assigned to edges. %s mins elapsed.",
                        round(end_time - start_time, digits = 3)))
        
        return(cell_edge)
}

# Given a matrix of cell-to-cell edges (cell_edges), the cell-backbone edge
# assignment matrix  (cell_edge_mat), the edges of the backbone graph (edges)
# and an initial edge (in_edge), recursively prunes the cell_edges that don't
# respect the connectivity imposed by the backbone graph.
prune_cell_edges <- function(cell_edges, cell_edge_mat, edges, in_edge) {
        # Ensure in_edge is a matrix
        if (!is.matrix(in_edge)) {
                in_edge <- matrix(in_edge, nrow = 1)
        }
        for (i in 1:nrow(in_edge)) {
                in_vec <- in_edge[i, ]
                if (!(in_vec[2] %in% edges[, 1])) {
                        return(cell_edges)
                } else {
                        out_edge <- edges[edges[, 1] == in_vec[2], , drop = FALSE]
                        # cells connected to incoming edge
                        cells_in <- cell_edge_mat$cell[
                                cell_edge_mat$from == in_vec[1] & cell_edge_mat$to == in_vec[2]
                        ]
                        # cells connected to outgoing edges
                        cells_out <- c()
                        for (j in 1:nrow(out_edge)) {
                                o <- out_edge[j, ]
                                co <- cell_edge_mat$cell[
                                        cell_edge_mat$from == o[1] & cell_edge_mat$to == o[2]
                                ]
                                cells_out <- c(cells_out, co)
                        }
                        # Discard cell edges that start in incoming edge
                        # but connect to cells that belong to edges
                        # other than incoming or outcoming.
                        bool_keep <- !(cell_edges[, 1] %in% cells_in &
                                               !(cell_edges[, 2] %in% c(cells_out, cells_in)))
                        cell_edges <- cell_edges[bool_keep, , drop = FALSE]
                        cell_edges <- prune_cell_edges(cell_edges,
                                                       cell_edge_mat,
                                                       edges, out_edge)
                }
        }
        return(cell_edges)
}

# Given the edges of the backbone graph (edges), its root (root), the
# cell-to-backbone edge assignments matrix (cell_edge_mat), the cell expression
# matrix (cell_mat), and a k value, computes a cell-to-cell edge matrix via KNN,
# prunes it based on the connectivity imposed by edges, and generates a
# cell-to-cell graph where weights are inversely proportional to Euclidean
# distance.
get_cell_graph <- function(edges,
                           root,
                           cell_edge_mat,
                           cell_mat,
                           k = (ncol(cell_mat) - 1),
                           cell_metdat,
                           cell_root,
                           algo = "mst"){
        #edges
        #root <- "youngest_cell"
        #cell_edge_mat <- cell_edge
        #k <- 100
        #cell_root <- youngest_cell
        
        # Build KNN graph of all the cells
        nn <- FNN::get.knn(t(cell_mat), k = k)
        n <- nrow(nn$nn.dist)
        cell_edges <- cbind(colnames(cell_mat)[rep(1:n, k)],
                            colnames(cell_mat)[as.vector(nn$nn.index)])
        # Obtain edge weights
        w <- 1 / (1 + as.vector(nn$nn.dist))
        names(w) <- paste(cell_edges[,1], cell_edges[,2], sep = "_")
        
        
        # Prune edges that do no satisfy backbone graph progression
        # based on edge assignments.
        in_edge <- edges[edges[, 1] == root, ]
        
        cell_edges <- prune_cell_edges(cell_edges, cell_edge_mat,
                                       edges, in_edge)
        w <- w[names(w) %in% paste(cell_edges[,1], cell_edges[,2], sep = "_")]
        cell_g <- graph_from_edgelist(cell_edges, directed = T)
        E(cell_g)$weight <- w
        if (algo == "mst"){
                cell_g <- graph_from_edgelist(cell_edges, directed = F)
                E(cell_g)$weight <- w
                tree <- mst(cell_g)
        }else if (algo == "edmonds"){
                cell_g <- graph_from_edgelist(cell_edges, directed = T)
                E(cell_g)$weight <- w
                nodes <- unique(c(cell_edges[,1], cell_edges[,2]))
                #gNEL <- new("graphNEL", nodes = nodes, edgemode = "directed")
                #edge_list <- split(cell_edges[,2], cell_edges[,1])
                #message("Creating graphNEL object...")
                gNEL <- ftM2graphNEL(cell_edges, W = w,
                                     V = nodes,
                                     edgemode="directed")
                
                
                #gNEL_edges <- lapply(nodes, function(n) {
                #        to_nodes <- edge_list[[n]]
                #        if (is.null(to_nodes)) to_nodes <- character(0)
                #        # each neighbor gets a weight
                #        list(edges = match(to_nodes, nodes), weights = w[cell_edges[,1] == n])
                #})
                #gNEL@edgeL <- gNEL_edges
                #names(gNEL@edgeL) <- nodes
                
                
                
                #pb <- txtProgressBar(max = nrow(cell_edges), style = 3)
                #for (i in seq_len(nrow(cell_edges))) {
                #        setTxtProgressBar(pb, value = i)
                #        from <- cell_edges[i,1]
                #        to <- cell_edges[i,2]
                #        gNEL <- addEdge(from, to, gNEL, weight = w[i])
                #}
                #close(pb)
                
                # Remove incoming edges to the youngest cell to force it as the root
                incoming <- inEdges(cell_root, gNEL)
                for (from in incoming) {
                        gNEL <- removeEdge(from, cell_root, gNEL)
                }
                
                arbo <- edmondsOptimumBranching(gNEL)
                edge_list <- t(arbo$edgeList)  # returns from-to as list
                tree <- subgraph.edges(cell_g, E(cell_g, P = t(edge_list)))
                E(tree)$weight <- w[match(paste(edge_list[, 1],
                                                edge_list[, 2],
                                                sep = "_"),
                                          names(w))]
        }
        V(tree)$age <- cell_metdat$age[match(names(V(tree)),
                                             cell_metdat$cell_name)]
        return(tree)
}


map_cells_to_backbone <- function(cell_mat,
                                  centroid_mat,
                                  backbone_graph,
                                  cell_metdat,
                                  k = 100,
                                  algo = "edmonds"){
        
        cell_mat <- X_cells
        backbone_graph <- g_pruned
        centroid_mat <- X
        k <- 100
        algo <- "edmonds"
        
        # Find leaves and root by finding endpoints of graph
        leaves <- names(which(igraph::degree(backbone_graph,
                                             mode = "out") == 0))
        root <- names(which(igraph::degree(backbone_graph,
                                           mode = "in") == 0))
        # Find all the paths from the root to the leaves
        all_paths <- list()
        for (leaf in leaves) {
                paths <- igraph::all_simple_paths(backbone_graph, from = root, to = leaf)
                all_paths <- c(all_paths, paths)
        }
        all_paths <- lapply(all_paths, names)
        
        # Find the trunk by finding the common nodes from across the paths.
        root_branch <- intersect_recursive(all_paths)
        
        # Find the last node in the trunk.
        last_root_branch <- last_in_branch(root_branch, g_pruned, root)
        
        # Normalice matrices, to the mean and sd of the centroids, as it's what
        # was used to build the backbone graph
        mean_centroids <- apply(centroid_mat, 1, mean)
        sd_centroids <- apply(centroid_mat, 1, sd)
        centroid_mat <- (centroid_mat - mean_centroids)/sd_centroids
        cell_mat <- (cell_mat - mean_centroids)/sd_centroids
        
        # Find the youngest cell by identifying the direction of change
        # across the trunk, and identifying the cell that goes most in the 
        # oposite direction.
        dir_root <- centroid_mat[, last_root_branch] - centroid_mat[, root]
        
        root_donors <- strsplit(root, split = ".", fixed = T)[[1]]
        root_cells <- cell_mat[, gsub("\\_.*",
                                      "",
                                      colnames(cell_mat)) %in% root_donors]
        
        # Also works, its the same. But crossprod is faster
        # youngest_cell <- colnames(root_cells)[which.min(as.vector(t(dir_root) %*% root_cells))]
        youngest_cell <- colnames(root_cells)[which.min(crossprod(root_cells,
                                                                  dir_root))]
        youngest_centroid <- root_cells[, youngest_cell]
        
        centroid_mat <- cbind(youngest_centroid,
                              centroid_mat)
        
        colnames(centroid_mat)[1] <- "youngest_cell"
        
        backbone_graph <- add_vertices(backbone_graph,
                                       1,
                                       name = "youngest_cell")
        
        dist_yng_cell <- sqrt(sum(youngest_centroid - centroid_mat[, root])^2)
        #1 / (1 + as.vector(nn$nn.dist))
        backbone_graph <- add_edges(backbone_graph, c("youngest_cell", root),
                                    attr = list(weight = 1 / (1 + dist_yng_cell)))
        
        # Assign cells to an edge
        edges <- as_edgelist(backbone_graph)
        
        # Compute projection of the vector delimited by the cell to the vector 
        # delimited by A and B (pair of centroids constituting an edge), and
        # then compute the distance between the cell and the projection
        # (perpendicular, or shortest path, from the cell to the edge). Does
        # this for each cell, and every edge, and then assigns the cell to the
        # centroid with the minimum distance
        cell_edge <- assign_cell_to_edge(cell_mat, centroid_mat, edges)
        
        # Prune the cell edges that don't respect the connectivity imposed by
        # the backbone graph.
        
        cell_g <- get_cell_graph(edges = edges,
                                 root = "youngest_cell",
                                 cell_edge_mat = cell_edge,
                                 cell_mat = cell_mat,
                                 k = k,
                                 cell_metdat = cell_metdat,
                                 cell_root = youngest_cell,
                                 algo = algo)
        return(cell_g)
}

cell_graph <- map_cells_to_backbone(cell_mat = X_cells,
                                    centroid_mat = X,
                                    backbone_graph = g_pruned,
                                    cell_metdat = cell_metdat,
                                    k = 20,
                                    algo = "edmonds")

pdf(sprintf("%scell_graph_fr.pdf", out_dir), width = 35, height = 35)
plot_graph(cell_graph,
           color = "age",
           plot_labels = F,
           node_size = 3,
           graph_layout = "fr")
dev.off()

pdf(sprintf("%scell_graph_tree.pdf", out_dir), width = 35, height = 35)
plot_graph(cell_graph,
           color = "age",
           plot_labels = F,
           node_size = 3,
           graph_layout = "tree",
           root = names(which(igraph::degree(cell_graph,
                                             mode = "in") == 0)))
dev.off()

coords <- apply(X_cells, 1, function(x) (x - mean(x))/sd(x))

pdf(sprintf("%scell_graph_coords.pdf", out_dir), width = 35, height = 35)
plot_graph(cell_graph,
           color = "age",
           plot_labels = F,
           node_size = 3,
           graph_layout = "coords",
           coords = coords)
dev.off()

# Smooth graph

igraph::get.adjacency(cell_graph)

W <- (as_adjacency_matrix(cell_graph, sparse = FALSE) +
              t(as_adjacency_matrix(cell_graph, sparse = FALSE))) / 2

D <- diag(rowSums(W))
L <- D - A

Y <- t(X_cells)

lambda <- 0.1   # smoothness
mu <- 0.05      # sparsity on edges

smoothness <- lambda * sum(diag(t(Y) %*% L %*% Y))
sparsity <- mu * sum(abs(W))

objective <- smoothness + sparsity



lambda <- 0.5   # smoothness strength
mu <- 0.1       # sparsity threshold for edges

# Simple iterative smoothing
for (iter in 1:20) {
        # Update Y: move each node toward weighted average of neighbors
        for (i in 1:nrow(Y)) {
                neighbors <- which(W[i,] > 0)
                if(length(neighbors) > 0){
                        Y[i,] <- (1 - lambda) * Y[i,] + lambda * colSums(W[i,neighbors] * Y[neighbors,]) / sum(W[i,neighbors])
                }
        }
        
        # Sparsity: remove very weak edges
        W[W < mu] <- 0
        
        # Recompute Laplacian
        D <- diag(rowSums(W))
        L <- D - W
}

# Plot smoothed tree
plot(graph_from_adjacency_matrix(W, mode="undirected", weighted=TRUE),
     layout=Y, edge.width=W*5, vertex.size=20)




# Compute pseudotime from a root node in a weighted, directed graph ---
compute_pseudotime <- function(g, root_node) {
        #g <- g_pruned
        #root_node <- "S1.S4.S2.S3"
        # Ensure weights reflect distances (if higher weight = closer, invert)
        w <- E(g)$weight
        if (max(w) > 1 || min(w) < 0) {
                w <- 1 / (w + 1e-6)   # avoid division by zero
        }
        
        # Compute shortest paths from root_node
        sp <- igraph::distances(g, 
                                v = root_node, 
                                to = V(g), 
                                weights = w, 
                                mode = "out") # follow directionality
        
        # Extract pseudotime vector
        pseudo <- as.numeric(sp[1, ])
        names(pseudo) <- names(V(g))
        
        # Normalize to 0–1
        pseudo <- (pseudo - min(pseudo, na.rm = TRUE)) / 
                (max(pseudo, na.rm = TRUE) - min(pseudo, na.rm = TRUE))
        
        return(pseudo)
}

# Example usage
root <- names(V(g_pruned))[which.min(V(g_pruned)$age)]
pseudotime <- compute_pseudotime(g_pruned, root)

pseudotime


A <- as_adjacency_matrix(g_pruned, attr = "weight", sparse = TRUE) # weighted

# Normalize rows to sum to 1 to get transition probabilities
P <- Diagonal(x = 1 / rowSums(A)) %*% A
rownames(P) <- colnames(P)

# Number of steps
t <- 1

# Compute t-step diffusion operator
P_t <- as.matrix(P) %^% t   # sparse matrix power (Matrix package)
#P_t <- as(P, "dgeMatrix") %^% t


library(RSpectra) # efficient for sparse matrices

k <- 2  # number of diffusion components

# Compute DCs with SVD
svd_res <- svds(P_t, k = k)
# Take U (diffusion components of each centroid--> Similarity of coinnectivity
# to the rest otf the graph over multiple steps)
DCs <- svd_res$u
rownames(DCs) <- rownames(P_t)
colnames(DCs) <- paste0("DC", 1:k)


plot(DCs[, 1], DCs[, 2])
text(DCs[, 1], DCs[, 2], labels = rownames(DCs), pos = 3, cex = 0.7)




library(igraph)
library(Matrix)

# Map cells onto graph edges and interpolate diffusion coordinates
map_cells_to_graph <- function(g, DCs, cell_mat, centroid_expr) {
        # Output matrix for cell-level diffusion coordinates
        g <- g_pruned
        DCs <- DCs
        centroid_expr <- X
        cell_mat <- X_cells
        
        DCs_cells <- matrix(NA, nrow = ncol(cell_mat), ncol = ncol(DCs),
                            dimnames = list(colnames(cell_mat),
                                            colnames(DCs)))
        
        edges <- as_edgelist(g, names = TRUE)
        
        # normalize expression for cosine similarity
        normalize <- function(mat){
                magns <- sqrt(colSums(mat^2))
                mat <- sweep(mat, 2, magns, FUN = "/")
                return(mat)
        }
        #normalize <- function(mat) mat / sqrt(colSums(mat^2))
        centroid_norm <- normalize(centroid_expr)
        cell_norm <- normalize(cell_mat)
        
        edge_vec <- c()
        
        for(cell_name in colnames(cell_mat)) {
                cell_vec <- cell_norm[, cell_name]
                
                best_edge <- NULL
                best_score <- -Inf
                w_A <- w_B <- NA
                
                # Find the edge (centroid pair) that best fits this cell
                for(e in 1:nrow(edges)) {
                        A <- edges[e,1]
                        B <- edges[e,2]
                        sim_A <- sum(cell_vec * centroid_norm[,A])
                        sim_B <- sum(cell_vec * centroid_norm[,B])
                        
                        total_sim <- sim_A + sim_B
                        if(total_sim == 0) total_sim <- 1e-6
                        wA <- sim_A / total_sim
                        wB <- sim_B / total_sim
                        
                        score <- wA*sim_A + wB*sim_B
                        if(score > best_score) {
                                best_score <- score
                                best_edge <- c(A,B)
                                w_A <- wA
                                w_B <- wB
                        }
                }
                edge_vec <- c(edge_vec,
                              paste(best_edge, collapse = "_"))
                # Interpolate the cell's DC along that edge
                DCs_cells[cell_name, ] <- w_A*DCs[best_edge[1], ] + w_B*DCs[best_edge[2], ]
        }
        
        cell_edge_df <- data.frame(cell = colnames(cell_mat),
                                   edge = edge_vec)
        
        return(DCs_cells)
}


map_cells_soft <- function(DCs, X, centroid_expr, metric = c("cosine", "euclidean")) {
        # DCs: diffusion coordinates of centroids (rows = centroids)
        # X: expression matrix of cells (rows = genes, cols = cells)
        # centroid_expr: expression matrix of centroids (rows = genes, cols = centroids)
        # metric: similarity metric ("cosine" or "euclidean")
        X <- X_cells
        g <- g_pruned
        DCs <- DCs
        centroid_expr <- X
        cell_mat <- X_cells
        
        metric <- match.arg(metric)
        
        # Output: diffusion coordinates for cells
        DCs_cells <- matrix(NA, nrow = ncol(X), ncol = ncol(DCs))
        rownames(DCs_cells) <- colnames(X)
        colnames(DCs_cells) <- colnames(DCs)
        
        # Optionally normalize for cosine similarity
        if(metric == "cosine") {
                normalize <- function(mat) mat / sqrt(colSums(mat^2))
                centroid_norm <- normalize(centroid_expr)
                cell_norm <- normalize(X)
        } else {
                centroid_norm <- centroid_expr
                cell_norm <- X
        }
        
        for(cell_name in colnames(X)) {
                cell_vec <- cell_norm[, cell_name]
                
                if(metric == "cosine") {
                        sim <- colSums(cell_vec * centroid_norm)  # cosine similarity
                } else {
                        # Euclidean distance → convert to similarity
                        dists <- sqrt(colSums((centroid_norm - cell_vec)^2))
                        sim <- exp(-dists^2 / 2)  # Gaussian kernel
                }
                
                # Normalize similarities to sum to 1
                sim[sim < 0] <- 0  # optional: ignore negative cosine
                weights <- sim / sum(sim)
                
                # Weighted average of centroid DCs
                DCs_cells[cell_name, ] <- t(weights) %*% DCs
        }
        
        return(DCs_cells)
}



# DCs_centroids = diffusion components of your centroids
# X_cells = expression matrix of individual cells
# centroid_expr = expression matrix of centroids used to build g_pruned

DCs_cells <- map_cells_to_graph(
        g_pruned = g_pruned,
        DCs = DCs_centroids,
        X = X_cells,
        centroid_expr = centroid_expr
)

# Compute pseudotime from root
root <- "centroid_1"  # earliest-age centroid
pseudotime_cells <- DCs_cells[, 1] - DCs_cells[root, 1]
pseudotime_cells <- pseudotime_cells - min(pseudotime_cells)
