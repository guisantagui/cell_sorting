# KNN
library(ggplot2)
library(RBGL)
library(parallel)

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
                theme_minimal()
        return(plt)
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

min_avg_weight_path_dp <- function(g, source, target,
                                   donor_info
                                   ) {
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
                        #if (donor_info$age[v] < donor_info$age[u]) next
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
        names(path) <- path
        #return(V(g)[path])
        return(path)
}

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

get_all_max_weight_g <- function(g, parallelize = T)
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

get_all_min_avg_weight_g <- function(g, donor_info, parallelize = F,
                                     mode = "shortest_path"){
        #g <- g_base
        #donor_info <- samp_info
        #mode <- "min_avg_weight_path"
        #parallelize <- T
        pairs <- t(combn(V(g)$name,2))
        pairs <- rbind(pairs,pairs[,c(2,1)])
        pairs <- data.frame(pairs)
        
        
        
        if (mode == "shortest_path"){
                func <- function(g, source, target){
                        out <- names(shortest_paths(g,
                                                    source, target,
                                                    mode = "out",
                                                    output = "both")$vpath[[1]])
                        names(out) <- out
                        return(out)
                }
        }else if (mode == "min_avg_weight_path"){
                func <- function(g, source, target){
                        out <- min_avg_weight_path_dp(g, source, target, donor_info = donor_info)
                        return(out)
                }
                #func <- min_avg_weight_path_dp
        }else if (mode == "min_avg_weight_path_slow"){
                func <- function(s, source, target){
                        out <- names(min_avg_weight_path(g, source, target))
                        names(out) <- out
                        return(out)
                }
        }
        
        # Keep only pairs where age of the first element is less or equal than age of
        # second.
        pairs$age_1 <- donor_info$age[match(pairs[, 1], donor_info$donor)]
        pairs$age_2 <- donor_info$age[match(pairs[, 2], donor_info$donor)]
        pairs <- as.matrix(pairs[pairs$age_1 <= pairs$age_2, c(1, 2)])
        
        
        pairs <- split(pairs,
                       1:nrow(pairs))
        
        if (!parallelize){
                pb <- txtProgressBar(min = 0, max = length(pairs), style = 3)
                all_paths <- list()
                for (i in seq_along(pairs)) {
                        pair <- pairs[[i]]
                        setTxtProgressBar(pb = pb, value = i)
                        path <- func(g, pair[1], pair[2])
                        if (!is.null(path)) {
                                all_paths <- append(all_paths, list(path))
                        }
                }
                close(pb)
        }else{
                all_paths <- mclapply(pairs, function(pair)
                {
                        func(g, pair[1], pair[2])
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

get_opt_branching <- function(g)
{
        #g <- sameage_g
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


plot_graph <- function(g, node_info, color){
        
        #g <- g_pruned_dp
        #node_info <- samp_info
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
        plot(g,
             layout=layout.fruchterman.reingold,
             vertex.color=V(g)$color,
             edge.arrow.size=.5,
             vertex.size = 8)
}

jitter_age <- function(node_info, seed = 1){
        set.seed(1)
        node_info$age <- node_info$age + rnorm(length(node_info$age),
                                               mean = 0,
                                               sd = 0.000001)
        return(node_info)
}

# This function computes KNN and builds a directed graph based on it, with the
# distances (euclidean) between the points as weights.
build_knn_graph <- function(X,
                            k = 6,
                            samp_info,
                            standardize = T#,
                            #jitterize = T
                            ){
        #X <- X
        #samp_info <- collapsed_sampinfo
        #standardize <- T
        if (standardize){
                X <- plotUtils::stand(t(X))
        }else{
                X <- t(X)
        }
        #if (jitterize){
        #        set.seed(1)
        #        samp_info$age <- samp_info$age + rnorm(length(samp_info$age),
        #                                               mean = 0,
        #                                               sd = 0.000001)
        #}
        
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
        V(g)
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

# Remove cycles caused by individuals of the same age by finding the shortest
# path that connects all of them using MST.
remove_same_age_cycles <- function(g, samp_info){
        #g <- g
        #samp_info <- samp_info
        dup_ages <- unique(samp_info$age[duplicated(samp_info$age)])
        all_edges <- as_edgelist(g)
        for(age in dup_ages){
                #age <- 1
                sameage_g <- subgraph(g, samp_info$donor[samp_info$age == age])
                # Remove edges of the subgraph
                all_edges <- all_edges[!(all_edges[, 1] %in% names(V(sameage_g)) & all_edges[, 2] %in% names(V(sameage_g))), ]
                sameage_g <- mst(sameage_g)
                # Add pruned edges
                all_edges <- rbind.data.frame(all_edges, as_edgelist(sameage_g)) 
        }
        edge_ids <- paste(all_edges[, 1],
                          all_edges[, 2],
                          sep = "|")
        g_nocycs <- subgraph.edges(g, edge_ids)
        return(g_nocycs)
}

library(igraph)

orient_same_age <- function(g){
        g <- g_base
        g_new <- g
        ages <- sort(unique(V(g)$age))
        name_to_idx <- function(names_vec) { match(names_vec, V(g)$name) }
        for (a in ages){
                a <- 1
                nodes_a <- V(g)$name[V(g)$age == a]
                if (length(nodes_a) <= 1) next
                
                # induced subgraph for these nodes (undirected for MST)
                subg <- induced_subgraph(g, vids = nodes_a)
                subg_undir <- as.undirected(subg, mode = "collapse", edge.attr.comb = "mean")
                subg_undir <- mst(subg_undir)
                mst_edges <- as_edgelist(subg_undir)
                mst_edges <- rbind(mst_edges,
                                   mst_edges[, c(2, 1)])
                # If subg has no edges (isolated points), nothing to do
                if (ecount(subg_undir) == 0) {
                        next
                }
                neighs <- neighborhood(g,
                                       order = 1,
                                       nodes = V(g)[nodes_a],
                                       mode = "all")
                subg <- induced_subgraph(g, vids = unlist(lapply(neighs, names)))
                
                min_age <- min(V(subg)$age)
                max_age <- max(V(subg)$age)
                
                min_vs <- names(V(subg))[V(subg)$age == min_age]
                max_vs <- names(V(subg))[V(subg)$age == max_age]
                
                subg_el <- as_edgelist(subg)
                
                # Remove edges that connect sink age nodes to each other
                subg_el <- subg_el[!(subg_el[, 1] %in% max_vs & subg_el[, 2] %in% max_vs), ]
                
                
                # If min_age is different than a (meaning that we are not in
                # the root), remove edges that connect younger age nodes to
                # each other
                if (min_age != a){
                        subg_el <- subg_el[!(subg_el[, 1] %in% min_vs & subg_el[, 2] %in% min_vs), ]
                }
                
                # Keep only edges from MST
                subg_el <- subg_el[!(subg_el[, 1] %in% nodes_a & subg_el[, 2] %in% nodes_a), ]
                subg_el <- rbind(subg_el, mst_edges)
                
                # Prune the subgraph
                subg <- subgraph.edges(subg, E(subg, P = t(subg_el)))
                
                pairs <- list()
                for (i in seq_along(min_vs)){
                        min_v <- min_vs[i]
                        for(j in seq_along(max_vs)){
                                max_v <- max_vs[j]
                                pairs[[(i - 1) * length(min_vs) + j]] <- c(min_v, max_v)
                        }
                }
                
                all_paths <- list()
                for(i in seq_along(pairs)){
                        pair <- pairs[[i]]
                        min_path <- names(min_avg_weight_path(subg,
                                                              from = pair[1],
                                                              to = pair[2]))
                        all_paths[[i]] <- min_path
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
                subg_pruned <- subgraph.edges(subg, E(subg, P = t(edge_list)))
                plot(subg_pruned)
                
                subg_el <- as_edgelist(subg)
                subg_el <- subg_el[!(subg_el[, 1] %in% nodes_a & subg_el[, 2] %in% nodes_a), ]
                subg_el <- rbind(mst_edges,
                                 subg_el)
                #edge_idx <- get.edge.ids(g, paste(subg_el[,1], subg_el[,2], sep="|"), directed=TRUE)
                edge_ids <- paste(subg_el[, 1],
                                  subg_el[, 2],
                                  sep = "|")
                subg <- subgraph.edges(g, edge_ids, delete.vertices = T)
                
                min_age <- 
                E(g)$weight
        }
}


orient_same_age_msts <- function(g) {
        # copy to modify
        #g <- g_base
        g_new <- g
        
        ages <- sort(unique(V(g)$age))
        name_to_idx <- function(names_vec) { match(names_vec, V(g)$name) } # helper
        
        for (a in ages) {
                # nodes in this age
                a <- 1
                nodes_a <- V(g)$name[V(g)$age == a]
                if (length(nodes_a) <= 1) next
                
                # induced subgraph for these nodes (undirected for MST)
                subg <- induced_subgraph(g, vids = nodes_a)
                subg_undir <- as.undirected(subg, mode = "collapse", edge.attr.comb = "mean")
                
                # If subg has no edges (isolated points), nothing to do
                if (ecount(subg_undir) == 0) {
                        next
                }
                
                # compute MST on this undirected subgraph (weights should be distances)
                # igraph::mst returns an undirected tree with length = vcount - 1
                mst_sub <- mst(subg_undir, weights = E(subg_undir)$weight)
                
                # Find nodes in this age that connect to younger / older nodes in the original graph g
                # We check edges *between* S and outside S:
                outside_edges <- E(g)[.from(nodes_a) & !.to(nodes_a)] # but inc() doesn't work this way
                # Simpler: find neighbors in g for each node
                neighs <- neighborhood(g, order = 1, nodes = V(g)[nodes_a], mode = "all") # list per node
                names(neighs) <- nodes_a
                younger_connected <- c()
                older_connected <- c()
                for (vname in nodes_a) {
                        nb <- setdiff(V(g)$name[unlist(neighs[[vname]])], vname) # neighbor names
                        if (length(nb) == 0) next
                        ages_nb <- V(g)$age[match(nb, V(g)$name)]
                        if (any(ages_nb < a, na.rm=TRUE)) younger_connected <- c(younger_connected, vname)
                        if (any(ages_nb > a, na.rm=TRUE)) older_connected <- c(older_connected, vname)
                }
                younger_connected <- unique(younger_connected)
                older_connected <- unique(older_connected)
                
                # compute distances on the MST from the Y set (multi-source)
                # need node indices relative to mst_sub
                mst_names <- V(mst_sub)$name
                if (length(younger_connected) > 0) {
                        src_idx <- which(mst_names %in% younger_connected)
                        dist_to_y <- distances(mst_sub, v = src_idx, to = V(mst_sub), weights = E(mst_sub)$weight)
                        # distances returns matrix with rows = sources; take min across sources
                        min_dist_to_y <- apply(dist_to_y, 2, min)
                        # orient edges from smaller min_dist -> larger min_dist
                        el <- as_edgelist(mst_sub, names = TRUE)
                        # prepare directed edges
                        oriented_edges <- matrix(ncol=2, nrow=0)
                        oriented_weights <- c()
                        for (i in seq_len(nrow(el))) {
                                u <- el[i,1]; v <- el[i,2]
                                du <- min_dist_to_y[which(mst_names == u)]
                                dv <- min_dist_to_y[which(mst_names == v)]
                                if (is.infinite(du) & is.infinite(dv)) {
                                        # no info, keep arbitrary orientation u->v
                                        oriented_edges <- rbind(oriented_edges, c(u,v))
                                } else if (du <= dv) {
                                        oriented_edges <- rbind(oriented_edges, c(u,v))
                                } else {
                                        oriented_edges <- rbind(oriented_edges, c(v,u))
                                }
                                oriented_weights <- c(oriented_weights, E(mst_sub)$weight[i])
                        }
                } else if (length(older_connected) > 0) {
                        # no younger connected nodes, orient toward older_connected
                        sink_idx <- which(mst_names %in% older_connected)
                        dist_to_o <- distances(mst_sub, v = V(mst_sub),
                                               to = sink_idx, weights = E(mst_sub)$weight)
                        min_dist_to_o <- apply(dist_to_o, 1, min)
                        el <- as_edgelist(mst_sub, names = TRUE)
                        oriented_edges <- matrix(ncol=2, nrow=0)
                        oriented_weights <- c()
                        for (i in seq_len(nrow(el))) {
                                u <- el[i,1]; v <- el[i,2]
                                du <- min_dist_to_o[which(mst_names == u)]
                                dv <- min_dist_to_o[which(mst_names == v)]
                                if (is.infinite(du) & is.infinite(dv)) {
                                        oriented_edges <- rbind(oriented_edges, c(u,v))
                                } else if (du >= dv) {
                                        # orient from larger -> smaller to point toward older sinks
                                        oriented_edges <- rbind(oriented_edges, c(u,v))
                                } else {
                                        oriented_edges <- rbind(oriented_edges, c(v,u))
                                }
                                oriented_weights <- c(oriented_weights, E(mst_sub)$weight[i])
                        }
                } else {
                        # neither younger nor older connected: keep tree as undirected edges; choose arbitrary orientation
                        el <- as_edgelist(mst_sub, names = TRUE)
                        oriented_edges <- el
                        oriented_weights <- E(mst_sub)$weight
                }
                
                # Now: remove all intra-age edges from g_new, then add the oriented MST edges with weights
                # Identify intra-age edges in original g by endpoints: both endpoints in nodes_a
                intra_edges_idx <- which(sapply(E(g_new), function(ei) {
                        ends(g_new, ei, names=TRUE)[1] %in% nodes_a && ends(g_new, ei, names=TRUE)[2] %in% nodes_a
                }))
                if (length(intra_edges_idx) > 0) g_new <- delete_edges(g_new, E(g_new)[intra_edges_idx])
                
                # Add oriented MST edges (directed)
                for (i in seq_len(nrow(oriented_edges))) {
                        u <- oriented_edges[i,1]; v <- oriented_edges[i,2]
                        w <- oriented_weights[i]
                        # avoid adding duplicate edges
                        if (!are.connected(g_new, u, v)) {
                                g_new <- add_edges(g_new, c(u, v))
                                # set weight on the last added edge
                                E(g_new)$weight[length(E(g_new))] <- w
                        } else {
                                # if edge exists, you may want to keep minimum weight
                                eid <- get.edge.ids(g_new, paste(u, v, sep="|"), directed=TRUE)
                                # get.edge.ids uses numeric ids, but simpler:
                                # find edge by endpoints
                                eids <- E(g_new)[from(V(g_new)[u]) & to(V(g_new)[v])]
                                if (length(eids) > 0) {
                                        E(g_new)[eids]$weight <- pmin(E(g_new)[eids]$weight, w)
                                } else {
                                        E(g_new)$weight[length(E(g_new))] <- w
                                }
                        }
                }
                
                # Finally, remove any edges that go from older -> younger across ages (optional)
                # This step can be done after processing all ages; we'll keep it per-age for safety
                inter_edges_idx <- which(sapply(E(g_new), function(ei) {
                        ep <- ends(g_new, ei, names=TRUE)
                        a1 <- V(g_new)$age[match(ep[1], V(g_new)$name)]
                        a2 <- V(g_new)$age[match(ep[2], V(g_new)$name)]
                        !is.na(a1) && !is.na(a2) && (a1 > a2)
                }))
                if (length(inter_edges_idx) > 0) g_new <- delete_edges(g_new, E(g_new)[inter_edges_idx])
        } # end per age loop
        
        return(g_new)
}



# Create toy dataset
toy_samps <- data.frame(donor = paste0("S", 1:length(rep(seq(1, 100, by = 10), each = 4))),
                        age = rep(seq(1, 100, by = 10), each = 4))
rownames(toy_samps) <- toy_samps$donor
toy_centroids <- matrix(data = NA, nrow = 2, ncol = nrow(toy_samps),
                        dimnames = list(paste0("G", 1:2),
                                        toy_samps$donor))

toy_centroids[1, ] <- gene_func1(toy_samps$age, slope = 2, intercept = 5, sd = 1)
toy_centroids[2, ] <- gene_func1(toy_samps$age, slope = 3, intercept = 10, sd = 1)

activtline_3 <- toy_centroids

activtline_3[2, ] <- rnorm(length(activtline_3[2, ]), sd = 1)

activtline_3[2, 1:ncol(activtline_3) %% 2 == 0] <- gene_func3(toy_samps$age[1:ncol(activtline_3) %% 2 == 0],
                                                              slope = 2, age_start = 40, sd = 1)

# KNN workflow
################################################################################
plot_centroids(activtline_3, 1, 2)


k <- 7

X <- activtline_3
#X <- plotUtils::stand(t(X))

samp_info <- toy_samps

#samp_info <- jitter_age(samp_info)
rownames(samp_info) <- samp_info$donor

g_base <- build_knn_graph(X = activtline_3,
                          samp_info = samp_info,
                          k = k,
                          standardize = T)

plot_graph(g_base, node_info = samp_info, color = "age")

g_pruned_dp <- get_all_min_avg_weight_g(g_base, samp_info, parallelize = T,
                                        mode = "min_avg_weight_path")
plot_graph(g_pruned_dp, node_info = samp_info, color = "age")

plot_graph(get_opt_branching(g_pruned_dp), node_info = samp_info, color = "age")


close_sameage_samps <- get_close_sameage_samps(g_base, resolution = 1, alg = "louvain")
collapsed_sampinfo <- build_collapsed_sample_info(samp_info, close_samps = close_sameage_samps)

X <- get_mean_close_samps(X, close_sameage_samps)

g <- build_knn_graph(X = X,
                     samp_info = collapsed_sampinfo, k = 6,
                     standardize = T)

#plot_graph(get_opt_branching(g), node_info = collapsed_sampinfo, color = "age")

plot_graph(g,
           node_info = collapsed_sampinfo, color = "age")

g_pruned <- get_all_max_avg_weight_g(g) 

g_pruned <- get_all_min_avg_weight_g(g,
                                     collapsed_sampinfo,
                                     parallelize = F,
                                     mode = "min_avg_weight_path") 

plot_graph(g_pruned,
           node_info = collapsed_sampinfo, color = "age")

plot_centroids(activtline_3, x = 1, y = 2)