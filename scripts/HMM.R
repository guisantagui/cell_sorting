#HMM

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


plot_centroids <- function(centroid_mat, x, y){
        #x <- 1
        #y <- 2
        centroid_mat <- as.matrix(centroid_mat)
        plot_df <- data.frame(x = centroid_mat[x, ],
                              y = centroid_mat[y, ],
                              labels = colnames(centroid_mat))
        plt <- ggplot(plot_df,
                      mapping = aes(x = x, y = y, label = labels)) +
                geom_text() +
                theme_minimal()
        return(plt)
}

# Functions for building graphs
create_base_graph <- function(donor_info, jitter = T, centroid_mat,
                              distance = "euclidean"){
        donors <- donor_info$donor
        ages_vec <- donor_info$age
        names(ages_vec) <- donors
        
        
        # Compute distance matrix
        if (distance == "euclidean"){
                dists <- as.matrix(dist(t(centroid_mat)))
        }else if (distance == "wasserstein"){
                dists <- get_wasserstein_dist(centroid_mat)
        }else if (distance == "maxnorm"){
                dists <- get_maxnorm_dist(centroid_mat)
        }
        
        
        # Add a bit of noise to the age to break potential cycles
        if (jitter){
                set.seed(1)
                ages_vec <- ages_vec + rnorm(length(ages_vec),
                                             mean = 0,
                                             sd = 0.000001)
        }
        g <- make_empty_graph(n = length(donors), directed = TRUE)
        for (i in 1:length(donors)) {
                for (j in 1:length(donors)) {
                        if (ages_vec[donors[j]] >= ages_vec[donors[i]] & i != j) {
                                g <- add_edges(g, c(i, j),
                                               attr = list(weight = dists[donors[i], donors[j]]))
                        }
                }
        }
        V(g)$name <- donors
        return(g)
}

min_avg_weight_path_dp <- function(g, source, target, donor_info) {
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
                        if (donor_info$age[v] < donor_info$age[u]) next
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

get_all_min_avg_weight_g <- function(g, donor_info, parallelize = F,
                                     mode = "shortest_path"){
        #g <- g_nocycs
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
                func <- function(g, souirce, target){
                        out <- min_avg_weight_path_dp(g, source, target, donor_info = donor_info)
                }
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
        #g <- g_all_const
        #node_info <- collapsed_sampinfo
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

get_transmat <- function(fitmod, age){
        #age <- 21
        n_states <- length(fitmod@transition)  # e.g., 5
        T_m <- matrix(0, n_states, n_states)
        for(s in 1:n_states){
                coefs_matrix <- fitmod@transition[[s]]@parameters$coefficients
                #age_c <- 0  # or set to the sample’s age if you want age-dependent transitions
                
                # Compute the logits for all next states (baseline = first state)
                logits <- rep(0, n_states)  # baseline state = 0
                for(t in 2:n_states){       # coefs for states 2..K
                        logits[t] <- coefs_matrix["(Intercept)", t] + coefs_matrix["ages", t] * age
                }
                
                # Convert to probabilities via softmax
                T_m[s, ] <- exp(logits - max(logits)) / sum(exp(logits - max(logits)))
        }
        return(T_m)
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

#activtline_3 <- apply(activtline_3, 1, function(x) (x - mean(x))/sd(x))

plot_centroids(activtline_3, x = 1, y = 2)

library(class)
#knn(train = t(stand(t(activtline_3))),
#    test =  t(stand(t(activtline_3))))

stand <- function (m, axis = 2, scale = T, center = T) 
{
        m <- activtline_3
        m_class <- class(m)[1]
        if ((!m_class %in% c("data.frame", "matrix")) | any(dim(m) == 
                                                            0)) {
                stop("Invalid input object.", call. = F)
        }
        if (!axis %in% c(1, 2)) {
                stop("Invalid axis.", call. = F)
        }
        if (!scale & !center) {
                warning("No processing will be done", call. = F)
        }
        na_mask <- is.na(m)
        if (axis == 1) {
                m <- t(m)
                na_mask <- t(na_mask)
        }
        if (any(na_mask)) {
                warning("The matrix contains NAs. They will be excluded during standardization.", 
                        call. = F)
                not_valid_entries <- colSums(!na_mask) <= 1
                if (any(not_valid_entries)) {
                        warning(sprintf("%s %s with one or less not NAs observations. Dropped from matrix.", 
                                        sum(not_valid_entries), c("rows", "columns")[axis]))
                }
                m <- m[, !not_valid_entries]
                na_mask <- na_mask[, !not_valid_entries]
        }
        format_m_mask <- function(m, na_mask, fun) {
                m_list <- m_lst <- apply(m, 2, fun)
                m <- matrix(NA, nrow = nrow(na_mask), ncol = ncol(na_mask), 
                            dimnames = dimnames(na_mask))
                m[!na_mask] <- unlist(m_list)
                return(m)
        }
        if (center) {
                m <- format_m_mask(m, na_mask, function(x) x[!is.na(x)] - 
                                           mean(x[!is.na(x)]))
        }
        if (scale) {
                keep <- apply(m, 2, function(x) sd(x[!is.na(x)])) > 0
                if (any(!keep)) {
                        warning(sprintf("%s %s with SD = 0. Dropped from matrix.", 
                                        sum(keep), c("rows", "columns")[axis]))
                }
                m <- m[, keep]
                na_mask <- na_mask[, keep]
                m <- format_m_mask(m, na_mask, function(x) x[!is.na(x)]/sd(x[!is.na(x)]))
        }
        if (axis == 1) {
                m <- t(m)
        }
        if (m_class == "data.frame") {
                m <- as.data.frame(m)
        }
        return(m)
}

plot_centroids(activtline_3, x = 1, y = 2)

if(!require("depmixS4", quietly = T)) install.packages("depmixS4")
library(depmixS4)

X <- activtline_3
ages <- toy_samps$age

data_hmm <- data.frame(sample = colnames(X),
                       ages = toy_samps$age,
                       G1 = activtline_3[1, ],
                       G2 = activtline_3[2, ])

nstates <- 5

mod <- depmix(list(G1 ~ 1, G2 ~ 1), 
              data = data_hmm, 
              nstates = nstates, 
              family = list(gaussian(), gaussian()), 
              transition = ~ages)

set.seed(123)
fitmod <- fit(mod)
getpars(fitmod)
summary(fitmod)

posterior(fitmod, type = "viterbi")


getpars(fitmod)


# Compute emission P(Xi | Z = s, agei)

ngenes <- nrow(X)
nsamples <- ncol(X)
emission_loglik <- matrix(0, nsamples, nstates)

for(s in 1:nstates){
        mu <- sapply(1:ngenes, function(g) fitmod@response[[s]][[g]]@parameters$coefficients)
        sd <- sapply(1:ngenes, function(g) fitmod@response[[s]][[g]]@parameters$sd)
        for(c in 1:nsamples){
                emission_loglik[c,s] <- sum(dnorm(X[,c], mean=mu, sd=sd, log=TRUE))
        }
}


ngenes <- nrow(X)
nsamples <- ncol(X)
emission_lik <- matrix(0, nsamples, nstates)

for(s in 1:nstates){
        mu <- sapply(1:ngenes, function(g) fitmod@response[[s]][[g]]@parameters$coefficients)
        sd <- sapply(1:ngenes, function(g) fitmod@response[[s]][[g]]@parameters$sd)
        for(c in 1:nsamples){
                emission_lik[c,s] <- prod(dnorm(X[,c], mean=mu, sd=sd, log=F))
        }
}

# Convert to probabilities with priors
state_prior <- rep(1/nstates, nstates) # or your HMM stationary distribution
emission_prob <- t(apply(emission_loglik, 1, function(loglik) {
        lik <- exp(loglik - max(loglik))  # avoid overflow
        prob <- lik * state_prior
        prob / sum(prob)
}))

rownames(emission_prob) <- colnames(X)


dist(emission_prob)

#create_base_graph(toy_samps, jitter = T, emission_prob)

g <- create_base_graph(toy_samps, jitter = T, t(emission_prob))

g <- create_base_graph(toy_samps, jitter = T, t(scale(t(activtline_3))), distance = "euclidean")

g <- get_all_min_avg_weight_g(g,
                              donor_info = toy_samps,
                              parallelize = T,
                              mode = "min_avg_weight_path")
plot_graph(get_opt_branching(g), toy_samps, color = "age")

plot(g)
plot(get_opt_branching(g))



nn_list <- list()
 
for(i in seq_along(V(g))) {
        v <- V(g)[i]
        neighbors <- neighbors(g, v, mode = "out") # outgoing edges
        edge_weights <- E(g)[from(v) & to(neighbors)]$weight
        nn <- neighbors[order(edge_weights)][1:k] # pick k closest by weight
        nn_list[[v$name]] <- names(nn)
}



k <- 6

X <- activtline_3
#X <- plotUtils::stand(t(X))

samp_info <- toy_samps
rownames(samp_info) <- samp_info$donor

build_knn_graph <- function(X, k = 6, samp_info, standardize = T){
        #X <- X
        #samp_info <- collapsed_sampinfo
        #standardize <- T
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
        edges$w <- as.vector(nn$nn.dist)
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

g_base <- build_knn_graph(X = activtline_3,
                          samp_info = samp_info,
                          k = 6)

plot_graph(g_base, node_info = samp_info, color = "age")

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

close_sameage_samps <- get_close_sameage_samps(g_base, resolution = 1, alg = "leiden")

collapsed_sampinfo <- build_collapsed_sample_info(samp_info, close_samps = close_sameage_samps)

X

# This function takes an exprtession matrix (genes in columns, samples in rows)
# and a vector of samples that are close in terms of KNN distance, and collapses
# similar samples into a single point by computing the average. Intended to be
# used to aggregate samples that are the same age and similar in terms of 
# transcriptome into a single point, to avoid cycles in the graph.

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

X <- get_mean_close_samps(X, close_sameage_samps)

g <- build_knn_graph(X = X,
                     samp_info = collapsed_sampinfo, k = 7)

#plot_graph(get_opt_branching(g), node_info = collapsed_sampinfo, color = "age")

plot_graph(g,
           node_info = collapsed_sampinfo, color = "age")

plot_centroids(t(X), 1, 2)

plot_graph(get_all_min_avg_weight_g(g,
                                    donor_info = collapsed_sampinfo,
                                    mode = "min_avg_weight_path"), node_info = collapsed_sampinfo, color = "age")

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


g_nocycs <- remove_same_age_cycles(g, samp_info = samp_info)

g_nocycs_pruned <- get_all_min_avg_weight_g(g_nocycs,
                                            donor_info = samp_info,
                                            parallelize = T,
                                            mode = "min_avg_weight_path")

plot_graph(g_nocycs_pruned,
           node_info = samp_info,
           color = "age")

plot(g)


