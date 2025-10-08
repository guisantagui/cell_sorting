if (!require(parallel, quietly = T)) install.packages("parallel")
if (!require(igraph, quietly = T)) install.packages("igraph")

toy_samps <- data.frame(donor = paste0("S", 1:length(rep(seq(1, 60, by = 10), each = 3))),
                        age = rep(seq(1, 60, by = 10), each = 3))
rownames(toy_samps) <- toy_samps$donor
toy_centroids <- matrix(data = NA, nrow = 2, ncol = nrow(toy_samps),
                        dimnames = list(paste0("G", 1:2),
                                        toy_samps$donor))


# Obtain Wasserstein distance matrix
get_wasserstein_dist <- function(centroid_mat){
        dists <- matrix(0,
                        nrow = ncol(centroid_mat),
                        ncol = ncol(centroid_mat),
                        dimnames = list(colnames(centroid_mat),
                                        colnames(centroid_mat)))
        for(i in 1:nrow(dists)){
                for(j in 1:ncol(dists)){
                        if(i == j){
                                dists[i,j] <- 0
                        }else{
                                dists[i,j] <- transport::wasserstein1d(centroid_mat[,i],
                                                                       centroid_mat[,j])
                        }
                }
        }
        return(dists)
}

get_maxnorm_dist <- function(centroid_mat){
        dists <- matrix(0,
                        nrow = ncol(centroid_mat),
                        ncol = ncol(centroid_mat),
                        dimnames = list(colnames(centroid_mat),
                                        colnames(centroid_mat)))
        for (i in 1:nrow(dists)){
                for (j in 1:ncol(dists)){
                        if (i == j){
                                dists[i, j] <- 0
                        }else{
                                dists[i, j] <- max(abs(centroid_mat[,i] - centroid_mat[,j]))
                        }
                }
        }
        return(dists)
}

# Create base centroid graph that 
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

# Find minimum average weight path based on dynamic programming.
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
        names(path) <- path
        #return(V(g)[path])
        return(path)
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

# Obtain the subgraph that includes all minimum average weight paths.
get_all_min_avg_weight_g <- function(g, donor_info, parallelize = F,
                                     mode = "shortest_path"){
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
                func <- min_avg_weight_path_dp
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

# Obtain optimal branching with Edmonds algorithm
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


plot_graph <- function(g, node_info, color){
        #g <- g_all_const
        #node_info <- toy_samps
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
        #centroid_mat <- accelerate_3
        plot_df <- data.frame(x = centroid_mat[x, ],
                              y = centroid_mat[y, ],
                              labels = colnames(centroid_mat))
        plt <- ggplot(plot_df,
                      mapping = aes(x = x, y = y, label = labels)) +
                geom_text() +
                theme_minimal()
        return(plt)
}

toy_centroids[1, ] <- gene_func1(toy_samps$age, slope = 2, intercept = 5, sd = 1)
toy_centroids[2, ] <- gene_func1(toy_samps$age, slope = 3, intercept = 10, sd = 1)
#toy_centroids[3, ] <- gene_func1(toy_samps$age, slope = 1, intercept = 10, sd = 10)
#toy_centroids[4, ] <- gene_func1(toy_samps$age, slope = 2, intercept = 10, sd = 10)
#toy_centroids[5, ] <- gene_func1(toy_samps$age, slope = 4.5, intercept = 10, sd = 10)

all_constant <- toy_centroids
accelerate_3 <- toy_centroids
activtline_3 <- toy_centroids



accelerate_3[1, 1:ncol(accelerate_3) %% 3 == 0] <- gene_func2(toy_samps$age[1:ncol(accelerate_3) %% 3 == 0],
                                                              slope = 2, intercept = 5, sd = 1,
                                                              thrshld = 20, slopemult = 8)
accelerate_3[2, 1:ncol(accelerate_3) %% 3 == 0] <- gene_func2(toy_samps$age[1:ncol(accelerate_3) %% 3 == 0],
                                                              slope = 3, intercept = 10, sd = 1,
                                                              thrshld = 20, slopemult = 8)
#accelerate_3[3, 1:ncol(accelerate_3) %% 3 == 0] <- gene_func2(toy_samps$age[1:ncol(accelerate_3) %% 3 == 0],
#                                                              slope = 1, intercept = 10, sd = 10,
#                                                              thrshld = 40, slopemult = 5)
#accelerate_3[4, 1:ncol(accelerate_3) %% 3 == 0] <- gene_func2(toy_samps$age[1:ncol(accelerate_3) %% 3 == 0],
#                                                              slope = 2, intercept = 10, sd = 10,
#                                                              thrshld = 40, slopemult = 5)
#accelerate_3[5, 1:ncol(accelerate_3) %% 3 == 0] <- gene_func2(toy_samps$age[1:ncol(accelerate_3) %% 3 == 0],
#                                                              slope = 4.5, intercept = 10, sd = 10,
#                                                              thrshld = 40, slopemult = 5)


activtline_3[2, ] <- rnorm(length(activtline_3[2, ]), sd = 1)

#activtline_3[4, ] <- rnorm(length(activtline_3[4, ]), sd = 10)
#activtline_3[5, ] <- rnorm(length(activtline_3[5, ]), sd = 10)

activtline_3[2, 1:ncol(activtline_3) %% 3 == 0] <- gene_func3(toy_samps$age[1:ncol(activtline_3) %% 3 == 0],
                                                              slope = 2, age_start = 40, sd = 1)
#activtline_3[5, 1:ncol(activtline_3) %% 3 == 0] <- gene_func3(toy_samps$age[1:ncol(activtline_3) %% 3 == 0],
#                                                              slope = 4.5, age_start = 50, sd = 10)

plot_centroids(all_constant, x = 1, y = 2)
plot_centroids(accelerate_3, x = 1, y = 2)
plot_centroids(activtline_3, x = 1, y = 2)


g_all_const <- create_base_graph(toy_samps, centroid_mat = all_constant, distance = "maxnorm")

g_acc_three <- create_base_graph(toy_samps, centroid_mat = accelerate_3, distance = "maxnorm")

g_actline_3 <- create_base_graph(toy_samps, centroid_mat = activtline_3, distance = "maxnorm")

is_dag(g_all_const)
is_dag(g_acc_three)
is_dag(g_actline_3)


plot(g_all_const)
plot(g_acc_three)
plot(g_actline_3)

g_all_const <- get_all_min_avg_weight_g(g_all_const,
                                        donor_info = toy_samps,
                                        parallelize = T,
                                        mode = "min_avg_weight_path")

g_acc_three <- get_all_min_avg_weight_g(g_acc_three,
                                        donor_info = toy_samps,
                                        parallelize = T,
                                        mode = "min_avg_weight_path")

g_actline_3 <- get_all_min_avg_weight_g(g_actline_3,
                                        donor_info = toy_samps,
                                        parallelize = T,
                                        mode = "min_avg_weight_path")

plot_graph(g_all_const, toy_samps, color = "age")

plot_graph(g_acc_three, toy_samps, color = "age")

plot_graph(g_actline_3, toy_samps, color = "age")


plot_graph(get_opt_branching(g_all_const), toy_samps, color = "age")

plot_graph(get_opt_branching(g_acc_three), toy_samps, color = "age")

plot_graph(get_opt_branching(g_actline_3), toy_samps, color = "age")




# Try HMM

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

plot_centroids(activtline_3, x = 1, y = 2)

if(!require("depmixS4", quietly = T)) install.packages("depmixS4")
library(depmixS4)

X <- activtline_3
ages <- toy_samps$age

data_hmm <- data.frame(sample = colnames(X),
                       ages = toy_samps$age,
                       G1 = activtline_3[1, ],
                       G2 = activtline_3[2, ])

# Toy example of transition matrix: 3 states only
# - State 1: early aging
# - State 2: Late aging. G2-off
# - State 3: Late aging. G2-on

#trMat <- matrix(c(
#        1, 0, 0,      # State1 → only State1 allowed
#        0, NA, NA,    # State2 → unknown, will be learned
#        0, 0, 1       # State3 → absorbin
#        ), nrow=3, byrow=TRUE)


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

ngenes <- nrow(X)
nsamples <- ncol(X)
emission_lik <- matrix(0, nsamples, nstates)

for(s in 1:nstates){
        for(c in 1:nsamples){
                mu <- sapply(1:ngenes,
                             function(g) fitmod@response[[s]][[g]]@parameters$coefficients)
                sd <- sapply(1:ngenes,
                             function(g) fitmod@response[[s]][[g]]@parameters$sd)
                emission_lik[c,s] <- prod(dnorm(X[, c], mean = mu, sd=sd))
        }
}

#T_m <- matrix(0, nstates, nstates)
#for(s in 1:n_states){
#        coefs <- coef(fitmod@transition[[s]]@parameters)
#        # Compute probabilities using multinomial logit formula
#        T_m[s,] <- softmax(c(0, coefs))  # first state is baseline
#}

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

edge_weights <- matrix(0, nsamples, nsamples)
for(i in 1:(nsamples-1)){
        T_i <- get_transmat(fitmod, ages[i])
        for(j in (i+1):nsamples){  # enforce age ordering
                w <- 0
                for(s in 1:n_states){
                        for(t in 1:n_states){
                                w <- w + emission_lik[i,s] * T_i[s,t] * emission_lik[j,t]
                        }
                }
                edge_weights[i,j] <- w
        }
}

edge_weights_log <- edge_weights#log(edge_weights + 1e-300)
edge_weights_log <- log(edge_weights + 1e-300)
edge_weights_norm <- (edge_weights_log - min(edge_weights_log, na.rm=TRUE)) /
        (max(edge_weights_log, na.rm=TRUE) - min(edge_weights_log, na.rm=TRUE))

cost_mat <- matrix(Inf,
                   nrow = nrow(edge_weights_norm),
                   ncol = ncol(edge_weights_norm),
                   dimnames = list(colnames(X),
                                   colnames(X)))

cost_mat[upper.tri(cost_mat)] <- 1 - edge_weights_norm[upper.tri(edge_weights_norm)]
diag(cost_mat) <- 0




# If get_transmat setting age to 0 the output is the same as the probabilities 
# at zero values given by fitmod@transition


softmax <- function(x){
        exp_x <- exp(x - max(x))
        return(exp_x / sum(exp_x))
}

fitmod@transition
