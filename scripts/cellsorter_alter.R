library(Seurat)
library(transformGamPoi)
library(FNN)
library(igraph)
library(parallel)
library(graph)
library(RBGL)
library(ggplot2)
library(RANN)
library(irlba)
library(matrixStats)

data_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/GSE254569/counts/extracted_data/"
outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/results/cellsorter/"
plotUtils::create_dir_if_not(outDir)
mdat_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/data/counts_data/GSE254569/metadata/metadata/metadata.csv"
target_cell <- "Exc_L2-3"
cell_ids_in <- "celltype_final"


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

# cellsorter funcs
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
                                dim_red = NULL,
                                n_var_features = NULL,
                                shiftlog = T,
                                center = T,
                                scale = F){
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
                        message(sprintf("Identifying the %s most variable genes...",
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
        # Scale and center
        means <- rowMeans(mat)
        sds <- sqrt(rowMeans(mat^2) - (rowMeans(mat))^2)
        if (any(sds == 0)){
                warning("There are %s genes with SD == 0. They will be removed.",
                        sum(sds == 0))
                mat <- mat[sds != 0, ]
                means <- means[sds != 0]
                sds <- sds[sds != 0]
        }
        if (center){
                cell_expr_stand <- (mat - means)
        }
        if (scale){
                cell_expr_stand <- cell_expr_stand / sds
        }
        if (scale | center){
                default_slot <- "cell_expr_stand"
        }else{
                default_slot <- "cell_expr"
        }
        meta <- seur@meta.data
        obj <- list(cell_expr = mat,
                    cell_expr_stand = cell_expr_stand,
                    pca_scores = NULL,
                    pca_loads = NULL,
                    cent_expr = NULL,
                    cent_expr_stand = NULL,
                    means = means,
                    sds = sds,
                    cell_metadata = meta,
                    centroid_metadata = NULL,
                    centroid_graph = NULL,
                    cell_graph = NULL,
                    default_slot = default_slot)
        if (!is.null(dim_red)){
                message(sprintf("Computing %s first PCs...", dim_red))
                pca <- prcomp_irlba(t(cell_expr_stand),
                                    n = dim_red,
                                    retx = T,
                                    center = F,
                                    scale. = F)
                obj$pca_scores <- t(pca$x)
                obj$pca_loads <- pca$rotation
                colnames(obj$pca_scores) <- colnames(cellsort$cell_expr)
                rownames(obj$pca_loads) <- rownames(cellsort$cell_expr)
                obj$default_slot <- "pca_scores"
        }
        class(obj) <- "cellsort"
        return(obj)
}

get_centroids <- function(cellsort, centroid_ids_in){
        #centroids_ids_in <- "donor"
        if (class(cellsort) != "cellsort"){
                stop("'cellsort' is not a cellsort object",
                     call. = F)
        }
        if (!centroid_ids_in %in% colnames(cellsort$cell_metadata)){
                stop(sprintf("%s is not a column in 'cellsort' cell_metadata slot",
                             centroid_ids_in),
                     call. = F)
        }
        if (!"age" %in% colnames(cellsort$cell_metadata)){
                stop("there is no 'age' column in cell_metadata")
        }
        centroid_ids <- unique(cellsort$cell_metadata[, centroid_ids_in])
        centroid_mat <- sapply(centroid_ids,
                               function(x) rowMeans(cellsort[[cellsort$default_slot]][, cellsort$cell_metadata[, centroid_ids_in] == x]))
        colnames(centroid_mat) <- centroid_ids
        cellsort$cent_expr <- centroid_mat
        cent_metdat <- data.frame(centroid_id = centroid_ids,
                                  age = cellsort$cell_metadata$age[match(centroid_ids,
                                                                         cellsort$cell_metadata[, centroid_ids_in])])
        cent_metdat <- cent_metdat[order(cent_metdat$age), ]
        rownames(cent_metdat) <- 1:nrow(cent_metdat)
        cellsort$cent_expr <- cellsort$cent_expr[, match(cent_metdat$centroid_id,
                                                         colnames(cellsort$cent_expr))]
        
        cellsort$centroid_metadata <- cent_metdat
        #means <- apply(cellsort$cent_expr, 1, mean)
        #sds <- apply(cellsort$cent_expr, 1, sd)
        #if (any(sds == 0)){
        #        warning("There are %s genes with SD == 0. They will be removed.",
        #                sum(sds == 0))
        #        means <- means[sds != 0]
        #        cellsort$cent_expr <- cellsort$cent_expr[sds != 0, ]
        #        cellsort$cell_expr <- cellsort$cell_expr[sds != 0, ]
        #        sds <- sds[sds != 0]
        #}
        #cellsort$cent_expr_stand <- sweep(sweep(cellsort$cent_expr,
        #                                        1,
        #                                        means,
        #                                        "-"),
        #                                  1,
        #                                  sds, "/")
        # Store means and sds to use them on cell data afterwards
        #cellsort$means <- means
        #cellsort$sds <- sds
        return(cellsort)
}

build_knn_graph <- function(cent_expr,
                            cent_mdat,
                            k = 10){
        # Obtain KNNs and create edge matrix
        nn <- get.knn(t(cent_expr),
                      k = k)
        n <- nrow(nn$nn.dist)
        edges <- cbind(as.character(cent_mdat$centroid_id[rep(1:n, k)]),
                       as.character(cent_mdat$centroid_id[as.vector(nn$nn.index)]))
        # Prune edges that don't satisfy age constraint
        edges <- as.data.frame(edges)
        edges$age_1 <- cent_mdat$age[match(edges[, 1],
                                           cent_mdat$centroid_id)]
        edges$age_2 <- cent_mdat$age[match(edges[, 2],
                                           cent_mdat$centroid_id)]
        edges$w <- 1 / (1 + as.vector(nn$nn.dist))
        edges <- edges[edges$age_2 >= edges$age_1, ]
        w <- edges$w
        edges <- as.matrix(edges[, 1:2])
        # Build the graph
        g <- graph_from_edgelist(edges, directed = T)
        #V(g)
        E(g)$weight <- w
        V(g)$age <- cent_mdat$age[match(names(V(g)),
                                        cent_mdat$centroid_id)]
        return(g)
}

group_close_sameage_centroids <- function(g, resolution = 1, alg = "louvain"){
        
        
        dup_ages <- unique(V(g)$age[duplicated(V(g)$age)])
        centr_info <- data.frame(centroid_id = names(V(g)),
                                 age = V(g)$age)
        #grpvec <- c()
        for (age in dup_ages){
                #age <- dup_ages[1]
                sameage_g <- as.undirected(subgraph(g,
                                                    names(V(g))[V(g)$age == age]),
                                           mode = "collapse",
                                           edge.attr.comb = "mean")
                subg_names <- names(V(sameage_g))
                centr_info <- centr_info[!centr_info$centroid_id %in% subg_names, ]
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
                to_bind <- data.frame(centroid_id = grp,
                                      age = rep(age, length(grp)))
                centr_info <- rbind.data.frame(centr_info, to_bind)
                #grpvec <- c(grpvec, grp)
        }
        centr_info <- centr_info[order(centr_info$age), ]
        rownames(centr_info) <- 1:nrow(centr_info)
        return(centr_info)
}

#build_collapsed_sample_info <- function(samp_info, close_samps){
#        
#        close_samps_split <- unique(unlist(sapply(close_samps,
#                                                  function(x) strsplit(x,
#                                                                       split = "-")[[1]])))
#        collapsed_sampinfo <- data.frame(donor = c(close_samps,
#                                                   samp_info$donor[!samp_info$donor %in% close_samps_split]))
#        collapsed_sampinfo$age <- samp_info$age[match(collapsed_sampinfo$donor, samp_info$donor)]
#        collapsed_sampinfo$age[is.na(collapsed_sampinfo$age)] <- as.numeric(sapply(close_samps,
#                                                                                   function(x) paste(unique(samp_info$age[match(strsplit(x, "-")[[1]], samp_info$donor)]),
#                                                                                                     collapse = "-")))
#        collapsed_sampinfo$donor <- make.names(collapsed_sampinfo$donor)
#        rownames(collapsed_sampinfo) <- collapsed_sampinfo$donor
#        return(collapsed_sampinfo)
#}

get_mean_close_samps <- function(cent_expr, cent_info_aggr,
                                 cell_expr, cell_info,
                                 centroid_info_in = "donor"){
        cent_expr_aggr <- matrix(nrow = nrow(cent_expr),
                                 ncol = 0,
                                 dimnames = list(rownames(cent_expr),
                                                 NULL))
        
        for(i in 1:nrow(cent_info_aggr)){
                s <- cent_info_aggr$centroid_id[i]
                if (s %in% colnames(cent_expr)){
                        exp_vec <- cent_expr[, s]
                        
                }else{
                        s_split <- strsplit(s, split = "-")[[1]]
                        # Compute a new centroid with all the merged cells
                        cent_expr_s <- cell_expr[, rownames(cell_info)[cell_info[, centroid_info_in] %in% s_split]]
                        exp_vec <- rowMeans(cent_expr_s)
                }
                cent_expr_aggr <- cbind(cent_expr_aggr,
                                        matrix(exp_vec,
                                               ncol = 1,
                                               dimnames = list(names(exp_vec),
                                                               s)))
        }
        return(cent_expr_aggr)
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
                cores <- detectCores() - 1
                if (cores < 1){
                        cores <- 1
                }
                all_paths <- mclapply(pairs, function(pair)
                {
                        max_avg_weight_path_dp(g, pair[1], pair[2])
                }, mc.cores = cores)
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

# Use Chiu-Edmonds algorithm to find optimum arborescence.
get_opt_branching <- function(g, root = NULL)
{
        edges <- as_edgelist(g)
        w <- E(g)$weight
        v <- names(V(g))
        gNEL <- ftM2graphNEL(edges, W = w,
                             V = v,
                             edgemode="directed")
        if (!is.null(root)){
                incoming <- inEdges(root, gNEL)
                for (from in incoming) {
                        gNEL <- removeEdge(from, root, gNEL)
                }
        }
        arbo <- edmondsOptimumBranching(gNEL)
        edge_list <- t(arbo$edgeList)  # returns from-to as list
        tree <- subgraph.edges(g, E(g, P = t(edge_list)))
        V(tree)$age <- V(g)$age[match(names(V(tree)),
                                      names(V(g)))]
        return(tree)
}

plot_graph <- function(g, node_info = NULL, color = NULL,
                       plot_labels = T,
                       node_size = 8,
                       edge_arrow_size = 0.5,
                       graph_layout = "fr",
                       root = NULL,
                       coords){
        #g <- cellsort_toy$cell_graph
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

build_centroid_graph <- function(cellsort,
                                 centroid_ids_in,
                                 k = 10,
                                 sameage_res = 1,
                                 sameage_clust_alg = "louvain",
                                 plot_intermed_graphs = T){
        message("Obtaining centroids...")
        #cellsort <- cellsort_toy
        #centroid_ids_in <- "orig.ident"
        cellsort <- get_centroids(cellsort = cellsort,
                                  centroid_ids_in = centroid_ids_in)
        #g <- build_knn_graph(cellsort$cent_expr_stand,
        #                     cellsort$centroid_metadata,
        #                     k = k)
        g <- build_knn_graph(cellsort$cent_expr,
                             cellsort$centroid_metadata,
                             k = k)
        print("jo")
        if(!is_dag(g)){
                cent_info_aggr <- group_close_sameage_centroids(g,
                                                                resolution = sameage_res,
                                                                alg = sameage_clust_alg)
                cent_expr_aggr <- get_mean_close_samps(cellsort$cent_expr,
                                                       cent_info_aggr = cent_info_aggr,
                                                       cell_expr = cellsort[[cellsort$default_slot]],
                                                       cell_info = cellsort$cell_metadata,
                                                       centroid_info_in = centroid_ids_in)
                cellsort$centroid_metadata <- cent_info_aggr
                cellsort$cent_expr <- cent_expr_aggr
                # Re-compute means and SDs for the aggregated centroids.
                #means <- apply(cellsort$cent_expr, 1, mean)
                #sds <- apply(cellsort$cent_expr, 1, sd)
                #if (any(sds == 0)){
                #        warning("There are %s genes with SD == 0. They will be removed.",
                #                sum(sds == 0))
                #        means <- means[sds != 0]
                #        cellsort$cent_expr <- cellsort$cent_expr[sds != 0, ]
                #        cellsort$cell_expr <- cellsort$cell_expr[sds != 0, ]
                #        sds <- sds[sds != 0]
                #}
                #cellsort$cent_expr_stand <- sweep(sweep(cellsort$cent_expr,
                #                                        1,
                #                                        means,
                #                                        "-"),
                #                                  1,
                #                                  sds, "/")
                # Store means and sds to use them on cell data afterwards
                #cellsort$means <- means
                #cellsort$sds <- sds
                #g <- build_knn_graph(cent_expr = cellsort$cent_expr_stand,
                #                     cent_mdat = cellsort$centroid_metadata,
                #                     k = k)
                g <- build_knn_graph(cent_expr = cellsort$cent_expr,
                                     cent_mdat = cellsort$centroid_metadata,
                                     k = k)
                if (plot_intermed_graphs){
                        plot_graph(g, color = "age")
                }
        }
        # Prune graph to keep only the maximum weight paths
        g <- get_all_max_avg_weight_g(g)
        if (plot_intermed_graphs){
                plot_graph(g, color = "age")
        }
        g <- get_opt_branching(g)
        if (plot_intermed_graphs){
                plot_graph(g, color = "age")
        }
        cellsort$centroid_graph <- g
        return(cellsort)
}

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
        message("Computing cell KNN graph...")
        start_time <- Sys.time()
        nn <- nn2(t(cell_mat), k = k)
        #nn <- FNN::get.knn(t(cell_mat), k = k)
        n <- nrow(nn$nn.dist)
        cell_edges <- cbind(colnames(cell_mat)[rep(1:n, k)],
                            colnames(cell_mat)[as.vector(nn$nn.idx)])
        # Obtain edge weights
        w <- 1 / (1 + as.vector(nn$nn.dist))
        names(w) <- paste(cell_edges[,1], cell_edges[,2], sep = "_")
        end_time <- Sys.time()
        elapsed <- round(as.numeric(end_time - start_time, units = "mins"),
                         digits = 3)
        message(sprintf("KNN graph is done. %s mins elapsed.",
                        elapsed, digits = 3))
        
        # Prune edges that do no satisfy backbone graph progression
        # based on edge assignments.
        in_edge <- edges[edges[, 1] == root, ]
        message("Pruning edges that don't satisfy backbone graph's constraints...")
        cell_edges <- prune_cell_edges(cell_edges, cell_edge_mat,
                                       edges, in_edge)
        w <- w[names(w) %in% paste(cell_edges[,1], cell_edges[,2], sep = "_")]
        #cell_g <- graph_from_edgelist(cell_edges, directed = T)
        E(cell_g)$weight <- w
        message(sprintf("Computing %s optimal graph...", algo))
        start_time <- Sys.time()
        if (algo == "mst"){
                cell_g <- graph_from_edgelist(cell_edges, directed = F)
                E(cell_g)$weight <- w
                tree <- mst(cell_g)
        }else if (algo == "edmonds"){
                cell_g <- graph_from_edgelist(cell_edges, directed = T)
                E(cell_g)$weight <- w
                nodes <- unique(c(cell_edges[,1], cell_edges[,2]))
                gNEL <- ftM2graphNEL(cell_edges, W = w,
                                     V = nodes,
                                     edgemode="directed")

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
                                             rownames(cell_metdat))]
        end_time <- Sys.time()
        elapsed <- round(as.numeric(end_time - start_time, units = "mins"),
                         digits = 3)
        message(sprintf("Optimal graph finished. %s mins elapsed.",
                        elapsed, digits = 3))
        return(tree)
}

# 
build_cell_graph <- function(cellsort,
                             centroid_ids_in,
                             k = 100,
                             algo = "edmonds",
                             parallelize = T,
                             cores = NULL){
        
        #cellsort <- cellsort_toy
        centroid_ids_in <- "donor"
        #centroid_ids_in <- "orig.ident"
        parallelize <- T
        k <- 20
        algo <- "edmonds"
        cores = 6
        
        start_time <- Sys.time()
        #cell_mat <- cellsort$cell_expr_stand
        cell_mat <- cellsort[[cellsort$default_slot]]
        backbone_graph <- cellsort$centroid_graph
        #centroid_mat <- cellsort$cent_expr_stand
        centroid_mat <- cellsort$cent_expr
        cell_metdat <- cellsort$cell_metadata
        
        #backbone_graph <- g_pruned
        
        
        # Find leaves and root by finding endpoints of graph
        leaves <- names(which(igraph::degree(backbone_graph,
                                             mode = "out") == 0))
        root <- names(which(igraph::degree(backbone_graph,
                                           mode = "in") == 0))
        # Find all the paths from the root to the leaves
        all_paths <- list()
        for (leaf in leaves) {
                paths <- igraph::all_simple_paths(backbone_graph,
                                                  from = root, to = leaf)
                all_paths <- c(all_paths, paths)
        }
        all_paths <- lapply(all_paths, names)
        
        if (length(all_paths) > 2){
                # Keep only paths longer than the median
                path_lens <- sapply(all_paths, length)
                all_paths <- all_paths[path_lens > mean(path_lens)]
        }
        
        
        # Find the trunk by finding the common nodes from across the paths.
        root_branch <- intersect_recursive(all_paths)
        
        trunk_g <- subgraph(backbone_graph, root_branch)
        
        # Find the last node in the trunk.
        last_root_branch <- last_in_branch(root_branch, trunk_g, root)
        
        # Normalize cell matrix to the mean and sd of the centroids, as it's what
        # was used to build the backbone graph
        #mean_centroids <- apply(centroid_mat, 1, mean)
        #sd_centroids <- apply(centroid_mat, 1, sd)
        #centroid_mat <- (centroid_mat - mean_centroids)/sd_centroids
        #message("Scaling cell expression...")
        #cell_expr_stand <- cellsort$cell_expr
        #cell_expr_stand <- t(cellsort$cell_expr)
        #cell_expr_stand <- sweep(cell_expr_stand, 2, cellsort$means, "-")
        #cell_expr_stand <- sweep(cell_expr_stand, 2, cellsort$sds, "/")
        #cellsort$cell_expr_stand <- t(cell_expr_stand)
        #cell_expr_stand <- (cell_mat - cellsort$means)/cellsort$sds
        #cellsort$cell_expr_stand <- cell_expr_stand
        
        # Find the youngest cell by identifying the direction of change
        # across the trunk, and identifying the cell that goes most in the 
        # oposite direction.
        dir_root <- centroid_mat[, last_root_branch] - centroid_mat[, root]
        
        root_donors <- strsplit(root, split = "-", fixed = T)[[1]]
        root_cells <- cellsort$cell_expr_stand[, cellsort$cell_metadata[, centroid_ids_in] %in% root_donors]
        
        # Also works, its the same. But crossprod is faster
        # youngest_cell <- colnames(root_cells)[which.min(as.vector(t(dir_root) %*% root_cells))]
        youngest_cell <- colnames(root_cells)[which.min(crossprod(as.matrix(root_cells),
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
        #cell_edge <- assign_cell_to_edge(cellsort$cell_expr_stand, centroid_mat, edges,
        #                                 parallelize = parallelize,
        #                                 cores = cores)
        cell_edge <- assign_cell_to_edge(cellsort$cell_expr, centroid_mat, edges,
                                         parallelize = parallelize,
                                         cores = cores)
        
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
        cellsort$cell_graph <- cell_g
        end_time <- Sys.time()
        elapsed <- round(as.numeric(end_time - start_time, units = "mins"),
                         digits = 3)
        message(sprintf("Total time: %s mins.",
                        elapsed))
        return(cellsort)
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
plot_centroids(all_constant, 1, 2)
plot_centroids(activtline_3, 1, 2)
plot_centroids(accelerate_3, 1, 2)

X <- activtline_3
#X <- accelerate_3

X_cells <- get_toy_cells(X, sd = 15, ncells = 10)
X_cells[X_cells < 0] <- 0
X_cells <- round(X_cells)

cell_metdat <- data.frame(cell_name = colnames(X_cells),
                          donor = gsub("\\_.*", "", colnames(X_cells)))

cell_metdat$age <- toy_samps$age[match(cell_metdat$donor,
                                       toy_samps$donor)]
X_cells <- as(X_cells, "dgCMatrix")
#rownames(X_cells) <- paste0("Gene", seq_len(nrow(X_cells)))
#colnames(X_cells) <- paste0("Cell", seq_len(ncol(X_cells)))

seur_toy <- CreateSeuratObject(counts = X_cells,
                               project = "MyProject",
                               min.cells = 3,
                               min.features = 1)
seur_toy@meta.data$age <- cell_metdat$age[match(rownames(seur_toy@meta.data),
                                                cell_metdat$cell_name)]

cellsort_toy <- create_cellsort_obj(seur_toy,
                                    cell_ids_in = NULL,
                                    target_cell = NULL,
                                    shiftlog = F)
cellsort_toy <- build_centroid_graph(cellsort = cellsort_toy,
                                     centroid_ids_in = "orig.ident",
                                     k = 7,
                                     sameage_res = 1,
                                     sameage_clust_alg = "louvain",
                                     plot_intermed_graphs = T)

cellsort_toy <- build_cell_graph(cellsort_toy,
                                 centroid_ids_in = "orig.ident",
                                 k = 100,
                                 algo = "edmonds",
                                 parallelize = T)

plot_graph(cellsort_toy$cell_graph,
           color = "age",
           plot_labels = F,
           node_size = 3,
           graph_layout = "tree",
           root = names(which(igraph::degree(cellsort_toy$cell_graph,
                                             mode = "in") == 0)))

################################################################################
# Real data                                                                    #
################################################################################

# Load data and create seurat object
################################################################################

seur <- Seurat::Read10X(data.dir = data_dir)
seur <- CreateSeuratObject(counts = seur)
metdat <- plotUtils::read_table_fast(mdat_file)

feats <- read.table(list.files(data_dir, full.names = T)[2])
rownames(feats) <- feats[, 2]

seur@meta.data$celltype_final <- metdat$celltypes_final[match(rownames(seur@meta.data),
                                                              metdat$index)]
seur@meta.data$celltype_major <- metdat$major_celltypes[match(rownames(seur@meta.data),
                                                              metdat$index)]
seur@meta.data$donor <- metdat$Donor[match(rownames(seur@meta.data),
                                           metdat$index)]
seur@meta.data$age <- metdat$Age[match(rownames(seur@meta.data),
                                       metdat$index)]
seur[["RNA"]]@meta.features <- feats

# Run workflow
################################################################################

cellsort <- create_cellsort_obj(seur,
                                cell_ids_in = cell_ids_in,
                                target_cell = target_cell,
                                n_var_features = 5000,
                                shiftlog = T,
                                dim_red = 100,
                                center = T,
                                scale = F)
saveRDS(cellsort, file = sprintf("%scellsort_alt.rds", outDir))
cellsort <- readRDS(sprintf("%scellsort_alt.rds", outDir))

cellsort <- build_centroid_graph(cellsort = cellsort,
                                 centroid_ids_in = "donor",
                                 k = 20,
                                 sameage_res = 1,
                                 sameage_clust_alg = "louvain",
                                 plot_intermed_graphs = F)

plot_graph(cellsort$centroid_graph, color = "age")


cellsort <- build_cell_graph(cellsort = cellsort,
                             centroid_ids_in = "donor",
                             k = 20,
                             algo = "edmonds",
                             parallelize = T,
                             cores = 8)

plot_graph(cellsort$cell_graph,
           color = "age",
           plot_labels = F,
           node_size = 3,
           graph_layout = "tree",
           root = names(which(igraph::degree(cellsort$cell_graph,
                                             mode = "in") == 0)))
