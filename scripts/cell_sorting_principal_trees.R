if(!require(DDRTree, quietly = T)) install.packages("DDRTree")
if(!require(Iso, quietly = T)) install.packages("Iso")
library(DDRTree)
library(ggplot2)
library(igraph)
library(Iso)


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
        X_cells#,
        #ncenter = length(unique(cell_metdat$donor))
)

plot(tree$Z[1, ], tree$Z[2, ])

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
