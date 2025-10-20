if(!require("devtools")){
        install.packages("devtools")
}
if (!require(distutils_agemod, quietly = T)){
        devtools::install_github("guisantagui/distutils_agemod", upgrade = "never")
}
library(distutils)

if (!require("ElPiGraph.R", quietly = T)){
        devtools::install_github("Albluca/ElPiGraph.R", upgrade = "never")
}
library(ElPiGraph.R)
library(ggplot2)
library(igraph)
library(Iso)
library(Seurat)
library(transformGamPoi)
library(igraph)

library(Rcpp)
#setwd("/Users/guillem.santamaria/Documents/postdoc/comput/distutils")

#Rcpp::compileAttributes()
#devtools::install()

# Dir stuff
################################################################################

outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/results/cellsort/"

seur_file <- sprintf("%sGSE254569_ctrls_seur.rds", outDir)
seur_toy_file <- sprintf("%sseur_toy.rds", outDir)

donor_ids_in <- "donor"
cell_ids_in <- "celltype_final"
target_cell <- "Exc_L2-3"

donor_ids_in_toy <- "orig.ident"
cell_ids_in_toy <- NULL
target_cell_toy <- NULL

# Load data
################################################################################

seur <- readRDS(seur_file)

seur_toy <- readRDS(seur_toy_file)

# cellsorter functions (only for normalization purposes, vestigial. Need to be
# removed)
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

# Mod ElPiGraph.R Funcs
################################################################################

# Edit computeElasticPrincipalTree to call custom
# computeElasticPrincipalGraphWithGrammars_edit.
computeElasticPrincipalTree_edit <- function(X,
                                             NumNodes,
                                             NumEdges = Inf,
                                             InitNodes = 2,
                                             Lambda = 0.01,
                                             Mu = 0.1,
                                             MaxNumberOfIterations = 10,
                                             TrimmingRadius = Inf,
                                             eps = .01,
                                             Do_PCA = TRUE,
                                             InitNodePositions = NULL,
                                             InitEdges = NULL,
                                             ElasticMatrix = NULL,
                                             AdjustVect = NULL,
                                             CenterData = TRUE,
                                             ComputeMSEP = TRUE,
                                             verbose = FALSE,
                                             ShowTimer = FALSE,
                                             ReduceDimension = NULL,
                                             drawAccuracyComplexity = TRUE,
                                             drawPCAView = TRUE,
                                             drawEnergy = TRUE,
                                             n.cores = 1,
                                             MinParOP = 20,
                                             ClusType = "Sock",
                                             nReps = 1,
                                             Subsets = list(),
                                             ProbPoint = 1,
                                             Mode = 1,
                                             FinalEnergy = "Base",
                                             alpha = 0,
                                             beta = 0,
                                             FastSolve = FALSE,
                                             ICOver = NULL,
                                             DensityRadius = NULL,
                                             AvoidSolitary = FALSE,
                                             EmbPointProb = 1,
                                             ParallelRep = FALSE,
                                             AvoidResampling = FALSE,
                                             SampleIC = TRUE,
                                             AdjustElasticMatrix = NULL,
                                             AdjustElasticMatrix.Initial = NULL,
                                             Lambda.Initial = NULL,
                                             Mu.Initial = NULL,
                                             age_vec = NULL,
                                             eta = 1e-03,
                                             ...) {
        
        
        # Difine the initial configuration
        
        if(is.null(ICOver)){
                Configuration <- "Line"
        } else {
                Configuration <- ICOver
        }
        
        return(
                computeElasticPrincipalGraphWithGrammars_edit(X = X,
                                                              NumNodes = NumNodes,
                                                              NumEdges = NumEdges,
                                                              InitNodes = InitNodes,
                                                              Lambda = Lambda,
                                                              Mu = Mu,
                                                              GrowGrammars = list(c('bisectedge','addnode2node'),c('bisectedge','addnode2node')),
                                                              ShrinkGrammars = list(c('shrinkedge','removenode')),
                                                              MaxNumberOfIterations = MaxNumberOfIterations,
                                                              TrimmingRadius = TrimmingRadius,
                                                              eps = eps,
                                                              Do_PCA = Do_PCA,
                                                              InitNodePositions = InitNodePositions,
                                                              InitEdges = InitEdges,
                                                              AdjustVect = AdjustVect,
                                                              Configuration = Configuration,
                                                              CenterData = CenterData,
                                                              ComputeMSEP = ComputeMSEP,
                                                              verbose = verbose,
                                                              ShowTimer = ShowTimer,
                                                              ReduceDimension = ReduceDimension,
                                                              drawAccuracyComplexity = drawAccuracyComplexity,
                                                              drawPCAView = drawPCAView,
                                                              drawEnergy = drawEnergy,
                                                              n.cores = n.cores,
                                                              MinParOP = MinParOP,
                                                              ClusType = ClusType,
                                                              nReps = nReps,
                                                              Subsets = list(),
                                                              ProbPoint = ProbPoint,
                                                              Mode = Mode,
                                                              FinalEnergy = FinalEnergy,
                                                              alpha = alpha,
                                                              beta = beta,
                                                              gamma = gamma,
                                                              FastSolve = FastSolve,
                                                              DensityRadius = DensityRadius,
                                                              AvoidSolitary = AvoidSolitary,
                                                              EmbPointProb = EmbPointProb,
                                                              SampleIC = SampleIC,
                                                              AvoidResampling = AvoidResampling,
                                                              ParallelRep = ParallelRep,
                                                              AdjustElasticMatrix = AdjustElasticMatrix,
                                                              AdjustElasticMatrix.Initial = AdjustElasticMatrix.Initial,
                                                              Lambda.Initial = Lambda.Initial, Mu.Initial = Mu.Initial,
                                                              age_vec = age_vec,
                                                              eta = eta,
                                                              ...)
        )
        
}

# Edit computeElasticPrincipalGraphWithGrammars to make it call customized
# functions that allocate age. Also allows bootstrapping with the age, by
# selecting nReps > 1.
computeElasticPrincipalGraphWithGrammars_edit <- function(X,
                                                          NumNodes,
                                                          NumEdges = Inf,
                                                          InitNodes = 2,
                                                          Lambda = 0.01,
                                                          Mu = 0.1,
                                                          GrowGrammars,
                                                          ShrinkGrammars,
                                                          GrammarOptimization = FALSE,
                                                          MaxSteps = Inf,
                                                          GrammarOrder = c("Grow", "Shrink"),
                                                          MaxNumberOfIterations = 10,
                                                          TrimmingRadius = Inf,
                                                          eps = .01,
                                                          Do_PCA = TRUE,
                                                          InitNodePositions = NULL,
                                                          AdjustVect = NULL,
                                                          ElasticMatrix = NULL,
                                                          InitEdges = NULL,
                                                          CenterData = TRUE,
                                                          ComputeMSEP = TRUE,
                                                          verbose = FALSE,
                                                          ShowTimer = FALSE,
                                                          ReduceDimension = NULL,
                                                          drawAccuracyComplexity = TRUE,
                                                          drawPCAView = TRUE,
                                                          drawEnergy = TRUE,
                                                          n.cores = 1,
                                                          ClusType = "Sock",
                                                          MinParOP = MinParOP,
                                                          nReps = 1,
                                                          ParallelRep = FALSE,
                                                          Subsets = list(),
                                                          ProbPoint = 1,
                                                          Mode = 1,
                                                          FinalEnergy = "Base",
                                                          alpha = 0,
                                                          beta = 0,
                                                          gamma = 0,
                                                          FastSolve = FALSE,
                                                          Configuration = NULL,
                                                          DensityRadius = NULL,
                                                          AvoidSolitary = FALSE,
                                                          EmbPointProb = 1,
                                                          SampleIC = TRUE,
                                                          AvoidResampling = TRUE,
                                                          AdjustElasticMatrix = NULL,
                                                          AdjustElasticMatrix.Initial = NULL,
                                                          Lambda.Initial = NULL, Mu.Initial = NULL,
                                                          age_vec = NULL,
                                                          eta = 1e-3,
                                                          ...) {
        
        # Create a cluster if requested
        
        if(n.cores > 1){
                if(ClusType == "Fork"){
                        print(paste("Creating a fork cluster with", n.cores, "nodes"))
                        cl <- parallel::makeCluster(n.cores, type="FORK")
                } else {
                        print(paste("Creating a sock cluster with", n.cores, "nodes"))
                        cl <- parallel::makeCluster(n.cores)
                        parallel::clusterExport(cl, varlist = c("X"), envir=environment())
                }
        } else {
                cl = 1
        }
        
        # Be default we are using a predefined initial configuration
        ComputeIC <- FALSE
        
        # Generate a dummy subset is not specified
        if(length(Subsets) == 0){
                Subsets[[1]] <- 1:ncol(X)
        }
        
        # Prepare the list to be returned
        ReturnList <- list()
        
        # Copy the original matrix, this is needed in case of subsetting
        Base_X <- X
        
        # For each subset
        for(j in 1:length(Subsets)){
                
                # Generate the appropriate matrix
                X <- Base_X[, Subsets[[j]]]
                
                # Export the data matrix to the cluster is needed
                if(n.cores > 1 & ClusType != "Fork"){
                        parallel::clusterExport(cl, varlist = c("X"), envir=environment())
                }
                
                # Define temporary variable to avoid excessing plotting
                Intermediate.drawPCAView <- drawPCAView
                Intermediate.drawAccuracyComplexity <- drawAccuracyComplexity 
                Intermediate.drawEnergy <- drawEnergy
                
                # print(Subsets[[j]])
                # print(Subsets[[j]] %in% colnames(Base_X))
                # print(dim(Base_X))
                
                if(ParallelRep & n.cores > 1){
                        
                        print("Using parallel sampling analysis. Limited output available")
                        
                        if(ProbPoint<1 & ProbPoint>0){
                                SelPoints <- lapply(as.list(1:nReps), function(i){return(runif(nrow(X)) <= ProbPoint)})
                        } else {
                                SelPoints <- lapply(as.list(1:nReps), function(i){return(rep(TRUE, nrow(X)))})
                        }
                        
                        # Do we need to compute the initial conditions?
                        if(is.null(InitNodePositions) | (is.null(InitEdges) & is.null(ElasticMatrix))){
                                
                                print("Generating the initial configurations")
                                
                                if(ClusType != "Fork"){
                                        print("Exporting initial configuration parameters")
                                        parallel::clusterExport(cl, varlist = c("InitNodes", "Configuration", "DensityRadius",
                                                                                "SelPoints", "Lambda", "Mu", "MaxNumberOfIterations",
                                                                                "TrimmingRadius", "eps", "Mode"), envir=environment())
                                }
                                
                                # We are computing the initial conditions. InitNodePositions need to be reset after each step!
                                ComputeIC <- TRUE
                                
                                if(SampleIC){
                                        
                                        if(AvoidResampling){
                                                
                                                Used <- rep(FALSE, nrow(X))
                                                InitialConf.List <- list()
                                                
                                                for(i in 1:nReps){
                                                        
                                                        print(paste("Rep", i, sum(Used), "out of", length(Used), "points have been used"))
                                                        
                                                        tryCatch(
                                                                InitialConf.List[[i]] <-
                                                                        ElPiGraph.R:::generateInitialConfiguration(X[SelPoints[[i]] & !Used, ],
                                                                                                                   Nodes = InitNodes,
                                                                                                                   Configuration = Configuration,
                                                                                                                   DensityRadius = DensityRadius),
                                                                error = function(e){
                                                                        print(e)
                                                                        print("Resetting Initial point set")
                                                                        Used <<- rep(FALSE, nrow(X))
                                                                        InitialConf.List[[i]] <<-
                                                                                ElPiGraph.R:::generateInitialConfiguration(X[SelPoints[[i]] & !Used, ],
                                                                                                                           Nodes = InitNodes,
                                                                                                                           Configuration = Configuration,
                                                                                                                           DensityRadius = DensityRadius)
                                                                },
                                                                warning = {})
                                                        
                                                        Dist <- apply(distutils::PartialDistance(InitialConf.List[[i]]$NodePositions, Br = X), 2, min)
                                                        
                                                        if(!is.null(DensityRadius)){
                                                                Used <- Used | (Dist < DensityRadius)
                                                        } else {
                                                                Used <- Used | (Dist < .Machine$double.xmin)
                                                        }
                                                        
                                                        
                                                }
                                                
                                                
                                        } else {
                                                
                                                # Construct the initial configuration
                                                InitialConf.List <- 
                                                        parallel::parLapply(cl, as.list(1:nReps), function(i){
                                                                ElPiGraph.R:::generateInitialConfiguration(X[SelPoints[[i]], ], Nodes = InitNodes,
                                                                                                           Configuration = Configuration, DensityRadius = DensityRadius)
                                                        })
                                                
                                        }
                                        
                                        
                                        
                                        
                                } else {
                                        
                                        if(AvoidResampling){
                                                
                                                Used <- rep(FALSE, nrow(X))
                                                InitialConf.List <- list()
                                                
                                                for(i in 1:nReps){
                                                        
                                                        print(paste("Rep", i, sum(Used), "out of", length(Used), "points have been used"))
                                                        
                                                        tryCatch(
                                                                InitialConf.List[[i]] <-
                                                                        ElPiGraph.R:::generateInitialConfiguration(X[!Used, ],
                                                                                                                   Nodes = InitNodes,
                                                                                                                   Configuration = Configuration,
                                                                                                                   DensityRadius = DensityRadius),
                                                                error = function(e){
                                                                        print(e)
                                                                        print("Resetting Initial point set")
                                                                        Used <<- rep(FALSE, nrow(X))
                                                                        InitialConf.List[[i]] <<-
                                                                                ElPiGraph.R:::generateInitialConfiguration(X[!Used, ],
                                                                                                                           Nodes = InitNodes,
                                                                                                                           Configuration = Configuration,
                                                                                                                           DensityRadius = DensityRadius)
                                                                },
                                                                warning = {})
                                                        
                                                        Dist <- apply(distutils::PartialDistance(InitialConf.List[[i]]$NodePositions, Br = X), 2, min)
                                                        
                                                        if(!is.null(DensityRadius)){
                                                                Used <- Used | (Dist < DensityRadius)
                                                        } else {
                                                                Used <- Used | (Dist <= .Machine$double.xmin)
                                                        }
                                                        
                                                }
                                                
                                                
                                        } else {
                                                
                                                # Construct the initial configuration
                                                InitialConf.List <- 
                                                        parallel::parLapply(cl, as.list(1:nReps), function(i){
                                                                ElPiGraph.R:::generateInitialConfiguration(X, Nodes = InitNodes,
                                                                                                           Configuration = Configuration, DensityRadius = DensityRadius)
                                                        })
                                                
                                        }
                                        
                                        
                                        
                                }
                                
                                # Set the initial edge configuration
                                InitEdges.List <- lapply(InitialConf.List, "[[", "Edges")
                                
                                # Compute the initial elastic matrices
                                ElasticMatrix.List <- parallel::parLapply(cl, InitEdges.List, function(Edges){
                                        ElPiGraph.R:::Encode2ElasticMatrix(Edges = Edges, Lambdas = Lambda, Mus = Mu)
                                }) 
                                
                                if(ClusType != "Fork"){
                                        print("Exporting initial configuration parameters")
                                        parallel::clusterExport(cl, varlist = c("InitialConf.List", "ElasticMatrix.List"),
                                                                envir=environment())
                                }
                                
                                # Compute the initial node position
                                InitNodePositions.List <- parallel::parLapply(cl, as.list(1:nReps), function(i){
                                        PrimitiveElasticGraphEmbedment_edit(
                                                X = X,
                                                NodePositions = InitialConf.List[[i]]$NodePositions,
                                                MaxNumberOfIterations = MaxNumberOfIterations,
                                                TrimmingRadius = TrimmingRadius, eps = eps,
                                                ElasticMatrix = ElasticMatrix.List[[i]],
                                                Mode = Mode,
                                                age_vec = age_vec,
                                                eta = eta)$EmbeddedNodePositions
                                })  
                        }
                        
                        
                        # Do we need to compute AdjustVect?
                        if(is.null(AdjustVect)){
                                AdjustVect <- rep(FALSE, nrow(InitNodePositions.List[[1]]))
                        }
                        
                        # Prevent plotting after a few examples
                        print("Graphical output will be suppressed")
                        Intermediate.drawPCAView <- FALSE
                        Intermediate.drawAccuracyComplexity <- FALSE 
                        Intermediate.drawEnergy <- FALSE
                        
                        # print(paste("Constructing tree", i, "of", nReps, "/ Subset", j, "of", length(Subsets)))
                        
                        
                        if(ClusType != "Fork"){
                                print("Exporting parameters")
                                parallel::clusterExport(cl, varlist = c("InitNodePositions.List", "InitEdges.List", "ElasticMatrix.List",
                                                                        "NumNodes", "NumEdges", "InitNodePositions.List", "InitEdges.List",
                                                                        "ElasticMatrix.List", "AdjustVect", "GrowGrammars", "ShrinkGrammars",
                                                                        "GrammarOptimization", "MaxSteps", "GrammarOrder", "MaxNumberOfIterations",
                                                                        "TrimmingRadius", "eps", "Lambda", "Mu", "Do_PCA", "CenterData", "ComputeMSEP",
                                                                        "ReduceDimension", "Mode", "FinalEnergy", "alpha", "beta", "gamma",
                                                                        "Intermediate.drawAccuracyComplexity", "Intermediate.drawPCAView", "Intermediate.drawEnergy",
                                                                        "FastSolve", "AvoidSolitary", "EmbPointProb", "AdjustElasticMatrix", "AdjustElasticMatrix.Initial",
                                                                        "Lambda.Initial", "Mu.Initial"),
                                                        envir=environment())
                        }
                        
                        
                        
                        print("Analysis is running ... no output will be shown")
                        print("Sit back and relax, this may take a long time ...")
                        
                        tictoc::tic()
                        
                        # Run the ElPiGraph algorithm
                        ReturnList <- parallel::parLapply(cl, as.list(1:nReps), function(i, ...){
                                computeElasticPrincipalGraph_edit(Data = X[SelPoints[[i]], ],
                                                                  NumNodes = NumNodes,
                                                                  NumEdges = NumEdges,
                                                                  InitNodePositions = InitNodePositions.List[[i]],
                                                                  InitEdges = InitEdges.List[[i]],
                                                                  ElasticMatrix = ElasticMatrix.List[[i]],
                                                                  AdjustVect = AdjustVect,
                                                                  GrowGrammars = GrowGrammars,
                                                                  ShrinkGrammars = ShrinkGrammars,
                                                                  GrammarOptimization = GrammarOptimization,
                                                                  MaxSteps = MaxSteps,
                                                                  GrammarOrder = GrammarOrder,
                                                                  MaxNumberOfIterations = MaxNumberOfIterations,
                                                                  TrimmingRadius = TrimmingRadius,
                                                                  eps = eps,
                                                                  Lambda = Lambda,
                                                                  Mu = Mu,
                                                                  Do_PCA = Do_PCA,
                                                                  CenterData = CenterData,
                                                                  ComputeMSEP = ComputeMSEP,
                                                                  verbose = FALSE,
                                                                  ShowTimer = FALSE,
                                                                  ReduceDimension = ReduceDimension,
                                                                  Mode = Mode,
                                                                  FinalEnergy = FinalEnergy,
                                                                  alpha = alpha,
                                                                  beta = beta,
                                                                  gamma = gamma,
                                                                  drawAccuracyComplexity = Intermediate.drawAccuracyComplexity,
                                                                  drawPCAView = Intermediate.drawPCAView,
                                                                  drawEnergy = Intermediate.drawEnergy,
                                                                  n.cores = 1,
                                                                  FastSolve = FastSolve,
                                                                  AvoidSolitary = AvoidSolitary,
                                                                  EmbPointProb = EmbPointProb,
                                                                  AdjustElasticMatrix = AdjustElasticMatrix,
                                                                  AdjustElasticMatrix.Initial = AdjustElasticMatrix.Initial,
                                                                  Lambda.Initial = Lambda.Initial,
                                                                  Mu.Initial = Mu.Initial,
                                                                  age_vec = age_vec[SelPoints[[i]]],
                                                                  eta = eta,
                                                                  ...)
                                
                                
                        })
                        
                        tictoc::toc()
                        
                        # Save extra information
                        for(i in 1:length(ReturnList)){
                                ReturnList[[i]]$SubSetID <- j
                                ReturnList[[i]]$ReplicaID <- i
                                ReturnList[[i]]$ProbPoint <- ProbPoint
                        }
                        
                        # Reset InitNodePositions for the next iteration
                        if(ComputeIC){
                                InitNodePositions <- NULL
                        }
                        
                        
                } else {
                        
                        Used <- rep(FALSE, nrow(X))
                        
                        for(i in 1:nReps){
                                
                                # Select the poits to be used
                                if(ProbPoint<1 & ProbPoint>0){
                                        SelPoints <- runif(nrow(X)) <= ProbPoint
                                } else {
                                        SelPoints <- rep(TRUE, nrow(X))
                                }
                                
                                # Do we need to compute the initial conditions?
                                if(is.null(InitNodePositions) | (is.null(InitEdges) & is.null(ElasticMatrix))){
                                        
                                        print("Generating the initial configuration")
                                        
                                        # We are computing the initial conditions. InitNodePositions need to be reset after each step!
                                        ComputeIC <- TRUE
                                        
                                        if(SampleIC){
                                                
                                                if(AvoidResampling){
                                                        
                                                        InitialConf <- 
                                                                ElPiGraph.R:::generateInitialConfiguration(X[SelPoints & !Used, ],
                                                                                                           Nodes = InitNodes,
                                                                                                           Configuration = Configuration,
                                                                                                           DensityRadius = DensityRadius)
                                                        
                                                        Dist <- apply(distutils::PartialDistance(InitialConf$NodePositions, Br = X), 2, min)
                                                        
                                                        if(!is.null(DensityRadius)){
                                                                Used <- Used | (Dist < DensityRadius)
                                                        } else {
                                                                Used <- Used | (Dist <= .Machine$double.xmin)
                                                        }
                                                        
                                                        if(sum(Used) < nrow(X)*.9){
                                                                print("90% of the points have been used as initial conditions. Resetting.")
                                                        }
                                                        
                                                } else {
                                                        # Construct the initial configuration
                                                        InitialConf <- 
                                                                ElPiGraph.R:::generateInitialConfiguration(X[SelPoints, ],
                                                                                                           Nodes = InitNodes,
                                                                                                           Configuration = Configuration,
                                                                                                           DensityRadius = DensityRadius)
                                                }
                                                
                                        } else {
                                                
                                                if(AvoidResampling){
                                                        
                                                        InitialConf <- 
                                                                ElPiGraph.R:::generateInitialConfiguration(X[!Used, ],
                                                                                                           Nodes = InitNodes,
                                                                                                           Configuration = Configuration,
                                                                                                           DensityRadius = DensityRadius)
                                                        
                                                        Dist <- apply(distutils::PartialDistance(InitialConf$NodePositions, Br = X), 2, min)
                                                        
                                                        if(!is.null(DensityRadius)){
                                                                Used <- Used | (Dist < DensityRadius)
                                                        } else {
                                                                Used <- Used | (Dist < .Machine$double.xmin)
                                                        }
                                                        
                                                        if(sum(Used) > nrow(X)*.9){
                                                                print("90% of the points have been used as initial conditions. Resetting.")
                                                        }
                                                        
                                                        
                                                } else {
                                                        
                                                        # Construct the initial configuration
                                                        InitialConf <- 
                                                                ElPiGraph.R:::generateInitialConfiguration(X[, ],
                                                                                                           Nodes = InitNodes,
                                                                                                           Configuration = Configuration,
                                                                                                           DensityRadius = DensityRadius)
                                                }
                                                
                                        }
                                        
                                        # Set the initial edge configuration
                                        InitEdges <- InitialConf$Edges
                                        
                                        # Compute the initial elastic matrix
                                        ElasticMatrix <- ElPiGraph.R:::Encode2ElasticMatrix(Edges = InitialConf$Edges, Lambdas = Lambda, Mus = Mu)
                                        
                                        # Compute the initial node position
                                        InitNodePositions <- PrimitiveElasticGraphEmbedment_edit(
                                                X = X, NodePositions = InitialConf$NodePositions,
                                                MaxNumberOfIterations = MaxNumberOfIterations, TrimmingRadius = TrimmingRadius, eps = eps,
                                                ElasticMatrix = ElasticMatrix, Mode = Mode,
                                                age_vec = age_vec,
                                                eta = eta)$EmbeddedNodePositions
                                }
                                
                                # Do we need to compute AdjustVect?
                                if(is.null(AdjustVect)){
                                        AdjustVect <- rep(FALSE, nrow(InitNodePositions))
                                }
                                
                                # Limit plotting after a few examples
                                if(length(ReturnList) == 3){
                                        print("Graphical output will be suppressed for the remaining replicas")
                                        Intermediate.drawPCAView <- FALSE
                                        Intermediate.drawAccuracyComplexity <- FALSE 
                                        Intermediate.drawEnergy <- FALSE
                                }
                                
                                print(paste("Constructing tree", i, "of", nReps, "/ Subset", j, "of", length(Subsets)))
                                
                                # Run the ElPiGraph algorithm
                                ReturnList[[length(ReturnList)+1]] <- computeElasticPrincipalGraph_edit(Data = X[SelPoints, ], NumNodes = NumNodes, NumEdges = NumEdges,
                                                                                                        InitNodePositions = InitNodePositions, InitEdges = InitEdges, ElasticMatrix = ElasticMatrix,
                                                                                                        AdjustVect = AdjustVect,
                                                                                                        GrowGrammars = GrowGrammars,
                                                                                                        ShrinkGrammars = ShrinkGrammars,
                                                                                                        GrammarOptimization = GrammarOptimization,
                                                                                                        MaxSteps = MaxSteps,
                                                                                                        GrammarOrder = GrammarOrder,
                                                                                                        MaxNumberOfIterations = MaxNumberOfIterations, TrimmingRadius = TrimmingRadius, eps = eps,
                                                                                                        Lambda = Lambda, Mu = Mu, Do_PCA = Do_PCA,
                                                                                                        CenterData = CenterData, ComputeMSEP = ComputeMSEP,
                                                                                                        verbose = verbose, ShowTimer = ShowTimer,
                                                                                                        ReduceDimension = ReduceDimension, Mode = Mode,
                                                                                                        FinalEnergy = FinalEnergy, alpha = alpha, beta = beta, gamma = gamma,
                                                                                                        drawAccuracyComplexity = Intermediate.drawAccuracyComplexity,
                                                                                                        drawPCAView = Intermediate.drawPCAView,
                                                                                                        drawEnergy = Intermediate.drawEnergy,
                                                                                                        n.cores = cl, ClusType = ClusType, MinParOP = MinParOP,
                                                                                                        FastSolve = FastSolve, AvoidSolitary = AvoidSolitary,
                                                                                                        EmbPointProb = EmbPointProb, AdjustElasticMatrix = AdjustElasticMatrix,
                                                                                                        AdjustElasticMatrix.Initial = AdjustElasticMatrix.Initial,
                                                                                                        Lambda.Initial = Lambda.Initial, Mu.Initial = Mu.Initial,
                                                                                                        age_vec = age_vec,
                                                                                                        eta = eta,
                                                                                                        ...)
                                
                                # Save extra information
                                ReturnList[[length(ReturnList)]]$SubSetID <- j
                                ReturnList[[length(ReturnList)]]$ReplicaID <- i
                                ReturnList[[length(ReturnList)]]$ProbPoint <- ProbPoint
                                
                                # Reset InitNodePositions for the next iteration
                                if(ComputeIC){
                                        InitNodePositions <- NULL
                                }
                                
                        }
                        
                }
                
                
                
                # Are we using bootstrapping (nRep > 1). If yes we compute the consensus tree
                if(nReps>1){
                        
                        print("Constructing average tree")
                        
                        # The nodes of the principal trees will be used as
                        # points to compute the consensus tree. The average
                        # ages, if age_vec selected, for these nodes will be
                        # used as age_vec
                        sel <- sapply(ReturnList, "[[", "SubSetID") == j
                        AllPoints <- do.call(rbind, lapply(ReturnList[sel], "[[", "NodePositions"))
                        AllAges   <- unlist(lapply(ReturnList[sel], function(x) x$age_tau$average_age))
                        
                        # De we need to compute the initial conditions?
                        if(is.null(InitNodePositions) | (is.null(InitEdges) & is.null(ElasticMatrix))){
                                
                                # areconstruct the initial configuration
                                InitialConf <- ElPiGraph.R:::generateInitialConfiguration(AllPoints, Nodes = InitNodes, Configuration = Configuration, DensityRadius = DensityRadius)
                                
                                # print(InitialConf)
                                
                                # Set the initial edge configuration
                                InitEdges <- InitialConf$Edges
                                
                                # Compute the initial elastic matrix
                                EM <- Encode2ElasticMatrix(Edges = InitialConf$Edges, Lambdas = Lambda, Mus = Mu)
                                
                                # Compute the initial node position
                                InitNodePositions <- PrimitiveElasticGraphEmbedment(
                                        X = X, NodePositions = InitialConf$NodePositions,
                                        MaxNumberOfIterations = MaxNumberOfIterations, TrimmingRadius = TrimmingRadius, eps = eps,
                                        ElasticMatrix = EM, Mode = Mode,
                                        age_vec = age_vec,
                                        eta = eta)$EmbeddedNodePositions
                        }
                        
                        
                        ReturnList[[length(ReturnList)+1]] <- computeElasticPrincipalGraph_edit(Data = AllPoints, NumNodes = NumNodes, NumEdges = NumEdges,
                                                                                                InitNodePositions = InitNodePositions, InitEdges = InitEdges, ElasticMatrix = ElasticMatrix,
                                                                                                AdjustVect = AdjustVect,
                                                                                                GrowGrammars = GrowGrammars,
                                                                                                ShrinkGrammars = ShrinkGrammars,
                                                                                                MaxNumberOfIterations = MaxNumberOfIterations, TrimmingRadius = TrimmingRadius, eps = eps,
                                                                                                Lambda = Lambda, Mu = Mu, Do_PCA = Do_PCA,
                                                                                                CenterData = CenterData, ComputeMSEP = ComputeMSEP,
                                                                                                verbose = verbose, ShowTimer = ShowTimer,
                                                                                                ReduceDimension = NULL, Mode = Mode,
                                                                                                FinalEnergy = FinalEnergy, alpha = alpha, beta = beta, gamma = gamma,
                                                                                                drawAccuracyComplexity = drawAccuracyComplexity,
                                                                                                drawPCAView = drawPCAView, drawEnergy = drawEnergy,
                                                                                                n.cores = cl, ClusType = ClusType, MinParOP = MinParOP,
                                                                                                FastSolve = FastSolve, AvoidSolitary = AvoidSolitary,
                                                                                                EmbPointProb = EmbPointProb, AdjustElasticMatrix = AdjustElasticMatrix,
                                                                                                AdjustElasticMatrix.Initial = AdjustElasticMatrix.Initial,
                                                                                                Lambda.Initial = Lambda.Initial, Mu.Initial = Lambda.Initial,
                                                                                                age_vec = AllAges,
                                                                                                eta = eta,
                                                                                                ...)
                        
                        # Run the ElPiGraph algorithm
                        ReturnList[[length(ReturnList)]]$SubSetID <- j
                        ReturnList[[length(ReturnList)]]$ReplicaID <- 0
                        ReturnList[[length(ReturnList)]]$ProbPoint <- 1
                        
                }
                
        }
        
        # if we created a cluster we need to stop it
        if(n.cores > 1){
                parallel::stopCluster(cl)
        }
        
        return(ReturnList)
        
}

# Edit computeElasticPrincipalGrah to make it call ElPrincGraph_edit function,
# instead of the native one.
computeElasticPrincipalGraph_edit <- function(Data,
                                              NumNodes,
                                              NumEdges = Inf,
                                              InitNodePositions,
                                              AdjustVect,
                                              InitEdges,
                                              ElasticMatrix = NULL,
                                              Lambda = 0.01,
                                              Mu = 0.1,
                                              MaxNumberOfIterations = 100,
                                              eps = .01,
                                              TrimmingRadius = Inf,
                                              Do_PCA = TRUE,
                                              CenterData = TRUE,
                                              ComputeMSEP = TRUE,
                                              verbose = FALSE,
                                              ShowTimer = FALSE,
                                              ReduceDimension = NULL,
                                              drawAccuracyComplexity = TRUE,
                                              drawPCAView = TRUE,
                                              drawEnergy = TRUE,
                                              n.cores = 1,
                                              ClusType = "Sock",
                                              MinParOP = 20,
                                              Mode = 1,
                                              FinalEnergy = "Base",
                                              alpha = 0,
                                              beta = 0,
                                              gamma = 0,
                                              GrowGrammars = list(),
                                              ShrinkGrammars = list(),
                                              GrammarOptimization = FALSE,
                                              MaxSteps = Inf,
                                              GrammarOrder = c("Grow", "Shrink"),
                                              FastSolve = FALSE,
                                              AvoidSolitary = FALSE,
                                              EmbPointProb = 1,
                                              AdjustElasticMatrix = NULL,
                                              AdjustElasticMatrix.Initial = NULL, 
                                              Lambda.Initial = NULL, Mu.Initial = NULL,
                                              age_vec = NULL,
                                              eta = 1e-3,
                                              ...) {
        ST <- date()
        tictoc::tic()
        
        if(is.null(ReduceDimension)){
                ReduceDimension <- 1:min(dim(Data))
        } else {
                
                if(!Do_PCA){
                        print("Cannot reduce dimensionality witout doing PCA.")
                        print("Dimensionality reduction will be ignored.")
                        ReduceDimension <- 1:min(dim(Data))
                }
                
        }
        
        if(CenterData){
                DataCenters <- colMeans(Data)
                
                Data <- t(t(Data) - DataCenters)
                InitNodePositions <- t(t(InitNodePositions) - DataCenters)
        }
        
        if(Do_PCA){
                
                print("Performing PCA on the data")
                
                ExpVariance <- sum(apply(Data, 2, var))
                
                if(length(ReduceDimension) == 1){
                        if(ReduceDimension < 1){
                                print("Dimensionality reduction via ratio of explained variance (full PCA will be computed)")
                                PCAData <- prcomp(Data, retx = TRUE, center = FALSE, scale. = FALSE)
                                
                                ReduceDimension <- 1:min(which(cumsum(PCAData$sdev^2)/ExpVariance >= ReduceDimension))
                        } else {
                                stop("if ReduceDimension is a single value it must be < 1")
                        }
                        
                } else {
                        
                        if(max(ReduceDimension) > min(dim(Data))){
                                print("Selected dimensions are outside of the available range. ReduceDimension will be updated")
                                ReduceDimension <- intersect(ReduceDimension, 1:min(dim(Data)))
                        }
                        
                        if(max(ReduceDimension) > min(dim(Data))*.75){
                                print("Using standard PCA")
                                PCAData <- prcomp(Data, retx = TRUE, center = FALSE, scale. = FALSE)
                        } else {
                                print("Using irlba PCA")
                                PCAData <- suppressWarnings(irlba::prcomp_irlba(Data, max(ReduceDimension), retx = TRUE, center = FALSE, scale. = FALSE))
                        }
                }
                
                ReduceDimension <- ReduceDimension[ReduceDimension <= ncol(PCAData$x)]
                
                perc <- 100*sum(PCAData$sdev[ReduceDimension]^2)/ExpVariance
                print(paste(length(ReduceDimension), "dimensions are being used"))
                print(paste0(signif(perc), "% of the original variance has been retained"))
                
                X <- PCAData$x[,ReduceDimension]
                InitNodePositions <- (InitNodePositions %*% PCAData$rotation)[,ReduceDimension]
                
        } else {
                X = Data
        }
        
        if(Do_PCA | CenterData){
                if(all(c("SOCKcluster", "cluster") %in% class(n.cores)) & ClusType != "Fork"){
                        print("Using a user supplied cluster. Updating the value of X")
                        parallel::clusterExport(n.cores, varlist = c("X"), envir=environment())
                }
        }
        
        if(is.null(Lambda.Initial)){
                Lambda.Initial <- Lambda
        }
        
        if(is.null(Mu.Initial)){
                Mu.Initial <- Mu
        } 
        
        if(is.null(ElasticMatrix)){
                InitElasticMatrix = ElPiGraph.R:::Encode2ElasticMatrix(Edges = InitEdges,
                                                                       Lambdas = Lambda.Initial,
                                                                       Mus = Mu.Initial)
        } else {
                print("The elastic matrix is being used. Edge configuration will be ignored")
                InitElasticMatrix = ElasticMatrix
        }
        
        if(nrow(InitElasticMatrix) != nrow(InitNodePositions) | ncol(InitElasticMatrix) != nrow(InitNodePositions)){
                stop("Elastic matrix incompatible with the node number. Impossible to proceed.")
        }
        
        
        
        # Computing the graph
        
        print(paste("Computing EPG with", NumNodes, "nodes on", nrow(X), "points and", ncol(X), "dimensions"))
        
        ElData <- ElPrincGraph_edit(X = X, NumNodes = NumNodes, NumEdges = NumEdges, Lambda = Lambda, Mu = Mu,
                                    MaxNumberOfIterations = MaxNumberOfIterations, eps = eps, TrimmingRadius = TrimmingRadius,
                                    NodesPositions = InitNodePositions, ElasticMatrix = InitElasticMatrix, AdjustVect = AdjustVect,
                                    CompileReport = TRUE, ShowTimer = ShowTimer,
                                    FinalEnergy = FinalEnergy, alpha = alpha, beta = beta, gamma = gamma, Mode = Mode,
                                    GrowGrammars = GrowGrammars, ShrinkGrammars = ShrinkGrammars,
                                    GrammarOptimization = GrammarOptimization, MaxSteps = MaxSteps, GrammarOrder = GrammarOrder,
                                    ComputeMSEP = ComputeMSEP, n.cores = n.cores, ClusType = ClusType, MinParOP = MinParOP,
                                    verbose = verbose, FastSolve = FastSolve, AvoidSolitary = AvoidSolitary,
                                    EmbPointProb = EmbPointProb, AdjustElasticMatrix = AdjustElasticMatrix,
                                    AdjustElasticMatrix.Initial = AdjustElasticMatrix.Initial, 
                                    age_vec = age_vec,
                                    eta = eta, ...)
        
        NodePositions <- ElData$NodePositions
        Edges <- DecodeElasticMatrix(ElData$ElasticMatrix)
        
        if(drawEnergy & !is.null(dim(ElData$ReportTable))){
                print(plotMSDEnergyPlot(ReportTable = ElData$ReportTable))
        }
        
        if(drawAccuracyComplexity & !is.null(dim(ElData$ReportTable))){
                print(accuracyComplexityPlot(ReportTable = ElData$ReportTable))
        }
        
        if(Do_PCA){
                NodePositions <- NodePositions %*% t(PCAData$rotation[, ReduceDimension])
        }
        
        EndTimer <- tictoc::toc()
        
        FinalPG <- list(NodePositions = NodePositions,
                        Edges = Edges,
                        ReportTable = ElData$ReportTable,
                        FinalReport = ElData$FinalReport,
                        ElasticMatrix = ElData$ElasticMatrix,
                        Lambda = ElData$Lambda,
                        Mu = ElData$Mu,
                        TrimmingRadius = TrimmingRadius,
                        FastSolve = ElData$FastSolve,
                        Mode = ElData$Mode,
                        MaxNumberOfIterations = ElData$MaxNumberOfIterations,
                        eps = ElData$eps,
                        Date = ST,
                        TicToc = EndTimer,
                        age_tau = ElData$age_tau,
                        g = ElData$g,
                        partitionData = ElData$partData)
        
        if(drawPCAView){
                p <- PlotPG(X = Data, TargetPG = FinalPG)
                print(p)
                
        }
        
        if(CenterData){
                NodePositions <- t(t(NodePositions) + DataCenters)
        }
        
        FinalPG$NodePositions <- NodePositions
        
        return(FinalPG)
        
}

# Edit calls to PrimitiveElasticGraphEmbedment to use the age-constraint one.
ApplyOptimalGraphGrammarOpeation_edit <- function(X,
                                                  NodePositions,
                                                  ElasticMatrix,
                                                  AdjustVect = NULL,
                                                  operationtypes,
                                                  SquaredX = NULL,
                                                  verbose = FALSE,
                                                  n.cores = 1,
                                                  EnvCl = NULL,
                                                  MinParOP = 20,
                                                  MaxNumberOfIterations = 100,
                                                  eps = .01,
                                                  TrimmingRadius = Inf,
                                                  Mode = 1,
                                                  FinalEnergy = "Base",
                                                  alpha = 1,
                                                  beta = 1,
                                                  gamma = 1,
                                                  FastSolve = FALSE,
                                                  EmbPointProb = 1,
                                                  AvoidSolitary = FALSE,
                                                  AdjustElasticMatrix = NULL,
                                                  age_vec = NULL,
                                                  eta = 1e-3,
                                                  ...) {
        
        if(is.null(AdjustVect)){
                AdjustVect <- rep(FALSE, nrow(NodePositions))
        }
        
        # % this function applies the most optimal graph grammar operation of operationtype
        # % the embedment of an elastic graph described by ElasticMatrix
        
        k=1
        
        if(is.null(SquaredX)){
                SquaredX = rowSums(X^2)
        }
        
        NodePositionArrayAll <- list()
        ElasticMatricesAll <- list()
        AdjustVectAll <- list()
        
        Partition = ElPiGraph.R::PartitionData(X = X,
                                               NodePositions = NodePositions,
                                               SquaredX = SquaredX,
                                               TrimmingRadius = TrimmingRadius)$Partition
        
        for(i in 1:length(operationtypes)){
                if(verbose){
                        tictoc::tic()
                        print(paste(i, 'Operation type =', operationtypes[i]))
                }
                
                NewMatrices <- ElPiGraph.R:::GraphGrammarOperation(X = X, NodePositions = NodePositions,
                                                                   ElasticMatrix = ElasticMatrix, type = operationtypes[i],
                                                                   Partition = Partition, AdjustVect = AdjustVect)
                
                NodePositionArrayAll <- c(NodePositionArrayAll,
                                          NewMatrices$NodePositionArray)
                ElasticMatricesAll <- c(ElasticMatricesAll,
                                        NewMatrices$ElasticMatrices)
                AdjustVectAll <- c(AdjustVectAll,
                                   NewMatrices$AdjustVect)
                
                if(verbose){
                        tictoc::toc()
                }
                
        }
        
        # minEnergy = .Machine$double.xmax
        # k = -1
        
        if(verbose){
                tictoc::tic()
                print("Optimizing graphs")
        }
        
        Valid <- as.list(1:length(NodePositionArrayAll))
        
        # Check that each point is associated with at least one point. Otherwise we exclude the configuration
        if(AvoidSolitary){
                Valid <- sapply(Valid, function(i){
                        Partition <- ElPiGraph.R::PartitionData(X = X,
                                                                NodePositions = NodePositionArrayAll[[i]],
                                                                SquaredX = SquaredX,
                                                                TrimmingRadius = TrimmingRadius)$Partition
                        if(all(1:nrow(NodePositionArrayAll[[i]]) %in% Partition)){
                                return(i)
                        }
                        return(0)
                })
                
                if(verbose){
                        tictoc::tic()
                        print(paste0(sum(Valid > 0), "configurations out of", length(Valid), "used"))
                }
                
                Valid <- Valid[Valid > 0]
        }
        
        CombinedInfo <- lapply(Valid, function(i){
                list(NodePositions = NodePositionArrayAll[[i]],
                     ElasticMatrix = ElasticMatricesAll[[i]],
                     AdjustVect = AdjustVectAll[[i]])
        })
        
        # We are adjusting the elastic matrix
        if(!is.null(AdjustElasticMatrix)){
                CombinedInfo <- lapply(CombinedInfo, AdjustElasticMatrix, ...)
        }
        
        # print(paste("DEBUG:", TrimmingRadius))
        
        DynamicProcess <- length(CombinedInfo) %/% MinParOP + 1
        if(DynamicProcess > n.cores){
                DynamicProcess <- n.cores
        }
        
        if((n.cores > 1) & (DynamicProcess > 1) ){
                
                if(is.null(EnvCl)){
                        cl <- parallel::makeCluster(n.cores)
                        parallel::clusterExport(cl, varlist = c("X"), envir = environment())
                } else {
                        cl <- EnvCl
                }
                
                Embed <- parallel::parLapply(cl[1:DynamicProcess], CombinedInfo, function(input){
                        PrimitiveElasticGraphEmbedment_edit(X = X,
                                                            NodePositions = input$NodePositions,
                                                            ElasticMatrix = input$ElasticMatrix,
                                                            SquaredX = SquaredX,
                                                            verbose = FALSE,
                                                            MaxNumberOfIterations = MaxNumberOfIterations,
                                                            eps = eps,
                                                            FinalEnergy = FinalEnergy,
                                                            alpha = alpha,
                                                            beta = beta,
                                                            gamma = gamma,
                                                            Mode = Mode,
                                                            TrimmingRadius = TrimmingRadius,
                                                            FastSolve = FastSolve,
                                                            prob = EmbPointProb,
                                                            age_vec = age_vec,
                                                            eta = eta)
                })
                
                if(is.null(EnvCl)){
                        parallel::stopCluster(cl)
                }
                
        } else {
                Embed <- lapply(CombinedInfo, function(input){
                        PrimitiveElasticGraphEmbedment_edit(X = X,
                                                            NodePositions = input$NodePositions,
                                                            ElasticMatrix = input$ElasticMatrix,
                                                            SquaredX = SquaredX,
                                                            verbose = FALSE,
                                                            MaxNumberOfIterations = MaxNumberOfIterations,
                                                            eps = eps,
                                                            FinalEnergy = FinalEnergy,
                                                            alpha = alpha,
                                                            beta = beta,
                                                            Mode = Mode,
                                                            TrimmingRadius = TrimmingRadius,
                                                            FastSolve = FastSolve,
                                                            prob = EmbPointProb,
                                                            age_vec = age_vec,
                                                            eta = eta)
                })
        }
        
        if(length(Embed)==0){
                return(NA)
        }
        
        Best <- which.min(sapply(Embed, "[[", "ElasticEnergy"))
        
        minEnergy <- Embed[[Best]]$ElasticEnergy
        NodePositions2 <- Embed[[Best]]$EmbeddedNodePositions
        ElasticMatrix2 <- CombinedInfo[[Best]]$ElasticMatrix
        AdjustVect2 <- CombinedInfo[[Best]]$AdjustVect
        
        if(verbose){
                tictoc::toc()
        }
        
        return(
                list(
                        NodePositions = NodePositions2,
                        ElasticMatrix = ElasticMatrix2,
                        ElasticEnergy = minEnergy,
                        MSE = Embed[[Best]]$MSE,
                        EP = Embed[[Best]]$EP,
                        RP = Embed[[Best]]$RP,
                        AdjustVect = AdjustVect2
                )
        )
        
}

# Add age_vec as an argument and remove references to package calls
ElPrincGraph_edit <- function(X,
                              NumNodes = 100,
                              NumEdges = Inf,
                              Lambda,
                              Mu,
                              ElasticMatrix,
                              NodesPositions,
                              AdjustVect,
                              verbose = FALSE,
                              n.cores = 1,
                              ClusType = "Sock",
                              MinParOP = 20,
                              CompileReport = FALSE,
                              ShowTimer = FALSE,
                              ComputeMSEP = TRUE,
                              FinalEnergy = "Base",
                              alpha = 0,
                              beta = 0,
                              gamma = 0,
                              Mode = 1,
                              MaxNumberOfIterations = 10,
                              MaxFailedOperations = Inf,
                              MaxSteps = Inf,
                              GrammarOptimization = FALSE,
                              eps = .01,
                              TrimmingRadius = Inf,
                              GrowGrammars = list(),
                              ShrinkGrammars = list(),
                              GrammarOrder = c("Grow", "Shrink"),
                              FastSolve = FALSE,
                              AvoidSolitary = FALSE,
                              EmbPointProb = 1,
                              AdjustElasticMatrix = NULL,
                              AdjustElasticMatrix.Initial = NULL,
                              age_vec = NULL,
                              eta = 1e-3,
                              ...) {
        if (!is.null(age_vec)){
                print(sprintf("Imposing aging constraint with Eta = %s", eta))
        }
        if(GrammarOptimization){
                print("Using grammar optimization")
                if(is.infinite(MaxSteps)){
                        warning("When setting GrammarOptimization to TRUE, MaxSteps must be finite. Using MaxSteps = 1")
                        MaxSteps = 1
                }
        }
        
        if(is.list(X)){
                warning("Data matrix must be a numeric matrix. It will be converted automatically. This can introduce inconsistencies")
                X <- data.matrix(X)
        }
        
        if(!CompileReport){
                verbose = FALSE
        }
        
        if(!is.null(ElasticMatrix)){
                if(any(ElasticMatrix != t(ElasticMatrix))){
                        stop('ERROR: Elastic matrix must be square and symmetric')
                }
        }
        
        if(verbose){
                cat('BARCODE\tENERGY\tNNODES\tNEDGES\tNRIBS\tNSTARS\tNRAYS\tNRAYS2\tMSE\tMSEP\tFVE\tFVEP\tUE\tUR\tURN\tURN2\tURSD\n')
        }
        
        if(!is.null(AdjustElasticMatrix.Initial)){
                tGraphInfo <- list(ElasticMatrix = ElasticMatrix, AdjustVect = AdjustVect)
                ElasticMatrix <- AdjustElasticMatrix.Initial(tGraphInfo, ...)$ElasticMatrix
                
                print(paste(sum(ElasticMatrix != tGraphInfo$ElasticMatrix), "values of the elastic matrix have been updated"))
        }
        
        InitNodePositions <- PrimitiveElasticGraphEmbedment_edit(
                X = X, NodePositions = NodesPositions,
                MaxNumberOfIterations = MaxNumberOfIterations, TrimmingRadius = TrimmingRadius, eps = eps,
                ElasticMatrix = ElasticMatrix, Mode = Mode,
                age_vec = age_vec,
                eta = eta)$EmbeddedNodePositions
        
        UpdatedPG <- list(ElasticMatrix = ElasticMatrix,
                          NodePositions = InitNodePositions,
                          AdjustVect = AdjustVect)
        
        ReportTable <- NULL
        ToSrink <- c(2, 9, 10, 11, 12, 13, 14, 15, 16, 17)
        
        SquaredX = rowSums(X^2)
        
        # now we grow the graph up to NumNodes
        
        GlobalCluster <- TRUE
        
        if(all(class(n.cores) %in% c("numeric", "integer"))){
                if(n.cores > 1){
                        if(ClusType == "Fork"){
                                print(paste("Creating a fork cluster with",
                                            n.cores,
                                            "nodes"))
                                cl <- parallel::makeCluster(n.cores,
                                                            type="FORK")
                                GlobalCluster <- FALSE
                        } else {
                                print(paste("Creating a sock cluster with", n.cores, "nodes"))
                                cl <- parallel::makeCluster(n.cores)
                                GlobalCluster <- FALSE
                                parallel::clusterExport(cl,
                                                        varlist = c("X",
                                                                    "SquaredX",
                                                                    "MaxNumberOfIterations",
                                                                    "TrimmingRadius",
                                                                    "eps",
                                                                    "verbose",
                                                                    "EmbPointProb",
                                                                    "alpha",
                                                                    "beta",
                                                                    "FinalEnergy"), envir=environment())
                        }
                } else {
                        print("Using a single core")
                        cl <- NULL
                }
        } else {
                if(all(c("SOCKcluster", "cluster") %in% class(n.cores))){
                        print("Using a user supplied cluster. It must contains the data points in a matrix X")
                        cl <- n.cores
                        CheckX <- unlist(parallel::clusterCall(cl, function(){exists("X")}))
                        if(all(CheckX)){
                                GlobalCluster <- TRUE
                                if(ClusType != "Fork"){
                                        print("Exporting the additional variables to the cluster")
                                        parallel::clusterExport(cl, varlist = c("SquaredX", "MaxNumberOfIterations", "TrimmingRadius", "eps", "verbose",
                                                                                "EmbPointProb", "alpha", "beta", "FinalEnergy"), envir=environment())
                                }
                                n.cores = length(CheckX)
                                
                        } else {
                                print("Unable to find X on the cluster. Single processor computation will be used")
                                n.cores = 1
                        }
                }
        }
        
        if(nrow(UpdatedPG$NodePositions) >= NumNodes & !GrammarOptimization){
                
                FinalReport <- ReportOnPrimitiveGraphEmbedment_edit(X = X, NodePositions = UpdatedPG$NodePositions,
                                                                    ElasticMatrix = UpdatedPG$ElasticMatrix,
                                                                    PartData = PartitionData(X = X,
                                                                                             NodePositions = UpdatedPG$NodePositions,
                                                                                             SquaredX = SquaredX,
                                                                                             TrimmingRadius = TrimmingRadius,
                                                                                             nCores = 1),
                                                                    ComputeMSEP = ComputeMSEP,
                                                                    age_vec = age_vec,
                                                                    eta = eta)
                
                return(
                        list(NodePositions = UpdatedPG$NodePositions, ElasticMatrix = UpdatedPG$ElasticMatrix,
                             ReportTable = unlist(FinalReport), FinalReport = FinalReport, Lambda = Lambda, Mu = Mu,
                             FastSolve = FastSolve)
                )
        }
        
        StartNodes <- nrow(UpdatedPG$NodePositions)
        
        # print(FinalReport)
        
        FailedOperations <- 0
        Steps <- 0
        FirstPrint <- TRUE
        
        while((nrow(UpdatedPG$NodePositions) < NumNodes) | GrammarOptimization){
                
                nEdges <- sum(UpdatedPG$ElasticMatrix[lower.tri(UpdatedPG$ElasticMatrix, diag = FALSE)] > 0)
                
                if((nrow(UpdatedPG$NodePositions) >= NumNodes | nEdges >= NumEdges) & !GrammarOptimization){
                        break()
                }
                
                if(!verbose & ShowTimer){
                        print(paste("Nodes = ", nrow(UpdatedPG$NodePositions)))
                }
                
                if(!verbose & !ShowTimer){
                        if(FirstPrint){
                                cat("Nodes = ")
                                FirstPrint <- FALSE
                        }
                        cat(nrow(UpdatedPG$NodePositions))
                        cat(" ")
                }
                
                OldPG <- UpdatedPG
                
                for(OpType in GrammarOrder){
                        
                        if(OpType == "Grow" & length(GrowGrammars)>0){
                                for(k in 1:length(GrowGrammars)){
                                        if(ShowTimer){
                                                print("Growing")
                                                tictoc::tic()
                                        }
                                        
                                        UpdatedPG <- ApplyOptimalGraphGrammarOpeation_edit(X = X,
                                                                                           NodePositions = UpdatedPG$NodePositions,
                                                                                           ElasticMatrix = UpdatedPG$ElasticMatrix,
                                                                                           AdjustVect = UpdatedPG$AdjustVect,
                                                                                           operationtypes = GrowGrammars[[k]],
                                                                                           SquaredX = SquaredX,
                                                                                           FinalEnergy = FinalEnergy,
                                                                                           alpha = alpha,
                                                                                           beta = beta,
                                                                                           gamma = gamma,
                                                                                           Mode = Mode,
                                                                                           MaxNumberOfIterations = MaxNumberOfIterations,
                                                                                           eps = eps,
                                                                                           TrimmingRadius = TrimmingRadius,
                                                                                           verbose = FALSE,
                                                                                           n.cores = n.cores,
                                                                                           EnvCl = cl,
                                                                                           MinParOP = MinParOP,
                                                                                           FastSolve = FastSolve,
                                                                                           AvoidSolitary = AvoidSolitary,
                                                                                           EmbPointProb = EmbPointProb,
                                                                                           AdjustElasticMatrix = AdjustElasticMatrix,
                                                                                           age_vec = age_vec,
                                                                                           eta = eta,
                                                                                           ...)
                                        
                                        
                                        if(!is.list(UpdatedPG)){
                                                
                                                FailedOperations <- FailedOperations + 1
                                                UpdatedPG <- OldPG
                                                break()
                                                
                                        } else {
                                                
                                                FailedOperations <- 0
                                                
                                                if(nrow(UpdatedPG$NodePositions) == 3){
                                                        # this is needed to erase the star elasticity coefficient which was initially assigned to both leaf nodes,
                                                        # one can erase this information after the number of nodes in the graph is > 2
                                                        
                                                        inds = which(colSums(UpdatedPG$ElasticMatrix-diag(diag(UpdatedPG$ElasticMatrix))>0)==1)
                                                        
                                                        UpdatedPG$ElasticMatrix[inds, inds] <- 0
                                                }
                                                
                                        }
                                        
                                        if(ShowTimer){
                                                tictoc::toc()
                                        }
                                        
                                }
                        }
                        
                        
                        if(OpType == "Shrink" & length(ShrinkGrammars)>0){
                                for(k in 1:length(ShrinkGrammars)){
                                        
                                        if(ShowTimer){
                                                print("Shrinking")
                                                tictoc::tic()
                                        }
                                        
                                        UpdatedPG <- ApplyOptimalGraphGrammarOpeation_edit(X = X,
                                                                                           NodePositions = UpdatedPG$NodePositions,
                                                                                           ElasticMatrix = UpdatedPG$ElasticMatrix,
                                                                                           AdjustVect = UpdatedPG$AdjustVect,
                                                                                           operationtypes = ShrinkGrammars[[k]],
                                                                                           SquaredX = SquaredX,
                                                                                           Mode = Mode,
                                                                                           FinalEnergy = FinalEnergy,
                                                                                           alpha = alpha,
                                                                                           beta = beta,
                                                                                           gamma = gamma,
                                                                                           MaxNumberOfIterations = MaxNumberOfIterations,
                                                                                           eps = eps,
                                                                                           TrimmingRadius = TrimmingRadius,
                                                                                           verbose = FALSE,
                                                                                           n.cores = n.cores,
                                                                                           MinParOP = MinParOP,
                                                                                           EnvCl = cl,
                                                                                           FastSolve = FastSolve,
                                                                                           AvoidSolitary = AvoidSolitary,
                                                                                           EmbPointProb = EmbPointProb,
                                                                                           AdjustElasticMatrix = AdjustElasticMatrix,
                                                                                           age_vec = age_vec,
                                                                                           eta = eta,
                                                                                           ...)
                                        
                                        
                                        if(!is.list(UpdatedPG)){
                                                
                                                FailedOperations <- FailedOperations + 1
                                                UpdatedPG <- OldPG
                                                break()
                                                
                                        } else {
                                                
                                                FailedOperations <- 0
                                                
                                        }
                                        
                                        
                                        if(ShowTimer){
                                                tictoc::toc()
                                        }
                                        
                                }
                        }
                        
                }
                
                if(CompileReport){
                        tReport <- ReportOnPrimitiveGraphEmbedment_edit(X = X, NodePositions = UpdatedPG$NodePositions,
                                                                        ElasticMatrix = UpdatedPG$ElasticMatrix,
                                                                        PartData = PartitionData(X = X,
                                                                                                 NodePositions = UpdatedPG$NodePositions,
                                                                                                 SquaredX = SquaredX,
                                                                                                 TrimmingRadius = TrimmingRadius,
                                                                                                 nCores = 1),
                                                                        ComputeMSEP = ComputeMSEP,
                                                                        age_vec,
                                                                        eta)
                        FinalReport <- tReport
                        tReport <- unlist(tReport)
                        tReport[ToSrink] <- sapply(tReport[ToSrink], function(x) {
                                signif(as.numeric(x), 4)
                        })
                        
                        ReportTable <- rbind(ReportTable, tReport)
                        
                        if(verbose){
                                cat(ReportTable[nrow(ReportTable), ], sep = "\t")
                                cat("\n")
                        }
                }
                
                # Count the execution steps
                Steps <- Steps + 1
                
                # If the number of execution steps is larger than MaxSteps stop the algorithm
                if(Steps > MaxSteps | FailedOperations > MaxFailedOperations){
                        break()
                }
                
        }
        
        # FinalReport <- NULL
        
        if(!verbose){
                if(!CompileReport){
                        tReport <- ReportOnPrimitiveGraphEmbedment_edit(X = X, NodePositions = UpdatedPG$NodePositions,
                                                                        ElasticMatrix = UpdatedPG$ElasticMatrix,
                                                                        PartData = PartitionData(X = X,
                                                                                                 NodePositions = UpdatedPG$NodePositions,
                                                                                                 SquaredX = SquaredX,
                                                                                                 TrimmingRadius = TrimmingRadius,
                                                                                                 nCores = 1),
                                                                        ComputeMSEP = ComputeMSEP,
                                                                        age_vec,
                                                                        eta)
                        FinalReport <- tReport
                        tReport <- unlist(tReport)
                        tReport[ToSrink] <- sapply(tReport[ToSrink], function(x) {
                                signif(as.numeric(x), 4)
                        })
                } else {
                        tReport <- ReportTable[nrow(ReportTable),]
                        tReport <- unlist(tReport)
                }
                
                cat("\n")
                cat('BARCODE\tENERGY\tNNODES\tNEDGES\tNRIBS\tNSTARS\tNRAYS\tNRAYS2\tMSE\tMSEP\tFVE\tFVEP\tUE\tUR\tURN\tURN2\tURSD\n')
                cat(tReport, sep = "\t")
                cat("\n")
        }
        
        # ReportTable <- rbind(ReportTable, tReport)
        
        if(!GlobalCluster){
                print("Stopping the cluster")
                parallel::stopCluster(cl)
        }
        
        # print(ReportTable)
        
        if(is.list(ReportTable)){
                ReportTable <- unlist(ReportTable)
        }
        
        # print(ReportTable)
        
        if(is.null(dim(ReportTable)) & !is.null(ReportTable)){
                RPNames <- names(ReportTable)
                ReportTable <- matrix(ReportTable, nrow = 1)
                colnames(ReportTable) <- RPNames
        }
        
        # print(ReportTable)
        if (!is.null(age_vec)){
                # Obtain partition data
                K <- nrow(UpdatedPG$NodePositions)
                partData <- PartitionData(X = X,
                                          NodePositions = UpdatedPG$NodePositions,
                                          SquaredX = SquaredX,
                                          TrimmingRadius = TrimmingRadius,
                                          nCores = 1)
                assignment <- partData$Partition
                mean_age_j <- get_node_meanage(UpdatedPG$NodePositions,
                                               assignment,
                                               age_vec)
                
                # Obtain final graph
                g <- get_node_graph(UpdatedPG$NodePositions,
                                    UpdatedPG$ElasticMatrix)
                
                # Obtain final pseudotimes
                tau_j <- get_node_pseudotime(UpdatedPG$NodePositions,
                                             mean_age_j,
                                             UpdatedPG$ElasticMatrix)
                
                
                #print(names(mean_list))
                #print(as.numeric(mean_list))
                #print(tau_j)
                out <- list(NodePositions = UpdatedPG$NodePositions, ElasticMatrix = UpdatedPG$ElasticMatrix,
                            ReportTable = ReportTable, FinalReport = FinalReport, Lambda = Lambda, Mu = Mu,
                            FastSolve = FastSolve, Mode = Mode, MaxNumberOfIterations = MaxNumberOfIterations,
                            eps = eps,
                            partData = partData,
                            g = g,
                            age_tau = data.frame(node = 1:K,
                                                 average_age = mean_age_j,
                                                 pseudotime = tau_j))
        }else{
                out <- list(NodePositions = UpdatedPG$NodePositions, ElasticMatrix = UpdatedPG$ElasticMatrix,
                            ReportTable = ReportTable, FinalReport = FinalReport, Lambda = Lambda, Mu = Mu,
                            FastSolve = FastSolve, Mode = Mode, MaxNumberOfIterations = MaxNumberOfIterations,
                            eps = eps,
                            partData = PartitionData(X = X,
                                                     NodePositions = UpdatedPG$NodePositions,
                                                     SquaredX = SquaredX,
                                                     TrimmingRadius = TrimmingRadius,
                                                     nCores = 1))
        }
        
        return(out)
        
}

ReportOnPrimitiveGraphEmbedment_edit <- function(X,
                                                 NodePositions,
                                                 ElasticMatrix,
                                                 PartData=NULL,
                                                 ComputeMSEP = FALSE,
                                                 age_vec,
                                                 eta) {
        
        # %   This function computes various measurements concerning a primitive
        # %   graph embedment
        # %
        # %           BARCODE is barcode in form ...S4|S3||N, where N is number of
        # %               nodes, S3 is number of 3-stars, S4 (S5,...) is number of
        # %               four (five,...) stars, etc.
        # %           ENERGY is total elastic energy of graph embedment (ENERGY = MSE + UE +
        #                                                                  %               UR)
        # %           NNODES is number of nodes.
        # %           NEDGES is number of edges
        # %           NRIBS is number of two stars (nodes with two otherr connected
        #                                           %               nodes).
        # %           NSTARS is number of stars with 3 and more leaves (nodes
        #                                                               %               connected with central node).
        # %           NRAYS2 is sum of rays minus doubled number of nodes.
        # %           MSE is mean square error or assessment of data approximation
        # %               quality.
        # %           MSEP is mean square error after piece-wise linear projection on the edges
        # %           FVE is fraction of explained variance. This value always
        # %               between 0 and 1. Greater value means higher quality of
        # %               data approximation.
        # %           FVEP is same as FVE but computed after piece-wise linear projection on the edges
        # %           UE is total sum of squared edge lengths.
        # %           UR is total sum of star deviations from harmonicity.
        # %           URN is UR * nodes
        # %           URN2 is UR * nodes^2
        # %           URSD is standard deviation of UR???
        
        Mu = diag(ElasticMatrix)
        L = ElasticMatrix - diag(Mu)
        Connectivities <- colSums(L>0)
        N <- table(factor(Connectivities, levels = 1:max(Connectivities)))
        DecodedMat <- DecodeElasticMatrix(ElasticMatrix)
        
        TotalVariance = sum(apply(X, 2, var))
        
        BARCODE = getPrimitiveGraphStructureBarCode(ElasticMatrix)
        
        if(is.null(PartData)){
                PartData <- PartitionData(X = X, NodePositions = NodePositions, rowSums(X^2))
        }
        
        target_position <- get_target_positions(NodePositions,
                                                PartData$Partition,
                                                age_vec,
                                                ElasticMatrix)
        Energies <- distutils::ElasticEnergy(X = X, NodePositions = NodePositions,
                                             ElasticMatrix = ElasticMatrix,
                                             Dists = PartData$Dists,
                                             partition = PartData$Partition,
                                             NodeAgeTargets = target_position,
                                             eta = eta)
        
        NNODES = nrow(NodePositions)
        NEDGES = nrow(DecodedMat$Edges)
        
        
        if(length(N)>1){
                NRIBS = N[2]
        } else{
                NRIBS = 0
        }
        
        if(length(N)>2){
                NSTARS = N[3]
        } else {
                NSTARS = 0
        }
        
        NRAYS = 0
        NRAYS2 = 0
        
        
        if(ComputeMSEP){
                NodeProj <- project_point_onto_graph(X, NodePositions = NodePositions,
                                                     Edges = DecodedMat$Edges, Partition = PartData$Partition)
                MSEP = NodeProj$MSEP
                FVEP = (TotalVariance-MSEP)/TotalVariance
        } else {
                MSEP = NA
                FVEP = NA
        }
        
        FVE = (TotalVariance-Energies$MSE)/TotalVariance
        URN = Energies$RP*NNODES
        URN2 = Energies$RP*NNODES*NNODES
        URSD = 0
        
        return(list(BARCODE = BARCODE, ENERGY = Energies$ElasticEnergy, NNODES = NNODES, NEDGES = NEDGES,
                    NRIBS = NRIBS, NSTARS = NSTARS, NRAYS = NRAYS, NRAYS2 = NRAYS2,
                    MSE = Energies$MSE, MSEP = MSEP, FVE = FVE, FVEP = FVEP, UE = Energies$EP, UR = Energies$RP,
                    URN = URN, URN2 = URN2, URSD = URSD))
        
        
}

# Given a node matrix, a partition vector and an age vector, compute mean
# age per node.
get_node_meanage <- function(nodes_mat, assignment, age_vec){
        K <- nrow(nodes_mat)
        
        ## counts and means
        fac <- factor(as.integer(assignment), levels = seq_len(K))
        n_j_tab <- table(fac)
        n_j <- integer(K)
        n_j[as.integer(names(n_j_tab))] <- as.integer(n_j_tab)
        mean_age_j <- rep(NA_real_, K)
        if(length(n_j_tab)>0){
                mean_list <- tapply(age_vec, assignment, mean)
                if(length(mean_list)>0){
                        idxs <- as.integer(names(mean_list))
                        mean_age_j[idxs] <- as.numeric(mean_list)
                }
        }
        return(mean_age_j)
}

# Given the elastic matrix and the nodes matrix, obtain an igraph graph object.
get_node_graph <- function(nodes_mat, ElasticMatrix){
        # Obtain decoded elastic matrix and
        # extract edges.
        dec_elMat <- ElPiGraph.R::DecodeElasticMatrix(ElasticMatrix)
        edgs <- dec_elMat$Edges
        # Obtain weights based on Euclidean
        # dists
        edge_dists <- sqrt(rowSums((nodes_mat[edgs[, 2], , drop = F] - nodes_mat[edgs[, 1], , drop = F])^2))
        edge_weights <- 1 / (1 + edge_dists)
        g <- igraph::graph_from_edgelist(edgs, directed = F)
        E(g)$weight <- edge_weights
        return (g)
}

# Computes pseudotime of nodes by extracting the graph, starting from the
# leaf with the lowest mean age.
get_node_pseudotime <- function(nodes_mat, mean_age, ElasticMatrix){
        K <- nrow(nodes_mat)
        tau_j <- rep(NA_real_, K)
        if(K > 1){
                g <- get_node_graph(nodes_mat, ElasticMatrix)
                # identify leaves and assign root to 
                # leaf with minimum age
                degs <- igraph::degree(g)
                leaves <- which(degs == 1)
                root <- NA_integer_
                if(length(leaves) > 0 && any(!is.na(mean_age[leaves]))){
                        root <- leaves[which.min(mean_age[leaves])]
                } else if(any(!is.na(mean_age))){
                        root <- which.min(mean_age)
                } else {
                        root <- 1L
                }
                dist_mat_nodes <- as.numeric(igraph::distances(g, v = root, to = V(g)))
                tau_j <- dist_mat_nodes
        } else {
                tau_j <- 0
        }
        return(tau_j)
}

# Get target node positions based on an isotonic regression on age vs pseudotime
get_target_positions <- function(nodes_mat, assignment, age_vec, ElasticMatrix){
        K <- nrow(nodes_mat)
        D <- ncol(nodes_mat)
        
        mean_age_j <- get_node_meanage(nodes_mat, assignment, age_vec)
        
        tau_j <- get_node_pseudotime(nodes_mat, mean_age_j, ElasticMatrix)
        
        
        ok <- !is.na(mean_age_j) & !is.na(tau_j)
        if(sum(ok) >= 2){
                iso_fit <- stats::isoreg(mean_age_j[ok], tau_j[ok])
                r_ok <- iso_fit$yf
                r_j <- rep(NA_real_, K)
                r_j[which(ok)] <- r_ok
                r_j[is.na(r_j)] <- tau_j[is.na(r_j)]
        } else {
                r_j <- tau_j
        }
        
        # interpolate per dimension
        ord <- order(tau_j)
        tau_ord <- tau_j[ord]
        nodes_ord <- nodes_mat[ord, , drop = FALSE]
        pos_target <- matrix(NA_real_, nrow = K, ncol = D)
        for(dimc in seq_len(D)){
                col_vals <- nodes_ord[, dimc]
                interp_vals <- stats::approx(x = tau_ord, y = col_vals, xout = r_j, rule = 2)$y
                pos_target[, dimc] <- interp_vals[order(order(tau_j))]  # keep original index order
        }
        return (pos_target)
}


PrimitiveElasticGraphEmbedment_edit <- function(X,
                                                NodePositions,
                                                ElasticMatrix,
                                                MaxNumberOfIterations,
                                                TrimmingRadius,
                                                eps,
                                                Mode = 1,
                                                FinalEnergy = "Base",
                                                SquaredX = NULL,
                                                verbose = FALSE,
                                                FastSolve = FALSE,
                                                DisplayWarnings = FALSE,
                                                alpha = 0,
                                                beta = 0,
                                                gamma = 0,
                                                prob = 1,
                                                ## NEW:
                                                age_vec = NULL,
                                                eta = 1e-3) {
        
        # Minimal argument checks
        if(!is.matrix(X)) X <- data.matrix(X)
        if(!is.matrix(NodePositions)) NodePositions <- data.matrix(NodePositions)
        N <- nrow(X)
        PointWeights <- rep(1, N)
        K_init <- nrow(NodePositions)
        
        # Auxiliary computations
        SpringLaplacianMatrix <- ComputeSpringLaplacianMatrix(ElasticMatrix)
        
        # Prepare squared norms
        if(is.null(SquaredX)){
                SquaredX <- rowSums(X^2)
        }
        
        # Initial partition
        PartDataStruct <- PartitionData(X = X, NodePositions = NodePositions,
                                        SquaredX = SquaredX,
                                        TrimmingRadius = TrimmingRadius)
        assignment <- PartDataStruct$Partition
        
        target_positions <- get_target_positions(NodePositions,
                                                 assignment,
                                                 age_vec,
                                                 ElasticMatrix)
        
        if(verbose | Mode == 2){
                OldPriGrElEn <-
                        distutils::ElasticEnergy(X = X,
                                                 NodePositions = NodePositions,
                                                 ElasticMatrix = ElasticMatrix,
                                                 Dists = PartDataStruct$Dists,
                                                 partition = assignment,
                                                 NodeAgeTargets = target_positions,
                                                 eta = eta)
        } else {
                OldPriGrElEn <- list(ElasticEnergy = NA, MSE = NA, EP = NA, RP = NA)
        }
        
        Converged <- FALSE
        NewNodePositions <- NodePositions
        
        for(i in 1:MaxNumberOfIterations){
                
                # if MaxNumberOfIterations == 0: do one fit and exit
                if(MaxNumberOfIterations == 0){
                        Converged <- TRUE
                        PriGrElEn <- OldPriGrElEn
                        difference <- NA
                        
                        NewNodePositions <- distutils::FitGraph2DataGivenPartition_Age(X = X,
                                                                                   PointWeights = PointWeights,
                                                                                   NodePositions = NodePositions,
                                                                                   SpringLaplacianMatrix = SpringLaplacianMatrix,
                                                                                   partition = PartDataStruct$Partition,
                                                                                   FastSolve = FastSolve,
                                                                                   NodeAgeTargets = target_positions,
                                                                                   eta = eta)
                        break()
                } # end MaxNumberOfIterations==0 case
                
                # Updated positions (possibly subsampled)
                if(prob < 1){
                        SelPoint <- runif(nrow(X)) < prob
                        while(sum(SelPoint) < min(3, nrow(X))){
                                SelPoint[ceiling(runif(1, min = 1, max = nrow(X)))] <- TRUE
                        }
                        NewNodePositions <- distutils::FitGraph2DataGivenPartition_Age(X = X[SelPoint, ],
                                                                                       PointWeights = PointWeights,
                                                                                       NodePositions = NodePositions,
                                                                                       SpringLaplacianMatrix = SpringLaplacianMatrix,
                                                                                       partition = PartDataStruct$Partition[SelPoint],
                                                                                       FastSolve = FastSolve,
                                                                                       NodeAgeTargets = target_positions,
                                                                                       eta = eta)
                } else {
                        NewNodePositions <- distutils::FitGraph2DataGivenPartition_Age(X = X,
                                                                                       PointWeights = PointWeights,
                                                                                       NodePositions = NodePositions,
                                                                                       SpringLaplacianMatrix = SpringLaplacianMatrix,
                                                                                       partition = PartDataStruct$Partition,
                                                                                       FastSolve = FastSolve,
                                                                                       NodeAgeTargets = target_positions,
                                                                                       eta = eta)
                }
                
                PartDataStruct <- PartitionData(X = X, NodePositions = NewNodePositions, SquaredX = SquaredX,
                                                TrimmingRadius = TrimmingRadius)
                assignment <- PartDataStruct$Partition
                
                target_positions <- get_target_positions(NewNodePositions,
                                                         assignment,
                                                         age_vec,
                                                         ElasticMatrix)
                # Compute energy if requested
                if(verbose | Mode == 2){
                        PriGrElEn <- distutils::ElasticEnergy(X = X,
                                                              NodePositions = NewNodePositions,
                                                              ElasticMatrix =  ElasticMatrix,
                                                              Dists = PartDataStruct$Dists,
                                                              partition = assignment,
                                                              NodeAgeTargets = target_positions,
                                                              eta = eta)
                } else {
                        PriGrElEn <- list(ElasticEnergy = NA, MSE = NA, EP = NA, RP = NA)
                }
                
                # compute difference
                if(Mode == 1){
                        difference <- ComputeRelativeChangeOfNodePositions(NodePositions, NewNodePositions)
                } else {
                        difference <- (OldPriGrElEn$ElasticEnergy - PriGrElEn$ElasticEnergy)/PriGrElEn$ElasticEnergy
                }
                
                if(verbose){
                        message(sprintf("Iteration %d diff=%g E=%g MSE=%g EP=%g RP=%g",
                                        i, signif(difference,5), signif(PriGrElEn$ElasticEnergy,5),
                                        signif(PriGrElEn$MSE,5), signif(PriGrElEn$EP,5), signif(PriGrElEn$RP,5)))
                }
                
                if(!is.finite(difference)) difference <- 0
                
                if(difference < eps){
                        Converged <- TRUE
                        #print(sprintf("Converged after %s iterations.", i))
                        break()
                } else {
                        if(i < MaxNumberOfIterations){
                                #PartDataStruct <- PartitionData(X = X, NodePositions = NewNodePositions, SquaredX = SquaredX,
                                #                                TrimmingRadius = TrimmingRadius)
                                NodePositions <- NewNodePositions
                                OldPriGrElEn <- PriGrElEn
                        }
                }
                
        } # end for loop
        #print(sprintf("N. iters: %s", i))
        
        if(DisplayWarnings & !Converged){
                warning(paste0("Maximum number of iterations (", MaxNumberOfIterations, ") has been reached. diff = ", difference))
        }
        
        # final energy computation if needed
        if( (FinalEnergy != "Base") | (!verbose & Mode != 2) ){
                if(FinalEnergy == "Base"){
                        PriGrElEn <-
                                distutils::ElasticEnergy(X = X,
                                                         NodePositions = NewNodePositions,
                                                         ElasticMatrix =  ElasticMatrix,
                                                         Dists = PartDataStruct$Dists,
                                                         partition = assignment,
                                                         NodeAgeTargets = target_positions,
                                                         eta = eta)
                }
                if(FinalEnergy == "Penalized"){
                        PriGrElEn <-
                                distutils::PenalizedElasticEnergy(X = X,
                                                                  NodePositions =  NewNodePositions,
                                                                  ElasticMatrix = ElasticMatrix,
                                                                  Dists = PartDataStruct$Dists,
                                                                  partition = assignment,
                                                                  NodeAgeTargets = target_positions,
                                                                  alpha = alpha,
                                                                  beta = beta,
                                                                  eta = eta)
                }
        }
        
        return(list(EmbeddedNodePositions = NewNodePositions,
                    ElasticEnergy = PriGrElEn$ElasticEnergy,
                    partition = PartDataStruct$Partition,
                    MSE = PriGrElEn$MSE,
                    EP = PriGrElEn$EP,
                    RP = PriGrElEn$RP))
}


# Create input matrices
################################################################################

cellsort <- create_cellsort_obj(seur = seur, cell_ids_in = cell_ids_in, target_cell = target_cell, n_var_features = 5000)
cellsort_toy <- create_cellsort_obj(seur = seur_toy, cell_ids_in = cell_ids_in_toy, target_cell = target_cell_toy, shiftlog = F)

start_time <- Sys.time()
x_pca <- irlba::prcomp_irlba(Matrix::t(cellsort$cell_expr), n = 100)
end_time <- Sys.time()
elapsed <- round(as.numeric(end_time - start_time, units = "mins"), digits = 3)
print(sprintf("PCA done. %s minutes elapsed", elapsed))

X <- x_pca$x

X <- X[, order(abs(apply(x_pca$x,
                         2,
                         function(x) cor(x,
                                         cellsort$cell_metadata$age))), decreasing = T)[1:10]]

X_toy <- as.matrix(Matrix::t(cellsort_toy$cell_expr))
age_vec_toy <- cellsort_toy$cell_metadata$age



# Do the tree
################################################################################

tree_toy <- computeElasticPrincipalTree_edit(X_toy,
                                             NumNodes = 19,
                                             Lambda = 0.03, Mu = 0.01,
                                             age_vec = cellsort_toy$cell_metadata$age,
                                             Do_PCA = F,
                                             eta = 0.5,
                                             FastSolve = T,
                                             MaxNumberOfIterations = 100)

plot(tree_toy[[1]]$age_tau$pseudotime, tree_toy[[1]]$age_tau$average_age)

cor(tree_toy[[1]]$age_tau$pseudotime,
    tree_toy[[1]]$age_tau$average_age)

plot(tree_toy[[1]]$NodePositions[, 1],
     tree_toy[[1]]$NodePositions[, 2])

tree_real <- computeElasticPrincipalTree_edit(X,
                                              NumNodes = 15,
                                              Lambda = 0.03, Mu = 0.01,
                                              age_vec = cellsort$cell_metadata$age,
                                              Do_PCA = F,
                                              eta = 0.5,
                                              FastSolve = T,
                                              MaxNumberOfIterations = 100)

plot(tree_real[[1]]$age_tau$pseudotime, tree_real[[1]]$age_tau$average_age)
cor(tree_real[[1]]$age_tau$pseudotime,
    tree_real[[1]]$age_tau$average_age)
