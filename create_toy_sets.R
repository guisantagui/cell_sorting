library(Seurat)

# Create synthetic seurat objects to test cell sorting algorithm
outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/results/cellsort/"

n_donors <- 109
ncells_per_donor <- 200

n_rep_age <- 6

modulo <- n_donors %% n_rep_age

if (modulo != 0){
        n_donors <- n_donors - modulo
}

toy_donors <- data.frame(donor = paste0("d", 1:n_donors),
                         age = rep(round(seq(20, 100,
                                             length.out = n_donors/n_rep_age)),
                                   each = n_rep_age))

set.seed(123)
toy_donors$age <- toy_donors$age + round(rnorm(nrow(toy_donors), mean = 0, sd = 1))

set.seed(234)
toy_donors$body_weight <- rnorm(n = nrow(toy_donors),
                                mean = 70,
                                sd = 20)

toy_donors$calorie_intake <- 1000 + 35 * toy_donors$body_weight - 0.15 * toy_donors$body_weight^2 + rnorm(nrow(toy_donors), mean = 0, sd = 3)


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

gene_func3 <- function(age, sd = 2, slope = 2, age_start = 50, baseline = 0){
        #age_start_lev <- age_start * slope
        expr <- (-age_start * slope + age * slope)# + baseline 
        expr[age < age_start] <- 0
        expr <- expr + baseline
        expr <- expr + rnorm(length(age), sd = sd)
        return(expr)
}

# No exprtession until an age_start, then exponential.
gene_func4 <- function(age, sd = 2, pow = 2, age_start, baseline = 10){
        #age <- toy_donors$age
        
        expr <- (age/(age_start/(baseline)^(1/pow)))^pow
        
        expr[age < age_start] <- rep(baseline, length(expr[age < age_start]))
        expr <- expr + rnorm(n = length(expr), mean = 0, sd = sd)
        return(expr)
}

get_toy_cells <- function(centroid_mat, sd = "adaptative", ncells = 100){
        #cell_mat <- matrix(nrow = 0, ncol = ncol(centroid_mat),
        #                   dimnames = list(NULL, colnames(centroid_mat)))
        cell_mat <- list()
        for(i in 1:nrow(centroid_mat)){
                cell_names <- paste(rownames(centroid_mat)[i],
                                    paste0("C", 1:ncells),
                                    sep = "_")
                don_mat <- list()
                #don_mat <- matrix(nrow = ncells, ncol = 0,
                #                  dimnames = list(paste(rownames(centroid_mat)[i],
                #                                        paste0("C", 1:ncells),
                #                                        sep = "_"),
                #                                  NULL))
                for(j in 1:ncol(centroid_mat)){
                        if (sd == "adaptative"){
                                sd <- sqrt(var(centroid_mat[, j]))
                        }
                        cell_vec <- rnorm(ncells, mean = centroid_mat[i, j],
                                          sd = sd)
                        to_bind <- matrix(cell_vec, nrow = length(cell_vec),
                                          ncol = 1,
                                          dimnames = list(cell_names,
                                                          colnames(centroid_mat)[j]))
                        don_mat[[j]] <- to_bind
                        #don_mat <- rbind(don_mat, to_bind)
                }
                don_mat <- do.call(cbind, don_mat)
                cell_mat[[i]] <- don_mat
                #cell_mat <- cbind(cell_mat, don_mat)
        }
        cell_mat <- do.call(rbind, cell_mat)
        return(cell_mat)
}

# Genes linearly associated to age
################################################################################

# Let's create a 50 gene dataset
# 10% of the genes (idxs 1:5) advance linearly with age, with different slopes

G1_to_G5 <- data.frame(G1 = gene_func1(toy_donors$age,
                                       sd = 1,
                                       intercept = 500,
                                       slope = 1.5),
                       G2 = gene_func1(toy_donors$age,
                                       sd = 2,
                                       intercept = 1000,
                                       slope = -.5),
                       G3 = gene_func1(toy_donors$age,
                                       sd = .4,
                                       intercept = 520,
                                       slope = .5),
                       G4 = gene_func1(toy_donors$age,
                                       sd = 3,
                                       intercept = 270,
                                       slope = .75),
                       G5 = gene_func1(toy_donors$age,
                                       sd = 2,
                                       intercept = 300,
                                       slope = -.2))


# First branch: genes exponentially activated after a given. Branching point: 35
################################################################################
first_branch_splitage <- 35

G6_baseline <- 300
G6_sd <- 10
G7_baseline <- 100
G7_sd <- 2
G8_baseline <- 200
G8_sd <- 6
G9_baseline <- 500
G9_sd <- 10
G10_baseline <- 800
G10_sd <- 2


# Now 10% of the genes are expressed at a baseline in half of the samples
# throughout their lifetime, but on the other half are activated 
G6_to_G10 <- data.frame(G6 = rnorm(n = nrow(toy_donors),
                                   mean = G6_baseline,
                                   sd = G6_sd),
                        G7 = rnorm(n = nrow(toy_donors),
                                   mean = G7_baseline,
                                   sd = G7_sd),
                        G8 = rnorm(n = nrow(toy_donors),
                                   mean = G8_baseline,
                                   sd = G8_sd),
                        G9 = rnorm(n = nrow(toy_donors),
                                   mean = G9_baseline,
                                   sd = G9_sd),
                        G10 = rnorm(n = nrow(toy_donors),
                                    mean = G10_baseline,
                                    sd = G10_sd))

# Even index samples are going to branch due to exponential activation of genes
# G6-G10 after age = 35
G6_to_G10[1:nrow(G6_to_G10) %% 2 == 0, 1] <- gene_func4(toy_donors$age[1:nrow(G6_to_G10) %% 2 == 0],
                                                        sd = G6_sd,
                                                        pow = 1.5,
                                                        age_start = first_branch_splitage,
                                                        baseline = G6_baseline)
G6_to_G10[1:nrow(G6_to_G10) %% 2 == 0, 2] <- gene_func4(toy_donors$age[1:nrow(G6_to_G10) %% 2 == 0],
                                                        sd = G7_sd,
                                                        pow = 1.3,
                                                        age_start = first_branch_splitage,
                                                        baseline = G7_baseline)
G6_to_G10[1:nrow(G6_to_G10) %% 2 == 0, 3] <- gene_func4(toy_donors$age[1:nrow(G6_to_G10) %% 2 == 0],
                                                        sd = G8_sd,
                                                        pow = 1.2,
                                                        age_start = first_branch_splitage,
                                                        baseline = G8_baseline)
G6_to_G10[1:nrow(G6_to_G10) %% 2 == 0, 4] <- gene_func4(toy_donors$age[1:nrow(G6_to_G10) %% 2 == 0],
                                                        sd = G9_sd,
                                                        pow = 1.5,
                                                        age_start = first_branch_splitage,
                                                        baseline = G9_baseline)
G6_to_G10[1:nrow(G6_to_G10) %% 2 == 0, 5] <- gene_func4(toy_donors$age[1:nrow(G6_to_G10) %% 2 == 0],
                                                        sd = G10_sd,
                                                        pow = 1.4,
                                                        age_start = first_branch_splitage,
                                                        baseline = G10_baseline)

plot(toy_donors$age, G6_to_G10[, 1])

# Branch 2.1. Take samples with odd indices and create a split in age = 55
################################################################################

second_branch_splitage <- 55

G11_baseline <- 200
G11_sd <- 3
G12_baseline <- 300
G12_sd <- 4
G13_baseline <- 500
G13_sd <- 6
G14_baseline <- 400
G14_sd <- 12
G15_baseline <- 250
G15_sd <- 5


# Now in 10% of the genes there's baseline expression. In the samples with
# indexes divisible by 3, the expression of these genes will deviate linearly
# from the baseline from age 60 onwards positively. In the odd ones that are
# not divisible by 3 they will deviate linearly from the baseline in a negative
# manner. In the even ones they will be constant
G11_to_G15 <- data.frame(G11 = rnorm(n = nrow(toy_donors),
                                    mean = G11_baseline,
                                    sd = G11_sd),
                         G12 = rnorm(n = nrow(toy_donors),
                                     mean = G12_baseline,
                                     sd = G12_sd),
                         G13 = rnorm(n = nrow(toy_donors),
                                     mean = G13_baseline,
                                     sd = G13_sd),
                         G14 = rnorm(n = nrow(toy_donors),
                                     mean = G14_baseline,
                                     sd = G14_sd),
                         G15 = rnorm(n = nrow(toy_donors),
                                      mean = G15_baseline,
                                      sd = G15_sd))


G11_to_G15[1:nrow(G11_to_G15) %% 3 == 0 & 1:nrow(G11_to_G15) %% 2 != 0, 1] <- gene_func3(toy_donors$age[1:nrow(G11_to_G15) %% 3 == 0 & 1:nrow(G11_to_G15) %% 2 != 0],
                                                                                         sd = G11_sd,
                                                                                         slope = 2,
                                                                                         age_start = second_branch_splitage,
                                                                                         baseline = G11_baseline)
G11_to_G15[1:nrow(G11_to_G15) %% 3 == 0 & 1:nrow(G11_to_G15) %% 2 != 0, 2] <- gene_func3(toy_donors$age[1:nrow(G11_to_G15) %% 3 == 0 & 1:nrow(G11_to_G15) %% 2 != 0],
                                                                                         sd = G12_sd,
                                                                                         slope = 1.5,
                                                                                         age_start = second_branch_splitage,
                                                                                         baseline = G12_baseline)
G11_to_G15[1:nrow(G11_to_G15) %% 3 == 0 & 1:nrow(G11_to_G15) %% 2 != 0, 3] <- gene_func3(toy_donors$age[1:nrow(G11_to_G15) %% 3 == 0 & 1:nrow(G11_to_G15) %% 2 != 0],
                                                                                         sd = G13_sd,
                                                                                         slope = 1,
                                                                                         age_start = second_branch_splitage,
                                                                                         baseline = G13_baseline)
G11_to_G15[1:nrow(G11_to_G15) %% 3 == 0 & 1:nrow(G11_to_G15) %% 2 != 0, 4] <- gene_func3(toy_donors$age[1:nrow(G11_to_G15) %% 3 == 0 & 1:nrow(G11_to_G15) %% 2 != 0],
                                                                                         sd = G14_sd,
                                                                                         slope = 2.5,
                                                                                         age_start = second_branch_splitage,
                                                                                         baseline = G14_baseline)
G11_to_G15[1:nrow(G11_to_G15) %% 3 == 0 & 1:nrow(G11_to_G15) %% 2 != 0, 5] <- gene_func3(toy_donors$age[1:nrow(G11_to_G15) %% 3 == 0 & 1:nrow(G11_to_G15) %% 2 != 0],
                                                                                         sd = G15_sd,
                                                                                         slope = 3,
                                                                                         age_start = second_branch_splitage,
                                                                                         baseline = G15_baseline)



G11_to_G15[1:nrow(G11_to_G15) %% 3 != 0 & 1:nrow(G11_to_G15) %% 2 != 0, 1] <- gene_func3(toy_donors$age[1:nrow(G11_to_G15) %% 3 != 0 & 1:nrow(G11_to_G15) %% 2 != 0],
                                                                                         sd = G11_sd,
                                                                                         slope = -2,
                                                                                         age_start = second_branch_splitage,
                                                                                         baseline = G11_baseline)
G11_to_G15[1:nrow(G11_to_G15) %% 3 != 0 & 1:nrow(G11_to_G15) %% 2 != 0, 2] <- gene_func3(toy_donors$age[1:nrow(G11_to_G15) %% 3 != 0 & 1:nrow(G11_to_G15) %% 2 != 0],
                                                                                         sd = G12_sd,
                                                                                         slope = -1.5,
                                                                                         age_start = second_branch_splitage,
                                                                                         baseline = G12_baseline)
G11_to_G15[1:nrow(G11_to_G15) %% 3 != 0 & 1:nrow(G11_to_G15) %% 2 != 0, 3] <- gene_func3(toy_donors$age[1:nrow(G11_to_G15) %% 3 != 0 & 1:nrow(G11_to_G15) %% 2 != 0],
                                                                                         sd = G13_sd,
                                                                                         slope = -1,
                                                                                         age_start = second_branch_splitage,
                                                                                         baseline = G13_baseline)
G11_to_G15[1:nrow(G11_to_G15) %% 3 != 0 & 1:nrow(G11_to_G15) %% 2 != 0, 4] <- gene_func3(toy_donors$age[1:nrow(G11_to_G15) %% 3 != 0 & 1:nrow(G11_to_G15) %% 2 != 0],
                                                                                         sd = G14_sd,
                                                                                         slope = -2.5,
                                                                                         age_start = second_branch_splitage,
                                                                                         baseline = G14_baseline)
G11_to_G15[1:nrow(G11_to_G15) %% 3 != 0 & 1:nrow(G11_to_G15) %% 2 != 0, 5] <- gene_func3(toy_donors$age[1:nrow(G11_to_G15) %% 3 != 0 & 1:nrow(G11_to_G15) %% 2 != 0],
                                                                                         sd = G15_sd,
                                                                                         slope = -3,
                                                                                         age_start = second_branch_splitage,
                                                                                         baseline = G15_baseline)
plot(toy_donors$age, G11_to_G15[, 1])


# Branch 2.2. Take samples with even indices and create a split in age = 70
################################################################################

third_branch_splitage <- 70

G16_baseline <- 100
G16_sd <- 4
G17_baseline <- 150
G17_sd <- 10
G18_baseline <- 240
G18_sd <- 3
G19_baseline <- 130
G19_sd <- 4
G20_baseline <- 210
G20_sd <- 6


G16_to_G20 <- data.frame(G16 = rnorm(n = nrow(toy_donors),
                                     mean = G16_baseline,
                                     sd = G16_sd),
                         G17 = rnorm(n = nrow(toy_donors),
                                     mean = G17_baseline,
                                     sd = G17_sd),
                         G18 = rnorm(n = nrow(toy_donors),
                                     mean = G18_baseline,
                                     sd = G18_sd),
                         G19 = rnorm(n = nrow(toy_donors),
                                     mean = G19_baseline,
                                     sd = G19_sd),
                         G20 = rnorm(n = nrow(toy_donors),
                                     mean = G20_baseline,
                                     sd = G20_sd))


G16_to_G20[1:nrow(G6_to_G10) %% 2 == 0 & 1:nrow(G6_to_G10) %% 4 == 0, 1] <- gene_func4(toy_donors$age[1:nrow(G16_to_G20) %% 2 == 0 & 1:nrow(G6_to_G10) %% 4 == 0],
                                                                                       sd = G16_sd,
                                                                                       pow = 1.3,
                                                                                       age_start = third_branch_splitage,
                                                                                       baseline = G16_baseline)
G16_to_G20[1:nrow(G6_to_G10) %% 2 == 0 & 1:nrow(G6_to_G10) %% 4 == 0, 2] <- gene_func4(toy_donors$age[1:nrow(G16_to_G20) %% 2 == 0 & 1:nrow(G6_to_G10) %% 4 == 0],
                                                                                       sd = G17_sd,
                                                                                       pow = 1.5,
                                                                                       age_start = third_branch_splitage,
                                                                                       baseline = G17_baseline)
G16_to_G20[1:nrow(G6_to_G10) %% 2 == 0 & 1:nrow(G6_to_G10) %% 4 == 0, 3] <- gene_func4(toy_donors$age[1:nrow(G16_to_G20) %% 2 == 0 & 1:nrow(G6_to_G10) %% 4 == 0],
                                                                                       sd = G18_sd,
                                                                                       pow = 1.6,
                                                                                       age_start = third_branch_splitage,
                                                                                       baseline = G18_baseline)
G16_to_G20[1:nrow(G6_to_G10) %% 2 == 0 & 1:nrow(G6_to_G10) %% 4 == 0, 4] <- gene_func4(toy_donors$age[1:nrow(G16_to_G20) %% 2 == 0 & 1:nrow(G6_to_G10) %% 4 == 0],
                                                                                       sd = G19_sd,
                                                                                       pow = 1.7,
                                                                                       age_start = third_branch_splitage,
                                                                                       baseline = G19_baseline)
G16_to_G20[1:nrow(G6_to_G10) %% 2 == 0 & 1:nrow(G6_to_G10) %% 4 == 0, 5] <- gene_func4(toy_donors$age[1:nrow(G16_to_G20) %% 2 == 0 & 1:nrow(G6_to_G10) %% 4 == 0],
                                                                                       sd = G20_sd,
                                                                                       pow = 1.4,
                                                                                       age_start = third_branch_splitage,
                                                                                       baseline = G20_baseline)

# Genes that change a lot randomly (noisy). 20%
################################################################################

G21_to_G30 <- data.frame(G21 = rnorm(n = nrow(toy_donors),
                                     mean = 2000,
                                     sd = 600),
                         G22 = rnorm(n = nrow(toy_donors),
                                     mean = 3000,
                                     sd = 400),
                         G23 = rnorm(n = nrow(toy_donors),
                                     mean = 5000,
                                     sd = 700),
                         G24 = rnorm(n = nrow(toy_donors),
                                     mean = 3500,
                                     sd = 500),
                         G25 = rnorm(n = nrow(toy_donors),
                                     mean = 2500,
                                     sd = 750),
                         G26 = rnorm(n = nrow(toy_donors),
                                     mean = 6000,
                                     sd = 650),
                         G27 = rnorm(n = nrow(toy_donors),
                                     mean = 2600,
                                     sd = 300),
                         G28 = rnorm(n = nrow(toy_donors),
                                     mean = 3100,
                                     sd = 600),
                         G29 = rnorm(n = nrow(toy_donors),
                                     mean = 1000,
                                     sd = 300),
                         G30 = rnorm(n = nrow(toy_donors),
                                     mean = 10000,
                                     sd = 1000))

# Genes that are associated to body weight
################################################################################

G31_to_G35 <- data.frame(G31 = gene_func4(toy_donors$body_weight,
                                          sd = 6,
                                          pow = 1.4, min(toy_donors$body_weight),
                                          baseline = 10),
                         G32 = gene_func4(toy_donors$body_weight,
                                          sd = 4,
                                          pow = 1.3, min(toy_donors$body_weight),
                                          baseline = 20),
                         G33 = gene_func4(toy_donors$body_weight,
                                          sd = 10,
                                          pow = 1.5, min(toy_donors$body_weight),
                                          baseline = 30),
                         G34 = gene_func4(toy_donors$body_weight,
                                          sd = 5,
                                          pow = 1.6, min(toy_donors$body_weight),
                                          baseline = 15),
                         G35 = gene_func4(toy_donors$body_weight,
                                          sd = 3,
                                          pow = 1.4, min(toy_donors$body_weight),
                                          baseline = 25))

# Create a split based on body weight (higher than mean trends diverge)
################################################################################

G36_to_G45 <- data.frame(G36 = gene_func3(toy_donors$body_weight,
                                          sd = 3.5,
                                          slope = 1.5,
                                          age_start = mean(toy_donors$body_weight),
                                          baseline = 300),
                         G37 = gene_func3(toy_donors$body_weight,
                                          sd = 3,
                                          slope = 2,
                                          age_start = mean(toy_donors$body_weight),
                                          baseline = 200),
                         G38 = gene_func3(toy_donors$body_weight,
                                          sd = 2.5,
                                          slope = 2.3,
                                          age_start = mean(toy_donors$body_weight),
                                          baseline = 250),
                         G39 = gene_func3(toy_donors$body_weight,
                                          sd = 4,
                                          slope = 0.5,
                                          age_start = mean(toy_donors$body_weight),
                                          baseline = 220),
                         G40 = gene_func3(toy_donors$body_weight,
                                          sd = 2,
                                          slope = 0.75,
                                          age_start = mean(toy_donors$body_weight),
                                          baseline = 1000),
                         G41 = gene_func3(toy_donors$body_weight,
                                          sd = 3.5,
                                          slope = 1.5,
                                          age_start = mean(toy_donors$body_weight),
                                          baseline = 600),
                         G42 = gene_func3(toy_donors$body_weight,
                                          sd = 3,
                                          slope = 2,
                                          age_start = mean(toy_donors$body_weight),
                                          baseline = 500),
                         G43 = gene_func3(toy_donors$body_weight,
                                          sd = 2.5,
                                          slope = 2.3,
                                          age_start = mean(toy_donors$body_weight),
                                          baseline = 320),
                         G44 = gene_func3(toy_donors$body_weight,
                                          sd = 4,
                                          slope = 0.5,
                                          age_start = mean(toy_donors$body_weight),
                                          baseline = 400),
                         G45 = gene_func3(toy_donors$body_weight,
                                          sd = 2,
                                          slope = 1,
                                          age_start = mean(toy_donors$body_weight),
                                          baseline = 300))

G36_to_G45[1:nrow(toy_donors) %% 2 == 0, 1] <- gene_func3(age = toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = 3.5,
                                                          slope = -1.5,
                                                          age_start = mean(toy_donors$body_weight),
                                                          baseline = 300)
G36_to_G45[1:nrow(toy_donors) %% 2 == 0, 2] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = 3,
                                                          slope = -2,
                                                          age_start = mean(toy_donors$body_weight),
                                                          baseline = 200)
G36_to_G45[1:nrow(toy_donors) %% 2 == 0, 3] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = 2.5,
                                                          slope = -2.3,
                                                          age_start = mean(toy_donors$body_weight),
                                                          baseline = 250)
G36_to_G45[1:nrow(toy_donors) %% 2 == 0, 4] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = 4,
                                                          slope = -0.5,
                                                          age_start = mean(toy_donors$body_weight),
                                                          baseline = 220)
G36_to_G45[1:nrow(toy_donors) %% 2 == 0, 5] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = 2,
                                                          slope = -0.75,
                                                          age_start = mean(toy_donors$body_weight),
                                                          baseline = 1000)
G36_to_G45[1:nrow(toy_donors) %% 2 == 0, 6] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = 3.5,
                                                          slope = -1.5,
                                                          age_start = mean(toy_donors$body_weight),
                                                          baseline = 600)
G36_to_G45[1:nrow(toy_donors) %% 2 == 0, 7] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = 3,
                                                          slope = -2,
                                                          age_start = mean(toy_donors$body_weight),
                                                          baseline = 500)
G36_to_G45[1:nrow(toy_donors) %% 2 == 0, 8] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = 2.5,
                                                          slope = -2.3,
                                                          age_start = mean(toy_donors$body_weight),
                                                          baseline = 320)
G36_to_G45[1:nrow(toy_donors) %% 2 == 0, 9] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = 4,
                                                          slope = -0.5,
                                                          age_start = mean(toy_donors$body_weight),
                                                          baseline = 400)
G36_to_G45[1:nrow(toy_donors) %% 2 == 0, 10] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                           sd = 2,
                                                           slope = -1,
                                                           age_start = mean(toy_donors$body_weight),
                                                           baseline = 300)

# Add another split of BW, from first quartile onwards
################################################################################
G46_baseline <- 800
G46_sd <- 5
G47_baseline <- 210
G47_sd <- 6
G48_baseline <- 300
G48_sd <- 1
G49_baseline <- 100
G49_sd <- 3
G50_baseline <- 200
G50_sd <- 2

G46_to_G50 <- data.frame(G46 = rnorm(nrow(toy_donors),
                                     mean = G46_baseline,
                                     sd = G46_sd),
                         G47 = rnorm(nrow(toy_donors),
                                     mean = G47_baseline,
                                     sd = G47_sd),
                         G48 = rnorm(nrow(toy_donors),
                                     mean = G48_baseline,
                                     sd = G48_sd),
                         G49 = rnorm(nrow(toy_donors),
                                     mean = G49_baseline,
                                     sd = G49_sd),
                         G50 = rnorm(nrow(toy_donors),
                                     mean = G50_baseline,
                                     sd = G50_sd))

G46_to_G50[1:nrow(toy_donors) %% 2 == 0, 1] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = G46_sd,
                                                          slope = -1,
                                                          age_start = quantile(toy_donors$body_weight)[2],
                                                          baseline = G46_baseline)

G46_to_G50[1:nrow(toy_donors) %% 2 == 0, 2] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = G47_sd,
                                                          slope = -1,
                                                          age_start = quantile(toy_donors$body_weight)[2],
                                                          baseline = G47_baseline)

G46_to_G50[1:nrow(toy_donors) %% 2 == 0, 3] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = G48_sd,
                                                          slope = -1,
                                                          age_start = quantile(toy_donors$body_weight)[2],
                                                          baseline = G48_baseline)

G46_to_G50[1:nrow(toy_donors) %% 2 == 0, 4] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = G49_sd,
                                                          slope = -1,
                                                          age_start = quantile(toy_donors$body_weight)[2],
                                                          baseline = G49_baseline)

G46_to_G50[1:nrow(toy_donors) %% 2 == 0, 5] <- gene_func3(toy_donors$body_weight[1:nrow(toy_donors) %% 2 == 0],
                                                          sd = G50_sd,
                                                          slope = -1,
                                                          age_start = quantile(toy_donors$body_weight)[2],
                                                          baseline = G50_baseline)

dataset <- cbind(G1_to_G5,
                 G6_to_G10,
                 G11_to_G15,
                 G16_to_G20,
                 G21_to_G30,
                 G31_to_G35,
                 G36_to_G45,
                 G46_to_G50)

rownames(dataset) <- toy_donors$donor



cell_dataset <- get_toy_cells(dataset, sd = 15, ncells = ncells_per_donor)

metadata <- data.frame(cell = rownames(cell_dataset),
                       donor = gsub("\\_.*", "", rownames(cell_dataset)),
                       age = rep(toy_donors$age, each = ncells_per_donor),
                       body_weight = rep(toy_donors$body_weight,
                                         each = ncells_per_donor))

cell_dataset_stand <- apply(cell_dataset, 2, function(x) (x - mean(x))/sd(x))

plot(dataset[, 1], dataset[, 6])

plot(cell_dataset[, 1], cell_dataset[, 6], xlim = c(500, 650), ylim = c(400, 650))

plot(cell_dataset[, 1], cell_dataset[, 10])

plot(cell_dataset[, 1], cell_dataset[, 11])

plot(cell_dataset_stand[, 1], cell_dataset_stand[, 10])

cell_pca <- prcomp(cell_dataset_stand)

metadata_4_pca <- metadata

metadata_4_pca$sample <- metadata$cell

plotUtils::plotPCA(cell_pca, samp_info = metadata_4_pca, col = "age")

plotUtils::plotPCA(cell_pca, samp_info = metadata_4_pca, col = "body_weight")


cell_dataset <- round(cell_dataset)

cell_dataset[cell_dataset < 0] <- 0



# Create the seurat object

cells_sparse <- as(t(cell_dataset), "dgCMatrix")

seur <- CreateSeuratObject(counts = cells_sparse,
                           project = "MyProject",
                           min.cells = 3,
                           min.features = 1)
seur@meta.data$donor <- metadata$donor[match(rownames(seur@meta.data),
                                             metadata$cell)]
seur@meta.data$age <- metadata$age[match(rownames(seur@meta.data),
                                         metadata$cell)]
seur@meta.data$body_weight <- metadata$body_weight[match(rownames(seur@meta.data),
                                                         metadata$cell)]

saveRDS(seur, file = sprintf("%sseur_synth_complex.rds", outDir))


# Create a simple synthetic dataset, with one bifurcation at 50 years old,
# caused by a linearly divergent gene. There will be another confounder,
# associated to body weight
################################################################################

# Create the age associated bifurcation
dataset_one_bifurcation <- data.frame(G1 = gene_func1(toy_donors$age,
                                                      sd = 2,
                                                      intercept = 100,
                                                      slope = 1),
                                      G2 = rnorm(nrow(toy_donors),
                                                 mean = 200,
                                                 sd = 3))

dataset_one_bifurcation$G2[1:nrow(dataset_one_bifurcation) %% 2 == 0] <- gene_func3(toy_donors$age[1:nrow(dataset_one_bifurcation) %% 2 == 0],
                                                                                    age_start = 50,
                                                                                    baseline = 200,
                                                                                    slope = 2)


# Create the body weight associated bifurcation
dataset_one_bifurcation$G3 <- gene_func1(toy_donors$body_weight,
                                         sd = 2,
                                         intercept = 350,
                                         slope = 1.3)

dataset_one_bifurcation$G4 <- rnorm(nrow(toy_donors),
                                    mean = 300,
                                    sd = 3)

dataset_one_bifurcation$G4[1:nrow(dataset_one_bifurcation) %% 2 != 0] <- gene_func3(toy_donors$body_weight[1:nrow(dataset_one_bifurcation) %% 2 != 0],
                                                                                    age_start = mean(toy_donors$body_weight),
                                                                                    baseline = 300,
                                                                                    slope = 2,
                                                                                    sd = 3)

# Add two noise variables, static
dataset_one_bifurcation$G5 <- rnorm(nrow(toy_donors),
                                    mean = 600,
                                    sd = 10)

dataset_one_bifurcation$G6 <- rnorm(nrow(toy_donors),
                                    mean = 500,
                                    sd = 7)

# Add two genes related to calorie intake
dataset_one_bifurcation$G7 <- gene_func1(toy_donors$calorie_intake,
                                         sd = 4,
                                         intercept = 5000,
                                         slope = -0.1)

dataset_one_bifurcation$G8 <- gene_func1(toy_donors$calorie_intake,
                                         sd = 3,
                                         intercept = 7000,
                                         slope = -0.12)

# Add a bifurcation related to calorie intake
dataset_one_bifurcation$G9 <- rnorm(nrow(toy_donors),
                                    mean = 500,
                                    sd = 2)

dataset_one_bifurcation$G9[1:nrow(dataset_one_bifurcation) %% 3 == 0] <- gene_func3(toy_donors$calorie_intake[1:nrow(dataset_one_bifurcation) %% 3 == 0],
                                                                                    age_start = mean(toy_donors$body_weight),
                                                                                    baseline = 500,
                                                                                    slope = -0.022,
                                                                                    sd = 2)

# Add another gene linearly related to age
dataset_one_bifurcation$G10 <- gene_func1(toy_donors$age,
                                          sd = 2.5,
                                          intercept = 200,
                                          slope = 1.2)

# Add another bifurcation related to body weight, but this time taking random
# samples

dataset_one_bifurcation$G11 <- rnorm(nrow(toy_donors),
                                     mean = 500,
                                     sd = 4)

set.seed(321)
samp_idxs <- sample(1:nrow(dataset_one_bifurcation), nrow(dataset_one_bifurcation), replace = F)

dataset_one_bifurcation$G11[samp_idxs %% 2 != 0] <- gene_func3(toy_donors$calorie_intake[samp_idxs %% 2 != 0],
                                                               age_start = mean(toy_donors$calorie_intake),
                                                               baseline = 500,
                                                               slope = -.5,
                                                               sd = 4)

dataset_one_bifurcation$G12 <-  gene_func3(toy_donors$calorie_intake,
                                           age_start = mean(toy_donors$calorie_intake),
                                           baseline = 300,
                                           slope = -.25,
                                           sd = 3)

dataset_one_bifurcation$G12[samp_idxs %% 2 != 0] <-  gene_func3(toy_donors$calorie_intake[samp_idxs %% 2 != 0],
                                                                age_start = mean(toy_donors$calorie_intake),
                                                                baseline = 300,
                                                                slope = .25,
                                                                sd = 3)


plot(toy_donors$calorie_intake, dataset_one_bifurcation$G11)
plot(toy_donors$calorie_intake, dataset_one_bifurcation$G12)

apply(dataset_one_bifurcation, 2, var)

rownames(dataset_one_bifurcation) <- toy_donors$donor

cell_dataset_one_bifurcation <- get_toy_cells(dataset_one_bifurcation,
                                              sd = 5,
                                              ncells = ncells_per_donor)

metadata_one_bifurcation <- data.frame(cell = rownames(cell_dataset_one_bifurcation),
                                       donor = gsub("\\_.*", "", rownames(cell_dataset_one_bifurcation)),
                                       age = rep(toy_donors$age, each = ncells_per_donor),
                                       body_weight = rep(toy_donors$body_weight,
                                                         each = ncells_per_donor),
                                       calorie_intake = rep(toy_donors$calorie_intake,
                                                            each = ncells_per_donor))
metadata_one_bifurcation_4pca_plot <- metadata_one_bifurcation
metadata_one_bifurcation_4pca_plot$sample <- metadata_one_bifurcation_4pca_plot$cell

# Round and remove potential negative values
cell_dataset_one_bifurcation <- round(cell_dataset_one_bifurcation)

cell_dataset_one_bifurcation[cell_dataset_one_bifurcation < 0] <- 0

pca_one_bifurcation <- prcomp(cell_dataset_one_bifurcation, scale. = T, center = T)

plotUtils::plotPCA(pca_one_bifurcation, samp_info = metadata_one_bifurcation_4pca_plot,
                   col = "age")

plotUtils::doPCAMultiPlot(pca_one_bifurcation, samp_info = metadata_one_bifurcation_4pca_plot,
                          col = "age",
                          nComps = 5)

plot(metadata_one_bifurcation$age, cell_dataset_one_bifurcation[, 2])

one_bifurcation_umap <- uwot::umap(cell_dataset_one_bifurcation, ret_model = T,
                                   min_dist = .3)

df <- as.data.frame(one_bifurcation_umap$embedding)

colnames(df) <- paste0("UMAP", 1:ncol(df))
df$age <- metadata_one_bifurcation$age
df$body_weight <- metadata_one_bifurcation$body_weight
df$calorie_intake <- metadata_one_bifurcation$calorie_intake
df$is_even <- as.numeric(gsub("d", "", metadata_one_bifurcation$donor)) %% 2 == 0


ggplot(data = df, mapping = aes(x = UMAP1, y = UMAP2, col = age)) +
        geom_point()

ggplot(data = df, mapping = aes(x = UMAP1, y = UMAP2, col = body_weight)) +
        geom_point()

ggplot(data = df, mapping = aes(x = UMAP1, y = UMAP2, col = calorie_intake)) +
        geom_point()

ggplot(data = df, mapping = aes(x = UMAP1, y = UMAP2, col = is_even)) +
        geom_point()

# Create seurat object
cells_sparse_one_bifurcation <- as(t(cell_dataset_one_bifurcation), "dgCMatrix")

seur_one_bifurcation <- CreateSeuratObject(counts = cells_sparse_one_bifurcation,
                                           project = "MyProject",
                                           min.cells = 3,
                                           min.features = 1)
seur_one_bifurcation@meta.data$donor <- metadata_one_bifurcation$donor[match(rownames(seur_one_bifurcation@meta.data),
                                                                             metadata_one_bifurcation$cell)]
seur_one_bifurcation@meta.data$age <- metadata_one_bifurcation$age[match(rownames(seur_one_bifurcation@meta.data),
                                                                         metadata_one_bifurcation$cell)]
seur_one_bifurcation@meta.data$body_weight <- metadata_one_bifurcation$body_weight[match(rownames(seur_one_bifurcation@meta.data),
                                                                                         metadata_one_bifurcation$cell)]

seur_one_bifurcation@meta.data$calorie_intake <- metadata_one_bifurcation$calorie_intake[match(rownames(seur_one_bifurcation@meta.data),
                                                                                               metadata_one_bifurcation$cell)]
saveRDS(seur_one_bifurcation, file = sprintf("%sseur_synth_one_bifurcation.rds", outDir))
