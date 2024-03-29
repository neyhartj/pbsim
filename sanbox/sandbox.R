## Testing sim.cross capability

# Let's try to speed this up using rQTL and simcross
# library(qtl)
library(simcross)
library(pbsim)
library(dplyr)
library(stringr)



n.mar = c(500, 500, 500)
len = c(120, 130, 140)
map <- NULL





# Population size per cross
n.ind = 30
# Number of crosses
n.crosses = 50
# Number of SNPs
n.snp <- 3072
# Selfing generations
n.self.gen = 2

# Chromosome lengths
chr.len <- seq(140, 170, by = 5)
# SNPs per chromosome
set.seed(313)
snp.per.chr <- sample(x = 1:7, size = n.snp, replace = T, prob = (chr.len / sum(chr.len))) %>%
  sort() %>%
  table()


# Genetic map
map <- mapply(chr.len, snp.per.chr, FUN = function(cM, snps) 
  runif(n = snps, min = 0, max = cM) %>% sort() )


# Create parents
founders <- replicate(n = n.crosses * 2, expr = {
  sample(x = c(1,3), size = n.snp, replace = T)
})

# Create crossing pairs
cross.pairs <- split(x = seq_len(n.crosses * 2), rep(seq_len(n.crosses), each = 2))

# Create a population
ped <- make.pedigree(n.self.gen = n.self.gen, n.ind = n.ind)
# ped <- data.frame(
#   id = 1:5,
#   mom = c(0,0,1,3,4),
#   dad = c(0,0,2,3,4),
#   sex = c(0,1,0,0,0),
#   gen = c(0,0,1,2,3)
# )


# IDs to extract
prog.ids <- which(ped$gen == (n.self.gen + 1))


genos <- lapply(X = cross.pairs, FUN = function(cross) {
    
    # Extract the parent genos
    parents <- founders[,cross]
    
    # Simulate cross-over data
    xo.data <- sim_from_pedigree_allchr(pedigree = ped, map = map, m = 0)
    
    # Generate genotypes
    convert2geno_allchr(xodat = xo.data, map = map,
                        founder_geno = parents, return.matrix = T,
                        id = prog.ids)
  })




# Attempt to force a cross object
library(qtl)
library(pbsim)


# Simulate a genome and genetic model
n.mar  <- c(505, 505, 505)
len <- c(120, 130, 140)

genome <- sim_genome(len, n.mar)

na_mat <- matrix(nrow = 15, ncol = 4)

qtl.model <- list(na_mat, na_mat)
add.dist <- "geometric"
prob.corr <- c(0.5, 0.1, 0.1, 0.1, 0.2)

genome <- sim_gen_model(genome, qtl.model, add.dist = add.dist, prob.corr = prob.corr)


cross <- sim.cross(map = genome$map, model = genome$gen_model[[1]], n.ind = 100,
                   type = "bcsft", cross.scheme = c(0, 4))







# Reorder based on chromosome, then position
qtl_specs1 <- lapply(qtl_specs1, function(mat) {
  mat1 <- mat[order(mat[,1], mat[,2]),]
  # Remove row names
  row.names(mat1) <- NULL
  # Return the matrix
  return(mat1) })

# Determine if any QTL line up exactly on marker positions
qtl_specs2 <- lapply(X = qtl_specs1, FUN = function(qtlmod) {
  perf.snp <- apply(X = qtlmod, MARGIN = 1, FUN = function(qtl) qtl[2] %in% genome$map[qtl[1]])
  cbind(qtlmod, perf.snp) })

## Check two or more qtl models for common QTL
plei <- vector("list", length(qtl_specs2))

# Iterate over all qtl models
for (i in seq_along(qtl_specs2)) {
  
  x <- qtl_specs2[[i]]
  
  # Iterate over remaining qtl models
  plei[[i]] <- sapply(X = qtl_specs2[-i], FUN = function(y)
    apply(X = x, MARGIN = 1, FUN = function(qtl1) 
      any(apply(X = y, MARGIN = 1, FUN = function(qtl2) 
        qtl1[1] == qtl2[1] & qtl1[2] == qtl2[2] ))))
  
}

# Add peiotropy info to the genetic models
qtl_specs3 <- mapply(qtl_specs2, plei, FUN = function(mod, pl) {
  plei.qtl <- rowSums(pl)
  cbind(mod, plei.qtl) }, SIMPLIFY = FALSE)



## Testing genotype filling
data("s2_genos")
data("s2_snp_info")

library(GSSimTPUpdate)

data("CAP.haploids")
data("CAP.markers")

line_names <- row.names(CAP.haploids) %>%
  str_replace(".[1-2]$", "") %>%
  unique()

split_vec <- split(x = seq(nrow(CAP.haploids)), f = cut(seq(nrow(CAP.haploids)), breaks = nrow(CAP.haploids) / 2))

CAP_geno <- lapply(split_vec, function(ind) CAP.haploids[ind,]) %>%
  lapply(colSums) %>%
  do.call("rbind", .)

row.names(CAP_geno) <- line_names

s6_cap_genos <- CAP_geno
s6_snp_info <- CAP.markers %>% 
  mutate(chrom = str_c(chrom, "H"), cM_pos = pos * 100, rs = as.character(rs)) %>% 
  select(-pos) %>%
  tbl_df

# Grab map positions for markers in s6 from the s2 data.frame
s6_snp_info <- s6_snp_info %>%
  mutate(cM_pos = ifelse(rs %in% s2_snp_info$rs, s2_snp_info$cM_pos[match(rs, s2_snp_info$rs)], cM_pos))

map <- s6_snp_info %>% 
  split(.$chrom) %>% 
  lapply(FUN = function(chr)
    structure(chr$cM_pos, names = as.matrix(chr$rs), class = "A") )


s6_snp_info$cM_pos <- do.call("c", new_map)

    
s2_cap_genos <- s2_genos



## Common markers





# geno <- s2_genos + 1
# 
# # Create a genome with genetic architecture
# len <- tapply(s2_snp_info$cM_pos, s2_snp_info$chrom, max)
# n_mar <- tapply(s2_snp_info$cM_pos, s2_snp_info$chrom, length)
# map <- lapply(split(s2_snp_info, s2_snp_info$chrom), function(chr) structure(chr$cM_pos, names = chr$rs) )


library(pbsim)

geno <- s6_cap_genos

len <- tapply(s6_snp_info$cM_pos, s6_snp_info$chrom, max)
n_mar <- tapply(s6_snp_info$cM_pos, s6_snp_info$chrom, length)
map <- lapply(split(s6_snp_info, s6_snp_info$chrom), function(chr) structure(chr$cM_pos, names = chr$rs) )

genome <- sim_genome(len = len, n.mar = n_mar, map = map)

na_mat <- matrix(nrow = 5, ncol = 4)

qtl.model <- list(na_mat, na_mat)
add.dist <- "geometric"
prob.corr <- cbind(c(1), c(1))

## Test various genetic architectures for genetic correlation between traits
test_cor <- replicate(n = 100, expr = {
  
  genome <- sim_gen_model(genome, qtl.model, add.dist = add.dist, max.qtl = 5,
                          prob.corr = prob.corr)
  
  # new_geno <- fill_qtl_geno(genome = genome, geno = geno)

  # edit the genome
  # genome <- adj_gen_model(genome = genome, geno = geno, pos.cor = TRUE)
  
  # Create pop
  pop <- create_pop(genome = genome, geno = geno)
  
  # # Try creating a sim.cross population from founders
  # founders <- sim_founders(object = genome)
  # 
  # # RIL pedigree
  # ril_ped <- sim_pedigree(n.ind = 100)
  # 
  # # Create a family
  # pop <- sim_family(genome = genome, pedigree = ril_ped, founder_geno = founders)
  
  c(par = cor(pop$geno_val[,"trait1"], pop$geno_val[,"trait2"]))
  
})



# Try creating a sim.cross population from founders
founders <- sim_founders(object = genome)

# RIL pedigree
ril_ped <- sim_pedigree(n.ind = 100)

# Create a family
family <- sim_family(genome = genome, pedigree = ril_ped, founder_geno = founders)
  
  
  # ## Create a population
  # cross <- qtl::sim.cross(map = genome$map, n.ind = 1000, type = "riself")
  # 
  # prog_geno <- lapply(cross$geno, FUN = "[[", "data") %>% do.call("cbind", .)
  # 
  # par_geno <- pop$geno %>% 
  #   do.call("cbind", .) %>%
  #   as.data.frame() %>% 
  #   sample_n(size = 2) %>% 
  #   as.matrix()
  # 
  # prog_geno_recode <- t(apply(X = prog_geno, MARGIN = 1, FUN = function(gen)
  #   par_geno[cbind(gen, seq(ncol(prog_geno)))] )) %>%
  #   structure(., dimnames = list(NULL, colnames(par_geno)))
  # 
  # prog_qtl <- pull_genotype(genome = genome, geno = prog_geno_recode, loci = qtl$qtl_name)
  # 
  # prog_val <- calc_genoval(genome = genome, geno = prog_geno_recode)
  # 
  # c(par = cor(pop$geno_val[,"trait1"], pop$geno_val[,"trait2"]),
  #   prog = cor(prog_val[,"trait1"], prog_val[,"trait2"]) )
  
  c(par = cor(pop$geno_val[,"trait1"], pop$geno_val[,"trait2"]))
  
  }, simplify = FALSE)

# For each QTL, find the correlation of that QTL's genotype with the closest QTL
# of the other trait
qtl_cor <- subset(genome$gen_model[[2]], select = c(qtl_name, qtl1_pair)) %>% 
  apply(MARGIN = 1, FUN = list) %>% 
  unlist(recursive = F) %>% 
  lapply(pull_genotype, genome = genome, geno = pop$geno) %>%
  lapply(cor)

qtl_temp <- pull_genotype(genome = genome, geno = pop$geno, 
                          loci = unlist(genome$gen_model[[2]][1,5:6]))






# What is the correlation between marker genotypes as a function of genetic distance?
## Pairwise correlation between markers on the same chromosome
marker_cor <- s6_cap_genos %>%
  split_geno(genome = genome) %>%
  lapply(cor) %>%
  lapply(FUN = function(corr) corr[lower.tri(corr)])

# Pairwise distance between markers on the same chromosome
marker_dist <- map %>% 
  lapply(dist) %>%
  lapply(as.numeric)

# Bind
cor_by_dist <- data.frame(marker_cor = unlist(marker_cor), 
                          marker_dist = unlist(marker_dist))


    

temp <- marker_cor$`1H` %>% .[lower.tri(.)]  
temp1 <- marker_dist$`1H` %>% as.numeric()

      
## Playing with Meiosis package
library(Meiosis)

n_marker <- pbsim::nmar(genome, by.chr = T)

xoparam <- create_xoparam(L = genome$len)
positions <- lapply(seq_len(nchr(genome)), function(i) sort(runif(n_marker[i], min = 0, max = genome$len[i])))


ind <- replicate(2L, lapply(nmar(genome), function(n) sample(c(0L, 1L), n, replace = TRUE)), 
                 simplify = FALSE)  ## simulate some genotypic data

p_geno <- Meiosis::cross_geno(father = ind, mother = ind, positions = positions, 
                              xoparam = xoparam)






library(pbsim)

map <- lapply(split(s6_snp_info, s6_snp_info$chrom), function(chr) 
  structure(chr$cM_pos, names = chr$rs) )

# CReate a hypred genome using the S6 data
genome <- sim_genome(map = map, type = "hypred")

# Simulate genetic architecture with two traits
na_mat <- matrix(nrow = 15, ncol = 4)

qtl.model <- list(na_mat, na_mat)
add.dist <- "geometric"
prob.corr <- cbind(c(0, 2, 15, 30, 50), c(0.5, 0.1, 0.1, 0.1, 0.2))

genome <- sim_gen_model(genome, qtl.model, de.novo = FALSE, add.dist = add.dist, 
                        prob.corr = prob.corr, max.qtl = 20)











## Sample two vectors given variance-covariance matrix

# Multivariate normal
library(mvtnorm)
library(Matrix)

n = c(100, 100)
corr <- 0.75
sigma <- rbind(c(1, corr),
               c(corr, 1))

# Cholesky composition of sigma
sigma_decomp <- chol(sigma, pivot = TRUE)
sigma_decomp <- sigma_decomp[, order(attr(sigma_decomp, "pivot"))]

# Geometric series
a <- (1 - n) / (1 + n)
a1 <- mapply(a, lapply(n, seq), FUN = `^`, SIMPLIFY = FALSE)

a_mat <- sapply(a1, `length<-`, max(n))
a_mat_sample <- apply(X = a_mat, MARGIN = 2, FUN = sample)
a_mat_sample1 <- ifelse(is.na(a_mat_sample), 0, a_mat_sample) %*% sigma_decomp
a_mat_sample1[is.na(a_mat_sample)] <- 0

cor(a_mat_sample1)


# Sample into matrix
a_sample <- replicate(n = ncol(sigma), sample(a))
# Multiple by decomposition
a_new <- a_sample %*% sigma_decomp



## Create a population with 3 sub-populations
library(pbsim)
library(rrBLUP)
library(tidyverse)

genome <- sim_genome(len = rep(150, 5), n.mar = rep(520, 5))
qtl.model <- matrix(NA, ncol = 4, nrow = 100)
genome1 <- sim_gen_model(genome = genome, qtl.model = qtl.model, add.dist = "geometric")

founders <- sim_founders(genome = genome1)

## Create three populations from the cross between founders
# 1. Backcross to the first founder and self
# 2. Self
# 3. Backcross to the second founder and self

pop1_ped <- sim_pedigree(n.ind = 1000, bcpar = 1, n.bcgen = 2, n.selfgen = 8)
pop2_ped <- sim_pedigree(n.ind = 1000, n.selfgen = 8)
pop3_ped <- sim_pedigree(n.ind = 1000, bcpar = 2, n.bcgen = 2, n.selfgen = 8)

# Combine the pedigrees

# Create the populations
pop1 <- sim_family(genome = genome1, pedigree = pop1_ped, founder.pop = founders, family.num = 1)
pop2 <- sim_family(genome = genome1, pedigree = pop2_ped, founder.pop = founders, family.num = 2)
pop3 <- sim_family(genome = genome1, pedigree = pop3_ped, founder.pop = founders, family.num = 3)

# Combine
pop <- combine_pop(list(pop1, pop2, pop3))

# Phenotype
pop_select <- sim_phenoval(pop = pop, h2 = 0.5, n.env = 1, n.rep = 1) %>%
  select_pop(intensity = 1000, index = 1)

# Genotype
marker_geno <- genotype(genome = genome1, pop = pop_select)

K <- A.mat(X = marker_geno, min.MAF = 0, max.missing = 1)
K_svd <- prcomp(x = K)

svd_df <- K_svd$rotation %>% 
  as.data.frame() %>% 
  rownames_to_column("line_name") %>% 
  mutate(family = str_extract(line_name, "[0-9]{4,}"))

qplot(x = PC2, y = PC3, color = family, data = svd_df)




## Test the expectation of genetic correlation in bi-parental populations
library(tidyverse)
library(pbsim)

example("calc_exp_genvar")

ped <- sim_pedigree(n.ind = 500, n.selfgen = Inf)

qtl.model <- replicate(2, matrix(NA, 50, 4), simplify = FALSE)
genome <- sim_multi_gen_model(genome = genome, qtl.model = qtl.model, corr = 0.5, 
                              prob.corr = cbind(0, 1), add.dist = "normal")

# Simulate the genotypes for 8 founders
founder.pop <- sim_pop(genome = genome, n.ind = 25)

cb <- sim_crossing_block(parents = indnames(founder.pop), n.crosses = 50)

expected <- calc_exp_genvar(genome = genome, pedigree = ped, founder.pop = founder.pop, 
                            crossing.block = cb) %>%
  distinct(parent1, parent2, exp_corG)

## Create bi-parental RIL families using this crossing block
families <- sim_family_cb(genome = genome, pedigree = ped, founder.pop = founder.pop, crossing.block = cb)

## Calculate the genetic correlation per family
realized <- families$geno_val %>% 
  mutate(family = str_extract(ind, "1[0-9]{3}")) %>% 
  group_by(family) %>% 
  summarize(gencor = cor(trait1, trait2))

## Plot
plot(expected$exp_corG, realized$gencor)
cor(expected$exp_corG, realized$gencor)





### Try to speed up the prediction of genetic correlation
library(pbsim)
library(tidyverse)

L <- 100
# rho <- 0.75
rho <- -0.75

# Simulate a genome
genome <- sim_genome(len = rep(150, 3), n.mar = rep(102, 3))
# Simulate a genetic model
qtl.model <- replicate(2, matrix(data = NA, ncol = 4, nrow = L), simplify = FALSE)


genome1 <- sim_multi_gen_model(genome = genome, qtl.model = qtl.model, corr = rho, prob.corr = cbind(30, 1),
                               add.dist = "geometric")



# Simulate a family from founders
founders <- sim_founders(genome = genome1, n.str = 2, c(0, 1))
# pop <- sim_family(genome = genome1, pedigree = sim_pedigree(n.ind = 200), founder.pop = founders)
pop <- sim_pop(genome = genome1, n.ind = 200)

# Measure the genetic correlation
cor(pop$geno_val[,-1])[1,2]

## What is the correlation between QTL genotype states
qtl_geno <- pull_genotype(genome = genome1, geno = pop$geno, loci = qtlnames(genome1)) - 1

D <- cor(qtl_geno)

## Generate new QTL effects
qtl_eff <- replicate(2, sample(((L - 1) / (L + 1)) ^ seq(L)), simplify = FALSE) %>%
  map2(.x = ., .y = lapply(X = genome1$gen_model, "[[", "qtl_name"), ~set_names(.x, .y))



# ###### Just look at 1 "pair" of QTL
# qtl_geno_test <- qtl_geno[,c(1,3)]
# qtl_eff_test <- sapply(X = qtl_eff, FUN = "[[", 1)
# # Muliply by the effects and correlate
# cor(qtl_geno_test * qtl_eff_test)
# 
# ## So if the effects are the same, the correlation among the genotype states is the same 
# ## as the correlation among breeding values




## Generate new QTL effects
qtl_eff <- replicate(2, sample(((L - 1) / (L + 1)) ^ seq(L)), simplify = FALSE) %>%
  map2(.x = ., .y = lapply(X = genome1$gen_model, "[[", "qtl_name"), ~set_names(.x, .y))

# Measure the correlation
bvs <- map(qtl_eff, ~qtl_geno[,names(.x)] %*% .x)
(cor1 <- cor(reduce(bvs, cbind))[1,2])

## Adjust the trait 2 effects by the desired correlation
# qtl_eff_adj <- reduce(map(qtl_eff, as.matrix), rbind)
qtl_eff_adj <- reduce(map(qtl_eff, as.matrix), cbind) %*% chol(rbind(c(1, rho), c(rho, 1)))


# Combine
qtl_eff1 <- rbind(qtl_eff_adj[,1,drop = F], qtl_eff_adj[,2,drop = F])

# Subset D for the two traits
D1 <- D[names(qtl_eff[[1]]), names(qtl_eff[[2]])]

# Revise the effects for trait2
qtl_eff2 <- qtl_eff
qtl_eff2[[2]] <- set_names(as.numeric(qtl_eff_adj[,2] %*% D1), names(qtl_eff[[2]])) 

# Measure the correlation
bvs2 <- map(qtl_eff2, ~qtl_geno[,names(.x)] %*% .x)
(cor2 <- cor(reduce(bvs2, cbind))[1,2])




## BGLR marker effects
library(pbsim)
library(tidyverse)

L <- 30

# Simulate a genome
genome <- sim_genome(len = rep(150, 3), n.mar = rep(100 + (L / 3), 3))
# Simulate a genetic model
qtl.model <- matrix(data = NA, ncol = 4, nrow = L)
genome1 <- sim_gen_model(genome = genome, qtl.model = qtl.model, add.dist = "geometric")

# Create a population
pop <- sim_pop(genome = genome1, n.ind = 600) %>%
  sim_phenoval(pop = ., h2 = 0.5, n.env = 3, n.rep = 1)

# Genotype
geno <- genotype(genome = genome1, pop = pop)


## Use BGLR to predict marker effects
library(BGLR)

y <- pop$geno_val$trait1
Z <- geno
X <- model.matrix(~ 1, pop$geno_val)

ETA <- list(list(X = Z, model = "BayesC"))
fit <- BGLR(y = y, ETA = ETA, verbose = FALSE, nIter = 6000, burnIn = 1000)

# Extract the marker effects
effs <- fit$ETA[[1]]$b
# Extract the pi hyperparameter
pi <- fit$ETA[[1]]$probIn

# Repeat 



# Simulate a genome
genome <- sim_genome(len = rep(150, 3), n.mar = rep(100 + L, 3))
# Simulate a genetic model
qtl.model <- replicate(2, matrix(data = NA, ncol = 4, nrow = L), simplify = FALSE)

genome1 <- sim_multi_gen_model(genome = genome, qtl.model = qtl.model, corr = rho, prob.corr = cbind(30, 1),
                               add.dist = "geometric")




### Ian's problems

library(pbsim)

# Simulate a genome
L <- 30

# Simulate a genome
chr_map <- qtl::sim.map(len = rep(150, 3), n.mar = 100, include.x = FALSE)
genome <- sim_genome(map = chr_map)



# Simulate a genetic model
qtl.model <- replicate(2, cbind(
  chr = c(1, 3),
  pos = sapply(chr_map[c(1,3)], sample, 1),
  add_eff = c(1, 1),
  dom_eff = 0
), simplify = FALSE)

qtl.model <- matrix(NA, nrow = L, ncol = 4)
  
  
genome1 <- sim_multi_gen_model(genome = genome, qtl.model = qtl.model)







library(pbsim)

# Simulate a genome
n.mar  <- c(505, 505, 505)
len <- c(120, 130, 140)

genome <- sim_genome(len, n.mar)


# Create a pedigree with 100 individuals selfed to the F_3 generation
ped <- sim_pedigree(n.ind = 100, n.selfgen = 2)

## If two traits are present, the genetic correlation is calculated
# Simulate two quantitative traits influenced by 50 pleiotropic QTL
qtl.model <- replicate(2, matrix(NA, 50, 4), simplify = FALSE)
genome <- sim_multi_gen_model(genome = genome, qtl.model = qtl.model, corr = 0.5, 
                              prob.corr = cbind(0, 1), add.dist = "normal")

# Simulate the genotypes for 8 founders
founder.pop <- sim_founders(genome = genome, n.str = 8)
training.pop <- sim_phenoval(founder.pop, h2 = 0.8)

# Generate a crossing block with 5 crosses
cb <- sim_crossing_block(parents = indnames(founder.pop), n.crosses = 5)

cb1 <- do.call("rbind", replicate(100, cb, simplify = FALSE))

pred_out <- pred_genvar(genome = genome, pedigree = ped, training.pop = training.pop, founder.pop = founder.pop, crossing.block = cb1)



genome = genome
pedigree = ped
training.pop = training.pop
founder.pop = founder.pop
crossing.block = cb
method = c("RRBLUP")
  

{

n_traits <- length(genome$gen_model)

# Predict marker effects - only if the TP does not have them
if (is.null(training.pop$mar_eff)) {
  
  marker_eff <- pred_mar_eff(genome = genome, training.pop = training.pop, method = method, n.iter = n.iter,
                             burn.in = burn.in, thin = thin, save.at = save.at)
  
} else {
  marker_eff <- training.pop
  
}

## Predict genotypic values in the founder population - this will use the marker effects in the tp
founder_pop1 <- pred_geno_val(genome = genome, training.pop = marker_eff, candidate.pop = founder.pop)
# Predict marker effects
marker_eff <- marker_eff$mar_eff

## Find the positions of these markers
marker_pos <- find_markerpos(genome = genome, marker = marker_eff$marker)
marker_pos$add_eff <- NA
marker_pos$dom_eff <- 0
marker_pos$qtl_name <- row.names(marker_pos)
marker_pos$qtl1_pair <- row.names(marker_pos)

# Duplicate by the number of traits
marker_pos_list <- replicate(n = n_traits, marker_pos, simplify = FALSE)
marker_pos_list[[1]]$qtl1_pair <- NA

# Add effects
for (i in seq_len(n_traits)) {
  marker_pos_list[[i]]$add_eff <- marker_eff[,-1][[i]]
}

## Create a new genome with markers as QTL
genome_use <- genome

# Add to the genome
genome_use$gen_model <- marker_pos_list

}


## Predict
predicted_genvar <- pred_genvar(genome = genome_use, pedigree = pedigree, founder.pop = founder.pop, 
                                    crossing.block = crossing.block)



genome <- genome_use





{

## How many traits
n_traits <- length(genome$gen_model)

# If it is more than 2, error out
stopifnot(n_traits <= 2)


## Calculate the expected genetic variance


## What are the expected allele frequencies in the population?
## Is there any backcrossing?
mom_ped <- pedigree[pedigree$mom == 1,]
dad_ped <- pedigree[pedigree$mom == 2,]

mom_dist_gen <- length(unique(mom_ped$gen))
dad_dist_gen <- length(unique(dad_ped$gen))

max_bc_gen <- pmax(mom_dist_gen, dad_dist_gen) - 1

# The expected frequency of the minor allele is 0.5 ^ n_bc_gen + 1
exp_q <- 0.5^(max_bc_gen + 1)
exp_p <- 1 - exp_q

# Get the QTL information - drop unused levels
qtl_info <- pull_qtl(genome, unique = FALSE)
# Filter out QTL with no additive effect
qtl_info <- droplevels(qtl_info[qtl_info$add_eff != 0,,drop = FALSE])
# Split by trait
qtl_info_split <- split(qtl_info, qtl_info$trait)


## Iterate over traits
qtl_covariance <- lapply(X = qtl_info_split, FUN = function(trait_qtl) {
  
  row.names(trait_qtl) <- trait_qtl[["qtl_name"]]
  
  ## Calculate the expected genetic variance and covariance of QTL
  qtl_info <- as.matrix(trait_qtl[,c("chr", "pos", "add_eff"), drop = FALSE])
  add_eff <- qtl_info[,"add_eff", drop = FALSE]
  pos <- qtl_info[,"pos", drop = FALSE]
  
  covar <- tcrossprod(add_eff)
  
  ## Create an empty matrix
  D <- matrix(0, nrow = nrow(pos), ncol = nrow(pos), dimnames = dimnames(covar))
  
  # Calculate separate distance matrices per chromosome
  chr_c <- lapply(X = split(trait_qtl, trait_qtl[,"chr",drop = FALSE]), FUN = function(x) as.matrix(dist(x[,"pos",drop = FALSE])))
  
  for (cr in chr_c) {
    cr2 <- qtl:::mf.h(cr)
    d <- ((1 - (2 * cr2)) / (1 + (2 * cr2)))
    D[row.names(cr), colnames(cr)] <- d
  }
  
  # The covariance is the QTL effect product multiplied by the expected D
  qtl_covar <- covar * D
  
})


if (n_traits > 1) {
  
  ## Calculate the genetic covariance between QTL for different traits
  # Split by chromosome
  qtl_chr_split <- split(qtl_info, qtl_info$chr)
  
  # Create an empty matrix of trait1 and trait2 QTL
  qtl_trait_covariance <- matrix(0, nrow = nrow(qtl_info_split[[1]]), ncol = nrow(qtl_info_split[[2]]),
                                 dimnames = list(qtl_info_split[[1]][["qtl_name"]], qtl_info_split[[2]][["qtl_name"]]))
  
  
  ## Iterate over chromosomes
  covar_list <- lapply(X = qtl_chr_split, FUN = function(chr_qtl) {
    
    # Split by trait
    trait_split <- split(chr_qtl, chr_qtl$trait)
    
    ## QTL names for each trait
    qtl_names <- lapply(X = trait_split, FUN = "[[", "qtl_name")
    qtl_pos <- lapply(X = trait_split, FUN = "[[", "pos")
    qtl_eff <- lapply(X = trait_split, FUN = function(q) as.matrix(q$add_eff))
    
    ## Calculate the pairwise distance
    d <- abs(outer(X = qtl_pos[[1]], Y = qtl_pos[[2]], FUN = `-`))
    # Calculate pairwise D (see Zhong and Jannink, 2007)
    # First convert cM to recombination fraction
    c <- qtl:::mf.h(d)
    D <- ((1 - (2 * c)) / (1 + (2 * c)))
    
    # Product of QTL effects
    qtl_crossprod <- tcrossprod(qtl_eff[[1]], qtl_eff[[2]])
    dimnames(qtl_crossprod) <- qtl_names
    
    # The covariance is the QTL effect product multiplied by the expected D
    qtl_crossprod * D
    
  })
  
  ## Add to the large matrix
  for (cov in covar_list) {
    qtl_trait_covariance[row.names(cov), colnames(cov)] <- cov
  }
  
} else {
  qtl_trait_covariance <- NULL
  
}


## Now we iterate over the parent pairs to determine the QTL that are segregating

# Replicate the crossing block 

## Add columns to the crossing.block for exp mu and exp varG
crossing_block <- crossing(crossing.block, trait = paste0("trait", seq(length(genome$gen_model))))
exp_mu <- list()
exp_varG <- list()
exp_corG <- list()



## Pull out the qtl genotypes for each trait
qtl_names <- lapply(X = qtl_info_split, FUN = "[[", "qtl_name")
qtl_geno <- lapply(X = qtl_names, function(q) pull_genotype(genome = genome, geno = founder.pop$geno, loci = q) - 1)

}

# Iterate over the crossing block
for (j in seq(nrow(crossing.block))) {
  
  pars <- as.character(crossing.block[j,1:2])
  
  ## Get a list of the polymorphic QTL
  poly_qtl_list <- lapply(X = qtl_geno, FUN = function(tr_qtl) {
    
    # Subset the parents
    par_qtl_geno <- tr_qtl[pars,,drop = FALSE]
    qtl_means <- colMeans(par_qtl_geno)
    par1_qtl <- par_qtl_geno[1,,drop = FALSE]
    
    par1_qtl[,qtl_means == 0, drop = FALSE]
    
  })
  
  
  # Iterate over the traits and calculate individual genetic variance
  trait_var <- mapply(poly_qtl_list, qtl_covariance, FUN = function(x, y) sum(crossprod(x) * y[colnames(x), colnames(x)]))
  
  
  if (!is.null(qtl_trait_covariance)) {
    
    ## Calculate the expected covariance
    trait_cov <- sum(qtl_trait_covariance[colnames(poly_qtl_list[[1]]), colnames(poly_qtl_list[[2]])] * crossprod(poly_qtl_list[[1]], poly_qtl_list[[2]]))
    
    # The expected correlation is calculated using the expected sd and expected cov
    exp_corG_j <- trait_cov / prod(sqrt(trait_var)) 
    exp_corG[[j]] <- rep(exp_corG_j, 2)
    
  }
  
  # The expected mu is simply the mean of the genotypic values of the two parents
  exp_mu_j <- colMeans(founder.pop$geno_val[founder.pop$geno_val$ind %in% pars,-1,drop = F])
  
  ## Add to the lists
  exp_mu[[j]] <- exp_mu_j
  exp_varG[[j]] <- trait_var
  
}

## Add the variances and means to the crossing block
crossing_block$exp_mu <- unlist(exp_mu)
crossing_block$exp_varG <- unlist(exp_varG)  
crossing_block$exp_corG <- unlist(exp_corG)

# Return the crossing block
return(crossing_block)

}













# PGVs
pgvs <- founder_pop1$pred_val

# Replace the expected mu with the predicted mu
for (i in seq_len(nrow(predicted_genvar))) {
  predicted_genvar$exp_mu[i] <- mean(pgvs[pgvs$ind %in% predicted_genvar[i,1:2,drop = T],predicted_genvar$trait[i]])
}


if (n_traits == 1) {
  names(predicted_genvar)[-1:-3] <- c("pred_mu", "pred_varG")
} else {
  names(predicted_genvar)[-1:-3] <- c("pred_mu", "pred_varG", "pred_corG")
}

# Return 
return(predicted_genvar)

}





### Increase marker density through HMM ###
###
###

library(pbsim)

# Load simulation data
load(file.path(gdrive_dir, "BarleyLab/Projects/SideProjects/Resources/s2_cap_simulation_data.RData"))

# Filter for breeding programs relevant to my data
s2_cap_genos <- s2_cap_genos[str_detect(string = row.names(s2_cap_genos), pattern = "AB|BA|WA|N2|MT"),]
s2_cap_genos <- s2_cap_genos[,!colMeans(s2_cap_genos) %in% c(0, 2)]
s2_snp_info <- subset(s2_snp_info, rs %in% colnames(s2_cap_genos))






















