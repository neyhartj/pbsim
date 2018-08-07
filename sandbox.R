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
