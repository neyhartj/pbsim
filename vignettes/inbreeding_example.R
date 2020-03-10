library(pbsim)
library(rrBLUP)
library(dplyr)

set.seed(118)

genome <- sim_gen_model(genome = sim_genome(len = rep(150, 7), n.mar = rep(505, 7)), qtl.model = matrix(NA, nrow = 10, ncol = 4), add.dist = "normal")

# F3 pedigree
pedF3 <- sim_pedigree(n.ind = 100, n.selfgen = 2)
# RIL pedigree
pedRIL <- sim_pedigree(n.ind = 100)

# Same founders
founders <- sim_founders(genome = genome)

## F3 family
familyF3 <- sim_family(genome = genome, pedigree = pedF3, founder.pop = founders)
## RIL family
familyRIL <- sim_family(genome = genome, pedigree = pedRIL, founder.pop = founders)

# Genotype both families
familyF3_geno <- genotype(genome = genome, pop = familyF3)
familyRIL_geno <- genotype(genome, familyRIL)

# A.mat for both matrices
familyF3_Amat <- A.mat(X = familyF3_geno, min.MAF = 0, max.missing = 1)
familyRIL_Amat <- A.mat(X = familyRIL_geno, min.MAF = 0, max.missing = 1)


# Mean of diagonal is 1 + f (inbreeding coefficient)
familyF3_f <- mean(diag(familyF3_Amat)) - 1
familyRIL_f <- mean(diag(familyRIL_Amat)) - 1

# Compare
c(F3 = familyF3_f, RIL = familyRIL_f)


## Measure allele frequencies in the RIL base pop
# Calculate the frequency of the 1 allele in the base training population
p_i <- apply(X = familyRIL_geno + 1, MARGIN = 2, FUN = mean) / 2
# Calculate the P matrix
P = matrix(2 * (p_i - 0.5))
# Calculate the normalization constant
c = 2 * sum(p_i * (1 - p_i))



## Select from the RIL family and mate
parents <- select_pop(pop = sim_phenoval(familyRIL, 0.2), intensity = 10, index = 1, type = "phenotypic")

## Produce 6 families
cycle2 <- sim_family_cb(genome = genome, pedigree = pedRIL, founder.pop = parents, cycle.num = 2,
                        sim_crossing_block(parents = indnames(parents), n.crosses = 5))

# Genotype
cycle2_geno <- genotype(genome, cycle2)

M <- cycle2_geno
# M <- familyRIL_geno

# Subtract P to make Z (need to convert P into a repeated matrix)
W = M - matrix(P, nrow(M), length(P), byrow = T)
# Calculate the relationship matrix
A = tcrossprod(W) / c
# Calculate average inbreeding
cycle2_f <- mean(diag(A)) - 1

cycles <- 5
popi <- cycle2
pop_f <- numeric(cycles); pop_f[1:2] <- c(familyRIL_f, cycle2_f)
pop_g <- numeric(cycles); pop_g[1:2] <- c(mean(familyRIL$geno_val$trait1), mean(cycle2$geno_val$trait1))

## Repeat
for (i in 3:cycles) {
  
  ## Select from the RIL family and mate
  parents <- select_pop(pop = sim_phenoval(popi, 0.2), intensity = 10, index = 1, type = "phenotypic")
  
  ## Produce 6 families
  popi <- sim_family_cb(genome = genome, pedigree = pedRIL, founder.pop = parents, cycle.num = i,
                          sim_crossing_block(parents = indnames(parents), n.crosses = 5))
  
  pop_g[i] <- mean(popi$geno_val$trait1)
  
  M <- genotype(genome, popi)

  # Subtract P to make Z (need to convert P into a repeated matrix)
  W = M - matrix(P, nrow(M), length(P), byrow = T)
  # Calculate the relationship matrix
  A = tcrossprod(W) / c
  # Calculate average inbreeding
  pop_f[i] <- mean(diag(A)) - 1
  
}




