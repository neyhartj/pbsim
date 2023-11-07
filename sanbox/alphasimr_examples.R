# AlphaSimR for plant breeding
# 
# Cases of using AlphaSimR to simulate a plant breeding program.
# These will be used to create helper functions for future simulations
# 

library(tidyverse)
library(AlphaSimR)
library(pbsim)
library(pbsimData)
library(qtl)

# Trait parameters
nQTL <- 300
# Narrow-sense heritability
h2 <- 0.5
# Additive genetic variance
varA <- 15
# Dominance variance
varD <- 3
# AxA epistatic variance
varAA <- 0.5
# GxE variance
varGE <- 10
# Environmental variance
varEnv <- (varA + varD + varAA) * 8

# Genome parameters
nMar <- 3000



# Simulate a genetic map
gen_map <- sim.map(len = seq(1.1, 1.7, length.out = 12), n.mar = rep((nQTL + nMar) / 12, 12), include.x = FALSE, eq.spacing = FALSE)
gen_map_df <- map2table(map = gen_map) %>%
  rownames_to_column("marker")

# Number of individuals
nInd <- 300

# simulate LE haplotypes
haplos <- t(replicate(n = nInd * 2, rbinom(n = nrow(gen_map_df), size = 1, prob = 0.5)))
colnames(haplos) <- gen_map_df$marker

# Load the genetic map and haplotypes into ASR
popHaplo <- importHaplo(haplo = haplos, genMap = gen_map_df, ploidy = 2)


# Create simulation parameters
SP <- SimParam$new(founderPop = popHaplo)


# Create a trait
SP$addTraitADEG(nQtlPerChr = nQTL / 12, mean = 0, var = varA, varEnv = varEnv,
                varGxE = varGE, meanDD = 0, varDD = varD, relAA = varAA/varA, 
                useVarA = TRUE, gamma = TRUE, name = "trait1")

# Set the residual variance via heritability
SP$setVarE(H2 = h2)

# Create a marker panel
SP$addSnpChip(nSnpPerChr = nMar / 12)


# Create a new phenotyped population from the founders
founders <- newPop(rawPop = popHaplo)






## Phenotype a population in multiple environments
# Range of p values for environments (this defines how far away the environments are from a target)
env_p_range <- c(0.2, 0.8)
# Number of environments
nEnv <- 5

# Generate phenotypes
multiEnvPheno <- sapply(X = runif(n = nEnv, min = min(env_p_range), max = max(env_p_range)), FUN = function(env_p) 
  setPheno(pop = founders, p = env_p, onlyPheno = TRUE) )

# Rename
colnames(multiEnvPheno) <- paste0(SP$traitNames, ".", paste0("env", seq_len(nEnv)))

# Tidy
multiEnvPheno_df <- `rownames<-`(multiEnvPheno, founders@id) %>%
  as.data.frame() %>%
  rownames_to_column("geno_id") %>%
  gather(trait_env, value, -geno_id) %>%
  separate(trait_env, c("trait", "env"), sep = "\\.")

# Model sanity check
fit <- lme4::lmer(value ~ env + (1|geno_id) + (1|geno_id:env), data = multiEnvPheno_df,
                  control = lme4::lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))



## Select individuals to use as parents



















