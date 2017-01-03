## Testing sim.cross capability

# Let's try to speed this up using rQTL and simcross
# library(qtl)
library(simcross)
library(qgsim)
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
























# Document and install
if(basename(getwd()) == "quantgen") setwd("..")
devtools::document("quantgen")
devtools::install("quantgen")


# quantgen sandbox
library(quantgen)

n.chr = 10
chr.len = seq(1.0, 1.45, by = 0.05)
chr.snps = sample(15:20, size = 10, replace = T)


# Make a genome
test.genome <- make.genome(n.chr = n.chr, chr.len = chr.len, chr.snps = chr.snps)

n.qtl = 50
genome.object <- test.genome
qtl.ind = NULL
qtl.add.eff = "geometric"
qtl.dom.eff = NULL
qtl.dom.ind = NULL

# Assign QTL
test.genome <- trait.architecture(genome.object = test.genome, n.qtl = 50)

genome.object <- test.genome

# Create two founder genotypes
founders <- make.founders(genome.object = test.genome)

hap.genome1 <- founders[1,]
hap.genome2 <- founders[2,]

mutate = T
mu.snp = 7e-8
mu.qtl = 7e-8

# Recombine
new.gamete <- recombine(genome.object = genome.object, hap.genome1 = hap.genome1,
                        hap.genome2 = hap.genome2, mutate = T, mu.snp = mu.snp, mu.qtl = mu.qtl)










#### Benchmarking ####

rec.prof <- lineprof(recombine(genome.object = genome.object, hap.genome1 = hap.genome1,
                               hap.genome2 = hap.genome2, mutate = T, mu.snp = mu.snp, mu.qtl = mu.qtl)
)








## 3 chromosome of lengths 1, 1.1, and 1.2 M and 100 SNP each
genome.hypred <- hypredGenome(num.chr = 3, len.chr = c(1.0, 1.1, 1.2), num.snp.chr = 100)

genome.quantgen <- make.genome(n.chr = 3, chr.len = c(1.0, 1.1, 1.2), chr.snps = c(100, 100, 100))

## produce two haploid founder line genomes
founder <- make.founders(genome.object = genome.quantgen)

# Create 1000 gametes
system.time(gametes.hypred <- replicate(10000, hypredRecombine(genome.hypred, genomeA = founder[1,], genomeB = founder[2,], mutate = T, mutation.rate.snp = 7e-8, mutation.rate.qtl = 7e-8, block = FALSE)))
system.time(gametes.quantgen <- replicate(10000, recombine(genome.quantgen, founder[1,], founder[2,], mutate = T, mu.snp = 7e-8, mu.qtl = 7e-8)))



print(gamete)



n.iter = 10000

system.time( 
  { empty <- vector("numeric", n.iter)
  for (i in 1:n.iter) {
    empty[[i]] <- rnorm(n = 1)
  }
  })











