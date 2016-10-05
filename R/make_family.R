#' Make a family
#' 
#' @param map A \code{list} of genetic positions of SNPs for each chromosome.
#' @param parent1.genome A vector of the parent 1 genotypes.
#' @param parent2.genome See \code{parent1.genome}
#' @param n.ind The number of F_1 individuals to generate.
#' @param n.gen The number of selfing generations of single-seed descent. 
#' For instance \code{generations = 5} would results in an F_6 population. 
#' For each generation, each of the F_1 individuals results in a single new 
#' individual that has undergone an additional generation of inbreeding.
#' @param m A cross-over interference parameter. See \code{\link[simcross]{
#' sim_from_pedigree}}.
#' @param cycle.number The integer number of the breeding cycle.
#' @param family.number The integer number of the family.
#' 
#' 
#' @import dplyr
#' @import stringr
#' @import simcross
#' 
make.family <- function(map, parent1.genome, parent2.genome, n.gen, n.ind,
                      m = 0, cycle.number = 1, family.number = 1) {
  
  # Create a pedigree
  family.pedigree <- make.pedigree(n.self.gen = n.gen, n.ind = n.ind)
  
  # Extract the index of the terminal taxa
  terminals <- family.pedigree %>% 
    filter(gen == n.gen + 1) %>% 
    select(id) %>%
    as.matrix() %>%
    as.numeric()
  
  # Create another vector of row.names to substitute
  new.names <- paste("C", cycle.number, "_", n.gen + 1, 
                     str_pad(family.number, width = 2, pad = 0), "-", 
                     str_pad(seq(n.ind), width = 3, pad = 0), sep = "")
  
  # Create a matrix of founder genotypes
  parents <- rbind(parent1.genome, parent2.genome) %>%
    t()
  
  # Generate cross-over data
  xo.data <- sim_from_pedigree_allchr(pedigree = family.pedigree, 
                                      map = map, m = m)
  
  # Simulate genotypic data
  prog.genos <- convert2geno_allchr(xodat = xo.data, map = map, id = terminals,
                                    founder_geno = parents)
  
  # Rename rows
  row.names(prog.genos) <- new.names
  
  # Return
  return(prog.genos)
  
} # Close the loop