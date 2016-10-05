#' Make a population
#' 
#' @description 
#' Takes genptype information, a crossing block, and other parameters to
#' create multiple inbred families into one population.
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
#' @import dplyr
#' @import stringr
#' @import simcross
#' 
#' 
make.population1 <- function(map, crossing.block, parent.genos, n.ind, n.self.gen, 
                            m = 0, p = 0, cycle.number = 1) {
  
  # Formatting
  m <- as.integer(m)
  p <- as.double(p)
  
  # Find the number of crosses
  n.crosses <- nrow(crossing.block)
  
  # Create the pedigree template
  # This is the same across all families in the population
  family.pedigree <- make.pedigree(n.self.gen = n.self.gen, n.ind = n.ind)
  
  # ID's to extract from the pedigree
  prog.ids <- family.pedigree %>% 
    filter(gen == n.self.gen + 1) %>% 
    select(id) %>% 
    as.matrix()
  
  # Iterate over the number of crosses
  prog.genos <- apply(X = crossing.block, seq_len(nrow(crossing.block)), 
                      MARGIN = 1, FUN = function(parents, i) {
    
    genos <- parent.genos[parents,] %>%
      t()
    
    # Generate cross-over data
    xodata <- sim_from_pedigree_allchr(pedigree = family.pedigree, map = map, 
                                       m = m, p = p, obligate_chiasma = FALSE)
    
    # Generate genotypes
    prog.genos <- convert2geno_allchr(xodat = xodata, map = map, id = prog.ids,
                                      founder_geno = genos, return.matrix = T)
    
    # Rename the progeny lines
    row.names(prog.genos) <- paste("C", cycle.number, "_", n.self.gen + 1, 
                                   str_pad(i, width = 2, pad = 0), "-", 
                                   str_pad(seq(n.ind), width = 3, pad = 0), sep = "")
    
    return(prog.genos)
  })
    
  #   
  #   
  # 
  # # Generate crossover data for each pedigree
  # xo.list <- lapply(X = seq(n.crosses), FUN = function(i)
  #   sim_from_pedigree_allchr(pedigree = family.pedigree, map = map, m = m, p = p) )
  #   
  # 
  # 
  # # Iterate over the number of crosses in the crossing block
  # pop <- lapply(X = seq(n.crosses), FUN = function(cross.ind) {
  #   
  #   # Get the cross
  #   cross <- crossing.block[cross.ind,]
  #   
  #   # Get the parent names
  #   par1 <- cross[1] %>% 
  #     as.matrix() %>% 
  #     as.character()
  #   par2 <- cross[2] %>% 
  #     as.matrix() %>% 
  #     as.character()
  #   
  #   # Extract the parental genotypes
  #   parent1.genos <- parent.genos[par1,]
  #   parent2.genos <- parent.genos[par2,]
  #   
  #   # Make the family and return it
  #   make.family(map = map, parent1.genome = parent1.genos,
  #               parent2.genome = parent2.genos, n.gen = n.gen,
  #               n.ind = n.ind, cycle.number = cycle.number,
  #               family.number = cross.ind) })
  #   
  # 
  # # Rename the list entries
  # names(pop) <- apply(X = crossing.block, MARGIN = 1, FUN = function(cross) 
  #   return(paste(cross[1],cross[2], sep = ".")))
  
  
  # Return the list
  return(prog.genos)
  
} # Close the function