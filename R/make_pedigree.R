#' Make a pedigree structure
#' 
#' @description 
#' Creates a pedigree for a bi-parental, recombinant inbred line family.
#' Selfing is done in the true sense. 
#' 
#' @param n.self.gen The number of selfing generations. For instance, 
#' \code{n.self.gen = 5} would produce a F_6 family
#' @param n.ind The number of initial F_1 individuals that are created. By
#' single-seed decent, this is also the number of final inbred individuals.
#' 
#' @export 
#' 
make.pedigree <- function(n.self.gen, n.ind) {
  
  # Determine the total number of individuals in the pedigree
  n.total <- 2 + (n.ind * (1 + n.self.gen))
  
  # The id vector
  id <- seq(n.total)
  
  # A vector of sex indicators
  sex <- rep(0, n.total)
  sex[2] <- 1
  
  # Vector of generation indicators
  n.gen <- n.self.gen + 1
  gen <- rep(0, n.total)
  gen[-c(1,2)] <- rep(seq(n.gen), each = n.ind)
  
  mom <- rep(0, n.total)
  dad <- mom
  
  for (g in seq(n.gen)) {
    
    # Fix for the first generation
    if (g == 1) {
      mom[gen == g] <- 1
      dad[gen == g] <- 2
    
    } else {
      
      mom[gen == g] <- id[gen == g-1]
      dad[gen == g] <- id[gen == g-1]
      
    }
  } # Close the loop
  
  # Turn the detailed pedigree into something compatable with qtl
  # abd return it
  data.frame(id = id, mom = mom, dad = dad, sex = sex, gen = gen)

} # Close the function