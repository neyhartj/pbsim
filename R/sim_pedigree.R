#' Create a pedigree for a cross 
#' 
#' @description 
#' Creates a pedigree for a bi-parental, recombinant inbred line family. Options
#' are available for partial selfing, backcrossing, and random mating.
#' 
#' @param n.ind The number of initial F_1 individuals that are created. By
#' single-seed decent, this is also the number of final inbred individuals.
#' @param bcpar The recurrant parent in backcrossing. If not \code{NULL} (default),
#' can be \code{1} or \code{2}.
#' @param n.bcgen The number of backcross generations to the parent specified
#' by \code{bc.par}.
#' @param n.selfgen The number of selfing generations. For instance, 
#' \code{n.self.gen = 5} would produce a F_6 family. If \code{Inf}, the pedigree
#' defines "complete" selfing. If not, selfing is "partial."
#' 
#' @details
#' 
#' @return 
#' A matrix of class \code{pedigree} with 5 columns. The first column contains
#' the individual id, the second column contains the individual ids for the mom,
#' the third column contains the individual ids for the dad, the fourth column
#' contains sex indicators, and the fifth column contains generations indicators.
#' 
#' @examples
#' # Create a BC2F5 family
#' ped <- sim_pedigree(n.ind = 50, bcpar = 1, n.bcgen = 2, n.selfgen = 4)
#' 
#' # Create a F3 family
#' ped <- sim_pedigree(n.ind = 50, n.bcgen = 0, n.selfgen = 2)
#' 
#' @import simcross
#' 
#' @export 
#' 
sim_pedigree <- function(n.ind, bcpar, n.bcgen = 0, n.selfgen = Inf) {
  
  # Error
  # n.ind must be integer
  if (!is.numeric(n.ind)) stop("The input 'n.ind' must be numeric.")
  
  # bcpar must be 1 or 2, if not NULL
  if (!missing(bcpar))
    if (!bcpar %in% c(1, 2))
      stop("The input 'bcpar' must be 1 or 2.")
  
  # Error if n.bcgen is greater than 0 and bcpar is NULL
  if (n.bcgen > 0 & missing(bcpar))
    stop("'bcpar' must be assigned if 'n.bcgen' is greater than 0.")
  
  # n.bcgen and n.selfgen must be numeric
  if (!is.numeric(n.bcgen)) stop("The input 'n.bcgen' must be numeric.")
  if (!is.numeric(n.selfgen)) stop("The input 'n.selfgen' must be numeric.")
  
  # n.bcgen and n.selfgen must not be negative
  if (n.bcgen < 0) stop("The input 'n.bcgen' must not be negative.")
  if (n.selfgen < 0) stop("The input 'n.selfgen' must not be negative.")
  
  # If n.selfgen is Inf, selfing is complete, otherwise selfing is partial
  if (n.selfgen == Inf) {
    selfing <- "complete"
    n.selfgen <- 0
    
  } else {
    selfing <- "partial"
    
  }
    
  # Determine the total number of entries in the pedigree
  n.par <- 2
  n.total <- n.par + (n.ind * (1 + n.bcgen + n.selfgen))
  
  # The id vector
  id <- seq_len(n.total)
  
  # A vector of sex indicators
  sex <- numeric(n.total)
  sex[2] <- 1
  
  n.gen <- 1 + n.bcgen + n.selfgen
  gen.ind <- seq_len(n.gen)
  gen <- numeric(n.total)
  gen[-c(1, 2)] <- rep(gen.ind, each = n.ind)
  
  mom <- dad <- numeric(n.total)
  
  # Assign parents for each entry in the pedigree
  # First assign the initial cross
  mom[gen == 1] <- 1
  dad[gen == 1] <- 2
  
  # Vector of bc or selfing generation indicators
  eff.gens <- setdiff(gen.ind, 1)
  
  # Extract the bc generation indicators
  bcind <- seq(2, 1 + n.bcgen, length.out = n.bcgen)
  # Extract the self generation indicators
  selfind <- setdiff(eff.gens, bcind)
  
  # Iterate over backcross generations
  for (bc in bcind) {
    
    # Assign the backcross parent - it's always the mom
    mom[gen == bc] <- bcpar
    # Assign the other parents
    dad[gen == bc] <- id[gen == (bc - 1)]

  }
  
  # Iterate over selfing generations
  for (s in selfind) {
    
    # Assign the selfing parents
    mom[gen == s] <- id[gen == (s - 1)]
    dad[gen == s] <- id[gen == (s - 1)]
    
  }
  
  # Combine into a matrix and return
  ped <- data.frame(id = id, mom = mom, dad = dad, sex = sex, gen = gen)
  
  # Set the selfing attribute
  attr(ped, "selfing") <- selfing
  
  return(ped)
  
} # Close the function