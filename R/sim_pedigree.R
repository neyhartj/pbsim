#' Create a pedigree for a cross 
#' 
#' @description 
#' Creates a pedigree for a bi-parental, recombinant inbred line family. Options
#' are available for partial selfing, backcrossing, and random mating.
#' 
#' @param n.par The number of founding parents. May be 2 (for a 2-way cross) or 
#' 4 (for a 4-way cross).
#' @param n.ind The number of initial F_1 individuals that are created. By
#' single-seed decent, this is also the number of final inbred individuals.
#' @param bcpar The recurrant parent in backcrossing. If not \code{NULL} (default),
#' must point to one of the two parents (for a 2-way cross) or four parents (for a
#' 4-way cross).
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
#' @import simcross
#' 
#' @export 
#' 
sim_pedigree <- function(n.par, n.ind, bcpar, n.bcgen = 0, n.selfgen = Inf) {
  
  # Error
  stopifnot(n.par %in% c(2, 4))
  
  # n.ind must be integer
  if (!is.numeric(n.ind)) stop("The input 'n.ind' must be numeric.")
  
  # If bcpar is not missing, it must be <= n.par
  if (!missing(bcpar)) {
    if (bcpar > n.par) stop("The bcpar must be one of the parents.")
  }
  
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
    selfing <- "infinite"
    n.selfgen <- 0
    
  } else {
    selfing <- "finite"
    
  }
    

  ## Calculate a parent quotient, for lack of better term
  n.par.quo <- n.par / 2
  n.par.quo1 <- ifelse(n.par == 4, 2, 0)
  # Determine the total number of entries in the pedigree
  n.total <- n.par + n.par.quo1 + (n.ind * (1 + n.bcgen + n.selfgen))

  # The id vector
  id <- seq_len(n.total)
  
  # A vector of sex indicators
  sex <- mom <- dad <- gen <- numeric(n.total)
  # Designate 1 (and 3, if applicable) as sex 1
  sex[(n.par %/% c(4, (4/3)))] <- 1
  
  ## Vector of generation numbers:
  ## Generation 0 = parents
  ## Generation 1 = crosses involving parents
  ## (the number of ids in this category is (n.par / 2)
  ## Generation 2:n = inbreeding
  ## 

  ## Create replication amounts for each generation
  par.rep <- n.par
  int.rep <- n.par.quo1
  f1.rep <- n.ind
  
  ## Generate generation reps
  i = 0
  gen1 <- rep(i, par.rep); i = i + 1
  if (int.rep > 0) {gen2 <- rep(i, int.rep); i = i + 1} else gen2 <- NULL
  gen3 <- rep(i, f1.rep); i = i + 1
  
  # Rep backrosses
  if (n.bcgen > 0) {gen4 <- rep(seq(i, i + n.bcgen - 1), each = n.ind); i = i + n.bcgen + 1} else gen4 <- NULL
  # Rep selfs
  if (n.selfgen > 0) {gen5 <- rep(seq(i, i + n.selfgen - 1), each = n.ind); i = i + n.selfgen + 1} else gen5 <- NULL
  
  # Combine gens
  gen <- c(gen1, gen2, gen3, gen4, gen5)
  

  # Assign parents for each entry in the pedigree
  # First assign the initial cross
  mom[gen == 1] <- seq(1, n.par, by = 2)
  dad[gen == 1] <- seq(2, n.par, by = 2)
  
  # Assign F1 generations, if not covered above
  if (!is.null(gen2)) {
    
    mom[gen == unique(gen3)] <- which(mom == 1)
    dad[gen == unique(gen3)] <- which(mom == 3)
    
  }
  
  # Assign backcross parents
  if (!is.null(gen4)) {
    # Vector of bc indices
    bcind <- unique(gen4)
    
    # Iterate over bc generations
    for (bc in bcind) {
      
      # Assign the backcross parent - it's always the mom
      mom[gen == bc] <- bcpar
      # Assign the other parents
      dad[gen == bc] <- id[gen == (bc - 1)]
      
    }
  }
  
  # Assign selfing parents
  if (!is.null(gen5)) {
    # Vector of selfing indices
    selfind <- unique(gen5)
    
    # Iterate over selfing generations
    for (s in selfind) {
      
      # Assign the selfing parents
      mom[gen == s] <- id[gen == (s - 1)]
      dad[gen == s] <- id[gen == (s - 1)]
      
    }
  }
  
  
  # Combine into a matrix and return
  ped <- data.frame(id = id, mom = mom, dad = dad, sex = sex, gen = gen,
                    observed = gen == max(gen))
  
  # Set the selfing attribute
  attr(ped, "selfing") <- selfing
  
  return(ped)
  
} # Close the function
