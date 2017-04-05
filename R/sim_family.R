#' Simulate a family from a pedigree
#' 
#' @param genome An object of class \code{genome}.
#' @param pedigree An object of class \code{pedigree} detailing the scheme to 
#' develop the family.
#' @param founder_geno A list of the same length as there are chromosomes in the 
#' \code{genome}. Each element should be a matrix with \code{n.founder} rows and \code{n.mar}
#' columns with founder genotypes. Must be coded as z = {0, 1, 2}.
#' @param ... Additional arguments. See \code{Details}. 
#' 
#' @details 
#' Other arguments can be passed to internal functions in \code{sim_family}:
#' \itemize{
#'   \item{\code{m}: The crossover interference parameter. See 
#'   \code{\link[simcross]{sim_from_pedigree}} for more information. }
#'   \item{\code{p}: The proportion of crosses from non-interference process.
#'   See \code{\link[simcross]{sim_from_pedigree}} for more information.}
#'   \item{\code{family.num}: A integer designator of the family.}
#'   \item{\code{cycle.num}: A integer designator of the breeding cycle.}
#' }
#' 
#' @examples 
#' # Simulate a genome
#' n.mar  <- c(505, 505, 505)
#' len <- c(120, 130, 140)
#' 
#' genome <- sim_genome(len, n.mar)
#' 
#' # Simulate a pedigree
#' ped <- sim_pedigree(n.ind = 50, n.bcgen = 0, n.selfgen = 2)
#' 
#' # Simulate the founder genotypes
#' founder_geno <- sim_founders(genome)
#' 
#' # Create the family
#' fam <- sim_family(genome, ped, founder_geno)
#'  
#' @import simcross
#' @importFrom stringr str_pad
#' 
#' @export
#' 
sim_family <- function(genome, pedigree, founder_geno, ...) {
  
  # Error
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Check the pedigree
  if (!check_pedigree(ped, ignore_sex = TRUE))
    stop("The pedigree is not formatted correctly.")
  
  # Verify the founder genotypes
  # How many founders?
  n.founders <- sapply(founder_geno, nrow)
  
  # Make sure this vector has the same elements
  if (length(unique(n.founders)) > 1) {
    stop("The 'founder_geno' input does not have a consistent number of founders.")
    
  } else {
    n.founders <- unique(n.founders)
    
  }
  
  # Are the founders coded correctly?
  if (!all(unlist(founder_geno) %in% c(0, 1, 2)))
    stop("The input 'founder_geno' must be encoded in z {0, 1, 2}.")
  
  # Are the number of founders correct vis a vis the pedigree?
  if (n.founders != sum(pedigree[,5] == 0))
    stop("The number of founders in the inpute 'fouders' is not equal to the
         number of founders in the 'pedigree.'")
  
  # Are the number of markers in the founders correct?
  if (!all.equal(target = unname(sapply(founder_geno, ncol)), current = genome$n.mar))
    stop("The number of markers in 'founder_geno' is not equal to the number
         of markers in the genome.")
  
  
  # Extract the individual ids of the finals
  final.id <- pedigree[,1][pedigree[,5] == max(pedigree[,5])]
  
  # Extract the map
  map <- genome$map
  
  
  # Parse other arguments
  other.args <- list(...)
  
  # For xodata
  m <- ifelse(is.null(other.args$m), 0, other.args$m)
  p <- ifelse(is.null(other.args$p), 1, other.args$p)
  
  # For naming
  family.num <- ifelse(is.null(other.args$family.num), 1, other.args$family.num)
  cycle.num <- ifelse(is.null(other.args$cycle.num), 1, other.args$cycle.num)
  
  # If selfing is partial, using simcross
  selfing <- attr(pedigree, "selfing")
  
  if (selfing == "partial") {
  
    # Generate cross-over data
    xo.data <- sim_from_pedigree_allchr(pedigree = pedigree, map = map, m = m, p = p)
    
    # Simulate genotypic data
    prog.genos <- convert2geno_allchr(xodat = xo.data, map = map, id = final.id)
    
  } else {
    
    # Otherwise use sim.cross
    cross_sim <- sim.cross(map = genome$map, n.ind = length(final.id), type = "riself", 
                           m = m, p = p)
    
    # Extract progeny genos
    prog.genos <- lapply(X = cross_sim$geno, FUN = function(geno_chr) 
      ifelse(geno_chr$data == 1, 1, ifelse(geno_chr$data == 2, 3, NA)) )
    
    # Bind
    prog.genos <- do.call("cbind", prog.genos)
    
    # Add row.names
    row.names(prog.genos) <- final.id
    
  }
  
  # Collapse the founder genotypes
  founder_geno <- do.call("cbind", founder_geno)
  
  # Convert the progeny genotypes to the parental states
  prog.genos.recode <- apply(X = rbind(founder_geno, prog.genos), MARGIN = 2, FUN = function(snp)
    ifelse(snp[-seq(n.founders)] == 1, snp[1], ifelse(snp[-seq(n.founders)] == 3, snp[2], 1)) )
  
  # Generate new progeny names
  n.ind <- sum(pedigree[,5] == max(pedigree[,5])) # Number of individuals
  gen <- max(pedigree[,5]) # Generation number
  
  # New names
  new.names <- str_c("C", cycle.num, "_", gen, str_pad(family.num, width = 3, pad = 0), 
                     "-", str_pad(seq_len(n.ind), width = 3, pad = 0))
  
  row.names(prog.genos.recode) <- new.names
  
  # Return
  return(prog.genos.recode)
  
} # Close the function
