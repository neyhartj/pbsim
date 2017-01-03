#' Simulate a family from a pedigree
#' 
#' @param genome An object of class \code{genome}.
#' @param pedigree An object of class \code{pedigree} detailing the scheme to 
#' develop the family.
#' @param founder_geno A matrix with \code{n.founder} rows and \code{n.mar}
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
#' founder_geno <- t(replicate(2, sample(x = c(0,2), size = sum(genome$n.mar), replace = T)))
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
  
  if (!inherits(pedigree, "pedigree"))
    stop("The input 'pedigree' must be of class 'pedigree.'")
  
  # Verify the founder genotypes
  # How many founders?
  n.founders <- nrow(founder_geno)
  
  # Are the founders coded correctly?
  if (!all(founder_geno %in% c(0, 1, 2)))
    stop("The input 'founder_geno' must be encoded in z {0, 1, 2}.")
  
  # Are the number of founders correct vis a vis the pedigree?
  if (n.founders != sum(pedigree[,5] == 0))
    stop("The number of founders in the inpute 'fouders' is not equal to the
         number of founders in the 'pedigree.'")
  
  # Are the number of markers in the founders correct?
  if (ncol(founder_geno) != sum(genome$n.mar))
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
  p <- ifelse(is.null(other.args$p), 0, other.args$p)
  
  # For naming
  family.num <- ifelse(is.null(other.args$family.num), 1, other.args$family.num)
  cycle.num <- ifelse(is.null(other.args$cycle.num), 1, other.args$cycle.num)
  
  
  # Generate cross-over data
  xo.data <- sim_from_pedigree_allchr(pedigree = pedigree, map = map, m = m, p = p)
  
  # Simulate genotypic data
  prog.genos <- convert2geno_allchr(xodat = xo.data, map = map, id = final.id)
  
  # Convert the progeny genotypes to the parental states
  prog.genos.recode <- apply(X = rbind(founder_geno, prog.genos), MARGIN = 2, FUN = function(snp)
    ifelse(snp[-1:-2] == 1, snp[1], ifelse(snp[-1:-2] == 3, snp[2], 1)) )
  
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
