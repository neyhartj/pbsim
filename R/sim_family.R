#' Simulate a family from a pedigree
#' 
#' @param genome An object of class \code{genome}.
#' @param pedigree An object of class \code{pedigree} detailing the scheme to 
#' develop the family.
#' @param founder_geno A matrix of dimension \code{n.founder} rows and \code{n.mar}
#' columns with founder genotypes. Must be coded as z = {0, 1, 2}. Can also be a list
#' of such matrices, where each matrix is the genotypic data from a chromosome.
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
#' @return 
#' An object of class \code{pop}.
#' 
#'  
#' @import simcross
#' @importFrom stringr str_pad
#' @importFrom stringr str_c
#' @importFrom qtl sim.cross
#' 
#' @export
#' 
sim_family <- function(genome, pedigree, founder_geno, ...) {
  
  # Error
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = founder_geno))
    stop("The geno did not pass. See warning for reason.")
  
  # Check the pedigree
  if (!check_pedigree(pedigree, ignore_sex = TRUE))
    stop("The pedigree is not formatted correctly.")
  
  
  # If the geno input is a list, recombine
  if (is.list(founder_geno))
    founder_geno <- do.call("cbind", founder_geno)
  
  # How many founders?
  n_founders <- nrow(founder_geno)
  
  # Are the founders coded correctly?
  if (!all(unlist(founder_geno) %in% c(0, 1, 2)))
    stop("The input 'founder_geno' must be encoded in z {0, 1, 2}.")
  
  # Are the number of founders correct vis a vis the pedigree?
  if (n_founders != sum(pedigree$gen == 0))
    stop("The number of founders in the inpute 'fouders' is not equal to the
         number of founders in the 'pedigree.'")
  
  # Are the number of markers in the founders correct?
  if (ncol(founder_geno) != nloci(genome))
    stop("The number of markers in 'founder_geno' is not equal to the number
         of markers in the genome.")
  
  
  # Extract the individual ids of the finals
  final_id <- subset(pedigree, gen == max(gen))$id
  
  # Extract the map
  blank_cross <- qtl::sim.cross(map = genome$map, n.ind = 1)
  map <- qtl::pull.map(cross = blank_cross)
  
  # Parse other arguments
  other.args <- list(...)
  
  # For xodata
  m <- ifelse(is.null(other.args$m), 0, other.args$m)
  p <- ifelse(is.null(other.args$p), 1, other.args$p)
  
  # For naming
  family_num <- ifelse(is.null(other.args$family.num), 1, other.args$family.num)
  cycle_num <- ifelse(is.null(other.args$cycle.num), 1, other.args$cycle.num)
  
  # If selfing is partial, using simcross
  selfing <- attr(pedigree, "selfing")
  
  if (selfing == "partial") {
  
    # Generate cross-over data
    xo_data <- sim_from_pedigree_allchr(pedigree = pedigree, map = map, m = m, p = p)
    
    # Simulate genotypic data
    prog_genos <- convert2geno_allchr(xodat = xo.data, map = map, id = final_id)
    
  } else {
    
    # Otherwise use sim.cross
    cross_sim <- qtl::sim.cross(map = genome$map, n.ind = length(final_id), type = "riself", 
                                m = m, p = p)
    
    # Extract progeny genos
    prog_genos <- lapply(X = cross_sim$geno, FUN = function(geno_chr) 
      ifelse(geno_chr$data == 1, 1, ifelse(geno_chr$data == 2, 3, NA)) )
    
    # Bind
    prog_genos <- do.call("cbind", prog_genos)
    
    # Add row.names
    row.names(prog_genos) <- final_id
    
  }
  
  # Convert the progeny genotypes to the parental states
  prog_genos_recode <- apply(X = rbind(founder_geno, prog_genos), MARGIN = 2, FUN = function(snp)
    ifelse(snp[-seq(n_founders)] == 1, snp[1], ifelse(snp[-seq(n_founders)] == 3, snp[2], 1)) )
  
  # Generate new progeny names
  n_ind <- sum(pedigree[,5] == max(pedigree[,5])) # Number of individuals
  gen <- max(pedigree[,5]) # Generation number
  
  # New names
  new_names <- str_c("C", cycle_num, "_", gen, str_pad(family_num, width = 3, pad = 0), 
                     "-", str_pad(seq_len(n_ind), width = 3, pad = 0))
  
  row.names(prog_genos_recode) <- new_names
  
  # Create the pop
  create_pop(genome = genome, geno = prog_genos_recode)
  
} # Close the function
