#' Simulate a family from a pedigree
#' 
#' @param genome An object of class \code{genome}.
#' @param pedigree A \code{pedigree} detailing the scheme to develop the family.
#' Use \code{\link{sim_pedigree}} to generate.
#' @param founder_pop An object of class \code{pop} with the geno information for
#' the founders. Use the \code{\link{subset.pop}} function to subset a \code{pop}
#' object.
#' 
#' @details 
#' Other arguments can be passed to internal functions in \code{sim_family}:
#' \describe{
#'   \item{\code{dh}}{A logical indicating if double-haploid lines should be created
#'   at the end of selfing. Doubled-haploids are generated at the last generation
#'   of selfing (e.g. if \emph{F_3} individuals are specified in the pedigree, DH
#'   lines are induced after the \emph{F_2}).}
#'   \item{\code{m}}{The crossover interference parameter. See 
#'   \code{\link[simcross]{sim_from_pedigree}} for more information. }
#'   \item{\code{p}}{The proportion of crosses from non-interference process.
#'   See \code{\link[simcross]{sim_from_pedigree}} for more information.}
#'   \item{\code{family.num}}{A integer designator of the family.}
#'   \item{\code{cycle.num}}{A integer designator of the breeding cycle.}
#' }
#' 
#' @return 
#' An object of class \code{pop}.
#' 
#' @examples 
#' 
#' # Simulate a genome
#' n.mar  <- c(505, 505, 505)
#' len <- c(120, 130, 140)
#' 
#' genome <- sim_genome(len, n.mar)
#' 
#' # Simulate a quantitative trait influenced by 50 QTL
#' qtl.model <- matrix(NA, 50, 4)
#' genome <- sim_gen_model(genome = genome, qtl.model = qtl.model, 
#'                         add.dist = "geometric", max.qtl = 50)
#' 
#' # Simulate the founder genotypes
#' founder_pop <- sim_founders(genome)
#' 
#' # Create a pedigree with 100 individuals selfed to the F_3 generation
#' ped <- sim_pedigree(n.ind = 100, n.selfgen = 2)
#' 
#' fam <- sim_family(genome = genome, pedigree = ped, founder_pop = founder_pop)
#' 
#' # Create a pedigree with 100 RIL individuals
#' ped <- sim_pedigree(n.ind = 100, n.selfgen = Inf)
#' 
#' fam <- sim_family(genome = genome, pedigree = ped, founder_pop = founder_pop)
#' 
#' # Create a pedigree with 100 doubled-haploid individuals induced after the F_2
#' # generation.
#' ped <- sim_pedigree(n.ind = 100, n.selfgen = 2)
#' 
#' fam <- sim_family(genome = genome, pedigree = ped, founder_pop = founder_pop,
#'                   dh = TRUE)
#' 
#'  
#' @import simcross
#' @importFrom stringr str_pad
#' @importFrom stringr str_c
#' @importFrom qtl sim.cross
#' 
#' @export
#' 
sim_family <- function(genome, pedigree, founder_pop, ...) {
  
  # Error
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Founder_pop needs to be a pop object
  if (!inherits(founder_pop, "pop"))
    stop("The input 'founder_pop' must be of class 'pop'")
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = founder_pop$geno))
    stop("The geno did not pass. See warning for reason.")
  
  # Check the pedigree
  if (!check_pedigree(pedigree, ignore_sex = TRUE))
    stop("The pedigree is not formatted correctly.")
  
  
  # Combine the founder geno input
  founder_geno <- do.call("cbind", founder_pop$geno)
  
  # How many founders?
  n_founders <- nrow(founder_geno)
  
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
  map <- genome$map
  
  # Parse other arguments
  other.args <- list(...)
  
  # For xodata
  m <- ifelse(is.null(other.args$m), 0, other.args$m)
  p <- ifelse(is.null(other.args$p), 1, other.args$p)
  
  # For naming
  family_num <- ifelse(is.null(other.args$family.num), 1, other.args$family.num)
  cycle_num <- ifelse(is.null(other.args$cycle.num), 1, other.args$cycle.num)
  
  # For doubled-haploids
  dh <- ifelse(is.null(other.args$dh), FALSE, other.args$dh)
  
  # If selfing is partial, using simcross
  selfing <- attr(pedigree, "selfing")
  
  if (selfing == "partial") {
  
    # Generate cross-over data
    xo_data <- sim_from_pedigree_allchr(pedigree = pedigree, map = map, m = m, p = p)
    
    # Simulate DH if called for
    if (dh) {
      xo_data <- induce_dh(xodat = xo_data, pedigree = pedigree)
    }
    
    # Simulate genotypic data
    prog_genos <- convert2geno_allchr(xodat = xo_data, map = map, id = final_id)
    
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
  
  # Create a matrix with the founder1, het, and founder2 genotypes
  founder_multipoint <- rbind(founder_geno[1,], colMeans(founder_geno), founder_geno[2,])
  
  # Convert the progeny genotypes to the parental states
  prog_genos_recode <- t(apply(X = prog_genos, MARGIN = 1, FUN = function(prog) {
    founder_multipoint[cbind(prog, seq(nrow(founder_multipoint)))] }))

  # Generate new progeny names
  n_ind <- sum(pedigree[,5] == max(pedigree[,5])) # Number of individuals
  gen <- max(pedigree[,5]) # Generation number
  
  # New names
  new_names <- str_c("C", cycle_num, "_", gen, str_pad(family_num, width = 3, pad = 0), 
                     "-", str_pad(seq_len(n_ind), width = 3, pad = 0))
  
  dimnames(prog_genos_recode) <- list(new_names, markernames(genome, include.qtl = TRUE))
  
  # Create the pop
  create_pop(genome = genome, geno = prog_genos_recode)
  
} # Close the function
