#' Generate founder inbreds
#' 
#' @description Generate a set of two completely inbred haploid parent genomes 
#' as founders. The alleles of the first founder are the complement of those of 
#' the second founder. Allele states can be generated randomly or alternating.
#' 
#' @param genome.object An object of class \code{genome}.
#' @param method Character describing how the allele states should be decided. 
#' If \code{"random"}, allele states are generating from a binomial distribution 
#' with success probability \code{theta}. If \code{"alternate"}, the first 
#' founder will have the "1" allele at odd-numbered loci, and the second founder 
#' will have the "1" allele at even-numbered loci.
#' @param theta Floating point probability of any one locus having the "1" 
#' allele in the first founder.
#' 
#' @details 
#' Since the genotypes generated from this function reflect completely inbred
#' founder lines, only the haploid genome is returned for each inbred.
#' 
#' @export
#' @return A \code{matrix} of dimensions 2 x m where m is the number of loci.
#' 
make.founders <- function(genome.object, method = "random", theta = 0.5) {
  
  # Error handling
  method <- as.character(method)
  theta <- as.numeric(theta)
  
  if (!method %in% c("random", "alternate")) stop("The input method must be either 'random' or 'alternate'.")
  
  # Pull out the number of loci across the genome
  n.snps <- slot(genome.object, "n.snps")
  
  # Generate haploid genomes based on the method
  if (method == "random") {
    founder.1 <- rbinom(n = n.snps, size = 1, prob = theta)
  }
  if (method == "alternate") {
    founder.1 <- rep(0, n.snps)
    founder.1[seq(1, n.snps, 2)] <- 1
  }
  
  # Generate founder2
  founder.2 <- 1 - founder.1
  
  # Generate the matrix
  return(rbind(founder.1, founder.2))
}