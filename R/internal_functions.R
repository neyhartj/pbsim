#' Simulate recombination using hypred
#' 
#' @description 
#' Internal function. Not generally to be called by the user.
#' 
#' @param genome An object of class \code{genome}.
#' @param haploids A list of length \code{n.chr}, the elements of which are matrices
#' of dimensions 2 x \code{n.loci} giving the haploid genotypes for a single individual
#' for a single chromosome.
#' @param mutate A logical as to generate random mutations when recombining.
#' @param mutation.rate.snp The per-base mutation rate of the SNPs.
#' @param mutation.rate.qtl The per-base mutation rate of the QTL.
#' 
#' @return 
#' A list of recombined haploid gametes (of class \code{matrix}) per chromosome.
#' 
#' @import hypred
#' @importFrom  purrr pmap
#' 
recombine_hypred <- function(genome, haploids, mutate = FALSE, 
                             mutation.rate.snp = 0, mutation.rate.qtl = 0) {

  # Sanity check
  if (any(sapply(haploids, nrow) != 2))
    stop("Each matrix in 'haploids' must have 2 rows.")
  
  if (!is.na(dim(haploids[[1]])[3]))
    stop("Only one individual must be represented in the haploids.")
  
  # Iterate over each chromosome
  pmap(list(genome$hypredGenomes, haploids), function(chr_genome, chr_haploids) {
    
    gamete <- hypredRecombine(object = chr_genome, genomeA = chr_haploids[1,],
                              genomeB = chr_haploids[2,], mutate = mutate, block = FALSE,
                              mutation.rate.snp = mutation.rate.snp, mutation.rate.qtl = mutation.rate.qtl)
    
    colnames(gamete) <- names(chr_genome@pos.snp)
    return(gamete)
    
  })
  
} # Close the function