#' Calculate the genotypic value of individuals
#' 
#' @param genome An object of class \code{genome}.
#' @param geno Genotype data on a population to phenotype. Can be a matrix of dimensions
#' \code{n.ind} x \code{n.loci}, the elements of which must be z {0, 1, 2}, or a list
#' of such matrices.
#' 
#' @details 
#' Calculates the genotypic value of an individual or one or more traits by adding
#' the allele effects and dominance deviations of the QTL carried by that individual.
#' 
#' @return
#' A matrix of dimensions \code{n.ind} x \code{n.trait} with genotypic values
#' of the individuals.
#' 
#' @import dplyr
#' 
#' @export 
#' 
calc_genoval <- function(genome, geno) {
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = geno))
    stop("The geno did not pass. See warning for reason.")
  
  # If geno is a list, recombine
  if (is.list(geno))
    geno <- do.call("cbind", geno)
  
  ## Iterate over traits in the genetic model
  geno_val <- lapply(X = genome$gen_model, FUN = function(qtlmod) {
    qtl_geno <- pull_genotype(genome = genome, geno = geno, loci = qtlmod$qtl_name)
    # Subtract 1
    qtl_geno1 <- qtl_geno - 1
    
    # Additive eff
    bv <- qtl_geno1 %*% as.matrix(qtlmod$add_eff)
    # Dominance - first convert matrix to dominance matrix, then calculate deviations
    qtl_dom_geno <- ifelse(qtl_geno1 != 0, 0, 1)
    dd <- qtl_dom_geno %*% as.matrix(qtlmod$dom_eff)
    
    # Sum
    bv + dd })
  
  # Bind columns
  geno_val1 <- do.call("cbind", geno_val)
  
  # Add trait names
  colnames(geno_val1) <- paste("trait", seq(ncol(geno_val1)), sep = "")
  
  # Convert the row.names to a column and return the data.frame
  data.frame(ind = row.names(geno_val1), geno_val1, row.names = NULL, stringsAsFactors = FALSE) %>%
    # Sort on individual name
    arrange(ind)
  
} # Close the function



  
  
  
  
  
  



