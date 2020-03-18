#' Check a matrix of genotypes for completeness
#' 
#' @param genome An object of class \code{genome}.
#' @param geno Genotype data on a population to phenotype. Must be a matrix of dimensions 
#' \code{n.ind} x \code{n.loci}, the elements of which must be z {0, 1, 2}, or a 
#' list of such matrices.
#' @param ignore.gen.model Logical - should the genetic model be ignored when
#' checking the matrix of genotypes?
#' 
#' @importFrom abind abind
#' 
#' @examples 
#' 
#' n.mar  <- c(505, 505, 505)
#' len <- c(120, 130, 140)
#' 
#' genome <- sim_genome(len, n.mar)
#' 
#' # Randomly generate 15 QTL with additive allelic effects following a
#' # genometric series
#' qtl.model <- matrix(nrow = 15, ncol = 4)
#' genome <- sim_gen_model(genome, qtl.model, add.dist = "geometric")
#' 
#' # Create a random population
#' pop <- sim_pop(genome = genome, n.ind = 200)
#' 
#' check_geno(genome, pop$geno)
#' 
#' @export
#' 
check_geno <- function(genome, geno, ignore.gen.model = FALSE) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  
  # Extract the number of markers
  n_marker <- nmar(genome)
  
  # The geno input may have n_marker + n_qtl - n_perf_mar columns
  tot_loci <- nloci(genome)
  
  
  # If the geno input is a list, recombine
  if (is.list(geno))
    geno <- do.call("cbind", geno)
  
  # # Make sure the genos are coded correctly
  # if (!all(geno %in% c(0, 1, 2))) {
  #   warning("The input 'geno' must be encoded in z {0, 1, 2}.") 
  #   return(FALSE)
  # }
  
  # Make sure the geno input has n_marker columns or n_marker + n_qtl columns
  if (ncol(geno) != tot_loci & ncol(geno) != n_marker) {
    warning("The number of loci in the geno input does not equal the number of markers or the number of markers plus the number of QTL in the genome.") 
    return(FALSE)
  }
  
  # Does the geno matrix have marker names
  markers <- colnames(geno)
  if (is.null(markers)) {
    warning("No marker names in the geno input.") 
    return(FALSE)
  }
  
  # Get the marker names, but include QTL if the gene model should be ignored.
  if (ignore.gen.model) {
    genome_markers <- markernames(genome, include.qtl = TRUE)
  } else {
    genome_markers <- markernames(genome, include.qtl = FALSE)
  }
  
  # Are the marker names consistent with the genome?
  if (!all(genome_markers %in% markers)) {
    warning("The marker names in the genome are not consistent with those in the geno input.") 
    return(FALSE)
  }
  
  # All clear
  return(TRUE)
  
} # Close the function



#' Split a genotype matrix into chromosomes
#' 
#' @param genome An object of class \code{genome}.
#' @param geno Genotype data on a population to phenotype. Must be a matrix of 
#' dimensions \code{n.ind} x \code{n.loci}, the elements of which must be z {0, 1, 2}, 
#' or a list of such matrices.
#' @param ignore.gen.model Logical - should the genetic model be ignored?
#' 
#' @return 
#' A list of geno matrices, split by chromosome.
#' 
#' @export
#' 
split_geno <- function(genome, geno, ignore.gen.model = FALSE) {
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = geno, ignore.gen.model = ignore.gen.model))
    stop("The geno did not pass. See warning for reason.")
  
  # If the geno input is a list, recombine
  if (is.list(geno))
    geno <- do.call("cbind", geno)
  
  # Total number of loci
  n_loci <- nloci(genome = genome, by.chr = TRUE)
  # Total number of markers
  n_marker <- nmar(genome = genome, by.chr = TRUE)
  
  # Create a splitting vector based on the number of columns in the geno
  if (ncol(geno) == sum(n_loci)) {
    
    # Index of all loci
    loci_ind <- seq(sum(n_loci))
    # Split vector
    chr_split <- mapply(seq(nchr(genome)), n_loci, FUN = rep)
    
    split_vec <- split(x = loci_ind, f = unlist(chr_split))
    names(split_vec) <- chrnames(genome)
  
  } else if (ncol(geno) == sum(n_marker)) {
    
    # Index of all markers
    loci_ind <- seq(sum(n_marker))
    # Split vector
    chr_split <- mapply(seq(nchr(genome)), n_marker, FUN = rep)
    
    split_vec <- split(x = loci_ind, f = unlist(chr_split))
    names(split_vec) <- chrnames(genome)
    
  } else {
    stop("The number of columns in the geno input is not equal to the number of total loci
         or the number of markers.")
  }
  
  geno_split <- lapply(X = split_vec, FUN = function(ind) geno[,ind])
  
  # Return
  return(geno_split)
  
}



#' Split a haploid array into chromosomes
#' 
#' @param genome An object of class \code{genome}.
#' @param geno Array of haploid genotypes. Must be an array of dimensions \code{2} 
#' x \code{n.loci} x \code{n.ind}, the elements of which must be z {0, 1}.
#' 
#' @return 
#' A list of haploid arrays, split by chromosome.
#' 
#' @export
#' 
split_haploid <- function(genome, geno) {
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = geno))
    stop("The geno did not pass. See warning for reason.")
  
  # If the geno input is a list, recombine
  if (is.list(geno))
    geno <- do.call("cbind", geno)
  
  # Total number of loci
  n_loci <- nloci(genome = genome, by.chr = TRUE)
  # Total number of markers
  n_marker <- nmar(genome = genome, by.chr = TRUE)
  
  # Create a splitting vector based on the number of columns in the geno
  if (ncol(geno) == sum(n_loci)) {
    
    # Index of all loci
    loci_ind <- seq(sum(n_loci))
    # Split vector
    chr_split <- mapply(seq(nchr(genome)), n_loci, FUN = rep)
    
    split_vec <- split(x = loci_ind, f = unlist(chr_split))
    names(split_vec) <- chrnames(genome)
    
  } else if (ncol(geno) == sum(n_marker)) {
    
    # Index of all markers
    loci_ind <- seq(sum(n_marker))
    # Split vector
    chr_split <- mapply(seq(nchr(genome)), n_marker, FUN = rep)
    
    split_vec <- split(x = loci_ind, f = unlist(chr_split))
    names(split_vec) <- chrnames(genome)
    
  } else {
    stop("The number of columns in the geno input is not equal to the number of total loci
         or the number of markers.")
  }
  
  geno_split <- lapply(X = split_vec, FUN = function(ind) geno[,ind,])
    
  
  # Return
  return(geno_split)
  
} # Close the function






#' Pull the genotype data for a named locus
#' 
#' @param genome An object of class \code{genome}.
#' @param geno Genotype data on a population to phenotype. Must be a matrix of dimensions
#' \code{n.ind} x \code{n.loci}, the elements of which must be z {0, 1, 2}.
#' @param loci A character vector of loci to subset from the geno matrix.
#' 
#' @return
#' A matrix of genotypes from the geno matrix, only for the loci requested.
#' 
#' @export 
#'
pull_genotype <- function(genome, geno, loci) {
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = geno))
    stop("The geno did not pass. See warning for reason.")
  
  # If the geno input is a list, recombine
  if (is.list(geno))
    geno <- do.call("cbind", geno)
  
  # If haploid, convert to geno
  if (is_haploid(geno))
    geno <- haploid_to_geno(geno = geno)
  
  # Make sure the loci are in the marker names
  marker_names <- markernames(genome, include.qtl = TRUE)
  
  if (!all(loci %in% marker_names))
    stop("Not all of the loci names are in the genome.")
  
  # Subset the geno data for those loci
  subset(x = geno, select = loci, drop = FALSE)
  
} # Close the function


  
  