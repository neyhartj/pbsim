#' Create a population object
#' 
#' @description 
#' Assembles genotype data and into a \code{pop} object.
#' 
#' @param genome An object of class \code{genome}.
#' @param geno Genotype data on a population to phenotype. Can be a matrix of dimensions
#' \code{n.ind} x \code{n.loci}, the elements of which must be z {0, 1, 2}, or a list
#' of such matrices.
#' 
#' 
#' @details 
#' A \code{pop} is similar to a \code{cross} object in \code{\link[qtl]{qtl-package}} 
#' (see \code{\link[qtl]{read.cross}}). The \code{pop} object stores information on the
#' genome, the genotypes at genetic markers, and phenotypes. The \code{pop} object is
#' meant to be a bit more flexible, without the pedigree or family structure required in
#' a \code{cross} object.
#' 
#' @return 
#' An object of class \code{pop} with genotype information for the individuals in
#' that population and the genotypic value of those individuals.
#' 
#' @examples 
#' 
#' # Load some historic data
#' data("s2_genos")
#' data("s2_snp_info")
#' 
#' geno <- s2_genos + 1
#' 
#' # Create a genome with genetic architecture
#' len <- tapply(s2_snp_info$cM_pos, s2_snp_info$chrom, max)
#' n_mar <- tapply(s2_snp_info$cM_pos, s2_snp_info$chrom, length)
#' map <- lapply(split(s2_snp_info, s2_snp_info$chrom), function(chr) structure(chr$cM_pos, names = chr$rs) )
#' 
#' genome <- sim_genome(len = len, n.mar = n_mar, map = map)
#' 
#' # Simulate a a trait with 15 QTL
#' qtl.model <- matrix(nrow = 15, ncol = 4)
#' 
#' genome <- sim_gen_model(genome, qtl.model, add.dist = "geometric")
#' 
#' # Add QTL to the geno matrix
#' new_geno <- fill_qtl_geno(genome = genome, geno = geno)
#' 
#' pop <- create_pop(genome = genome, geno = new_geno)
#' 
#' @import dplyr
#' 
#' @export
#' 
create_pop <- function(genome, geno) {
  
  ## Error handling
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = geno))
    stop("The geno did not pass. See warning for reason.")
  
  # Create empty pop list
  pop <- structure(vector("list"), class = "pop")
  
  # Split the geno matrix into chromosomes
  geno_split <- split_geno(genome = genome, geno = geno)
  
  # Calculate the genotypic value
  geno_val <- calc_genoval(genome = genome, geno = geno_split)
  
  
  # Add data to the pop
  pop[["geno"]] <- geno_split
  pop[["geno_val"]] <- geno_val
  
  # Return
  return(pop)
  
} # Close the function
