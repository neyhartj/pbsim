#' Create a population object
#' 
#' @description 
#' Assembles genotype data and into a \code{pop} object.
#' 
#' @param genome An object of class \code{genome}.
#' @param geno Genotype data on a population to phenotype. If the genome type is 
#' \code{"pbsim"}, must be a matrix of dimensions \code{n.ind} x \code{n.loci}, 
#' the elements of which must be z {0, 1, 2}, or a list of such matrices. If the 
#' genome type is \code{"hypred"}, must be an array of dimensions \code{2} x 
#' \code{n.loci} x \code{n.ind}, the elements of which must be z {0, 1}.
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
#' The genotypic value of individuals is calculcated as the sum of the QTL effects
#' carried by each individual. The genetic variance is calculated as the variance
#' of these genotypic values (\eqn{V_G = var(g)}).
#' 
#' @examples 
#' 
#' # Load historic data
#' data("s2_cap_genos")
#' data("s2_snp_info")
#' 
#' # Create a genome with genetic architecture
#' map <- lapply(split(s2_snp_info, s2_snp_info$chrom), function(chr) structure(chr$cM_pos, names = chr$rs) )
#' 
#' genome <- sim_genome(map = map)
#' 
#' # Simulate a a trait with 15 QTL
#' qtl.model <- matrix(nrow = 15, ncol = 4)
#' 
#' genome <- sim_gen_model(genome, qtl.model, add.dist = "geometric", max.qtl = 15)
#' 
#' pop <- create_pop(genome = genome, geno = s2_cap_genos)
#' 
#' 
#' ## Use haploid data and a 'hypred' genome
#' data("s2_cap_haploid")
#' 
#' genome <- sim_genome(map = map, type = "hypred")
#' 
#' genome <- sim_gen_model(genome, qtl.model, add.dist = "geometric", max.qtl = 15)
#' 
#' pop <- create_pop(genome = genome, geno = s2_cap_haploid)
#' 
#' 
#' @import dplyr
#' @importFrom abind abind
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
  
  # Genome type
  type <- attr(genome, "type")
  
  
  # Create empty pop list
  pop <- structure(vector("list", 2), class = "pop", names = c("geno", "geno_val"))
  
  # If the genome type is "hypred", recode the genos, but save the haploids
  if (type == "hypred") {
    
    # If the geno input is a list of arrays, recombine
    if (is.list(geno)) {
      
      # Calculate total dimensions
      n_row <- nrow(geno[[1]])
      n_col <- sum(sapply(geno, ncol))
      n_ind <- dim(geno[[1]])[3]
      
      # Column names
      loci_names <- unlist(sapply(geno, colnames, simplify = FALSE))
      
      geno <- abind::abind(geno, along = 2)
    }
    
    # Sort the haploids
    haploids <- geno[,,order(dimnames(geno)[[3]])]
    geno <- apply(X = haploids, MARGIN = 2, FUN = colSums)
    
    # Split the haploids
    haploid_split <- split_haploid(genome = genome, geno = haploids)
    
    # Add the haploids
    pop$haploids <- haploid_split
    
  } 
  
  # Sort the genos on individual names
  geno <- geno[order(row.names(geno)),]
  
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
