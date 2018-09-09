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
#' @param ignore.gen.model Logical - should any gene model be ignored when creating
#' the population? Use this to force a population without a gene model.
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
#' \dontrun{
#' 
#' # Use data from the PopVar package
#' library(PopVar)
#' data("think_barley")
#' 
#' # Format the map correctly and simulate a genome
#' map_in <- map.in_ex[,-1]
#' row.names(map_in) <- map.in_ex$mkr
#' genome <- sim_genome(map = table_to_map(map_in))
#' 
#' genos <- apply(X = G.in_ex[-1,-1], MARGIN = 2, FUN = as.numeric)
#' dimnames(genos) <- list(as.character(G.in_ex$V1[-1]), as.character(unlist(G.in_ex[1,-1])))
#' 
#' # Impute with the rounded mean
#' genos1 <- apply(X = genos, MARGIN = 2, FUN = function(snp) {
#'   mean <- ifelse(mean(snp, na.rm = T) < 0, -1, 1)
#'   snp[is.na(snp)] <- mean
#'   return(snp) })
#' 
#' ## Create a population without a genetic model
#' pop <- create_pop(genome = genome, geno = genos1 + 1, ignore.gen.model = T)
#' 
#' ## Create a genetic model with 15 QTL
#' qtl.model <- matrix(NA, ncol = 4, nrow = 15)
#' genome <- sim_gen_model(genome = genome, qtl.model = qtl.model, add.dist = "geometric")
#' 
#' pop <- create_pop(genome = genome, geno = genos1 + 1)
#' 
#' }
#' 
#' 
#' @importFrom abind abind
#' 
#' @export
#' 
create_pop <- function(genome, geno, ignore.gen.model = FALSE) {
  
  ## Error handling
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = geno, ignore.gen.model = ignore.gen.model))
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
  geno_split <- split_geno(genome = genome, geno = geno, ignore.gen.model = ignore.gen.model)
  
  # Calculate the genotypic value - ignore if told so
  if (ignore.gen.model) {
    geno_val <- NULL
  } else {
    geno_val <- calc_genoval(genome = genome, geno = geno_split)
  }
  
  
  # Add data to the pop
  pop[["geno"]] <- geno_split
  pop[["geno_val"]] <- geno_val
  
  # Return
  return(pop)
  
} # Close the function
