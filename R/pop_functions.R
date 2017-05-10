#' Number of individuals in a population
#' 
#' @param pop An object of class \code{pop}.
#' 
#' @return 
#' Scalar number of individuals.
#' 
#' @export
#' 
nind <- function(pop) {
  
  # Make sure pop inherits the class "pop"
  if (!inherits(pop, "pop"))
    stop("The input 'pop' must be of class 'pop'.")
  
  # Does the pop object have genotypic values
  if (is.null(pop$geno_val))
    stop("The 'pop' object must have the data.frame of genotypic values")
  
  # Number of rows in the genotypic value matrix
  nrow(pop$geno_val)
  
} # Close function


#' Names of individuals
#' 
#' @param pop An object of class \code{pop}.
#' 
#' @return 
#' Character vector of individual names.
#' 
#' @export
#' 
indnames <- function(pop) {
  
  # Make sure pop inherits the class "pop"
  if (!inherits(pop, "pop"))
    stop("The input 'pop' must be of class 'pop'.")
  
  # Does the pop object have genotypic values
  if (is.null(pop$geno_val))
    stop("The 'pop' object must have the data.frame of genotypic values")
  
  # Return the individual names from the genotypic value df
  pop$geno_val$ind
  
} # Close function




#' Genotype a population
#' 
#' @param genome An object of class \code{genome}.
#' @param pop An object of class \code{pop}.
#' @param error.rate The genotyping error rate. This argument is not yet operational.
#' 
#' @return 
#' A matrix of dimensions \code{nind} x \code{nmarkers} in format z {-1, 0, 1}.
#' 
#' @export
#' 
genotype <- function(genome, pop, error.rate = 0) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Make sure there is a genetic model
  if (is.null(genome$gen_model))
    stop("No genetic model has been declared for the genome")
  
  # Make sure pop inherits the class "pop"
  if (!inherits(pop, "pop"))
    stop("The input 'pop' must be of class 'pop'.")
  
  # Get the names of the markers
  marker_names <- markernames(genome)
  
  # Combine the genos per chromosome
  geno <- do.call("cbind", pop$geno)
  
  # Subset only the markers and subtract 1
  subset(x = geno, select = marker_names, drop = FALSE) - 1
  
} # Close the function 



#' Subset a population object for specific individuals
#' 
#' @param pop An object of class \code{pop}.
#' @param individual A character of individuals to subset from the \code{pop}.
#' 
#' @details 
#' If \code{pheno_val} is present in the \code{pop}, the variance components are
#' dropped.
#' 
#' @examples 
#' 
#' # Load some historic data
#' data("s2_cap_genos")
#' data("s2_snp_info")
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
#' genome <- sim_gen_model(genome, qtl.model, add.dist = "geometric", max.qtl = 15)
#' 
#' pop <- create_pop(genome = genome, geno = s2_cap_genos)
#' 
#' individual <- c("2ND27380", "06WA-406.6")
#' 
#' @import dplyr
#' 
#' @export
#' 
subset_pop <- function(pop, individual) {
  
  # Error handling
  # Make sure pop inherits the class "pop"
  if (!inherits(pop, "pop"))
    stop("The input 'pop' must be of class 'pop'.")

  # Convert individual to character
  individual <- as.character(individual)
    
  # Make sure the individuals specified are in the the pop
  if (!all(individual %in% indnames(pop)))
    stop("Not all of the individuals in 'individual' are in the 'pop' object.")

  
  # Empty pop object
  new_pop <- structure(vector("list", length(pop)), class = "pop", names = names(pop))
  
  # Subset various components
  new_pop$geno <- lapply(X = pop$geno, FUN = "[", individual, , drop = FALSE)
  new_pop$geno_val <- filter(pop$geno_val, ind %in% individual)

  if (!is.null(pop$pheno_val)) {
    # Get rid of the variance component estimate
    new_pop$pheno_val <- pop$pheno_val[-1]
    
    # Subset the phenotypic observations and pheno_mean
    new_pop$pheno_val <- lapply(X = new_pop$pheno_val, filter, ind %in% individual)
    
  }
  
  # Return the population
  return(new_pop)
  
} # Close the function



#' Combine a list of populations
#' 
#' @param pop_list A list of objects of class \code{pop}.
#' 
#' @details 
#' If \code{pheno_val} is present in the \code{pop}, the variance components are
#' dropped.
#' 
#' @importFrom purrr pmap
#' 
#' @export
#' 
combine_pop <- function(pop_list) {
  
  # Make sure each element of 'pop_list' is a pop
  if (!all(sapply(X = pop_list, FUN = inherits, "pop")))
    stop("One of more of the elements in 'pop_list' is not a 'pop' object.")
  
  # Create a new pop object
  new_pop <- structure(vector("list", length(pop_list[[1]])), class = "pop", 
                       names = names(pop_list[[1]]))
  
  # Combine genotypes
  # First extract the 'geno' element from each pop
  geno_list <- lapply(pop_list, "[[", "geno")
  # Combine
  new_pop$geno <- purrr::pmap(geno_list, rbind)
  
  # Combine genotypic values
  new_pop$geno_val <- do.call("rbind", lapply(pop_list, "[[", "geno_val"))

  # Combine phenotypic values if present
  if (!is.null(pop_list[[1]]$pheno_val)) {
    # Get rid of the variance component estimate
    new_pop$pheno_val <- structure(vector("list", 2), names = names(pop_list[[1]]$pheno_val)[-1])
    
    # Subset the 'pheno_val' element
    pheno_list <- lapply(X = pop_list, FUN = "[[", "pheno_val")
    
    # Combine
    new_pop$pheno_val$pheno_obs <- do.call("rbind", lapply(X = pheno_list, FUN = "[[", "pheno_obs"))
    new_pop$pheno_val$pheno_mean <- do.call("rbind", lapply(X = pheno_list, FUN = "[[", "pheno_mean"))
    
  }
  
  # Return the new pop
  return(new_pop)
  
} # Close the function
