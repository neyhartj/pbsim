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
  
  # Number of rows in the genotypic value matrix
  row.names(pop$geno_val)
  
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
#' @export
#' 
subset.pop <- function(pop, individual) {
  
  # Error handling
  # Make sure pop inherits the class "pop"
  if (!inherits(pop, "pop"))
    stop("The input 'pop' must be of class 'pop'.")
  
  # Make sure the individuals specified are in the the pop
  if (!all(individual) %in% indnames(pop))
    stop("Not all of the individuals in 'individual' are in the 'pop' object.")
  
  # Empty pop object
  new_pop <- structure(vector("list", length(pop)), class = "pop", names = names(pop))
  
  # Subset various components
  new_pop$geno <- lapply(X = pop$geno, FUN = "[", individual,, drop = FALSE)
  new_pop$geno_val <- pop$geno_val[individual, , drop = FALSE]

  if (!is.null(pop$pheno_val)) {
    # Get rid of the variance component estimate
    new_pop$pheno_val <- pop$pheno_val[-1]
    
    # Subset the phenotypic observations and pheno_mean
    new_pop$pheno_val <- lapply(X = new_pop$pheno_val, subset, subset = ind %in% individual, drop = FALSE)
    
  }
  
  # Return the population
  return(new_pop)
  
} # Close the function
    
    
  
  



