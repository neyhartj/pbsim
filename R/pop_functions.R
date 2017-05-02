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



