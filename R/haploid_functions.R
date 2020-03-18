#' Haploid genotype data
#' 
#' @param geno Population genotype data.
#' 
#' @details 
#' 
#' Haploid genotype data is stored as an array of dimensions \code{2} x 
#' \code{n.loci} x \code{n.ind}, the elements of which must be z {0, 1}.
#' 
#' Haploids can be checked by using the \code{\link{is_haploid}} function, and
#' can be converted to a standard genotype matrix by using the 
#' \code{\link{haploid_to_geno}} function.
#' 
#' 
#' @export
#' 
is_haploid <- function(geno) {
  
  # Does the object have 3 dimensions
  if (length(dim(geno)) != 3)
    return(FALSE)
  
  # Is the object encoded in z {0, 1}
  if (!all(geno %in% c(0, 1))) 
    return(FALSE)
  
  # If passing, return TRUE
  return(TRUE)
  
  
} # Close the function


#' @rdname is_haploid
#' @export
#' 
haploid_to_geno <- function(geno) {
  
  # Check if haploid
  stopifnot(is_haploid(geno))
  
  # colSums
  apply(X = geno, MARGIN = 2, FUN = colSums)
  
} # Close the function
