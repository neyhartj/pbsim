#' Generate doubled-hapoids from crossover data
#' 
#' @param xodat Crossover data generated from a pedigree. See \code{\link[simcross]{sim_from_pedigree}}.
#' @param pedigree A \code{pedigree} detailing the scheme to develop the family.
#' Use \code{\link{sim_pedigree}} to generate.
#' 
#' @return 
#' A new xodat object to generate doubled-haploids. Only information on the final
#' individuals in the pedigree is returned.
#' 
#' @importFrom simcross check_pedigree
#' 
#' @export
#' 
induce_dh <- function(xodat, pedigree) {
  
  # Check the pedigree
  if (!check_pedigree(pedigree, ignore_sex = TRUE))
    stop("The pedigree is not formatted correctly.")
  
  # Subset the pedigree for the finals
  final_id <- subset(pedigree, gen == max(gen))$id
  
  # Copy the xo_dat
  xodat_copy <- xodat
  
  # Iterate over chromosomes in the xodat
  xodat_dh <- lapply(X = xodat, FUN = function(chr) {
    # Iterate over the finals
    for (ind in final_id) {
      # Sample which of mat or pat
      which_par <- sample(length(chr[[ind]]), size = 1)
      chr[[ind]]$mat <- chr[[ind]]$pat <- chr[[ind]][[which_par]]
    }
    return(chr) })
  
  # Return
  return(xodat_dh)
  
  
} # Close the function
  
  