#' Create a simulated genome
#' 
#' @description
#' Creates a list containing information on a simulated genome
#' 
#' @param len A vector specifying the chromosome lengths (in cM).
#' @param n.mar A vector specifying the umber of markers per chromosome.
#' @param map A list of marker positions (in cM), where each list element is a
#' named vector of marker positions. If \code{NULL} (default), marker positions
#' are drawn from a uniform distribution. See \code{\link[qtl]{sim.map}} for more
#' information.
#' @param eq.spacing If TRUE, markers will be equally spaced. See \code{\link[qtl]{sim.map}}.
#' 
#' @return 
#' Object of class "genome" with the length of each chromosome, the number of
#' markers per chromosome, and the genetic map.
#' 
#' @examples
#' n.mar  <- c(505, 505, 505)
#' len <- c(120, 130, 140)
#' 
#' genome <- sim_genome(len, n.mar)
#' 
#' @importFrom qtl sim.map
#' 
#' @export 
#' 
#' 
sim_genome <- function(len, n.mar, map = NULL, eq.spacing = FALSE) {
  
  # Error handling
  len <- as.numeric(len)
  n.mar <- as.numeric(n.mar)
  
  # The length of the len vector must be the same as that of n.mar
  if (length(len) != length(n.mar)) 
    stop("The length of the input 'len' vector does not equal the length of
         'n.mar.'")
  
  # If map is not NULL, its length must be the same as len and n.mar
  if (!is.null(map)) {
    if (length(len) != length(map) | length(n.mar) != length(map))
      stop("The length of input 'map' must be the same as the length
           of inputs 'n.mar' and 'len.'")
    }

  # Create an empty list of length n.chr
  genome <- list()
  class(genome) <- "genome"
  
  # If genetic map is NULL, sample from uniform distribution
  if (is.null(map)) {
    
    # Use R/qtl
    map <- sim.map(len = len, n.mar = n.mar, include.x = FALSE, eq.spacing = eq.spacing)
    
    # If map is not null, check compatability
  } else {
    
    # First check that the markers are named
    if (any(sapply(map, function(m) is.null(names(m)))))
      stop("Marker positions in the input 'map' must be named.")
    
    # Next, check that the length of each chromosome in the map is within the
    # chromosome length set by 'len'
    if (!all(mapply(len, map, FUN = function(l, m) max(m) <= l)))
      stop("The maximum cM position in 'map' is not less than or equal to
           the chromosome length in 'len.'")
    
    # Now check that the number of markers is equal to n.mar
    if (!all(mapply(n.mar, map, FUN = function(n, m) length(m) == n)))
      stop("The number of loci in 'map' is not equal to the number of 
           markers specified in 'n.mar.'")
    
    # Make sure the map has the chromosome class attribute
    map <- lapply(map, structure, class = "A")
    # Reclass the map
    class(map) <- "map"
    
  }
  
  # Assemble the genome
  genome[["len"]] <- len
  genome[["n.mar"]] <- n.mar
  genome[["map"]] <- map
  
  # Return the genome
  return(genome)

}