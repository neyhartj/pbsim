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
#' @param type The type of genome output. If \code{"pbsim"}, the genome will include a 
#' map that is compatible with \code{\link[qtl]{qtl-package}}, and if \code{"hypred"},
#' the genome will be a list of genomes compatible with \code{\link[hypred]{hypred}}.
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
#' # Use a map instead
#' data("s2_snp_info")
#' map <- lapply(split(s2_snp_info, s2_snp_info$chrom), function(chr) structure(chr$cM_pos, names = chr$rs) )
#' 
#' genome <- sim_genome(map = map)
#' 
#' # Use 'hypred'
#' genome <- sim_genome(map = map, type = "hypred")
#' 
#' @importFrom qtl sim.map
#' @import hypred
#' 
#' @export 
#' 
#' 
sim_genome <- function(len, n.mar, map, eq.spacing = FALSE, type = c("pbsim", "hypred")) {
  
  # Error handling
  # len and n.mar must be present if map is not
  if (missing(map)) {
    len <- as.numeric(len)
    n.mar <- as.numeric(n.mar)
    
  } else {
    
    if (missing(len))
      len <- sapply(map, max)
    
    if (missing(n.mar))
      n.mar <- sapply(map, length)

  }
  
  # Match the type argument
  type <- match.arg(type)
  
  # The length of the len vector must be the same as that of n.mar
  if (length(len) != length(n.mar)) 
    stop("The length of the input 'len' vector does not equal the length of
         'n.mar.'")

  # Create an empty list of length n.chr
  genome <- structure(vector("list"), class = "genome", type = type)
  
  # If genetic map is NULL, sample from uniform distribution
  if (missing(map)) {
    
    # Use R/qtl
    map <- qtl::sim.map(len = len, n.mar = n.mar, include.x = FALSE, eq.spacing = eq.spacing)
    
    # If map is not null, check compatability
  } else {
    
    # First check that the markers are named
    if (any(sapply(map, function(m) is.null(names(m)))))
      stop("Marker positions in the input 'map' must be named.")
    
    # Make sure the map has the chromosome class attribute
    map <- lapply(map, structure, class = "A")
    # Reclass the map
    class(map) <- "map"
    
  }


  # What type of genome was requested?
  if (type == "pbsim") {
    
    # Assemble the genome
    genome[["len"]] <- len
    genome[["n.mar"]] <- n.mar
    genome[["map"]] <- map
    
  } else {
    
    # Remove the map classes
    map <- sapply(X = map, FUN = structure, class = NULL, simplify = FALSE)
    
    # Apply a function over the number of chromosomes
    hypred_list <- lapply(X = map, FUN = function(chr) {
      # Create a new base genome
      hp_genome <- hypredGenome(num.chr = 1, len.chr = max(chr) / 100, num.snp.chr = length(chr))
      # Add genetic map
      hypredNewMap(hp_genome, new.map = chr / 100) })
  
    # Add this information to the genome
    genome[["hypredGenomes"]] <- hypred_list
    
  }
  
  
  # Return the genome
  return(genome)
  
} # Close the function