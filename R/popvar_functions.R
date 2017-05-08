#' Convert data for use in PopVar
#' 
#' @param genome An object of class \code{genome}.
#' @param geno Genotype data on a population to phenotype. Must be a matrix of dimensions
#' \code{n.ind} x \code{n.loci}, the elements of which must be z {0, 1, 2}.
#' 
#' @return 
#' A \code{data.frame} of marker genotypes ready for use in \code{\link[PopVar]{pop.predict}}.
#' 
#' @export
#' 
geno_to_popvar <- function(genome, geno) {
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = geno))
    stop("The geno did not pass. See warning for reason.")
  
  # If the geno input is a list, recombine
  if (is.list(geno))
    geno <- do.call("cbind", geno)
  
  # Get the individual names from the geno
  ind_names <- row.names(geno)
  # Marker names
  marker_names <- markernames(genome, include.qtl = FALSE)
  
  # Is the length of the marker names the same as the columns in 'geno'?
  # If not, subset the geno for markers, not QTL
  if (length(marker_names) != ncol(geno)) {
    geno_to_convert <- pull_genotype(genome = genome, geno = geno, loci = marker_names)
  } else {
    geno_to_convert <- geno
  }

  # Recode
  geno_recode <- geno_to_convert - 1
  
  # Create data.frame for output
  as.data.frame(cbind( c("", ind_names), rbind(marker_names, geno_recode)) )
  
} # Close the function


#' Convert data for use in PopVar
#' 
#' @describeIn geno_to_popvar 
#' 
#' @import dplyr
#' @importFrom qtl sim.cross
#' 
#' @export
#' 
map_to_popvar <- function(genome) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Get the map by first calling a blank_cross
  blank_cross <- sim.cross(map = genome$map, n.ind = 1)
  map <- pull.map(cross = blank_cross)
  
  # Convert and export
  map_recode <- lapply(X = map, FUN = function(chr)
    data.frame(marker = names(chr), 
               pos = as.numeric(chr)) )
  
  # Add chromosome names
  map_recode1 <- mapply(names(map_recode), map_recode, FUN = function(chr_name, chr_map)
    data.frame(chr_map, chr = chr_name), SIMPLIFY = FALSE) %>%
    do.call("rbind", .)
  
  # Rearrange and return
  select(map_recode1, marker, chr, pos)
  
}
