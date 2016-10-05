#' Create a genome
#' 
#' @description Create an object of class "genome" with initial parameters of chromosome size, number of loci on each chromsome, and the position of those loci.
#' 
#' @param n.chr The integer number of chromosomes in the genome.
#' @param chr.len Vector specifying the Morgan length of each chromosome. Vector length must equal \code{chr.num}.
#' @param chr.snps Vector specifying the number of SNP loci per chromosome. Vector length must equal \code{chr.num}.
#' @param genetic.map List of Morgan positions of loci across the genome. List length must equal \code{chr.num} and list elements must be objects of class \code{numeric} with lengths equal to the elements in \code{chr.len}. If \code{NULL}, positions are sampled from a uniform distribution.
#' 
#' @return Object of class "genome" with the number of chromosomes, the number of snps,
#' the Morgan length of the whole genome, the number of QTL, and a list of \code{chromosome} objects.
#' 
#' @details 
#' Defining the genome is the first step in performing simulations using the \code{quantgen} package. Genomes are defined as
#' objects of the same name. See \code{\link[quantgen]{genome}} for information about the components of a \code{genome}
#' object and see \code{\link[quantgen]{chromosome}} for information about the components of a \code{chromosome} object.
#' 
#' @examples
#' n.chr = 10
#' chr.len = seq(1.0, 1.45, by = 0.05)
#' chr.snps = sample(15:20, size = 10, replace = T)
#' genome <- make.genome(n.chr = n.chr, chr.len = chr.len, chr.snps = chr.snps)
#' 
#' # Using the same parameters above, but defining a genetic map manually
#' genetic.map <- sapply(X = 1:n.chr, FUN = function(i) runif(n = chr.snps[i], min = 0, max = chr.len[i]))
#' genome <- make.genome(n.chr = n.chr, chr.len = chr.len, chr.snps = chr.snps, genetic.map = genetic.map)
#' 
#' @export 
#' 
#' 
make.genome <- function(n.chr, chr.len, chr.snps, genetic.map = NULL) {
  
  # Error handling
  n.chr <- as.integer(n.chr)
  chr.len <- as.numeric(chr.len)
  chr.snps <- as.numeric(chr.snps)
  
  # If the length of the chr.len vector is not the same as n.chr
  if (length(chr.len) != n.chr) stop("The length of the chr.len vector does not equal n.chr.")
  if (length(chr.snps) != n.chr) stop("The length of the chr.snps vector does not equal n.chr.")
  
  # If genetic map is NULL, sample from uniform distribution
  if (is.null(genetic.map)) {
    
    # Apply a function over the number of chromosomes
    chromosomes <- lapply(X = 1:n.chr, FUN = function(i) {
      
      # Pull out the chromosome length and number of loci
      chr.len.i <- chr.len[i]
      n.chr.snps.i <- chr.snps[i]
      # Sample a normal distribution
      snp.pos.i <- sort(runif(n = n.chr.snps.i, min = 0, max = chr.len.i))
      
      # Create list using positions
      snp.pos.i.list <- list(index = seq(n.chr.snps.i), M = snp.pos.i)
      
      # Create a new chromosome
      new(Class = "chromosome",
          chr.len = chr.len.i,
          n.snps = as.integer(n.chr.snps.i),
          pos.snps = snp.pos.i,
          n.markers = as.integer(n.chr.snps.i),
          pos.markers = snp.pos.i.list,
          n.add.qtl = as.integer(0),
          n.dom.qtl = as.integer(0),
          pos.add.qtl = list(index = as.numeric(0), M = as.numeric(0)),
          pos.dom.qtl = list(index = as.numeric(0), M = as.numeric(0)),
          qtl.effects = list(a = as.numeric(0),
                             d = as.numeric(0)) )
    })
    
  } else { # Otherwise use the genetic map provided
    genetic.map <- as.list(genetic.map)
    
    # Error handling
    if (length(genetic.map) != n.chr) stop("The length of genetic.map is not equal to chr.num")
    genetic.map.loci <- as.numeric(unlist(lapply(genetic.map, FUN = length)))
    if (!identical(genetic.map.loci, chr.snps)) stop("The length of the elements in the genetic.map list do not equal the elements in the chr.snps vector.")
    
    # Make sure that the range of loci positions for each chromosome is within the 
    ## length of the chromosome
    map.pos.interval <- sapply(X = 1:n.chr, FUN = function(i) {
      chr.len.i <- chr.len[i]
      map.range.i <- range(genetic.map[i])
      # Determine if the output of findInterval is c(1,1), indicating that
      ## the map positions are within the range.
      all(findInterval(vec = c(0, chr.len.i), x = map.range.i) == 1) })
    # Error
    if (!all(map.pos.interval)) stop("The genetic map positions for one or more of the chromosomes are beyond the range designated to the chromosome.")
    
    # Apply a function over the number of chromosomes
    chromosomes <- lapply(X = 1:n.chr, FUN = function(i) {
      
      # Pull out the chromosome length and number of loci
      chr.len.i <- chr.len[i]
      n.chr.snps.i <- chr.snps[i]
      # Sample a normal distribution
      snp.pos.i <- genetic.map[[i]]
      
      # Create list using positions
      snp.pos.i.list <- list(index = seq(n.chr.snps.i), M = snp.pos.i)
      
      # Create a new chromosome
      new(Class = "chromosome",
          chr.len = chr.len.i,
          n.snps = as.integer(n.chr.snps.i),
          pos.snps = snp.pos.i,
          n.markers = as.integer(n.chr.snps.i),
          pos.markers = snp.pos.i.list,
          n.add.qtl = as.integer(0),
          n.dom.qtl = as.integer(0),
          pos.add.qtl = list(index = as.numeric(0), M = as.numeric(0)),
          pos.dom.qtl = list(index = as.numeric(0), M = as.numeric(0)),
          qtl.effects = list(a = as.numeric(0),
                             d = as.numeric(0)) )
    })
    
  }
  
  # Find the number of snps
  n.snps <- as.integer(sum(chr.snps))
  # Find the total length of the genome
  genome.len <- sum(chr.len)
  
  # Return the genome
  genome <- new(Class = "genome", 
                n.chr = n.chr,
                n.snps = n.snps,
                genome.len = genome.len,
                n.qtl = as.integer(0),
                n.perfect.snps = as.integer(0),
                chromosomes = chromosomes)
                
  
  return(genome)
}