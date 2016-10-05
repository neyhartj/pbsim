#' Define the genetic architecture of a trait
#' 
#' @description 
#' Creates a list of genetic characteristics, including positon of QTL, additive
#' effect, dominance effect, and positon of perfect loci (i.e. SNP markers
#' that are also QTL).
#' 
#' @param genome The list of hypred genomes.
#' @param n.QTL The number of SNPs to become QTL.
#' @param qtl.index A list of SNP indices per chromosome to become QTL. If NULL, 
#' indices are randomly sampled.
#' @param qtl.dom.index A list of QTL indices per chromosome to have dominance 
#' effects.
#' @param qtl.perf.index A list of QTL indices per chromosome to be perfect
#' markers (i.e. the observed SNP is the QTL).
#' @param qtl.add.eff A vector of length \code{n.QTL} of additive effects of QTL.
#' The vector is interpreted as such: the ith element in the vector will be the
#' additive effect assigned to the ith QTL. Can also be \code{"normal"} to draw 
#' effects from a standard normal distribution or \code{"geometric"} to draw
#' effects from a geometric series.
#' @param qtl.dom.eff A vector of length \code{n.QTL} of dominance effects of QTL.
#' The vector is interpreted as such: the ith element in the vector will be the
#' dominance effect assigned to the ith QTL.
#' 
#' @return 
#' A list of hypred genomes with assigned trait architecture.
#' 
#' @import hypred
#' 
#' @export
#' 
genetic.architecture <- function(map, n.QTL, qtl.index = NULL, 
                                 qtl.dom.index = NULL, qtl.perf.index = NULL, 
                                 qtl.add.eff = "normal", qtl.dom.eff = NULL) {
  
  # Find the number of chromosomes from the length of the map list
  n.chr <- length(map)
  
  # If the qtl.ids are NULL, randomly sample the loci
  if (is.null(qtl.index)) {

    # First create a vector of chromosomes for each QTL
    qtl.per.chr <- sample(n.chr, size = n.QTL, replace = T) %>%
      table()
    
    # Iterate over chromosomes
    qtl.index.per.chr <- mapply(qtl.per.chr, map, FUN = function(n.qtl.p, chr) {
      
      # Pull out the number of loci on the chromsome and create an index
      loci.index.p <- length(chr) %>% seq()
      
      # Sample the index and add to the list
      sort(sample(x = loci.index.p, size = n.qtl.p)) })
    
    # If the list is specified, check it
  } else {
    
    if (length(qtl.index) != n.chr) stop("The length of the qtl.index list is not the same as the number of chromosomes.")
    # Rename it
    qtl.index.per.chr <- qtl.index
    # Figure out the QTL per chromosomes
    qtl.per.chr <- sapply(X = qtl.index.per.chr, FUN = length)
  }
  
  # QTL dominance IDs
  if (is.null(qtl.dom.index)) {
    # Create a list of NULLs
    qtl.dom.index.per.chr <- vector("list", n.chr)
    
  } else { # Otherwise check the list
    if (length(qtl.dom.index) != n.chr) stop("The length of the qtl.dom.index list is not the same as the number of chromsomes.")
    
    # Rename it
    qtl.dom.index.per.chr <- qtl.dom.index
  }
  
  # Perfect QTL IDs
  if (is.null(qtl.perf.index)) {
    # Create a list of NULLs
    qtl.perf.index.per.chr <- vector("list", n.chr)
    
  } else { # Otherwise check the list
    if (length(qtl.perf.index) != n.chr) stop("The length of the qtl.perf.index list is not the same as the number of chromosomes.")
    
    # Rename it
    qtl.perf.index.per.chr <- qtl.perf.index
  }
  
  # Assign qtl effects
  if (is.character(qtl.add.eff)) {
    
    # Draw from standard normal
    if (qtl.add.eff == "normal") {
      qtl.add.eff <- abs(rnorm(n = n.QTL, mean = 0, sd = 1))
      # Randomly assign negative values to the qtl effects - this corresponds to the value of the 1 allele
      qtl.add.eff <- qtl.add.eff * sample(c(1,-1), n.QTL, replace = T)
      
    } else { # Draw from geometric series
      
      if (qtl.add.eff == "geometric") {
        a = (n.QTL - 1) / (n.QTL + 1)
        qtl.add.eff <- sample(a^(1:n.QTL))
        # Randomly assign negative values to the qtl effects - this corresponds to the value of the 1 allele
        qtl.add.eff <- qtl.add.eff * sample(c(1,-1), n.QTL, replace = T)
        
        # Break up the effects into chromosomes with the same number of elements
        ## as QTL on those chromosomes
        qtl.add.eff.per.chr <- split(qtl.add.eff, rep(seq(n.chr), qtl.per.chr))
        
      } else {
        # Otherwise, make sure the length of the additive effects is the same as the number of chr
        if (length(qtl.add.eff) != n.chr) stop("The length of the QTL effects is not the same as the number of chromosomes")
      }}}
  
  # Dominance effects
  if (is.null(qtl.dom.eff)) {
    # Create a list of NULLs
    qtl.dom.eff.per.chr <- vector("list", n.chr)
    
  } else { # Otherwise check it
    if (length(qtl.dom.eff) != n.chr) stop("The length of the qtl.dom.eff list is not the same as the number of chromosomes.")
    
    # Make sure each element in the qtl.dom.eff list is the same as the qtl.dom.index.per.chr list
    dom.elements.same <- sapply(X = 1:n.chr, FUN = function(i) length(qtl.dom.eff[[i]]) == length(qtl.dom.index.per.chr[[i]]) )
    
    if (!all(dom.elements.same)) stop("One or more of the elements in the qtl.dom.eff list are not the same length as the corresponding qtl.dom.index list.")
  }
  
  # Iterate over the chromosomes and qtl indices, effects, etc
  genome1 <- mapply(genome, qtl.index.per.chr, qtl.dom.index.per.chr, 
                    qtl.perf.index.per.chr, qtl.add.eff.per.chr, 
                    qtl.dom.eff.per.chr, FUN = function(chr, a.i, d.i, p.i, a, d) {
                      
                      chr <- hypredNewQTL(chr,
                                          new.id.add = a.i,
                                          new.id.dom = d.i,
                                          new.id.per.mar = p.i,
                                          new.eff.add = a,
                                          new.eff.dom = d ) })
  
  # Return the new genome
  return(genome1)
  
} # Close the function
