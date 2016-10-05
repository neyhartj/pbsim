#' Create a crossing block
#' 
#' @description Creates a crossing block of parent names for use in creating a family or population.
#' 
#' @param genome.object An object of class \code{genome} with declared positions of
#' SNP markers and QTL.
#' @param hap.genome1 A numeric vector of a parent's 1st haploid genome that is 
#' the same length as the number of SNP loci in the genome. The elements of the 
#' vector should be 0 or 1.
#' @param hap.genome2 A numeric vector of a parent's 1st haploid genome that is 
#' the same length as the number of SNP loci in the genome. The elements of the 
#' vector should be 0 or 1.
#' @param mutate A logical switch whether or not to mutate haploids before recombination.
#' @param mu.snp A numeric vector giving the probability of a mutation at any given
#' SNP marker per meiotic event. Mutations are assumed completely independent.
#' @param mu.qtl A numeric vector giving the probability of a mutation at any given
#' QTL per meiotic event. Mutations are assumed completely independent.
#' 
#' @details TBD
#' 
#' @examples 
#' n.chr = 3
#' chr.len = c(1.0, 1.1, 1.2)
#' chr.snps = c(15, 16, 17)
#' genome <- make.genome(n.chr = n.chr, chr.len = chr.len, chr.snps = chr.snps)
#' genome <- trait.architecture(genome.object = genome, n.qtl = 10)
#' 
#' hap1 <- sample(c(1,0), sum(chr.snps), replace = T)
#' hap2 <- sample(c(1,0), sum(chr.snps), replace = T)
#' 
#' recombine(genome, hap1, hap2, mutate = T, mu.snp = 7e-8, mu.qtl = 7e-8)
#' 
#' @return A \code{numeric} vector representing the haploid genome of a gamete that
#' is the product of recombination between the two haploid genomes of a parent and 
#' following any mutation events.
#' 
#' 
#' @export
#' 
make.crossing.block <- function(parent1.lines, # Character vector of lines for the first parent
                                parent2.lines, # Character vector of lines for the second parent
                                n.crosses, # Number of crosses to make
                                method = "random", # Method to assign parents. Can be "random" for random pairs of the parent1 and parent2, "chain" for sequential crosses (Note: can only be used if the parent1.line and parent2.lines vectors are identical), or "all.pairwise" for all possible pairwise crosses.
                                use.parents.once = FALSE) {
  
  # Find the intersection of parent1.lines and parent2.lines
  par.intersect <- intersect(parent1.lines, parent2.lines)
  n.intersect <- length(par.intersect)
  # Find the number of parent1 and parent2 lines
  n.par1 <- length(parent1.lines)
  n.par2 <- length(parent2.lines)
  
  if ( all(n.intersect == n.par1, n.intersect == n.par2) ) {
    n.possible.crosses <- (n.par1 * (n.par2 - 1)) / 2
    same.lines = T
  } else {
    n.possible.crosses <- n.par1 * n.par2
    same.lines = F
  }
  
  # Create all pairwise crosses
  sample.crosses <- expand.grid(parent1.lines, parent2.lines)[,c(2,1)]
  # Remove selfs
  sample.crosses <- subset(x = sample.crosses, !apply(X = sample.crosses, MARGIN = 1, FUN = function(cross) any(duplicated(cross))) )
  
  # If statements for methods
  if (method == "all.pairwise") {
    return(sample.crosses)
  }
  if (method == "random") {
    # First see if the number of requested crosses is more than possible
    if (n.crosses > n.possible.crosses) stop("n.crosses is more than possible.")
    
    # If parents should only be used once, sample without replacement all lines and put them into a matrix
    if (use.parents.once) {
      # Sample into a matrix
      if (same.lines) {
        random.crosses <- as.data.frame(matrix(sample(parent1.lines), ncol = 2, byrow = T))
        crosses.ind <- sort(sample(1:nrow(random.crosses), n.crosses))
        random.crosses <- random.crosses[crosses.ind,]
      } else {
        # Are there any overlapping lines in the two parent vectors?
        if (n.intersect == 0) {
          # If so just sample each vector into a matrix
          random.crosses <- as.data.frame(cbind( sample(parent1.lines), sample(parent2.lines) ))
          crosses.ind <- sort(sample(1:nrow(random.crosses), n.crosses))
          random.crosses <- random.crosses[crosses.ind,]
        } else {
          stop("Function not available.")
        }
      }
      
    } else {
      
      # Remove reciprocal
      sample.crosses <- subset(x = sample.crosses, subset = !duplicated(t(apply(sample.crosses, MARGIN = 1, FUN = sort))))
      crosses.ind <- sort(sample(1:nrow(sample.crosses), n.crosses))
      random.crosses <- sample.crosses[crosses.ind,]
    }
    
    colnames(random.crosses) <- c("Parent1", "Parent2")
    row.names(random.crosses) <- NULL
    return(random.crosses)
  }
  
} # Close the function