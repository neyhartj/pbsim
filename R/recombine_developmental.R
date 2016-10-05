#' Simulate a recombined gamete from a two haploid parental genomes
#' 
#' @description Creates a single haploid gamete from two haploid parental genomes. Recombination is modeled on the Haldane mapping function, which assumes no cross-over interference. Number of cross-over events per chromosome is samples from a Poisson distribution with lambda = L where L is the Morgan length of the chromosome. Mutation rates on a per-SNP and per-QTL basis may be included.
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
#' 
recombine <- function(genome.object, hap.genome1, hap.genome2, mutate = TRUE, mu.snp, mu.qtl) {
  
  ### Error handling and validation ###
  hap.genome1 <- as.vector(hap.genome1)
  hap.genome2 <- as.vector(hap.genome2)
  mu.snp <- as.numeric(mu.snp)
  mu.qtl <- as.numeric(mu.qtl)
  
  # If the haploids are not the right length, error out
  if (length(hap.genome1) != genome.object@n.snps) stop("The length of hap.genome1 is not the same as the number of SNPs")
  if (length(hap.genome2) != genome.object@n.snps) stop("The length of hap.genome1 is not the same as the number of SNPs")
  
  
  # Pull out the number of chromosomes
  n.chr <- genome.object@n.chr
  
  # Pull out the number of snps per chromosome
  chr.snps <- sapply(X = genome.object@chromosomes, FUN = function(chromosome) chromosome@n.snps)
  # Pull out the number of markers per chromosome
  chr.markers <- sapply(X = genome.object@chromosomes, FUN = function(chromosome) chromosome@n.markers)
  # Pull out the number of qtl per chromosome
  chr.qtl <- sapply(X = genome.object@chromosomes, FUN = function(chromosome) chromosome@n.add.qtl)
  # Pull out the length of each chromosome
  chr.len <- sapply(X = genome.object@chromosomes, FUN = function(chromosome) chromosome@chr.len)
  # Pull out the map of each chromosome
  chr.maps <- lapply(genome.object@chromosomes, function(chromosome) chromosome@pos.snps)
  # Pull out the index of the markers and qtl on each chromosome
  marker.indices <- lapply(genome.object@chromosomes, function(chromosome) chromosome@pos.markers$chr.index)
  qtl.indicies <- lapply(genome.object@chromosomes, function(chromosome) chromosome@pos.add.qtl$chr.index)

  
  
  ### More error handling and valiation ###
  # 1. The haploid genomes must be the same length as the number of snps
  if (length(hap.genome1) != sum(chr.snps) | length(hap.genome2) != sum(chr.snps)) 
    stop ("The length of the haploid genome(s) provided is not equal to the number of SNPs in the genome.")
  
  
  # Split the vector of haploid genotypes into a list of chromosomes
  hap.genome1.split <- split(hap.genome1, rep(1:n.chr, chr.snps))
  hap.genome2.split <- split(hap.genome2, rep(1:n.chr, chr.snps))
  
  # Create a list of the list
  hap.genome.split.list <- list(hap.genome1.split = hap.genome1.split, hap.genome2.split = hap.genome2.split)
  
  
  ### Mutation ###
  
  if (mutate) {
    
    # Apply a function over the list of split haploid genome
    new.hap.genome.split.list <- lapply(X = hap.genome.split.list, FUN = function(hap.genome.split)
    
      # Whether a snp is mutated will be determined by drawing from a Binomial 
      ## distribtion with success probability theta = mu.snp
      lapply(X = seq(n.chr), FUN = function(i) {
        
        markers.to.mutate.i <- rbinom(n = chr.markers[i], size = 1, prob = mu.snp)
        qtl.to.mutate.i <- rbinom(n = chr.qtl[i], size = 1, prob = mu.qtl)
        
        # If a marker is to be mutated, change the 1 to 0 or 0 to 1
        if (sum(markers.to.mutate.i) > 0) {
          
          # Find the absolute value of the marker genotypes minus the output
          ## from sampling the binomial distribution
          ## The reason for the absolute value is that it makes it simple to convert
          ## 1 to 0 and 0 to 1. For instance, abs( 1 - 1 ) = 0, and abs( 0 - 1 ) = 1
          new.marker.geno.i <- abs( hap.genome.split[[i]][marker.indices[[i]]] - markers.to.mutate.i )
          # Do the same for qtl
          new.qtl.geno.i <- abs( hap.genome.split[[i]][qtl.indicies[[i]]] - qtl.to.mutate.i )
          
          # Create a combined vector of genotypes
          new.genotypes.i <- numeric()
          new.genotypes.i[marker.indices[[i]]] <- new.marker.geno.i
          new.genotypes.i[qtl.indicies[[i]]] <- new.qtl.geno.i
          
          # Return the new chromosome
          return(new.genotypes.i)
          
        } else { # With no mutations, just return the chromosome as is
          
          return(hap.genome.split[[i]])
          
        }
      }) # Close the sapply
    ) # Close the lapply
    
  } else { # If not mutating, simply move the unmutated hap.genome.split.list to
    ## a "new" version
    new.hap.genome.split.list <- hap.genome.split.list
    
  } # Close the mutate if statement
  
  
  ### Recombination ###
  
  # Apply across the number of chromosomes
  recom.hap.split <- lapply(X = seq(n.chr), FUN = function(i) {
    
    # The number of cross-overs is sampled from a Poisson distribution 
    ## with parameter lambda = length of chromosome
    n.cx <- rpois(n = 1, lambda = chr.len[i])
    
    # Randomly determine which homologue will be inherited to the gamete
    hap.order <- sample(c(1,2))
    
    # If the number of cross-overs is 0, simply sample a homologue to be intactly
    ## inherited
    if (n.cx == 0) return(new.hap.genome.split.list[[hap.order[1]]][[i]])
    
    # The position of the cross-overs is sampled from a uniform distribution with
    ## n = n.cx in the interval (0, length of chromosome)
    pos.cx <- sort(runif(n = n.cx, min = 0, max = chr.len[i]))
    
    # Determine the SNP after which a cross-over occurs, plus add the last snp index
    cx.snp.ind <- c(findInterval(x = pos.cx, chr.maps[[i]]), chr.snps[i])
    # If there is a zero, remove it
    cx.snp.ind <- cx.snp.ind[cx.snp.ind != 0]
    
    # Concatenate the different homologue sections that result
    do.call("c", lapply(X = 1:length(cx.snp.ind), FUN = function(m) {
      
      # Determine the index of snps
      homologue.snp.ind <- setdiff(seq(cx.snp.ind[m]), seq(cx.snp.ind[m-1]))
      
      # Determine the homologue that the section will come from
      which.homologue <- hap.order[(m %% 2) + 1]
      
      # Pull out the homologue for the snp indices
      new.hap.genome.split.list[[which.homologue]][[i]][homologue.snp.ind]
      
    }) ) # Close the inner lapply
    
  }) # Close the outer lapply 
  
  # Return a concatenated, mutated, recombined haploid gamete
  return(do.call("c", recom.hap.split))

} # Close the function
