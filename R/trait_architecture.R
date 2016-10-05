#' Define the genetic architecture of a quantitative trait
#' 
#' @description Take an object of class "genome" and define the genetic architecture of a quantitative trait by specifying the number of QTL, the positon of the QTL, and the allelic effects of the QTL.
#' 
#' @param genome.object An object of class \code{genome} created by the 
#' \code{make.genome} function.
#' @param n.qtl Integer of the total number of QTL.
#' @param qtl.ind Numeric specifying the index of SNPs to become QTL. For 
#' instance, if the 1st, 25th, and 50th SNPs (which may be on different 
#' chromosomes) are to become QTL, then \code{qtl.ind} would take the form 
#' \code{c(1, 25, 50)}. If \code{NULL}, QTL are randomly sampled from all 
#' possible SNPs. 
#' Default is \code{NULL}.
#' @param qtl.dom.ind Numeric specifying the index of QTL exhibiting dominance. 
#' For instance, if 100 QTL are assigned positions by \code{qtl.pos} and the 
#' 10th, 50th, and 70th QTL (which may be on different chromsomes) will exhibit 
#' dominance, then \code{qtl.dom.ind} would take the form \code{c(10, 50, 70)}. 
#' If \code{NULL}, no QTL will be assigned dominance effects. Default is \code{NULL}.
#' @param qtl.add.eff Character or list of numeric vectors specifying the additive 
#' allelic effects (\code{a}) of QTL. Character options include "normal" or 
#' "geometric." See Details for a description of these options. This argument 
#' may also take the form of a numeric vector of length \code{n.qtl} with the 
#' values of \code{a} for each QTL. See Details for more information.
#' @param qtl.dom.eff Numeric specifying the dominance effect (\code{d}) of QTL. 
#' If specified, must be of length \code{n.qtl} with QTL with no desired 
#' dominance effect being assigned a value of 0. If \code{NULL}, no QTL will 
#' exhibit dominance effects. Default is \code{NULL}.
#' 
#' @details 
#' Once a genome is created, a next step is to define the architecture of a 
#' simulated quantitative trait. Empirical evidence suggests that different quantitative
#' traits possess different genetic characteristics, including number and effect 
#' of quantitative trait loci (QTL). This function provides flexibility as to defining
#' those characteristics of a simulated trait. 
#' 
#' Besides providing the \code{genome} object, the only other requirement of this
#' function is the number of QTL. The only restriction to this input is that
#' it cannot be greater than the total number of SNP loci in the genome.
#' 
#' The remaining arguments of this function allow fine-tuning of trait architecture.
#' The \code{qtl.ind} arguments takes a numeric vector of indices of SNPs that will be converted to QTL. The 
#' remaining SNPs become markers. For instance, if the genome contains 500 SNP 
#' loci, and the user wants the trait to be controlled by 50 QTL, 450 SNPs will remain 
#' as markers. If this argument is \code{NULL}, SNPs across the genome are randomly
#' sampled to become QTL, regardless of the number of SNPs on each chromosome.
#' 
#' The \code{qtl.dom.ind} is similar to \code{qtl.ind}, but it instead defines the
#' index of \emph{QTL} that have a non-zero dominance (\eqn{d}) value. For instance,
#' if 50 QTL control your simulated trait, and you want the 10th, 50th, and 70th
#' QTL to have some dominance effect, the input for \code{qtl.dom.ind} would be
#' \code{c(10, 50, 70)}.
#' 
#' The genotypic values of QTL are defined by the \code{qtl.add.eff} and 
#' \code{qtl.dom.eff} arguments. The \code{qtl.add.eff} argument determines the
#' genotypic values of the two homozygotes at a particular QTL. These values are
#' defined relative to the midparent value (\eqn{P.bar} = 0), where the value
#' of the unfavorable homozygote is \eqn{-a} and the values of the favorable 
#' homozygote is \eqn{a}, assuming that higher trait values are favorable. The
#' \code{qtl.dom.eff} determines the genotypic value of the heterzygote (\eqn{d}).
#' 
#' Two presets are included for defining the genotypic values of QTL. The \code{"geometric"}
#' option sets genotypic values from a geometric series:
#' \deqn{a = \frac{{L - 1}{L + 1}}}{a = (L - 1) / (L + 1)}
#' where \eqn{L} is the number of QTL and the value of the favorable homozygote 
#' for the \emph{k}th QTL is \eqn{a^k} for \eqn{k = 1, 2, \ldots, L}.
#' The justification behind this distribution of genotypic values is many quantitative
#' traits seem to be controlled by few loci of major effect and many loci of 
#' minor effect (Lande and Thompson, 1990). Note that when choosing this option,
#' QTL are randomly selected to have the value \eqn{-a} or \eqn{a} for homozygotes
#' of the "1" allele.
#' 
#' The other preset is \code{"normal"}, where values are sampled from a normal
#' distribution such that \eqn{a_k ~ N(0,1)}.
#' 
#' Besides presets, a user may specify the values \eqn{a} of QTL manually.
#' 
#' Finally, when determining breeding values, \eqn{a} is considered the value of
#' homozygotes for the "1" allele at a locus.
#' 
#' @examples
#' n.chr = 3
#' chr.len = c(1.0, 1.1, 1.2)
#' chr.snps = c(15, 16, 17)
#' genome <- make.genome(n.chr = n.chr, chr.len = chr.len, chr.snps = chr.snps)
#' genome <- trait.architecture(genome.object = genome, n.qtl = 10)
#' 
#' @export
#' @return An object of class \code{genome} with the specified trait architecture attributes.
#' 
trait.architecture <- function(genome.object, n.qtl, qtl.ind = NULL, qtl.dom.ind = NULL, qtl.add.eff = "geometric", qtl.dom.eff = NULL) {
  
  # Pull out genomic information
  n.chr <- genome.object@n.chr
  chr.len <- sapply(X = genome.object@chromosomes, FUN = function(chr) chr@chr.len)
  n.snps <- genome.object@n.snps
  n.snps.chr <- sapply(X = genome.object@chromosomes, FUN = function(chr) chr@n.snps)
  
  # Vector of loci indicies from 1 to n.snps
  snp.ind <- seq(n.snps)
  # Label the vector with the chromosome
  names(snp.ind) <- rep(1:n.chr, n.snps.chr)
  # Split the index into chromosomes
  snp.ind.chr <- split(x = snp.ind, rep(1:n.chr, n.snps.chr))
  
  # Error handling
  n.qtl <- as.integer(n.qtl)
  
  # If the positions of the QTL are not provided, randomly sample all possible
  ## snp positions
  if (is.null(qtl.ind)) {
    
    # Randomly sample the total number of snps
    qtl.ind <- sort(sample(n.snps, n.qtl))
    
    # Find the intersection of qtl.ind in snp.ind
    qtl.ind <- snp.ind[snp.ind %in% qtl.ind]
    
  } else {
    
    # Error handling
    qtl.ind <- as.numeric(qtl.ind)
    # If the length of qtl.ind is not the number of qtl, error out
    if (length(qtl.ind) != n.qtl) stop("The length of qtl.ind does not equal the number of QTL")
    # If the index is not within the total snp index, error out
    if (!all(findInterval(x = qtl.ind, c(0, n.snps+1)) == 1)) stop("The indices called by qtl.ind are not in the range of all SNP indices.")

    # Determine the chromosomes by intersecting to snp.ind
    qtl.ind <- snp.ind[snp.ind %in% qtl.ind]
    
  }
  
  # Verify the index of dominance QTL
  # If qtl.dom.ind is null, assign a value of 0 to every qtl
  if (is.null(qtl.dom.ind)) {
    qtl.dom.ind <- NULL
    
  } else {
    
    # Error checking
    qtl.dom.ind <- as.numeric(qtl.dom.ind)
    # If the index is not within the total qtl index, error out
    if (!all(findInterval(x = qtl.dom.ind, c(0, n.qtl+1)) == 1)) stop("The indices called by qtl.ind are not in the range of all SNP indices.")
    
    # Determine the chromosomes by intersecting to snp.ind
    names(qtl.dom.ind) <- names(qtl.ind)[qtl.dom.ind]
  }
  
  
  # Determine the QTL additive effects
  # If qtl.add.eff is a character...
  if (is.character(qtl.add.eff)) {
    
    # Make sure it is within the acceptable bounds
    acceptable.inputs <- c("normal", "geometric")
    if (!qtl.add.eff %in% acceptable.inputs) stop("The input for qtl.add.eff is not acceptable. Choices are normal or geometric.")
    
    # If normal, sample a normal distribution
    if (qtl.add.eff == "normal") {
      
      qtl.add.eff <- rnorm(n.qtl, 0, 1)
      
      # Assign chromosomes
      names(qtl.add.eff) <- names(qtl.ind)
      
    }
    
    # If geometric, sample a geometric distribution
    if (qtl.add.eff == "geometric") {
      
      # Set the base a
      a = (n.qtl - 1) / (n.qtl + 1)
      # Determine additive effects
      qtl.add.eff <- a ^ (1:n.qtl)
      # Randomly sample QTL to have either a positive or negative effect of allele 1
      qtl.add.eff <- qtl.add.eff * sample(c(-1,1), n.qtl, replace = T)
      
      # Assign chromosomes
      names(qtl.add.eff) <- names(qtl.ind)
      
    }
    
  } else { # If qtl.add.eff is not a character
    
    if (is.numeric(qtl.add.eff)) { # If numeric...
      
      # Error handling
      if (length(qtl.add.eff) != n.qtl) stop("The length of qtl.add.eff does not equal the number of QTL.")
      
      # Assign chromosomes
      names(qtl.add.eff) <- names(qtl.ind)
      
    } else { # If not numeric, exit
      stop("The class of qtl.add.eff is not character or numeric.")
    }
  }
  
  
  # Determine qtl dominance effects
  if (is.null(qtl.dom.eff)) {
    
    # If no input is given for qtl.dom.eff, but an input is given for qtl.dom.ind, error out
    if(!is.null(qtl.dom.ind)) stop("QTL with dominance were specified by qtl.dom.ind, but no dominance effects were assigned.")
    
    # Assign 0 as the dominance effect
    qtl.dom.eff <- rep(0, n.qtl)
    # Rename
    names(qtl.dom.eff) <- names(qtl.ind)
    
  } else {
    
    # Error handling
    qtl.dom.eff <- as.numeric(qtl.dom.eff)
    # If the length is not the same as that of qtl.dom.ind, error out
    if(length(qtl.dom.eff) != length(qtl.dom.ind)) stop("The length of qtl.dom.eff is not the same as the number of QTL designated to have dominance effects.")
    
    # Fill in the dominance effects of QTL without dominance with 0
    qtl.dom.eff.i <- rep(0, n.qtl)
    qtl.dom.eff.i[qtl.dom.ind] <- qtl.dom.eff
    qtl.dom.eff <- qtl.dom.eff.i
    
    # Rename
    names(qtl.dom.eff) <- names(qtl.ind)
  }
  
  
  # Split the qtl index by chromosomes
  qtl.chr.split <- split(seq(n.qtl), as.numeric(names(qtl.ind)))
  
  # Add the information to the chromosomes
  genome.object@chromosomes <- sapply(X = seq(n.chr), FUN = function(i) {
    
    # Pull out the chromosome
    chromosome <- genome.object@chromosomes[[i]]
    
    qtl.ind.i <- qtl.ind[qtl.chr.split[[i]]]
    qtl.chr.ind.i <- which(snp.ind.chr[[i]] %in% qtl.ind.i)
    qtl.add.eff.i <- qtl.add.eff[qtl.chr.split[[i]]]
    qtl.dom.ind.i <- split(as.numeric(qtl.dom.ind), as.numeric(names(qtl.dom.ind)))[[as.character(i)]]
    qtl.dom.chr.ind.i <- which(qtl.chr.split[[i]] %in% qtl.dom.ind.i)
    qtl.dom.eff.i <- qtl.dom.eff[qtl.chr.split[[i]]]
    
    # Create the content to add to chromosomes
    n.add.qtl <- as.integer(length(qtl.ind.i))
    n.dom.qtl <- as.integer(length(qtl.dom.ind.i))
    pos.add.qtl <- list( genome.index = as.numeric(qtl.ind.i), # Index of QTL across genome
                         chr.index = qtl.chr.ind.i, # Index of QTL on chromosome
                         M = as.numeric(chromosome@pos.snps[qtl.chr.ind.i]) )
    pos.dom.qtl <- list( genome.index = qtl.dom.ind.i, 
                         chr.index = qtl.dom.chr.ind.i,
                         M = as.numeric(chromosome@pos.add.qtl$M)[qtl.dom.ind.i] )
    
    n.markers <- chromosome@n.markers - n.add.qtl
    pos.markers <- list( genome.index = setdiff(snp.ind.chr[[i]], pos.add.qtl$genome.index),
                         chr.index = setdiff(chromosome@pos.markers$index, pos.add.qtl$chr.index),
                         M = setdiff(chromosome@pos.markers$M, pos.add.qtl$M) )
    
    
    
    # Add content to the chromosomes
    chromosome@n.add.qtl <- n.add.qtl
    chromosome@n.dom.qtl <- n.dom.qtl
    chromosome@pos.add.qtl <- pos.add.qtl
    chromosome@pos.dom.qtl <- pos.dom.qtl
    chromosome@qtl.effects$a <- as.numeric(qtl.add.eff.i)
    chromosome@qtl.effects$d <- as.numeric(qtl.dom.eff.i)
    
    # Adjust content of chromosomes
    chromosome@n.markers <- n.markers
    chromosome@pos.markers <- pos.markers
    
    # Return the chromosome
    return(chromosome) })
  
  # Adjust content for the whole genome
  genome.object@n.qtl <- as.integer(n.qtl)
  
  # Return the genome
  return(genome.object)
  
} # Close the function