#' Fill in genotype data for QTL
#' 
#' @description 
#' Generates genotype data for the QTL in a genome. 
#' 
#' @param genome An object of class \code{genome}.
#' @param geno Genotype data on a population to phenotype. Must be a matrix of dimensions
#' \code{n.ind} x \code{n.marker}, the elements of which must be z {0, 1, 2}. Genotype data for
#' the QTL should not be included.
#' 
#' @details 
#' For now, genotypes of the QTL are copied exactly from the genotypes of the nearest 
#' marker. Other options may be considered later.
#' 
#' @return 
#' A matrix of dimensions \code{n.ind} x \code{n.loci}, with genotypes for the QTL included.
#' 
#' @importFrom qtl sim.geno
#' @importFrom qtl find.flanking
#' 
#' @export
#' 
fill_qtl_geno <- function(genome, geno) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # The geno input should have n_marker columns
  if (ncol(geno) != nmar(genome) )
    stop("The number of markers in the 'geno' input is not equal to the number
         of markers in the genome.")
  
  ## Find the flanking markers of all of the QTL
  # First create a blank cross object
  blank_cross <- qtl::sim.cross(map = genome$map, n.ind = 1)
  
  # Pull the unique QTL
  unique_qtl <- pull_qtl(genome)
  
  # Vector of chromosomes of QTL
  chr <- chrnames(genome)[unique_qtl$chr]
  
  ## Find the closest marker to each QTL
  # First pull out the map without QTL
  map_noqtl <- qtl::pull.map(blank_cross)
  
  # Find the closest marker
  closest_mar <- apply(X = unique_qtl, MARGIN = 1, FUN = function(qtl) {
    ind <- which.min(abs(as.numeric(qtl[2]) - map_noqtl[[as.numeric(qtl[1])]]))
    names(map_noqtl[[as.numeric(qtl[1])]])[ind] })
  
  # Add to qtl info
  qtl_info <- unique_qtl %>% 
    mutate(close = closest_mar)
  
  
  # Fill in the qtl genotype by selecting the marker that is closest
  qtl_geno <- pull_genotype(genome = genome, geno = geno, loci = qtl_info$close)
  # Rename the QTL
  colnames(qtl_geno) <- paste("QTL", seq(nrow(qtl_info)), sep = "")
  
  
  # Names of all loci
  locinames <- unlist(lapply(X = genome$map, names))
  
  
  ## Add the QTL genotypes to the geno matrix
  # First make a blank matrix
  new_geno <- matrix(data = NA, nrow = nrow(geno), ncol = nloci(genome), 
                     dimnames = list(row.names(geno), locinames))
  
  # Use the new map to add genotypes to the new genotype matrix
  new_geno[,colnames(new_geno) %in% colnames(geno)] <- geno
  new_geno[,colnames(new_geno) %in% colnames(qtl_geno)] <- qtl_geno
  
  # Return the genotypes
  return(new_geno)
  
} # Close the function