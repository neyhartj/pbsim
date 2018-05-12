#' Simulate founders
#' 
#' @description
#' Simulates the SNP genotype states for a set of founders. This is a simple wrapper of 
#' \code{\link[qtl]{simFounderSnps}} that allows for bi-parental populations
#' 
#' @param genome An object of class \code{genome}.
#' @param n.str The number of founders (either 2, 4, or 8)
#' @param pat.freq A vector of length \code{n.str}/2 + 1. Frequency of SNP genotype 
#' patterns in the founder, where the first vector element is the frequency of monoallelic
#' genotypes (i.e. fixed), and the second element is the frequency of polymorphic genotypes. 
#' In the case of 4 or 8 founders, the second element is the frequency of SNPs unique to 1 founder, 
#' and the third element is the frequency of SNPs in 2 founders. In the case of 8 founders,
#' the fourth element is the frequency of SNPs in 3 out of 8 founders, and the fifth element
#' is the frequency of SNPs in 4 out of 8 founders.
#' 
#' @details 
#' See \code{\link[qtl]{simFounderSnps}}. SNPs are simulated to be linkage equilibrium.
#' 
#' @return 
#' See \code{\link[qtl]{simFounderSnps}}. A \code{list} of the same length as \code{map},
#' with each element being a matrix of 0's and 2's.
#' 
#' @examples 
#' # Simulate a genome
#' n.mar  <- c(505, 505, 505)
#' len <- c(120, 130, 140)
#' 
#' genome <- sim_genome(len, n.mar)
#' 
#' # Simulate a quantitative trait influenced by 50 QTL
#' qtl.model <- matrix(NA, 50, 4)
#' genome <- sim_gen_model(genome = genome, qtl.model = qtl.model, 
#'                         add.dist = "geometric", max.qtl = 50)
#' 
#' # Simulate the founder genotypes
#' founder_geno <- sim_founders(genome)
#' 
#' @importFrom qtl sim.cross
#' @importFrom qtl pull.map
#' 
#' @export
#' 
sim_founders <- function(genome, n.str = c("2", "4", "8"), pat.freq) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Get the genome type
  type <- attr(genome, "type")
  
  # Extract the map, depending on type
  if (type == "pbsim") {
    map <- genome$map
    
  } else {
    map <- lapply(X = genome$hypredGenomes, function(hyp_chr) {
      # Extract and convert to cM
      structure(slot(object = hyp_chr, "pos.snp") * 100, class = "A") })
    
    class(map) <- "map"
  
  }
    
    
  ## Borrowed from simFounderSnps
  if (is.numeric(n.str)) 
    n.str <- as.character(n.str)
  
  n.str <- as.numeric(match.arg(n.str))
  
  # Check if pat.freq is called
  if (missing(pat.freq)) {
    
    if (n.str == 8) 
      pat.freq <- c(0, 0.4, 0.3, 0.2, 0.1)
      
    if (n.str == 4) 
      pat.freq <- c(0, 0.7, 0.3)
    
    if (n.str == 2)
      pat.freq <- c(0.3, 0.7)
    
  }
  
  # Check the length of pat.freq
  if (length(pat.freq) < n.str/2 + 1) {
    # Append 0s
    pat.freq <- c(pat.freq, rep(0, n.str/2 + 1 - length(pat.freq)))
    
  } else {
    pat.freq <- pat.freq[seq(n.str/2 + 1)]
    
  }
  
  # If the number of founders is 2, the first element of pat.freq should be 1 - the second element
  if (n.str == 2)
    pat.freq <- c(1 - pat.freq[2], pat.freq[2])
  
  # Proportions
  pat.freq <- pat.freq / sum(pat.freq)
  
  # Find the number of markers
  n.mar <- nloci(genome, by.chr = TRUE)
  # Get the names of markers
  mar_names <- markernames(genome, include.qtl = TRUE)
  # Create founder names
  founder_names <- paste("founder", seq(n.str), sep = "")
  
  
  # Empty vector
  output <- structure(vector("list", nchr(genome)), names = chrnames(genome))

  # Iterate over the map
  for (i in seq_along(map)) {
    
    # Sample 0s and 1s 
    thepat <- sample(seq_along(pat.freq) - 1, n.mar[i], prob = pat.freq, replace = TRUE)
    
    # Empty matrix
    output[[i]] <- matrix(0, nrow = n.str, ncol = n.mar[i])
      
    for (j in seq_along(thepat)) {
      
      output[[i]][sample(seq(n.str), thepat[j]), j] <- 2
    }
  }
  
  # Bind the list into a matrix
  geno <- do.call("cbind", output)
  
  # Add row and column names
  dimnames(geno) <- list(founder_names, mar_names)
  
  # If the genome is hypred, convert to haploids
  if (type == "hypred") {
    
    # Split the genos by individual
    geno_split <- split(geno, 1:nrow(geno), drop = FALSE)
    geno_split <- lapply(geno_split, FUN = function(ind) rbind(ind / 2, ind / 2) )
    
    # Empty array
    haploid <- array(data = NA, dim = c(2, ncol(geno), length(geno_split)),
                     dimnames = list(NULL, colnames(geno), row.names(geno)))
    
    for (k in seq_along(geno_split)) {
      haploid[,,k] <- geno_split[[k]]
    }
    
    geno <- haploid
     
  }
  
  # Create the pop object and return
  create_pop(genome = genome, geno = geno)

} # Close the function
