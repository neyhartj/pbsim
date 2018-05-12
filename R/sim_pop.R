#' Simulate a population
#' 
#' @description
#' Simulates a random population of given population size. SNPs are simulated
#' to be in linkage equilibrium.
#' 
#' @param genome A genome object.
#' @param n.ind The number of individual in the population.
#' 
#' @return 
#' A \code{pop} object.
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
#' # Simulate the population
#' pop <- sim_pop(genome = genome, n.ind = 100)
#' 
#' @export
#' 
sim_pop <- function(genome, n.ind) {
  
  # Check inputs
  stopifnot(class(genome) == "genome")
  stopifnot(class(n.ind) == "numeric")
  
  n.ind <- as.integer(n.ind)
  
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
  
  # Create individual names
  ind_names <- paste("ind", formatC(x = seq(n.ind), width = nchar(n.ind), flag = 0, format = "d"), sep = "")
  
  # Simulate markers
  geno <- lapply(X = map, FUN = function(m) {
    mar_chr <- length(m)
    M <- replicate(n = mar_chr, ifelse(runif(n = n.ind) < 0.5, 0, 2))
    `dimnames<-`(M, list(ind_names, names(m))) })
  
  geno1 <- do.call("cbind", geno)
  
  # Create the population
  create_pop(genome = genome, geno = geno1)
  
}
    
  
  
  
  
  
  
  
  
  