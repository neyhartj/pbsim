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
  
  # Convert the map to a table
  map_table <- map_to_table(genome)
  
  map_table %>% 
    mutate(marker = row.names(.)) %>% 
    select(marker, chr, pos)
  
}



#' Quicker procedures for PopVar
#' 
#' @param genome An object of class \code{genome}.
#' @param map An object of class map to simulate the virtual bi-parental populations.
#' @param training.pop An object of class \code{pop} with the elements \code{geno} and 
#' \code{pheno_val}. This is used as the training population.
#' @param founder.pop An object of class \code{pop} with the geno information for
#' the parents. Additional individuals can be present in \code{parent_pop}. They
#' will be filtered according to the parents in the \code{crossing_block}.
#' @param crossing.block A crossing block detailing the crosses to make. Must be a
#' \code{data.frame} with 2 columns: the first gives the name of parent 1, and the 
#' second gives the name of parent 2. See \code{\link{sim_crossing_block}}.
#' @param n.rep The number of virtual bi-parental populations to simulate per cross.
#' @param n.ind The number of individual in each virtual bi-parental population.
#' @param tail.p The percentile to calculate the superior progeny mean. The 
#' \emph{tail.p}% of predicted genotypic values is used to predict the superior
#' progeny mean.
#' 
#' @examples 
#' 
#' # Simulate a genome
#' n.mar  <- c(505, 505, 505)
#' len <- c(120, 130, 140)
#' 
#' genome <- sim_genome(len, n.mar)
#' 
#' # Simulate another genome for use in creating the bi-parental populations
#' genome_use <- sim_genome(len, n.mar)
#' 
#' # Simulate a quantitative trait influenced by 50 QTL
#' qtl.model <- matrix(NA, 50, 4)
#' genome <- sim_gen_model(genome = genome, qtl.model = qtl.model, 
#'                         add.dist = "geometric", max.qtl = 50)
#'                         
#' # Remove QTL from the second genome map
#' map_use <- lapply(genome_map$map, function(chr) chr[setdiff(names(chr), qtlnames(genome))])
#' 
#' # Simulate the genotypes of eight founders
#' founder_pop <- sim_founders(genome, n.str = 8)
#' founder_pop <- sim_phenoval(pop = founder_pop, h2 = 0.5)
#' 
#' ped <- sim_pedigree(n.ind = 100, n.selfgen = 2)
#' 
#' # Extract the founder names
#' parents <- indnames(founder_pop)
#' 
#' # Generate a crossing block with 5 crosses
#' cb <- sim_crossing_block(parents = parents, n.crosses = 5)
#' 
#' # Use the map in the second genome to predict the genetic variance in the 
#' # bi-parental populations
#' pop_predict_quick(genome = genome, map = map_use, training.pop = founder_pop, founder.pop = founder_pop, 
#'                   crossing.block = cb, n.rep = 3, n.ind = 50, tail.p = 0.1)
#' 
#' 
#' 
#' @import dplyr
#' @importFrom qtl sim.cross
#' @importFrom purrr map
#' @importFrom purrr map_df
#' 
#' @export
#' 
pop_predict_quick <- function(genome, map, training.pop, founder.pop, crossing.block,
                              n.rep, n.ind, tail.p) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # First predict marker effects
  mar_eff <- pred_mar_eff(genome = genome, training.pop = training.pop)$mar_eff %>%
    select(-marker) %>%
    as.matrix()
  
  # Seq over the number of markers
  j <- seq(nmar(genome))
  
  # Iterate over the crossing block
  sim_fam <- crossing.block %>%
    group_by(parent1, parent2) %>%
    do({
      
      # Subset the founders
      parent_genos <- subset_pop(founder.pop, c(.$parent1, .$parent2)) %>%
        genotype(genome = genome, pop = .)
      
      # Simulate populations using qtl - replicate
      cross_preds <- replicate(n.rep, expr = {
        sim.cross(map = map, n.ind = n.ind, type = "riself") }, 
        simplify = FALSE) %>%
        
        map(function(cross)
          do.call("cbind", lapply(cross$geno, "[[", "data")) ) %>%
        map_df(function(geno) {
          
          # Recode to parental
          recode <- t(apply(geno, MARGIN = 1, FUN = function(ind) parent_genos[cbind(ind, j)]))

          # Calculate PGV
          pgv <- recode %*% mar_eff
          
          # Mean, variance, and superior progeny
          data.frame(pred_mu = mean(pgv), pred_mu_sp = mean(sort(pgv, decreasing = TRUE)[seq(n)]),
            pred_V_G = as.numeric(var(pgv))) })
      
      
      # Summarize
      cross_preds %>% 
        summarise_each(funs(mean, sd)) })
  
  # Output
  as.data.frame(sim_fam)
  
} # Close the function
        
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

