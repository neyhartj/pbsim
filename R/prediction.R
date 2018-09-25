#' Predict marker effects in a training population
#' 
#' 
#' @param genome An object of class \code{genome}.
#' @param training.pop An object of class \code{pop} with the elements \code{geno} and 
#' \code{pheno_val}. This is used as the training population.
#' @param method The statistical method to predict marker effects. If \code{"RRBLUP"}, the
#' \code{\link[qtl]{mixed.solve}} function is used. Otherwise, the \code{\link[BGLR]{BGLR}}
#' function is used.
#' @param n.iter,burn.in,thin Number of iterations, number of burn-ins, and thinning, respectively. See 
#' \code{\link[BGLR]{BGLR}}.
#' @param save.at See \code{\link[BGLR]{BGLR}}.
#' 
#' @details 
#' The \code{training.pop} must have phenotypic values associated with each entry.
#' The mean phenotype is used as training data in the model. Genotypic data (excluding
#' QTL) are used to predict marker effects.
#' 
#' 
#' @return 
#' The \code{training.pop} with predicted marker effects.
#' 
#' @examples 
#' 
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
#' # Simulate the populations according to the crossing block
#' pop <- sim_family_cb(genome = genome, pedigree = ped, founder.pop = founder_pop, 
#'                      crossing.block = cb)
#'                      
#' # Use the founders as a training population for the progeny
#' # Predict marker effects
#' training.pop <- pred_mar_eff(genome = genome, training.pop = founder_pop)
#'                      
#' @importFrom  rrBLUP mixed.solve
#' @importFrom BGLR BGLR
#' @import dplyr
#' 
#' @export
#' 
pred_mar_eff <- function(genome, training.pop, method = c("RRBLUP", "BRR", "BayesA", "BL", "BayesB", "BayesC"), n.iter = 1500, 
                         burn.in = 500, thin = 5, save.at = "") {

  # Check the populations
  if (!inherits(training.pop, "pop"))
    stop("The input 'training.pop' must be of class 'pop'.")
  
  # Make sure the training population has phenotypes
  if (is.null(training.pop$pheno_val))
    stop("The 'training.pop' must have phenotypic values.")
  
  # Match the method arg
  method <- match.arg(method)
  
  # Genotype the training set and candidates
  training_geno <- genotype(genome = genome, pop = training.pop)

  # Pull out the training phenos
  training_pheno <- training.pop$pheno_val$pheno_mean %>% 
    data.frame(row.names = .$ind) %>% 
    select(-ind) %>%
    as.matrix()
  
  # Run by method
  if (method == "RRBLUP") {
    
    # Iterate over traits
    mar_eff <- apply(X = training_pheno, MARGIN = 2, FUN = function(trait) {
      
      # Solve the mixed model
      solve_out <- mixed.solve(y = trait, Z = training_geno, method = "REML")
      # Pull out marker effects
      list(effects = solve_out$u)
      
    })
    
  } else {
    
    ## Use BGLR
    # Iterate over traits
    mar_eff <- apply(X = training_pheno, MARGIN = 2, FUN = function(trait) {
      
      # Solve the mixed model
      solve_out <- BGLR(y = trait, ETA = list(list(X = training_geno, model = method)), nIter = n.iter,
                        burnIn = burn.in, thin = thin, verbose = FALSE, saveAt = save.at)
      
      # Return the marker effects and the probIn parameter
      list(effects = solve_out$ETA[[1]]$b, pi = solve_out$ETA[[1]]$probIn)
      
    })
    
  }

  
  ## Reorganize the marker effects
  marker_eff <- sapply(X = mar_eff, "[[", "effects")
  hyperparam <- sapply(X = mar_eff, "[[", "pi")
  
  # Return the TP with the marker effects
  training.pop$mar_eff <- data.frame(marker = row.names(marker_eff), marker_eff, row.names = NULL, stringsAsFactors = FALSE)
  training.pop$mar_eff_meta <- data.frame(param = "pi", t(hyperparam), row.names = NULL, stringsAsFactors = FALSE)
  
  return(training.pop)
  
}
  
  

#' Predict genotypic values using genomewide markers
#' 
#' @param genome An object of class \code{genome}.
#' @param training.pop An object of class \code{pop} with the elements \code{geno} and 
#' \code{pheno_val}. This is used as the training population. If marker effects
#' are present, they are used.
#' @param candidate.pop An object of class \code{pop} with the element \code{geno}.
#' Genotypic values are predicted for individuals in this object.
#' @param method The statistical method to predict marker effects. If \code{"RRBLUP"}, the
#' \code{\link[qtl]{mixed.solve}} function is used. Otherwise, the \code{\link[BGLR]{BGLR}}
#' function is used.
#' @param n.iter,burn.in,thin Number of iterations, number of burn-ins, and thinning, respectively. See 
#' \code{\link[BGLR]{BGLR}}.
#' @param save.at See \code{\link[BGLR]{BGLR}}.
#' 
#' @details 
#' The \code{training.pop} must have phenotypic values associated with each entry.
#' The mean phenotype is used as training data in the model. Genotypic data (excluding
#' QTL) are used to predict marker effects, which are then used to predict the 
#' genotypic value of the individuals in the \code{candidate.pop}.
#' 
#' 
#' @return 
#' The \code{candidate.pop} with predicted genotypic values.
#' 
#' @examples 
#' 
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
#' # Simulate the populations according to the crossing block
#' pop <- sim_family_cb(genome = genome, pedigree = ped, founder_pop = founder_pop, 
#'                      crossing_block = cb)
#'                      
#' # Use the founders as a training population for the progeny
#' pop <- pred_geno_val(genome = genome, training.pop = founder_pop, candidate.pop = pop)
#' 
#' ## Alternatively, predict marker effects first, then predict genotypic values
#' training.pop <- pred_mar_eff(genome = genome, training.pop = founder_pop)
#' pop <- pred_geno_val(genome = genome, training.pop = founder_pop, candidate.pop = pop)
#'                      
#' @importFrom  rrBLUP mixed.solve
#' @importFrom BGLR BGLR
#' @import dplyr
#' 
#' @export
#' 
pred_geno_val <- function(genome, training.pop, candidate.pop, method = c("RRBLUP", "BRR", "BayesA", "BL", "BayesB", "BayesC"), 
                          n.iter = 1500, burn.in = 500, thin = 5, save.at = "") {
  
  # Check the populations
  if (!inherits(training.pop, "pop"))
    stop("The input 'training.pop' must be of class 'pop'.")
  
  if (!inherits(candidate.pop, "pop"))
    stop("The input 'candidate.pop' must be of class 'pop'.")
  
  # Make sure the training population has phenotypes
  if (is.null(training.pop$pheno_val))
    stop("The 'training.pop' must have phenotypic values.")
  
  # Match the method arg
  method <- match.arg(method)

  # Genotype the training set and candidates
  training_geno <- genotype(genome = genome, pop = training.pop)
  candidate_geno <- genotype(genome = genome, pop = candidate.pop)
  
  # Pull out the training phenos
  training_pheno <- training.pop$pheno_val$pheno_mean %>% 
    data.frame(row.names = .$ind) %>% 
    select(-ind) %>% 
    as.matrix()
  
  # Pull out marker effect
  mar_eff <- training.pop$mar_eff
  
  # If marker effects are present, use them
  if (is.null(mar_eff)) {
    
    # Predict the marker effects
    mar_eff <- pred_mar_eff(genome = genome, training.pop = training.pop, method = method, n.iter = n.iter, burn.in = burn.in, thin = thin)$mar_eff
    
    # Run by method
    if (method == "RRBLUP") {
      
      # Iterate over traits
      pgv <- apply(X = training_pheno, MARGIN = 2, FUN = function(trait) {
        
        # Solve the mixed model
        solve_out <- mixed.solve(y = trait, Z = training_geno, method = "REML")
        
        # Pull out marker effects
        marker_effects <- solve_out$u
        
        # Predict GVs
        candidate_geno %*% marker_effects })
      
    } else {
      
      ## Use BGLR
      # Iterate over traits
      mar_eff <- apply(X = training_pheno, MARGIN = 2, FUN = function(trait) {
        
        # Solve the mixed model
        solve_out <- BGLR(y = trait, ETA = list(list(X = training_geno, model = method)), nIter = n.iter,
                          burnIn = burn.in, thin = thin, verbose = FALSE, saveAt = save.at)
        
        # Return the marker effects and the probIn parameter
        list(effects = solve_out$ETA[[1]]$b, pi = solve_out$ETA[[1]]$probIn)
        
      })
      
    }
    
  }
  
  # Make sure all markers are in the candidate geno
  if (!all(mar_eff$marker %in% colnames(candidate_geno)))
    stop("Not all markers with predicted marker effects are in the 'candidate.pop'.")
  
  # Convert to matrix
  mar_eff_mat <- mar_eff %>% 
    data.frame(row.names = .$marker) %>% 
    select(-marker) %>% 
    as.matrix()
  
  # Iterate over traits
  pgv <- apply(X = mar_eff_mat, MARGIN = 2, FUN = function(trait) {
    
    # Predict GVs
    candidate_geno %*% trait })
  
  
  # Convert to data.frame
  pgv_df <- data.frame(ind = indnames(candidate.pop), pgv)
  
  # Add to the population and return
  candidate.pop$pred_val <- pgv_df
  
  return(candidate.pop)
  
} # Close the function