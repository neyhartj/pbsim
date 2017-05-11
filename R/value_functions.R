#' Calculate the genotypic value of individuals
#' 
#' @param genome An object of class \code{genome}.
#' @param geno Genotype data on a population to phenotype. Can be a matrix of dimensions
#' \code{n.ind} x \code{n.loci}, the elements of which must be z {0, 1, 2}, or a list
#' of such matrices.
#' 
#' @details 
#' Calculates the genotypic value of an individual or one or more traits by adding
#' the allele effects and dominance deviations of the QTL carried by that individual.
#' 
#' @return
#' A matrix of dimensions \code{n.ind} x \code{n.trait} with genotypic values
#' of the individuals.
#' 
#' @export 
#' 
calc_genoval <- function(genome, geno) {
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = geno))
    stop("The geno did not pass. See warning for reason.")
  
  # If geno is a list, recombine
  if (is.list(geno))
    geno <- do.call("cbind", geno)
  
  ## Iterate over traits in the genetic model
  geno_val <- lapply(X = genome$gen_model, FUN = function(qtlmod) {
    qtl_geno <- pull_genotype(genome = genome, geno = geno, loci = qtlmod$qtl_name)
    # Subtract 1
    qtl_geno1 <- qtl_geno - 1
    
    # Additive eff
    bv <- qtl_geno1 %*% as.matrix(qtlmod$add_eff)
    # Dominance - first convert matrix to dominance matrix, then calculate deviations
    qtl_dom_geno <- ifelse(qtl_geno1 != 0, 0, 1)
    dd <- qtl_dom_geno %*% as.matrix(qtlmod$dom_eff)
    
    # Sum
    bv + dd })
  
  # Bind columns
  geno_val1 <- do.call("cbind", geno_val)
  
  # Add trait names
  colnames(geno_val1) <- paste("trait", seq(ncol(geno_val1)), sep = "")
  
  # Convert the row.names to a column and return the data.frame
  data.frame(ind = row.names(geno_val1), geno_val1, row.names = NULL, stringsAsFactors = FALSE)
  
} # Close the function


#' Simulate phenotypic data
#' 
#' @description 
#' Simulates phenotypic observations given the genotypic value of individuals and the 
#' heritability of a quantitative trait.
#' 
#' @param pop An object of class \code{pop}.
#' @param h2 The heritability of the trait or traits. May be a numeric vector of 
#' length 1 (heritability is the same for all traits) or a numeric vector of 
#' length n_trait. 
#' @param n.env The number of environments in which to phenotype.
#' @param n.rep The number of replicates of each individual in each environment.
#' @param ... Other arguments. See \code{Details}.
#' 
#' @details
#' 
#' Other arguments that can be specified are :
#' \describe{
#'   \item{\code{V_E}}{The variance of environmental effects. May be a numeric vector of 
#' length 1 (V_E is the same for all traits) or a numeric vector of length n_trait. }
#'   \item{\code{V_R}}{The variance of the residual effects.May be a numeric vector of 
#' length 1 (V_R is the same for all traits) or a numeric vector of length n_trait. }
#' }
#' 
#' Environmental effects are drawn from a normal distribution such that \eqn{e ~ N(0, V_E)},
#' where \code{V_E} is the environmental variance. If \code{V_E} is not provided, it 
#' is 8 times the genotypic variance. Residual effects are drawn from a normal distribution 
#' such that \eqn{\epsilon ~ N(0, V_R)}, where \code{V_R} is the residual variance.
#' If \code{V_R} is not provided, it is calculated from \code{h2}, where
#' 
#' \eqn{V_R = n.rep * n.env * (\frac{Vg}{h2} - Vg)}.
#' 
#' @return 
#' An object of class \code{pop} with all information in the input \code{pop} object,
#' plus simulated phenotypes.
#' 
#' @examples 
#' 
#' # Load some historic data
#' data("s2_cap_genos")
#' data("s2_snp_info")
#' 
#' # Create a genome with genetic architecture
#' len <- tapply(s2_snp_info$cM_pos, s2_snp_info$chrom, max)
#' n_mar <- tapply(s2_snp_info$cM_pos, s2_snp_info$chrom, length)
#' map <- lapply(split(s2_snp_info, s2_snp_info$chrom), function(chr) structure(chr$cM_pos, names = chr$rs) )
#' 
#' genome <- sim_genome(len = len, n.mar = n_mar, map = map)
#' 
#' # Simulate a a trait with 15 QTL
#' qtl.model <- matrix(nrow = 15, ncol = 4)
#' 
#' genome <- sim_gen_model(genome, qtl.model, add.dist = "geometric", max.qtl = 15)
#' 
#' pop <- create_pop(genome = genome, geno = s2_cap_genos)
#' 
#' pop <- sim_phenoval(pop = pop, h2 = 0.5)
#' 
#' @import tidyr
#' @import dplyr
#' 
#' @export 
#' 
sim_phenoval <- function(pop, h2, n.env = 1, n.rep = 1, ...) {
  
  # Make sure pop inherits the class "pop"
  if (!inherits(pop, "pop"))
    stop("The input 'pop' must be of class 'pop'.")
  
  # Grab the genotypic values
  geno_val <- pop$geno_val
  
  # Does the pop object have genotypic values
  if (is.null(geno_val))
    stop("The 'pop' object must have the data.frame of genotypic values")
  
  # Number of traits
  n_trait <- ncol(geno_val) - 1
  
  # Capture the other arguments
  other.args <- list(...)
  
  # Other variances
  V_E <- other.args$V_E
  V_R <- other.args$V_R
  
  # If the h2 vector is longer than 1, the length must be the same as the number
  # of traits
  if (length(h2) > 1 & length(h2) != n_trait)
    stop("The length of h2, if not 1, must be the same as the number of traits.")
  
  # The same goes for V_E and V_R - this will pass if both or either are NULL.
  if (length(V_E) > 1 & length(V_E) != n_trait)
    stop("The length of V_E, if not 1, must be the same as the number of traits.")
  
  if (length(V_R) > 1 & length(V_R) != n_trait)
    stop("The length of h2, if not 1, must be the same as the number of traits.")
  
  
  # Heritability must be between 0 and 1
  if (!all(h2 >= 0, h2 <= 1))
    stop("Heritability must be between 0 and 1.")
  

  # If this list is not empty, make sure that the elements are correctly named
  if (length(other.args) != 0)
    # Warn if both are not provided
    if (!all(c("V_E", "V_R") %in% names(other.args)))
      warning("V_E and V_R might have been passed, but were not detected. Check 
              your arguments.")
  
  # Calculate genetic variance for each trait
  V_G <- geno_val %>% 
    summarize_at(vars(-ind), var)
  
  # Calculate environment variance if not provided
  if (is.null(V_E))
    V_E <- V_G * 8
  
  
  # Calculate residual variance if not provided
  if (is.null(V_R)) 
    V_R <- n.rep * n.env * ((V_G / h2) - V_G)

  
  # Number of individuals
  n_ind <- nind(pop)
  
  # List of variance components
  var_comp <- list(V_G = V_G, V_E = V_E, V_R = V_R)
  
  # Generate environment effects
  e <- lapply(V_E, sqrt) %>% 
    lapply(rnorm, n = n.env, mean = 0) %>%
    lapply(matrix, nrow = n_ind, ncol = n.env * n.rep, byrow = T)
  
  # Generate residual effects
  epsilon <- lapply(V_R, sqrt) %>% 
    lapply(rnorm, n = n.env * n.rep * n_ind, mean = 0) %>%
    lapply(matrix, nrow = n_ind, ncol = n.env * n.rep, byrow = T)
    
  g <- subset(pop$geno_val, select = -ind, drop = FALSE)
  
  # Apply over all of the list
  p <- mapply(g, e, epsilon, FUN = function(g1, e1, ep) {
    
    # Sum
    p <- structure(g1 + e1 + ep, dimnames = list(indnames(pop),
                                   paste( paste("env", seq(n.env), sep = ""), 
                                          rep(paste("rep", seq(n.rep), sep = ""), 
                                              each = n.env), sep = "_" )) )

    # return 
    return(p)
    
  }, SIMPLIFY = FALSE)
  
  
  # Tidy the phenotypes
  p_df <- mapply(names(p), p, FUN = function(tr, ph)
    data.frame(ind = row.names(ph), trait = tr, ph), SIMPLIFY = FALSE) %>%
    bind_rows()
  
  # Further tidying the phenotypes
  p_df1 <- p_df %>% 
    gather(obs, phenoval, -ind, -trait) %>%
    separate(col = obs, into = c("env", "rep"), sep = "_", remove = TRUE)
  
  # Calculate the mean phenotypic value
  mu_p <- p_df1 %>% 
    group_by(ind, trait) %>% 
    summarize(pheno_mean = mean(phenoval)) %>% 
    spread(trait, pheno_mean) %>%
    as.data.frame()
  
  # Add data to the pop object
  pheno_val <- list(
    var_comp = var_comp,
    pheno_obs = p_df1,
    pheno_mean = mu_p
  )
  
  pop[["pheno_val"]] <- pheno_val
  

  # Return the pop
  return(pop)
  
} # Close the function



#' Predict genotypic values using genomewide markers
#' 
#' 
#' @details 
#' The \code{training.pop} must have phenotypic values associated with each entry.
#' The mean phenotype is used as training data in the model. Genotypic data (excluding
#' QTL) are used to predict marker effects, which are then used to predict the 
#' genotypic value of the individuals in the \code{candidate.pop}.
#' 
#' When solving the mixed model, if \code{method = "RRBLUP"}, marker effects are
#' predicted by REML.
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
#' @import rrBLUP
#' 
#' @export
#' 
pred_geno_val <- function(genome, training.pop, candidate.pop, method = "RRBLUP") {
  
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
    
  }
  
  # Convert to data.frame
  pgv_df <- data.frame(ind = indnames(candidate.pop), pgv)
  
  # Add to the population and return
  candidate.pop$pred_val <- pgv_df
  
  return(candidate.pop)
  
} # Close the function
  
  
  
  
  
  



