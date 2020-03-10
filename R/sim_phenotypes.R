#' Simulate phenotypic values
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
#' @param return.eff Should the variance components and effects be returned?
#' @param ... Other arguments. See \code{Details}.
#' 
#' @details
#' 
#' Other arguments that can be specified are :
#' \describe{
#'   \item{\code{V_E}}{The variance of environmental effects. May be a scalar 
#'   (V_E is the same for all traits) or a numeric vector of length n_trait. }
#'   \item{\code{V_GE}}{The variance of genotype-environment interaction effects. May be a scalar 
#'   (V_GE is the same for all traits) or a numeric vector of length n_trait. }
#'   \item{\code{V_R}}{The variance of the residual effects. May be a scalar 
#'   (V_R is the same for all traits) or a numeric vector of length n_trait. }
#' }
#' 
#' Environmental effects are drawn from a normal distribution such that \eqn{e ~ N(0, V_E)},
#' where \code{V_E} is the environmental variance. If \code{V_E} is not provided, it 
#' is 8 times the genotypic variance. Genotype-environmental interaction effects are drawn
#' from a normal distribution such that \eqn{ge ~ N(0, V_GE)}, where \code{V_GE} is the genotype-environment 
#' interaction variance. If \code{V_GE} is not provided, it defaults to 2 times the genotypic variance. 
#' Residual effects are drawn from a normal distribution such that \eqn{\epsilon ~ N(0, V_R)}, where \code{V_R} is 
#' the residual variance. If \code{V_R} is not provided, it is calculated from \code{h2}, where 
#' \eqn{V_R = n.rep * n.env * (\frac{Vg}{h2} - Vg)}.
#' 
#' @return 
#' An object of class \code{pop} with all information in the input \code{pop} object,
#' plus simulated phenotypes.
#' 
#' @examples 
#' 
#' # Create a genome
#' n.mar  <- c(505, 505, 505)
#' len <- c(120, 130, 140)
#' 
#' genome <- sim_genome(len, n.mar)
#' 
#' # Simulate a a trait with 15 QTL
#' qtl.model <- matrix(nrow = 15, ncol = 4)
#' 
#' genome <- sim_gen_model(genome, qtl.model, add.dist = "geometric", max.qtl = 15)
#' 
#' pop <- sim_pop(genome = genome, n.ind = 200)
#' 
#' ## Default simulates environment and residual effects, but no GxE
#' pop1 <- sim_phenoval(pop = pop, h2 = 0.5, n.env = 2, n.rep = 2)
#' 
#' 
#' ## Simulate GxE
#' pop1 <- sim_phenoval(pop = pop, h2 = 0.5, n.env = 2, n.rep = 2, V_GE = var(pop$geno_val$trait1) * 3, V_E = NULL, V_R = NULL)
#' 
#' 
#' 
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
  n_trait <- sum(startsWith(x = names(geno_val), "trait"))
  
  # Capture the other arguments
  other.args <- list(...)
  
  # Other variances
  V_E <- other.args$V_E
  V_GE <- other.args$V_GE
  V_R <- other.args$V_R
  
  # The same goes for V_E and V_R - this will pass if both or either are NULL.
  if (length(V_E) > 1 & length(V_E) != n_trait)
    stop("The length of V_E, if not 1, must be the same as the number of traits.")
  
  # The same goes for V_E and V_R - this will pass if both or either are NULL.
  if (length(V_GE) > 1 & length(V_GE) != n_trait)
    stop("The length of V_GE, if not 1, must be the same as the number of traits.")
  
  if (length(V_R) > 1 & length(V_R) != n_trait)
    stop("The length of h2, if not 1, must be the same as the number of traits.")
  
  
  # If this list is not empty, make sure that the elements are correctly named
  if (length(other.args) != 0) {
    
    # Warn if all are not provided
    if (!all(c("V_E", "V_R", "V_GE") %in% names(other.args)))
      warning("V_E, V_GE, or V_R might have been passed, but were not detected. Check 
              your arguments.")
    
    # If none of the variance components are passed; you must provide the heritability
  } else {
    
    if (missing(h2)) stop("Variance components were not passed; you must provide a hertiability (h2).")
    
    # If the h2 vector is longer than 1, the length must be the same as the number
    # of traits
    if (length(h2) > 1 & length(h2) != n_trait)
      stop("The length of h2, if not 1, must be the same as the number of traits.")
    
    # Heritability must be between 0 and 1
    if (!all(h2 >= 0, h2 <= 1))
      stop("Heritability must be between 0 and 1.")
    
  }
  
  # Calculate genetic variance for each trait
  V_G <- sapply(geno_val[-1], var) 
    
  
  # Calculate environment variance if not provided
  if (is.null(V_E)) V_E <- V_G * 8
  
  ## V_GE is 0 if not passed
  if (is.null(V_GE)) V_GE <- 0
  
  # Calculate residual variance if not provided
  # if (is.null(V_R)) V_R <- n.rep * n.env * ( (V_G / h2) - V_G - (V_GE / n.env) )
  if (is.null(V_R)) V_R <- n.rep * n.env * ( (V_G / h2) - V_G )
  
  
  
  
  # Number of individuals
  n_ind <- nind(pop)
  
  # List of variance components
  var_comp <- list(V_G = V_G, V_E = V_E, V_GE = V_GE, V_R = V_R)
  
  # Generate environment effects
  t_eff <- lapply(X = V_E, function(varE) {
    sim <- rnorm(n = n.env, mean = 0, sd = sqrt(varE))
    structure(sim, names = paste("env", seq(n.env), sep = ""))
  })
    
  t_eff1 <- lapply(t_eff, matrix, nrow = n_ind, ncol = n.env * n.rep, byrow = TRUE)
  
  # Generate GxE effects
  gt_eff <- lapply(X = V_GE, function(varGE) {
    sim <- rnorm(n = n.env * n_ind, mean = 0, sd = sqrt(varGE))
    matrix(sim, nrow = n_ind, ncol = n.env * n.rep, byrow = FALSE)
  })
  
  # Generate residual effects
  epsilon <- lapply(X = V_R, function(varR) {
    sim <- rnorm(n = n.env * n.rep * n_ind, mean = 0, sd = sqrt(varR))
    matrix(sim, nrow = n_ind, ncol = n.env * n.rep, byrow = TRUE)
  })
  
  
  g <- subset(pop$geno_val, select = -ind, drop = FALSE)
  
  # Apply over all of the list
  p <- mapply(g, t_eff1, gt_eff, epsilon, FUN = function(g1, t1, gt1, ep) {
    
    # Reform the g1 matrix
    g1 <- matrix(g1, nrow = n_ind, ncol = n.env * n.rep)
    
    # Sum
    p <- g1 + t1 + gt1 + ep
    
    p1 <- structure(p, dimnames = list(indnames(pop), paste( paste("env", seq(n.env), sep = ""), rep(paste("rep", seq(n.rep), sep = ""), each = n.env), sep = "_" )) )
    
    p2 <- data.frame(ind = row.names(p1), p1, stringsAsFactors = FALSE, row.names = NULL)
    p2_reshape <- reshape(p2, direction = "long", varying = list(2:ncol(p2)), timevar = "env", v.names = "phenoval", idvar = "ind", new.row.names = NULL)
    p2_reshape$env <- names(p2)[-1][p2_reshape$env]
    p2_reshape$rep <- sapply(strsplit(p2_reshape$env, split = "_"), FUN = "[", 2)
    p2_reshape$env <- sapply(strsplit(p2_reshape$env, split = "_"), FUN = "[", 1)
 
    return(p2_reshape)
    
  }, SIMPLIFY = FALSE)
  
  
  # Tidy the phenotypes
  p_df <- mapply(names(p), p, FUN = function(tr, ph) data.frame(trait = tr, ph, stringsAsFactors = FALSE), SIMPLIFY = FALSE)
  p_df1 <- do.call("rbind", p_df)
  row.names(p_df1) <- NULL
    
  
  # Calculate the mean phenotypic value
  mu_p <- aggregate(formula = phenoval ~ ind + trait, data = p_df1, FUN = mean)
  mu_p1 <- reshape(data = mu_p, idvar = "ind", direction = "wide", v.names = "phenoval", timevar = "trait")
  names(mu_p1)[-1] <- names(p)
  
  pheno_val <- list(
      var_comp = var_comp,
      pheno_obs = p_df1,
      pheno_mean = mu_p1
  )
    

  pop[["pheno_val"]] <- pheno_val
  
  
  # Return the pop
  return(pop)
  
} # Close the function




#' Simulate an experimental field trial
#' 
#' @description 
#' Simulates a trial in which phenotypic data on a quantitative trait is gathered.
#' 
#' @param pop An object of class \code{pop}, the individuals in which are to be
#' phenotyped.
#' @param h2 The heritability of the trait or traits. May be a \code{double} 
#' (heritability is the same for all traits) or a numeric vector with length equal
#' to the number of traits.
#' @param n.env The number of environments in which to phenotype.
#' @param n.rep The number of replicates of each experimental (i.e. non-check)
#' individual in each environment.
#' @param check.pop A \code{pop} object with individuals to be used as repeated 
#' checks (e.g. for an augmented design). If missing (default), then no checks
#' are included.
#' @param check.rep A scalar specifying the number of times each check is replicated
#' in each environment.
#' @param ... Other arguments. See \code{Details}.
#' 
#' @return 
#' An object of class \code{pop} with all information in the input \code{pop} object,
#' plus simulated phenotypes.
#' 
#' @details
#' 
#' Other arguments that can be specified are :
#' \describe{
#'   \item{\code{V_E}}{The variance of environmental effects. May be a \code{double} 
#'   (V_E is the same for all traits) or a numeric vector with length equal
#'   to the number of traits. }
#'   \item{\code{V_R}}{The variance of the residual effects. May be a \code{double} 
#'   (V_R is the same for all traits) or a numeric vector with length equal
#'   to the number of traits. }
#' }
#' 
#' @examples 
#' 
#' # Load some historic data
#' data("s2_cap_genos")
#' data("s2_snp_info")
#' 
#' map <- s2_snp_info %>% 
#'   select(-alleles) %>% 
#'   data.frame(row.names = .$rs, stringsAsFactors = FALSE) %>% 
#'   select(-rs) %>% 
#'   table_to_map()
#' 
#' # Create a genome with genetic architecture
#' genome <- sim_genome(map = map)
#' 
#' # Simulate a a trait with 15 QTL
#' qtl.model <- matrix(nrow = 15, ncol = 4)
#' 
#' genome <- sim_gen_model(genome, qtl.model, add.dist = "geometric", max.qtl = 15)
#' 
#' pop <- create_pop(genome = genome, geno = s2_cap_genos)
#' 
#' # Sample the population for checks
#' set.seed(1035)
#' check_pop <- pop %>%
#'   subset_pop(individual = sample(indnames(.), 5))
#'   
#' # Remove the checks from the original pop
#' exp_pop <- pop %>% 
#'   subset_pop(setdiff(indnames(.), indnames(check_pop)))
#' 
#' # Simulate a trial
#' pop_pheno <- sim_trial(pop = exp_pop, h2 = 0.5, n.env = 3, n.rep = 1,
#'   check.pop = check_pop, check.rep = 3)
#'   
#' @import dplyr
#' @import tidyr
#' 
#' 
sim_trial <- function(pop, h2, n.env = 1, n.rep = 1, check.pop, check.rep, ...) {
  
  # Make sure pop inherits the class "pop"
  if (!inherits(pop, "pop"))
    stop("The input 'pop' must be of class 'pop'.")
  
  # n.rep and n.env must be integers
  n.rep <- as.integer(n.rep)
  n.env <- as.integer(n.env)
  
  # Check the checks, if present
  if (!missing(check.pop)) {
    if (!inherits(check.pop, "pop"))
      stop("The input 'pop' must be of class 'pop'.")
    
    if (missing(check.rep))
      stop("If 'check.pop' is passed, so too must 'check.rep'")
    
  }
  
  # Check the checks, if present
  if (!missing(check.rep)) {
    if (missing(check.pop))
      stop("If 'check.rep' is passed, so too must 'check.pop'")
    
    check.rep <- as.integer(check.rep)
    
  }
  
  # Capture the other arguments
  other.args <- list(...)
  
  # Other variances
  V_E <- other.args$V_E
  V_R <- other.args$V_R
  
  # Simulate phenotypes for the experimental individuals
  exp_pheno <- sim_phenoval(pop = pop, h2 = h2, n.env = n.env, n.rep = n.rep, 
                            V_E = V_E, V_R = V_R)
  
  # Extract the variance components from the exp_pheno for use with the checks
  V_R <- as.numeric(exp_pheno$pheno_val$var_comp$V_R)
  V_E <- 0
  
  # Note the environment variance is 0 so as to add the environmental effects
  # from the previous simulation
  check_pheno <- sim_phenoval(pop = check_pop, h2 = h2, n.env = n.env, n.rep = check.rep,
                              V_E = V_E, V_R = V_R)
  
  # Extract the phenotypic observations and add the environmental effects
  pheno_obs <- check_pheno$pheno_val$pheno_obs
  
  env_eff <- exp_pheno$pheno_val$effects$env %>% 
    lapply(function(eff) data.frame(env = names(eff), effect = eff, row.names = NULL))
  
  env_eff <- lapply(names(env_eff), function(trait) mutate(env_eff[[trait]], trait = trait)) %>%
    bind_rows()
  
  # Merge and add env eff
  pheno_obs1 <- merge(x = pheno_obs, y = env_eff) %>% 
    mutate(phenoval = phenoval + effect) %>% 
    select(ind, trait, env, rep, phenoval) %>%
    arrange(trait, env, rep)
  
  # Calculate the means
  pheno_mean <- pheno_obs1 %>% 
    group_by(ind, trait) %>% 
    summarize(phenomean = mean(phenoval)) %>% 
    spread(trait, phenomean) %>% 
    as.data.frame()
  
  # Adjust the check phenos
  check_pheno$pheno_val <- list(var_comp = exp_pheno$pheno_val$var_comp, 
                                effects = exp_pheno$pheno_val$effects, 
                                pheno_obs = pheno_obs1, 
                                pheno_mean = pheno_mean)
  
  # Merge populations
  final_pop <- combine_pop(list(exp_pheno, check_pheno))
  
  # Add variance components
  final_pop$pheno_val <- list(var_comp = exp_pheno$pheno_val$var_comp, 
                              effects = exp_pheno$pheno_val$effects, 
                              pheno_obs = final_pop$pheno_val$pheno_obs, 
                              pheno_mean = final_pop$pheno_val$pheno_mean)
  
  # Return
  return(final_pop)
  
}

  
  
  
  
  
  
  
  
  
  
  