#' Simulate a phenotypic trial
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
#' @export 
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
  
  
  
  
  
  
  
  
  
  
  
  
  