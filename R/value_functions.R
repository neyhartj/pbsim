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
#' @import dplyr
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
  data.frame(ind = row.names(geno_val1), geno_val1, row.names = NULL, stringsAsFactors = FALSE) %>%
    # Sort on individual name
    arrange(ind)
  
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
#' @param return.eff Should the variance components and effects be returned?
#' @param ... Other arguments. See \code{Details}.
#' 
#' @details
#' 
#' Other arguments that can be specified are :
#' \describe{
#'   \item{\code{V_E}}{The variance of environmental effects. May be a scalar 
#'   (V_E is the same for all traits) or a numeric vector of length n_trait. }
#'   \item{\code{V_R}}{The variance of the residual effects.May be a scalar 
#'   (V_R is the same for all traits) or a numeric vector of length n_trait. }
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
sim_phenoval <- function(pop, h2, n.env = 1, n.rep = 1, return.eff = TRUE, ...) {
  
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
  t_eff <- lapply(V_E, sqrt) %>% 
    lapply(rnorm, n = n.env, mean = 0) %>%
    lapply(structure, names = paste("env", seq(n.env), sep = ""))
  
  t <- t_eff %>%
    lapply(matrix, nrow = n_ind, ncol = n.env * n.rep, byrow = T)
  
  # Generate residual effects
  epsilon <- lapply(V_R, sqrt) %>% 
    lapply(rnorm, n = n.env * n.rep * n_ind, mean = 0) %>%
    lapply(matrix, nrow = n_ind, ncol = n.env * n.rep, byrow = T)
    
  g <- subset(pop$geno_val, select = -ind, drop = FALSE)
  
  # Apply over all of the list
  p <- mapply(g, t, epsilon, FUN = function(g1, t1, ep) {
    
    # Reform the g1 matrix
    g1 <- matrix(g1, nrow = n_ind, ncol = n.env * n.rep)
    
    # Sum
    p <- structure(g1 + t1 + ep, dimnames = list(indnames(pop),
                                   paste( paste("env", seq(n.env), sep = ""), 
                                          rep(paste("rep", seq(n.rep), sep = ""), 
                                              each = n.env), sep = "_" )) )

    # return 
    return(p)
    
  }, SIMPLIFY = FALSE)
  
  # List of environmental effects
  effects <- list(env = t_eff)
  
  
  # Tidy the phenotypes
  p_df <- mapply(names(p), p, FUN = function(tr, ph)
    data.frame(ind = row.names(ph), trait = tr, ph, stringsAsFactors = FALSE), SIMPLIFY = FALSE) %>%
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
  if (return.eff) {
    pheno_val <- list(
      var_comp = var_comp,
      effects = effects,
      pheno_obs = p_df1,
      pheno_mean = mu_p
    )
    
  } else {
    pheno_val <- list(
      pheno_obs = p_df1,
      pheno_mean = mu_p
    )
    
  }
  
  pop[["pheno_val"]] <- pheno_val
  

  # Return the pop
  return(pop)
  
} # Close the function


  
  
  
  
  
  



