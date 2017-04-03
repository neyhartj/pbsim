#' Simulate phenotypes given QTL genotypes
#' 
#' @description 
#' Simulates phenotypic observations given the genotypic value of individuals and the 
#' heritability of a quantitative trait.
#' 
#' @param genome A \code{genome} object.
#' @param geno Genotype data on a population to phenotype. Must be a matrix of dimensions
#' \code{n.ind} x \code{n.mar}, the elements of which must be z {0, 1, 2}.
#' @param h2 The heritability of the trait.
#' @param n.env The number of environments in which to phenotype.
#' @param n.rep The number of replicates of each individual in each environment.
#' 
#' @details
#' The genotypic value of individuals is calculcated as the sum of the QTL effects
#' carried by each individual. The genetic variance is calculated as the variance
#' of these genotypic values (\eqn{V_G = var(g)}).
#' 
#' Environmental effects are draw from a normal distribution such that \eqn{e ~ N(0, V_E)},
#' where \code{V_E} is the environmental variance and is 8 times the genotypic variance.
#' The residual variance (\code{V_R}) is calculated from \code{h2}, where
#' 
#' \eqn{V_R = n.rep * n.env * (\frac{Vg}{h2} - Vg)}.
#' 
#' @return 
#' A list with \code{geno_val} (the genotypic values), \code{pheno_val} (the phenotypic
#' observations of each rep in each environment), and \code{pheno_mean} (the mean phenotypic
#' value for each individual).
#' 
#' @export 
#' 
sim_pheno <- function(genome, geno, h2, n.env = 1, n.rep = 1) {
  
  # Check class
  stopifnot(inherits(x = genome, what = "genome"))
  
  # Make sure the genome has a genetic model
  if (is.null(genome$gen_model))
    stop("The 'genome' must have a genetic model. Use 'sim_gen_model' to create one.")
  
  # Verify the genotypes
  # How many individuals?
  n.ind <- nrow(geno)

  
  # Are the genos coded correctly?
  if (!all(unlist(geno) %in% c(0, 1, 2)))
    stop("The input 'geno' must be encoded in z {0, 1, 2}.")
  
  # Grab the genetic model
  gen_model <- genome$gen_model
  
  
  # Split the genos by chromsome
  split_ind <- rep(seq_along(genome$map), genome$n.mar)
  
  # Empty vector for genotypic values
  geno_val <- vector("list", length(genome$map))
  
  # Iterate over chromosomes
  for (i in seq_along(genome$map)) {
    
    # Split the genos
    geno_split <- geno[,split_ind == i]
    
    # Subset the model for that chromosome
    gen_model_chr <- subset(gen_model, chr == i)
    
    # Find the index of QTL in the chromosome
    qtl_ind <- gen_model_chr$qtl.ind
    
    # Find the effects of the QTL
    add_eff <- gen_model_chr$add.eff
    dom_eff <- gen_model_chr$dom.eff
    
    
    # If no QTL exist on the chromosome, add a 0 to the geno_val list and continue
    if (length(qtl_ind) == 0) {
      geno_val[[i]] <- rep(0, n.ind)
      
    } else {
      
      # Get the QTL genotypes and subtract 1
      qtl_geno <- geno_split[,qtl_ind, drop = FALSE] - 1
      
      # Empty matrix
      qtl_geno_val <- matrix(NA, nrow = nrow(qtl_geno), ncol = ncol(qtl_geno),
                             dimnames = dimnames(qtl_geno))
      
      # Iterate over QTL
      for (q in seq_along(qtl_ind)) {
        
        # Get the genotypes for that qtl
        q_geno <- qtl_geno[,q]
        
        # Get the genotypic values at that qtl
        qtl_geno_val[,q] <- ifelse(q_geno == 0, q_geno * dom_eff[q], q_geno * add_eff[q])
        
      }
      
      # Sum to get genotypic value
      geno_val[[i]] <- rowSums(qtl_geno_val)
      
    }
    
  }
  
  # Column bind and sum
  geno_val <- as.matrix(rowSums(do.call("cbind", geno_val)))
  g <- matrix(geno_val, nrow = n.ind, ncol = n.env * n.rep, byrow = F,
              dimnames = dimnames(geno_val))
  
  # Calculate genetic variance
  Vg <- var(geno_val)
  
  # Calculate environment variance
  Ve <- Vg * 8
  
  # Calculate residual variance
  Vr <- n.rep * n.env * ((Vg / h2) - Vg)
  
  
  # Generate environment effects
  e <- matrix(data = rnorm(n = n.env, mean = 0, sd = sqrt(Ve)), nrow = n.ind, 
              ncol = n.env * n.rep, byrow = T)
  
  # Generate residual effects
  epsilon <- matrix(data = rnorm(n = n.env * n.rep * n.ind, mean = 0, sd = sqrt(Ve)), 
                    nrow = n.ind, ncol = n.env * n.rep, byrow = T)
  
  # Sum
  p <- g + e + epsilon
  
  # Calculate the mean phenotypic value
  mu_p <- as.matrix(rowMeans(p))
  
  # Return a list of genotypic values, phenotypic values, and mean phenotypic values
  list(
    geno_val = geno_val,
    pheno_val = p,
    pheno_mean = mu_p
  )
  
} # Close the function

  