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
#' @param ... Other arguments. See \code{Details}.
#' 
#' @details
#' Other arguments that can be specified are :
#' \itemize{
#'   \item{\code{V_E}: The variance of environmental effects.}
#'   \item{\code{V_R}: The variance of the residual effects.}
#'   }
#' 
#' The genotypic value of individuals is calculcated as the sum of the QTL effects
#' carried by each individual. The genetic variance is calculated as the variance
#' of these genotypic values (\eqn{V_G = var(g)}).
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
#' A list with \code{geno_val} (the genotypic values), \code{pheno_val} (the phenotypic
#' observations of each rep in each environment), and \code{pheno_mean} (the mean phenotypic
#' value for each individual).
#' 
#' @examples 
#' # Simulate a genome
#' n.mar  <- c(505, 505, 505)
#' len <- c(120, 130, 140)
#' 
#' genome <- sim_genome(len, n.mar)
#' 
#' # Simulate a quantitative trait influenced by 50 QTL
#' genome <- sim_gen_model(genome = genome, qtl.model = matrix(NA, 50, 4), add.dist = "geometric")
#' 
#' # Simulate a pedigree
#' ped <- sim_pedigree(n.ind = 50, n.bcgen = 0, n.selfgen = 2)
#' 
#' # Simulate the founder genotypes
#' founder_geno <- sim_founders(genome)
#' 
#' # Simulate phenotypes for a trait with a heritability of 0.5
#' phenos <- sim_pheno(genome, founder_geno, h2 = 0.5)
#' 
#' @import tidyr
#' @import dplyr
#' 
#' @export 
#' 
sim_pheno <- function(genome, geno, h2, n.env = 1, n.rep = 1, ...) {
  
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
  
  # Capture the other arguments
  other.args <- list(...)
  # If this list is not empty, make sure that the elements are correctly named
  if (length(other.args) != 0)
    # Warn if both are not provided
    if (!all(c("V_E", "V_R") %in% names(other.args)))
      warning("Either V_E or V_R were not passed. This may be a mistake.")
  
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
  
  # Cbind the list
  geno_val <- as.matrix(rowSums(do.call("cbind", geno_val)))
  # Add column names
  colnames(geno_val) <- "geno_val"
  
  g <- matrix(geno_val, nrow = n.ind, ncol = n.env * n.rep, byrow = F,
              dimnames = list(row.names(geno_val), NULL))
  
  # Calculate genetic variance
  Vg <- as.numeric(var(geno_val))
  
  # Calculate environment variance if not provided
  if (is.null(other.args$V_E)) {
    Ve <- Vg * 8
    
  } else{
    Ve <- other.args$V_E
    
  }
  
  # Calculate residual variance if not provided
  if (is.null(other.args$V_R)) {
    Vr <- n.rep * n.env * ((Vg / h2) - Vg)
    
  } else{
    Vr <- other.args$V_R
    
  }
  
  # List of variance components
  var_comp <- list(V_G = Vg, V_E = Ve, V_R = Vr)
  
  # Generate environment effects
  e <- matrix(data = rnorm(n = n.env, mean = 0, sd = sqrt(Ve)), nrow = n.ind, 
              ncol = n.env * n.rep, byrow = T)
  
  # Generate residual effects
  epsilon <- matrix(data = rnorm(n = n.env * n.rep * n.ind, mean = 0, sd = sqrt(Ve)), 
                    nrow = n.ind, ncol = n.env * n.rep, byrow = T)
  
  # Sum
  p <- g + e + epsilon
  
  # Add environment/rep names
  colnames(p) <- paste( 
    paste("env", seq(n.env), sep = ""), 
    rep(paste("rep", seq(n.rep), sep = ""), each = n.env), sep = "_" )
    
  # Tidy the genotypes
  g_df <- data.frame(entry = row.names(geno_val), geno_val, row.names = NULL)
  
  # Tidy the phenotypes
  p_df <- data.frame(entry = row.names(p), p) %>% 
    gather(key = obs, value = pheno_val, -entry) %>%
    separate(col = obs, into = c("env", "rep"), sep = "_", remove = TRUE)
  
  # Calculate the mean phenotypic value
  mu_p <- p_df %>% 
    group_by(entry) %>% 
    summarize(pheno_mean = mean(pheno_val)) %>% 
    as.data.frame()
  
  # Return a list of variance components, genotypic values, phenotypic values, 
  # and mean phenotypic values as data.frames
  list(
    var_comp = var_comp,
    geno_val = g_df,
    pheno_val = p_df,
    pheno_mean = mu_p
  )
  
} # Close the function

  