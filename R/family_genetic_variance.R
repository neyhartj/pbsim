#' Calculate the expected genetic variance in simulated families
#' 
#' 
#' @description 
#' Calculates the expected genetic variance of a cross, assuming complete selfing.
#' 
#' @param genome An object of class \code{genome}.
#' @param pedigree A \code{pedigree} detailing the scheme to develop the family.
#' Use \code{\link{sim_pedigree}} to generate.
#' @param founder.pop An object of class \code{pop} with the geno information for
#' the parents. Additional individuals can be present in \code{parent_pop}. They
#' will be filtered according to the parents in the \code{crossing.block}.
#' @param crossing.block A crossing block detailing the crosses to make. Must be a
#' \code{data.frame} with 2 columns: the first gives the name of parent 1, and the 
#' second gives the name of parent 2. See \code{\link{sim_crossing.block}}.
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
#' # Simulate the genotypes for 8 founders
#' founder.pop <- sim_founders(genome = genome, n.str = 8)
#' 
#' # Generate a crossing block with 5 crosses
#' cb <- sim_crossing_block(parents = indnames(founder.pop), n.crosses = 5)
#' 
#' # Create a pedigree with 100 individuals selfed to the F_3 generation
#' ped <- sim_pedigree(n.ind = 100, n.selfgen = 2)
#' 
#' calc_exp_genvar(genome = genome, pedigree = ped, founder.pop = founder.pop, 
#'                 crossing.block = cb)
#'                 
#' 
#' ## If two traits are present, the genetic correlation is calculated
#' # Simulate two quantitative traits influenced by 50 pleiotropic QTL
#' qtl.model <- replicate(2, matrix(NA, 50, 4), simplify = FALSE)
#' genome <- sim_multi_gen_model(genome = genome, qtl.model = qtl.model, corr = 0.99, 
#'                               prob.corr = cbind(0, 1), add.dist = "normal")
#' 
#' # Simulate the genotypes for 8 founders
#' founder.pop <- sim_founders(genome = genome, n.str = 8)
#' 
#' calc_exp_genvar(genome = genome, pedigree = ped, founder.pop = founder.pop, 
#'                 crossing.block = cb)
#' 
#' 
#' 
#' @importFrom qtl mf.h
#' @importFrom simcross check_pedigree
#' @importFrom Matrix .bdiag
#' @importFrom tidyr crossing
#' @importFrom purrr pmap_dbl
#' 
#' @export
#' 
calc_exp_genvar <- function(genome, pedigree, founder.pop, crossing.block) {
  
  # Error handling
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Check the pedigree
  if (!check_pedigree(pedigree, ignore_sex = TRUE))
    stop("The pedigree is not formatted correctly.")
  
  # Check the crossing block
  if (ncol(crossing.block) != 2) {
    stop("The crossing block should have two columns.")
  } else {
    crossing.block <- as.data.frame(crossing.block)
  }
  
  # founder.pop needs to be a pop object
  if (!inherits(founder.pop, "pop"))
    stop("The input 'founder.pop' must be of class 'pop'")
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = founder.pop$geno))
    stop("The geno did not pass. See warning for reason.")
  
  ## How many traits
  n_traits <- length(genome$gen_model)
  
  # If it is more than 2, error out
  stopifnot(n_traits <= 2)
  
  
  ## Calculate the expected genetic variance
  
  
  ## What are the expected allele frequencies in the population?
  ## Is there any backcrossing?
  mom_ped <- pedigree[pedigree$mom == 1,]
  dad_ped <- pedigree[pedigree$mom == 2,]
  
  mom_dist_gen <- length(unique(mom_ped$gen))
  dad_dist_gen <- length(unique(dad_ped$gen))
  
  max_bc_gen <- pmax(mom_dist_gen, dad_dist_gen) - 1
  
  # The expected frequency of the minor allele is 0.5 ^ n_bc_gen + 1
  exp_q <- 0.5^(max_bc_gen + 1)
  exp_p <- 1 - exp_q
  
  # Get the QTL information - drop unused levels
  qtl_info <- pull_qtl(genome, unique = FALSE)
  # Filter out QTL with no additive effect
  qtl_info <- droplevels(qtl_info[qtl_info$add_eff != 0,,drop = FALSE])
  # Split by trait
  qtl_info_split <- split(qtl_info, qtl_info$trait)
  
  
  ## Iterate over traits
  qtl_covariance <- lapply(X = qtl_info_split, FUN = function(trait_qtl) {
    
    # Get the map and genotypes of only the QTL
    qtl_geno <- pull_genotype(genome = genome, geno = founder.pop$geno, loci = trait_qtl$qtl_name) - 1
    
    ## Calculate the expected genetic variance and covariance of QTL
    trait_qtl_split <- split(trait_qtl, trait_qtl$chr)
    
    
    ## Calculate the expected covariance between QTL
    qtl_covar <- lapply(X = trait_qtl_split, FUN = function(qtl_chr) {
      d <- as.matrix(dist(qtl_chr$pos))
      
      # Calculate pairwise D (see Zhong and Jannink, 2007)
      # First convert cM to recombination fraction
      c <- qtl:::mf.h(d)
      D <- ((1 - (2 * c)) / (1 + (2 * c)))
      # # The diagonals are 0
      # diag(D) <- 0
      
      # Calculate the pairwise product of all QTL effects
      qtl_crossprod <- tcrossprod(qtl_chr$add_eff)
      dimnames(qtl_crossprod) <- list(qtl_chr$qtl_name, qtl_chr$qtl_name)
      
      # The covariance is the QTL effect product multiplied by the expected D
      qtl_crossprod * D
      
    })
    
    # Combine into a block diagonal, since the covariance between QTL on different chromosomes
    # is expected to be 0
    Cov <- as.matrix(.bdiag(qtl_covar))
    dimnames(Cov) <- list(trait_qtl$qtl_name, trait_qtl$qtl_name)
    # Return
    return(Cov)
    
  })
  
  ## Calculate the genetic covariance between QTL for different traits
  # Split by chromosome
  qtl_chr_split <- split(qtl_info, qtl_info$chr)
  
  
  if (n_traits > 1) {
    
    ## Iterate over chromosomes
    qtl_trait_covariance <- lapply(X = qtl_chr_split, FUN = function(chr_qtl) {
      
      # Split by trait
      trait_split <- split(chr_qtl, chr_qtl$trait)
      
      ## QTL names for each trait
      qtl_names <- lapply(X = trait_split, FUN = "[[", "qtl_name")
      qtl_pos <- lapply(X = trait_split, FUN = "[[", "pos")
      qtl_eff <- lapply(X = trait_split, FUN = function(q) as.matrix(q$add_eff))
      
      ## Calculate the pairwise distance
      d <- abs(outer(X = qtl_pos[[1]], Y = qtl_pos[[2]], FUN = `-`))
      # Calculate pairwise D (see Zhong and Jannink, 2007)
      # First convert cM to recombination fraction
      c <- qtl:::mf.h(d)
      D <- ((1 - (2 * c)) / (1 + (2 * c)))
      
      # Product of QTL effects
      qtl_crossprod <- tcrossprod(qtl_eff[[1]], qtl_eff[[2]])
      dimnames(qtl_crossprod) <- qtl_names
      
      # The covariance is the QTL effect product multiplied by the expected D
      qtl_crossprod * D
      
    })
    
    # Combine into a block diagonal, since the covariance between QTL on different chromosomes
    # is expected to be 0
    qtl_trait_covariance <- as.matrix(.bdiag(qtl_trait_covariance))
    dimnames(qtl_trait_covariance) <- list(qtl_info_split[[1]]$qtl_name, qtl_info_split[[2]]$qtl_name)
    
  } else {
    qtl_trait_covariance <- NULL
    
  }
  
  
  
  ## Now we iterate over the parent pairs to determine the QTL that are segregating
  
  # Replicate the crossing block 
  
  ## Add columns to the crossing.block for exp mu and exp varG
  crossing_block <- crossing(crossing.block, trait = paste0("trait", seq(length(genome$gen_model))))
  exp_mu <- list()
  exp_varG <- list()
  exp_corG <- list()
  
  
  
  ## Pull out the qtl genotypes for each trait
  qtl_names <- lapply(X = qtl_info_split, FUN = "[[", "qtl_name")
  qtl_geno <- lapply(X = qtl_names, function(q) pull_genotype(genome = genome, geno = founder.pop$geno, loci = q) - 1)
  
  
  # Iterate over the crossing block
  for (j in seq(nrow(crossing.block))) {
    
    pars <- as.character(crossing.block[j,1:2])
    
    ## Get a list of the polymorphic QTL
    poly_qtl_list <- lapply(X = qtl_geno, FUN = function(tr_qtl) {
      
      # Subset the parents
      par_qtl_geno <- tr_qtl[pars,,drop = FALSE]
      qtl_means <- colMeans(par_qtl_geno)
      par1_qtl <- par_qtl_geno[1,]
      
      t(ifelse(qtl_means == 0, par1_qtl, 0))
      
      # poly_qtl <- names(qtl_means)[qtl_means == 0]
      # # Get the marker states for parent 1 - use this to determine phase
      # par_qtl_geno[1, poly_qtl, drop = F]
      
    })
    
    
    # Iterate over the traits and calculate individual genetic variance
    trait_var <- pmap_dbl(list(poly_qtl_list, qtl_covariance), ~sum(crossprod(.x) * .y))
    
    
    # trait_var <- mapply(poly_qtl_list, qtl_covariance, FUN = function(poly_mat, qtl_cov) {
    #   
    #   poly_qtl <- colnames(poly_mat)
    #   exp_covar1 <- qtl_cov[poly_qtl, poly_qtl] * crossprod(poly_mat)
    #   
    #   # The expected variance is the sum of the variances at the polymorphic QTL, plus 2 times
    #   # the expected covariance between all polymorphic QTL
    #   # One can simply sum up the elements in the covariance matrix
    #   sum(exp_covar1)
    #   
    # })
    
    if (!is.null(qtl_trait_covariance)) {
      
      
      ## Calculate the expected covariance
      trait_cov <- sum(qtl_trait_covariance * crossprod(poly_qtl_list[[1]], poly_qtl_list[[2]]))
      
      # The expected correlation is calculated using the expected sd and expected cov
      exp_corG_j <- trait_cov / prod(sqrt(trait_var)) 
      exp_corG[[j]] <- rep(exp_corG_j, 2)
      
    }
    
    # The expected mu is simply the mean of the genotypic values of the two parents
    exp_mu_j <- colMeans(founder.pop$geno_val[founder.pop$geno_val$ind %in% pars,-1,drop = F])
    
    ## Add to the lists
    exp_mu[[j]] <- exp_mu_j
    exp_varG[[j]] <- trait_var
    
  }
  
  ## Add the variances and means to the crossing block
  crossing_block$exp_mu <- unlist(exp_mu)
  crossing_block$exp_varG <- unlist(exp_varG)  
  crossing_block$exp_corG <- unlist(exp_corG)
  
  # Return the crossing block
  return(crossing_block)
  
}


#' Predict the genetic variance in prospective crosses
#' 
#' @description 
#' Uses the expected genetic variance formula and marker effects to predict the
#' genetic variance and correlation in potential crosses.
#' 
#' @param genome An object of class \code{genome}.
#' @param pedigree A \code{pedigree} detailing the scheme to develop the family.
#' Use \code{\link{sim_pedigree}} to generate.
#' @param founder.pop An object of class \code{pop} with the geno information for
#' the parents. Additional individuals can be present in \code{parent_pop}. They
#' will be filtered according to the parents in the \code{crossing.block}.
#' @param training.pop An object of class \code{pop} with the elements \code{geno} and 
#' \code{pheno_val}. This is used as the training population.
#' @param crossing.block A crossing block detailing the crosses to make. Must be a
#' \code{data.frame} with 2 columns: the first gives the name of parent 1, and the 
#' second gives the name of parent 2. See \code{\link{sim_crossing.block}}.
#' @param method The statistical method to predict marker effects. If \code{"RRBLUP"}, the
#' \code{\link[qtl]{mixed.solve}} function is used. Otherwise, the \code{\link[BGLR]{BGLR}}
#' function is used.
#' @param n.iter,burn.in,thin Number of iterations, number of burn-ins, and thinning, respectively. See 
#' \code{\link[BGLR]{BGLR}}.
#' @param save.at See \code{\link[BGLR]{BGLR}}.
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
#' # Simulate the genotypes for 8 founders
#' founder.pop <- sim_founders(genome = genome, n.str = 8)
#' training.pop <- sim_phenoval(founder.pop, h2 = 0.8)
#' 
#' # Generate a crossing block with 5 crosses
#' cb <- sim_crossing_block(parents = indnames(founder.pop), n.crosses = 5)
#' 
#' # Create a pedigree with 100 individuals selfed to the F_3 generation
#' ped <- sim_pedigree(n.ind = 100, n.selfgen = 2)
#' 
#' pred_genvar(genome = genome, pedigree = ped, training.pop = training.pop, 
#'             founder.pop = founder.pop, crossing.block = cb)
#'                 
#' 
#' ## If two traits are present, the genetic correlation is calculated
#' # Simulate two quantitative traits influenced by 50 pleiotropic QTL
#' qtl.model <- replicate(2, matrix(NA, 50, 4), simplify = FALSE)
#' genome <- sim_multi_gen_model(genome = genome, qtl.model = qtl.model, corr = 0.99, 
#'                               prob.corr = cbind(0, 1), add.dist = "normal")
#' 
#' # Simulate the genotypes for 8 founders
#' founder.pop <- sim_founders(genome = genome, n.str = 8)
#' training.pop <- sim_phenoval(founder.pop, h2 = 0.8)
#' 
#' pred_genvar(genome = genome, pedigree = ped, training.pop = training.pop, 
#'             founder.pop = founder.pop, crossing.block = cb)
#' 
#' @importFrom simcross check_pedigree
#' 
#' @export
#' 
pred_genvar <- function(genome, pedigree, training.pop, founder.pop, crossing.block, 
                        method = c("RRBLUP", "BRR", "BayesA", "BL", "BayesB", "BayesC"), 
                        n.iter = 1500, burn.in = 500, thin = 5, save.at = "") {
  
  # Error handling
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Check the pedigree
  if (!check_pedigree(pedigree, ignore_sex = TRUE))
    stop("The pedigree is not formatted correctly.")
  
  # Check the crossing block
  if (ncol(crossing.block) != 2) {
    stop("The crossing block should have two columns.")
  } else {
    crossing.block <- as.data.frame(crossing.block)
  }
  
  # founder.pop needs to be a pop object
  if (!inherits(founder.pop, "pop"))
    stop("The input 'founder.pop' must be of class 'pop'")
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = founder.pop$geno))
    stop("The geno did not pass. See warning for reason.")
  
  # Check the populations
  if (!inherits(training.pop, "pop"))
    stop("The input 'training.pop' must be of class 'pop'.")
  
  # Make sure the training population has phenotypes
  if (is.null(training.pop$pheno_val))
    stop("The 'training.pop' must have phenotypic values.")
  
  n_traits <- length(genome$gen_model)
  
  # Check the method
  method <- match.arg(method)
  
  # Predict marker effects - only if the TP does not have them
  if (is.null(training.pop$mar_eff)) {
    
    marker_eff <- pred_mar_eff(genome = genome, training.pop = training.pop, method = method, n.iter = n.iter,
                               burn.in = burn.in, thin = thin, save.at = save.at)
    
  } else {
    marker_eff <- training.pop
    
  }
  
  ## Predict genotypic values in the founder population - this will use the marker effects in the tp
  founder_pop1 <- pred_geno_val(genome = genome, training.pop = marker_eff, candidate.pop = founder.pop)
  # Predict marker effects
  marker_eff <- marker_eff$mar_eff
  
  ## Find the positions of these markers
  marker_pos <- find_markerpos(genome = genome, marker = marker_eff$marker)
  marker_pos$add_eff <- NA
  marker_pos$dom_eff <- 0
  marker_pos$qtl_name <- row.names(marker_pos)
  marker_pos$qtl1_pair <- row.names(marker_pos)
  
  # Duplicate by the number of traits
  marker_pos_list <- replicate(n = n_traits, marker_pos, simplify = FALSE)
  marker_pos_list[[1]]$qtl1_pair <- NA
  
  # Add effects
  for (i in seq_len(n_traits)) {
    marker_pos_list[[i]]$add_eff <- marker_eff[,-1][[i]]
  }
  
  ## Create a new genome with markers as QTL
  genome_use <- genome
  
  # Add to the genome
  genome_use$gen_model <- marker_pos_list
  
  
  ## Predict
  predicted_genvar <- calc_exp_genvar(genome = genome_use, pedigree = pedigree, founder.pop = founder.pop, 
                                      crossing.block = crossing.block)
  
  # PGVs
  pgvs <- founder_pop1$pred_val
  
  # Replace the expected mu with the predicted mu
  for (i in seq_len(nrow(predicted_genvar))) {
    predicted_genvar$exp_mu[i] <- mean(pgvs[pgvs$ind %in% predicted_genvar[i,1:2,drop = T],predicted_genvar$trait[i]])
  }
  
  
  if (n_traits == 1) {
    names(predicted_genvar)[-1:-3] <- c("pred_mu", "pred_varG")
  } else {
    names(predicted_genvar)[-1:-3] <- c("pred_mu", "pred_varG", "pred_corG")
  }
  
  # Return 
  return(predicted_genvar)

}
  



