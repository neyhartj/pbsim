#' Define the genetic model of one or more traits
#' 
#' @description 
#' Defines the genetic architecture of a trait.
#' 
#' @param genome An object of class \code{genome}.
#' @param qtl.model  A matrix specifying the QTL model. Each row corresponds to
#' a different QTL. The first column gives the chromosome number, the second 
#' column gives marker index, the third column gives the additive effect of
#' the favorable QTL allele (\code{a}) and the fourth column gives the dominance
#' effect at that QTL (\code{d}). If the matrix is one of NA, QTL will be
#' randomly assigned based on the number of rows in the matrix. If the genetic 
#' architecture of multiple traits is desired, a list of length \code{n_trait}
#' of \code{qtl.model} matrices can be provided. See \code{Details} for more 
#' information.
#' @param ... Other arguments. See \code{Details} for more information.
#' 
#' @details
#' To simulate QTL, genetic markers in the genome are assigned to become QTL. The 
#' \code{qtl.model} matrix specifies the information for this assignment. The 
#' first column in this matrix is the chromosome number. The second column is 
#' the marker index. For instance, to assign the 100th marker on a chromosome
#' to become a QTL, the marker index would be \code{100}. The third column is
#' the additive effect of the favorable QTL allele (\code{a}). Genotypes homozygous
#' for the favorable allele are assigned a genotypic value of \code{a} and
#' genotypes homozygous for the unfavorable allele are assigned a genotypic value
#' of \code{-a}. If dominance is absent, the genotypic value of heterozygotes is
#' \code{0}. The fourth column is the dominance effect at the QTL. This value can
#' be larger that \code{a} or smaller than \code{-a} (overdominance). The genotypic
#' value of heterozygotes at the QTL is \code{d}.
#' 
#' Other arguments include:
#' \itemize{
#'   \item{\code{add.dist}: The distribution of additive effects of QTL (if additive 
#'   effects are not provided in the \code{qtl.model} input). Can be 
#'   \code{"normal"} or \code{"geometric"}. For a 
#'   distribution of \code{"normal"}, additive effects are generated via the 
#'   \code{\link[stats]{rnorm}} function.
#'   For a distribution of \code{"geometric"}, additive effects are calculated for
#'   the k-th QTL as \eqn{a^k} where \eqn{a = (1 - L) / (1 + L)} and \eqn{L} is 
#'   the number of QTL (Lande and Thompson, 1990).}
#'   \item{\code{dom.dist}: The distribution of dominance effects of QTL (if dominance
#'   effects are not provided in the \code{qtl.model} input). Can be 
#'   \code{"normal"} for normally-distributed dominance effects.}
#'   \item{\code{prob.corr}: A matrix of two columns defining the probabilities
#'   that QTL for two or more traits are in some form that might confer genetic 
#'   correlation. The first column sets the maximum distance between a QTL from 
#'   a second trait and a QTL from the first trait, and the second column is the
#'   probability that QTL from a second trait have that maximum distance. Pleiotropic
#'   QTL can be simulated by providing a 0 in the first column, and no genetic 
#'   linkage is simulated if 50 is in the first column.}
#' }
#' 
#' Also note the following rules that apply when the \code{qtl.model} input is 
#' completely NA:
#' \itemize{
#'   \item{QTL positions are randomly drawn, with no regard to uniformity
#'   over chromosomes.}
#'   \item{The genetic architecture of multiple traits can be simulated by providing
#'   a list of \code{qtl.model} matrices. The length of this list will be the number
#'   of traits. The same rules apply for each \code{qtl.model} matrix. Pleitropy 
#'   can be simulated by designating QTL at the same location in the genome.}
#' }
#' 
#' @return 
#' A \code{genome} object with added information for the gentic model.
#' 
#' @examples 
#' n.mar  <- c(505, 505, 505)
#' len <- c(120, 130, 140)
#' 
#' genome <- sim_genome(len, n.mar)
#' 
#' chromosome <- c(1, 1, 2, 2, 3, 3)
#' pos <- as.numeric(sapply(X = genome$len, FUN = runif, n = 2, min = 0))
#' a <- c(1, 0.25, 0.5, 0.25, 0.25, 0.5)
#' d <- 0
#' 
#' qtl.model <- cbind(chromosome, pos, a, d)
#' 
#' genome <- sim_gen_model(genome, qtl.model)
#' 
#' # Randomly generate 15 QTL with additive allelic effects following a
#' # genometric series
#' qtl.model <- matrix(nrow = 15, ncol = 4)
#' 
#' genome <- sim_gen_model(genome, qtl.model, add.dist = "geometric")
#' 
#' # Randomly generate 15 QTL for each of two traits. 50% of the QTL are pleiotropic
#' na_mat <- matrix(nrow = 15, ncol = 4)
#' 
#' qtl.model <- list(na_mat, na_mat)
#' add.dist <- "geometric"
#' prob.corr <- cbind(c(0, 2, 15, 30, 50), c(0.5, 0.1, 0.1, 0.1, 0.2))
#' 
#' genome <- sim_gen_model(genome, qtl.model, add.dist = add.dist, prob.corr = prob.corr)
#'  
#' @import dplyr
#' @importFrom qtl sim.cross
#' @importFrom qtl pull.map
#' 
#' @export
#' 
sim_gen_model <- function(genome, qtl.model, ...) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Is there already genetic architecture?
  # If so clear it and replace the map
  if (!is.null(genome$gen_model)) {
    genome$gen_model <- NULL
    blank_cross <- qtl::sim.cross(map = genome$map, n.ind = 1)
    genome$map <- qtl::pull.map(cross = blank_cross)
  }
  
  # Is the qtl.model a list? If not make it one
  # Also determine the number of traits
  if (!is.list(qtl.model)) {
    qtl.model <- list(qtl.model)
    n_trait <- 1
  } else {
    n_trait <- length(qtl.model)
  }
  
  # Make sure each qtl.model matrix has four columns
  if (any(sapply(qtl.model, ncol) < 4))
    stop("The matrix or matrices in 'qtl.model' must have at least 4 columns.")
  
  # The QTL model must be all NA or no NA
  random_qtl <- sapply(qtl.model, function(mat) all(is.na(mat)))
  if (!sum(random_qtl) %in% c(0, length(random_qtl)))
    stop("For one or more traits, the qtl.model must be either totally complete (i.e.
         no NA) or totally missing (i.e. all NA).")
  
  
  # Extract other arguments
  other.args <- list(...)
  
  # Number of chromosomes
  n_chr <- nchr(genome)
  
  # If any should be randomly generated, continue
  if (all(random_qtl)) {
    
    # Empty list to store information
    qtl_specs <- vector("list", n_trait)
    
    ## Error handling
    # Determine the additive effect distribution
    add.dist <- other.args$add.dist
    dom.dist <- other.args$dom.dist
    
    # If the additive dist is NULL, error, otherwise check for consistency
    if (is.null(add.dist)) {
      stop("The qtl.model is NA but the distribution from which to draw additive
           effects was not provided. Options include 'normal' or 'geometric.'")
      
    } else {
      if (!add.dist %in% c("normal", "geometric"))
        stop("The argument 'add.dist' must be 'geometric' or 'normal.'")
    }
    
    # If dom.dist is not NULL, it must be "normal"
    if (!is.null(dom.dist)) 
      if(dom.dist != "normal")
        stop("If the argument 'dom.dist' is provided, it must equal 'normal.'")
    
    # If there is more than one trait, extract the correlation probability vector
    if (n_trait > 1) {
      
      prob.corr <- other.args$prob.corr
    
      # Error if null
      if (is.null(prob.corr))
        stop("If QTL for more than one trait are to be generated, the 'prob.corr'
             argument must be passed.")
      
      # Make sure it has two columns
      if (ncol(prob.corr) != 2)
        stop("The 'prob.corr' input must have 2 columns.")
      
      # Are all probabilities between 0 and 1?
      if (!all(prob.corr[,2] >= 0 & prob.corr[,2] <= 1))
        stop("The probabilities in 'prob.corr' are not all between 0 and 1.")
      
      # Do the probabilities sum to 1
      if (!sum(prob.corr[,2]) == 1)
        stop("The probabilities in 'prob.corr' do not sum to 1.")
      
      # Are any of the levels of p greater than 50 or less than 0
      if (!all(prob.corr[,1] >= 0 & prob.corr[,1] <= 50))
        stop("The distances in 'prob.corr' must be between 0 and 50.")
      
      # Sort the matrix on order of p
      prob.corr <- prob.corr[order(prob.corr[,1]),, drop = FALSE]
      
      
    }
    
    
    ## Simulate QTL for the first trait
    n_qtl <- nrow(qtl.model[[1]])
    
    # Simulate additive effects
    if (add.dist == "normal") 
      add.eff <- rnorm(n_qtl)
    
    if (add.dist == "geometric") {
      a <- (1 - n_qtl) / (1 + n_qtl)
      add.eff <- sample(abs( a ^ (seq(n_qtl)) ))
      # Randomly assign favorable or unfavorable for the first allele
      add.eff <- add.eff * sample(x = c(-1, 1), size = n_qtl, replace = TRUE)
    }
    
    # Simulate dominance effect 
    if (is.null(dom.dist)) {
      dom.eff <- rep(0, n_qtl)
    } else {
      dom.eff <- rnorm(n_qtl)
    }
    
    # Randomly draw chromosomes
    chr <- sort(sample(seq_len(n_chr), n_qtl, replace = TRUE))
    
    # Randomly designate genetic positions for QTL
    qtl.pos <- lapply(X = genome$len[chr], runif, n = 1, min = 0)
    pos <- unlist(qtl.pos)
      
    # Assemble the matrix
    qtl_specs[[1]] <- data.frame(chr = chr, pos = pos, add_eff = add.eff, dom_eff = dom.eff,
                                 qtl_name = paste("QTL", seq_along(chr), sep = ""),
                                 qtl1_pair = NA, stringsAsFactors = FALSE)
  
    
    
    ## Iterate over the remaining traits if more than 1
    if (n_trait > 1)
      for (t in seq(2, n_trait)) {
        
        # Number of QTL
        n_qtl <- nrow(qtl.model[[t]])
        
        # Simulate additive effects
        if (add.dist == "normal") 
          add.eff <- rnorm(n_qtl)
        
        if (add.dist == "geometric") {
          a <- (1 - n_qtl) / (1 + n_qtl)
          add.eff <- sample(abs( a ^ (seq(n_qtl)) ))
          # Randomly assign favorable or unfavorable for the first allele
          add.eff <- add.eff * sample(x = c(-1, 1), size = n_qtl, replace = TRUE)
        }
        
        # Simulate dominance effect 
        if (is.null(dom.dist)) {
          dom.eff <- rep(0, n_qtl)
        } else {
          dom.eff <- rnorm(n_qtl)
        }
        
        # Empty df with additive and dominance effects
        qtl_df <- data.frame(chr = NA, pos = NA, add_eff = add.eff, dom_eff = dom.eff, 
                             qtl_name = NA, qtl1_pair = NA, stringsAsFactors = FALSE)
        
        # If the length of prob.corr is 1, output a vector of that prob.corr
        if (nrow(prob.corr) == 1) {
          qtl_designator <- rep(x = prob.corr[,1], times = n_qtl)
        
        } else {
          # Randomly designated each QTL to share some correlation with QTL of the
          # first trait
          qtl_designator <- sample(prob.corr[,1], size = n_qtl, prob = prob.corr[,2], replace = TRUE)
        }
        
        # Pull out the correlation levels
        corr_level <- prob.corr[,1, drop = TRUE]
        
        ## Iterate over the correlation levels
        for (p in corr_level) {
          
          # If prob is 0, simulate pleiotropy by drawing chr and pos from the first trait
          if (p == 0) {
            
            pleiotropic <- sample(nrow(qtl_specs[[1]]), size = sum(qtl_designator == p))
            
            # Subset chromosomes and positions for these first QTL
            qtl_df[qtl_designator == p, 1:2] <- qtl_specs[[1]][pleiotropic, 1:2]
            
            # Extract the QTL names that are pleiotropic
            qtl_df[qtl_designator == p, 5:6] <- qtl_specs[[1]]$qtl_name[pleiotropic]
            
            
            # If p is 50, randomly draw positions and chromosomes, regardless of the first trait
          } else if (p == 50) {
            
            # Randomly draw chromosomes
            chr <- sort(sample(seq_len(n_chr), sum(qtl_designator == p), replace = TRUE))
            
            # Randomly designate genetic positions for QTL
            qtl.pos <- lapply(X = genome$len[chr], runif, n = 1, min = 0)
            pos <- unlist(qtl.pos)
            
            qtl_df[qtl_designator == p, 1:2] <- c(chr, pos)
            
            
            # Else sample QTL for linkage
          } else {
            
            # Sample QTL from trait 1
            qtl_one_sample <- sample(nrow(qtl_specs[[1]]), size = sum(qtl_designator == p))
            
            # Get those QTL positions
            qtl_one_pos <- qtl_specs[[1]][qtl_one_sample, 1:2, drop = FALSE]
            
            # What is the previous level of p? If it is the first, set the minimum distance to 0.000001
            prev_p <- match(p, prob.corr[,1]) - 1
            min_dist <- ifelse(prev_p == 0, 1e-6, corr_level[prev_p])
            
            # If min_dist is 0, set to 1e-6
            min_dist <- ifelse(min_dist == 0, 1e-6, min_dist)
            
            # Max dist is equal to p
            max_dist <- p
            
            
            # Sample a deviation based on p and then sample whether it is positive or negative
            sample_deviation <- runif(n = length(qtl_one_sample), min = min_dist, max = max_dist) %>%
            {. * sample(c(-1, 1), size = length(.), replace = TRUE)}
            
            # Designate the new position of this QTL
            new_pos <- qtl_one_pos[,2] + sample_deviation
            
            qtl_two_pos <- cbind(qtl_one_pos[,1, drop = FALSE], new_pos)
            
            # Reject if the new position is less than 0 or more than the chromosome length
            while (any(apply(qtl_two_pos, MARGIN = 1, function(qtl) 
              qtl[2] < 0 | qtl[2] > genome$len[qtl[1]]))) {
              
              # Sample a deviation based on p and then sample whether it is positive or negative
              sample_deviation <- runif(n = length(qtl_one_sample), min = min_dist, max = max_dist) %>%
              {. * sample(c(-1, 1), size = length(.), replace = TRUE)}
              
              # Designate the new position of this QTL
              new_pos <- qtl_one_pos[,2] + sample_deviation
              
              qtl_two_pos <- cbind(qtl_one_pos[,1, drop = FALSE], new_pos)
              
            }
            
            # Add new positions to the data.frame
            qtl_df[qtl_designator == p, 1:2] <- qtl_two_pos
            
            # Extract the QTL names that are paired
            qtl_df[qtl_designator == p, 6] <- qtl_specs[[1]]$qtl_name[qtl_one_sample]
            
          }
        } # Close the loop
        
        
        ## Any qtl_names that are NA get the next highest QTL name
        # number of QTL from the previous trait
        start_name <- nrow(qtl_specs[[t - 1]]) + 1
        
        # Number of QTL missing names
        qtl_missing <- sum(is.na(qtl_df$qtl_name))
        
        qtl_df$qtl_name[is.na(qtl_df$qtl_name)] <- 
          paste("QTL", seq(start_name, start_name + (qtl_missing - 1)), sep = "")
        
        # Add the mat to the specs
        qtl_specs[[t]] <- qtl_df
        
      } # Close the loop
    
  } else {
    # Else verify that the provided qtl.model matrices are sufficient
    
    for (mat in qtl.model) {
      
      # The qtl model must not have any NAs
      if (any(is.na(mat)))
        stop("Unless the 'qtl.model' is completely NA, there can be no NA elements.")
      
      if (!all(mat[,1] %in% seq(n_chr)))
        stop("The chromosome numbers in 'qtl.model' are not chromosomes in the 
             genome.")
      
      # Are the marker positions correct?
      qtl.pos <- split(mat[,2], mat[,1])
      
      if (!all(mapply(qtl.pos, genome$len, FUN = function(q, n) q <= n)))
        stop("The QTL positions in 'qtl.model' are not within the length of the 
             chromosomes.")
    }
    
    # Rename to qtl.specs - and convert to df
    qtl_specs <- lapply(qtl.model, as.data.frame) %>% 
      lapply(structure, names = c("chr", "pos", "add_eff", "dom_eff"))
    
  }
  
  # Add the genetic model to the genome
  genome[["gen_model"]] <- qtl_specs %>%
    lapply(arrange, chr, pos)
  
  ## Add names of QTL if not present
  # Pull out all QTL
  all_qtl <- pull_qtl(genome = genome, unique = FALSE)
  
  if (!"qtl_name" %in% names(all_qtl)) {
    
    # Subset the unique and add names
    unique_qtl <- all_qtl %>%
      distinct(chr, pos, .keep_all = TRUE) %>%
      mutate(qtl_name = paste("QTL", seq(n()), sep = ""),
             qtl1_pair = NA) %>%
      select(chr, pos, qtl_name)
    
    # Merge back with the total QTL
    all_qtl <- full_join(all_qtl, unique_qtl, by = c("chr" = "chr", "pos" = "pos")) %>% 
      split(.$trait) %>% 
      lapply(select, -trait)
    
  
    # Add back to genome
    genome$gen_model <- all_qtl
    
  }
  
  # Get unique QTL
  unique_qtl <- pull_qtl(genome = genome, unique = TRUE)
  
  # Add QTL to the map
  qtl_pos <- split(x = structure(unique_qtl$pos, names = unique_qtl$qtl_name), f = unique_qtl$chr)
  new_map <- mapply(genome$map, qtl_pos, FUN = c, SIMPLIFY = FALSE) %>%
    lapply(sort) %>%
    lapply(structure, class = "A")
  
  class(new_map) <- "map"
  
  genome$map <- new_map
  
  return(genome)
    
} # Close the function





#' Adjust the genetic model of one or more traits
#' 
#' @description 
#' Edits the additive effect of QTL of a second or more trait by assessing the
#' correlation among QTL for that trait with the designated pair in the first trait.
#' 
#' @param genome An object of class \code{genome}.
#' @param geno Genotype data on a population to phenotype. Can be a matrix of dimensions
#' \code{n.ind} x \code{n.loci}, the elements of which must be z {0, 1, 2}, or a list
#' of such matrices.  
#' @param pos.cor A logical indicating whether the two traits should be positively correlated.
#' 
#' @details 
#' Although QTL for two traits can be linked, depending on the population, they
#' may not be positively correlated. If a correlation between two traits is desired,
#' this function will edit the additive effects in all but the first trait to 
#' induce the desired correlation between traits.
#' 
#' @import dplyr
#' 
#' @export
#' 
adj_gen_model <- function(genome, geno, pos.cor = TRUE) {
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = geno))
    stop("The geno did not pass. See warning for reason.")
  
  # Make sure the genome has a genetic model
  if (is.null(genome$gen_model))
    stop("The 'genome' must have a genetic model. Use 'sim_gen_model' to create one.")
  
  # If the geno input is a list, recombine
  if (is.list(geno))
    geno <- do.call("cbind", geno)
  
  stopifnot(is.logical(pos.cor))
  
  # Are the QTL in the geno object?
  qtl_names <- pull_qtl(genome, unique = TRUE)$qtl_name
  
  if (!all(qtl_names %in% colnames(geno))) 
    stop("The QTL genotypes are not in the 'geno' object.")
  
  # Pull out the genotypic data for that QTL and the QTL1 pair
  qtl_geno <- pull_genotype(genome = genome, geno = geno, loci = qtl_names)
  
  # For each QTL in the second or more traits...
  for (t in seq(2, length(genome$gen_model))) {
    
    # Pull out the genetic model
    qtlmod_t <- genome$gen_model[[t]]
    
    # Empty vector to store modified qtl effects
    adj_add_eff <- vector("numeric", nrow(qtlmod_t))
    
    # Iterate over qtl
    for (i in seq(nrow(qtlmod_t))) {
      
      q <- qtlmod_t[i,]
      
      # If the qtl1_pair is NA, skip
      if (is.na(q$qtl1_pair)) {
        adj_add_eff[i] <- q$add_eff
        
      } else {
        
        # If the qtl_t name is the same as the qtl1 pair, the cor sign is 1
        if (q$qtl_name == q$qtl1_pair) {
          cor_sign <- 1
          
        } else {
          # Pull out the genotypic data for that QTL and the QTL1 pair
          qtl_geno_sub <- subset(qtl_geno, select = c(q$qtl_name, q$qtl1_pair))
        
          # What is the sign of the correlation?
          cor_sign <- sign(cor(qtl_geno[,1], qtl_geno[,2]))
          
        }
        
        # What is the sign of the qtl_t add_eff?
        qtl_t_sign <- sign(q$add_eff)
        # What is the sign of the qtl1 pair?
        qtl1_sign <- sign(subset(genome$gen_model[[1]], qtl_name == q$qtl1_pair, add_eff))
        
        # If the intended correlation is positive...
        if (pos.cor) {
          # If the geno corr sign is positive, the add_eff for qtl_t should be the same 
          # sign as that for qtl1
          if (cor_sign == 1) {
            adj_add_eff[i] <- ifelse(qtl_t_sign == qtl1_sign, q$add_eff, q$add_eff * -1)
          } else{
            adj_add_eff[i] <- ifelse(qtl_t_sign == qtl1_sign, q$add_eff * -1, q$add_eff)
          }
          
          # If the intended correlation is negative...
        } else {
          # If the geno corr sign is positive, the add_eff for qtl_t should be the opposite 
          # sign as that for qtl1
          if (cor_sign == 1) {
            adj_add_eff[i] <- ifelse(qtl_t_sign == qtl1_sign, q$add_eff * -1, q$add_eff)
          } else{
            adj_add_eff[i] <- ifelse(qtl_t_sign == qtl1_sign, q$add_eff, q$add_eff * -1)
          }
          
        }
        
      }} # Close the loop
    
    # Edit the qtl model
    genome$gen_model[[t]] <- qtlmod_t %>%
      mutate(add_eff = adj_add_eff)
    
  } # Close the loop
  
  # Return the genome
  return(genome)
  
} # Close the function