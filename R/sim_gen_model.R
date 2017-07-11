#' Define the genetic model of a trait
#' 
#' @description 
#' Defines the genetic architecture of a trait.
#' 
#' @param genome An object of class \code{genome}.
#' @param qtl.model  A matrix specifying the QTL model. Each row corresponds to
#' a different QTL. The first column gives the chromosome number, the second 
#' column gives the locus position (in cM), the third column gives the additive effect of
#' the favorable QTL allele (\code{a}) and the fourth column gives the dominance
#' effect at that QTL (\code{d}). If the matrix is one of NA, QTL will be
#' randomly assigned based on the number of rows in the matrix.
#' @param ... Other arguments. See \emph{Details} for more information.
#' 
#' @details
#' QTL are simulated by sampling or specifying existing markers, which become "hidden."
#' 
#' The \code{qtl.model} matrix specifies the information for this assignment. The 
#' first column in this matrix is the chromosome number. The second column is 
#' the QTL position. The third column is the additive effect of the "1" allele
#' at the QTL (\code{a}). Genotypes homozygous for the "1" allele are assigned a 
#' genotypic value of \code{a} and genotypes homozygous for the "-1" allele are 
#' assigned a genotypic value of \code{-a}. The value of \code{a} provided in 
#' \code{qtl.model} can be negative. The fourth column is the dominance effect at the 
#' QTL. If non-zero, this value can be larger that \code{a} or smaller than \code{-a} 
#' (overdominance). The genotypic value of heterozygotes at the QTL is \code{d}.
#' 
#' Other arguments include:
#' \describe{
#'   \item{\code{add.dist}}{The distribution of additive effects of QTL (if additive 
#'   effects are not provided in the \code{qtl.model} input). Can be 
#'   \code{"normal"} or \code{"geometric"}. For a distribution of \code{"normal"}, 
#'   additive effects are generated via the \code{\link[stats]{rnorm}} function.
#'   For a distribution of \code{"geometric"}, additive effects are calculated for
#'   the k-th QTL as \eqn{a^k} where \eqn{a = (1 - L) / (1 + L)} and \eqn{L} is 
#'   the number of QTL (Lande and Thompson, 1990).}
#'   \item{\code{dom.dist}}{The distribution of dominance effects of QTL (if dominance
#'   effects are not provided in the \code{qtl.model} input). Can be 
#'   \code{"normal"} for normally-distributed dominance effects.}
#'   \item{\code{max.qtl}}{The maximum number of QTL in the simulation experiment. Must
#'   be passed if the QTL are randomly sampled. If a trait is controlled by \emph{L} QTL
#'   and \code{max.qtl = M}, then \emph{M} - \emph{L} QTL are given NULL effects.
#'   This is useful if you want to simulate variable genetic architecture, but keep
#'   the number of SNP markers constant.}
#' }
#' 
#' Also note the following rules that apply when the \code{qtl.model} input is 
#' completely NA:
#' \itemize{
#'   \item{QTL positions are randomly drawn, with no regard to uniformity
#'   over chromosomes.}
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
#' @import dplyr
#' 
#' @export
#' 
sim_gen_model <- function(genome, qtl.model, ...) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Extract information based on the type
  type <- attr(genome, "type")
  
  # Is there already genetic architecture?
  # If so clear it
  if (!is.null(genome$gen_model)) {
    genome$gen_model <- NULL
  }
  
  
  # The qtl.model should be a matrix, not a list
  if (!is.matrix(qtl.model))
    stop("The 'qtl.model' input must be a matrix. Note that multi-trait genetic models
         are not accetpted by this function. Use the 'sim_multi_gen_model for multi-trait
         genetic models.")
  
  # Make sure each qtl.model matrix has four columns
  if (ncol(qtl.model) < 4)
    stop("The matrix or matrices in 'qtl.model' must have at least 4 columns.")
  
  # The QTL model must be all NA or no NA
  random_qtl <- all(is.na(qtl.model))
  
  # Extract other arguments
  other.args <- list(...)
  
  # Number of chromosomes
  n_chr <- nchr(genome)
  
  # length of chromosomes
  chr_len <- chrlen(genome)
  
  # If any should be randomly generated, continue
  if (random_qtl) {
    
    # Empty matrix to store information
    qtl_specs <- vector("list")
    
    ## Error handling
    # Determine the additive effect distribution
    add.dist <- other.args$add.dist
    dom.dist <- other.args$dom.dist
    
    # If the additive dist is NULL, error, otherwise check for consistency
    if (is.null(add.dist)) {
      stop("The qtl.model is NA but the distribution from which to draw additive
           effects was not provided. Options include 'normal' or 'geometric.'")
      
    } else if (!add.dist %in% c("normal", "geometric")) {
        stop("The argument 'add.dist' must be 'geometric' or 'normal.'")
    }
    
    # If dom.dist is not NULL, it must be "normal"
    if (!is.null(dom.dist)) 
      if(dom.dist != "normal")
        stop("If the argument 'dom.dist' is provided, it must equal 'normal.'")
  
    ## Simulate QTL for the first trait
    n_qtl <- nrow(qtl.model)
    
    # Simulate additive effects
    if (add.dist == "normal") {
      add.eff <- rnorm(n_qtl) 
      
    } else if (add.dist == "geometric") {
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
    
    

    # If the max.qtl argument is FALSE, set it to the maximum number of QTL
    # for any trait
    max.qtl <- other.args$max.qtl
    
    if (is.null(max.qtl))
      max.qtl <- nrow(qtl.model)
    

    
    # Sample marker names to become QTL
    marker_sample <- sample(x = markernames(genome, include.qtl = TRUE), size = max.qtl)
    
    # Get the positions of those markers
    marker_sample_pos <- find_markerpos(genome = genome, marker = marker_sample) %>%
      mutate(marker = row.names(.))
    
    # Set the new length of add.eff and dom.eff
    length(add.eff) <- max.qtl
    length(dom.eff) <- max.qtl
    
    # Pad with 0, not NA, and add to the df
    marker_sample_pos1 <- marker_sample_pos %>%
      mutate(add.eff = ifelse(is.na(add.eff), 0, add.eff),
             dom.eff = ifelse(is.na(dom.eff), 0, dom.eff)) %>%
      arrange(chr, pos)
    
    # Extract the chr and pos
    chr <- marker_sample_pos1$chr
    pos <- marker_sample_pos1$pos
    add.eff <- marker_sample_pos1$add.eff
    dom.eff <- marker_sample_pos1$dom.eff
    
    # Assign marker name for these QTL
    qtl_marker_name <- marker_sample_pos1$marker
      
    # Assemble the matrix
    qtl_specs[[1]] <- data.frame(chr = chr, pos = pos, add_eff = add.eff, dom_eff = dom.eff,
                                 qtl_name = qtl_marker_name, qtl1_pair = NA, 
                                 stringsAsFactors = FALSE)
    
  } else {
    # Else verify that the provided qtl.model matrices are sufficient

    # The qtl model must not have any NAs
    if (any(is.na(qtl.model)))
      stop("Unless the 'qtl.model' is completely NA, there can be no NA elements.")
    
    if (!all(qtl.model[,1] %in% seq(n_chr)))
      stop("The chromosome numbers in 'qtl.model' are not chromosomes in the 
           genome.")
    
    # Are the marker positions correct?
    qtl.pos <- split(qtl.model[,2], qtl.model[,1])
    
    if (!all(mapply(qtl.pos, genome$len, FUN = function(q, n) q <= n)))
      stop("The QTL positions in 'qtl.model' are not within the length of the 
           chromosomes.")
    
    # Rename to qtl.specs - and convert to df
    qtl_specs <- list(as.data.frame(qtl.model)) %>% 
      lapply(structure, names = c("chr", "pos", "add_eff", "dom_eff"))
    
  }
  
  # Add the genetic model to the genome
  genome[["gen_model"]] <- qtl_specs %>%
    lapply(arrange, chr, pos) %>%
    lapply(mutate, chr = factor(chr, levels = chrnames(genome)))
  
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
  
  
  return(genome)
    
} # Close the function



#' Define a genetic model for two or more traits
#' 
#' @description 
#' Generates a genetic model for two or more traits with desired genetic correlation.
#' 
#' @param genome An object of class \code{genome}.
#' @param qtl.model  A matrix specifying the QTL model. Each row corresponds to
#' a different QTL. The first column gives the chromosome number, the second 
#' column gives the locus position (in cM), the third column gives the additive effect of
#' the favorable QTL allele (\code{a}) and the fourth column gives the dominance
#' effect at that QTL (\code{d}). If the matrix is one of NA, QTL will be
#' randomly assigned based on the number of rows in the matrix. If the genetic 
#' architecture of multiple traits is desired, a list of length \code{n_trait}
#' of \code{qtl.model} matrices can be provided. See \emph{Details} for more 
#' information.
#' @param geno Genotype data on a population to phenotype. Can be a matrix of dimensions
#' \code{n.ind} x \code{n.loci}, the elements of which must be z {0, 1, 2}, or a list
#' of such matrices.
#' @param prob.corr A matrix of three columns defining the probabilities
#' that QTL for two or more traits are in some form that might confer genetic 
#' correlation. The first column sets the maximum distance between a QTL from 
#' a second trait and a QTL from the first trait, the second column is the
#' probability that QTL from a second trait have that maximum distance, and the 
#' third column specifies the minimum correlation (i.e. linkage disequilibrium) 
#' between QTL of one trait and QTL of another. Pleiotropic QTL can be simulated 
#' by providing a 0 in the first column, and no genetic linkage is simulated if 
#' 50 is in the first column.
#' @param pos.cor A logical indicating whether the two traits should be positively correlated.
#' If not passed, the correlation is not forced to be positive or negative.
#' 
#' @details
#' QTL are simulated by sampling or specifying existing markers, which become "hidden."
#' 
#' The \code{qtl.model} matrix specifies the information for this assignment. The 
#' first column in this matrix is the chromosome number. The second column is 
#' the QTL position. The third column is the additive effect of the "1" allele
#' at the QTL (\code{a}). Genotypes homozygous for the "1" allele are assigned a 
#' genotypic value of \code{a} and genotypes homozygous for the "-1" allele are 
#' assigned a genotypic value of \code{-a}. The value of \code{a} provided in 
#' \code{qtl.model} can be negative. The fourth column is the dominance effect at the 
#' QTL. If non-zero, this value can be larger that \code{a} or smaller than \code{-a} 
#' (overdominance). The genotypic value of heterozygotes at the QTL is \code{d}.
#'  
#' Although QTL for two traits can be linked, depending on the population, they
#' may not be highly correlated. If a correlation between two traits is desired,
#' this function will first resample markers to become QTL based on the desired
#' linkaged and correlations between the genotypes at those QTL. Next, the function
#' will edit the additive effects in all but the first trait to 
#' induce the desired correlation between traits.
#' 
#' Other arguments include:
#' \describe{
#'   \item{\code{add.dist}}{The distribution of additive effects of QTL (if additive 
#'   effects are not provided in the \code{qtl.model} input). Can be 
#'   \code{"normal"} or \code{"geometric"}. For a distribution of \code{"normal"}, 
#'   additive effects are generated via the \code{\link[stats]{rnorm}} function.
#'   For a distribution of \code{"geometric"}, additive effects are calculated for
#'   the k-th QTL as \eqn{a^k} where \eqn{a = (1 - L) / (1 + L)} and \eqn{L} is 
#'   the number of QTL (Lande and Thompson, 1990).}
#'   \item{\code{dom.dist}}{The distribution of dominance effects of QTL (if dominance
#'   effects are not provided in the \code{qtl.model} input). Can be 
#'   \code{"normal"} for normally-distributed dominance effects.}
#'   \item{\code{max.qtl}}{The maximum number of QTL in the simulation experiment. Must
#'   be passed if the QTL are randomly sampled. If a trait is controlled by \emph{L} QTL
#'   and \code{max.qtl = M}, then \emph{M} - \emph{L} QTL are given NULL effects.
#'   This is useful if you want to simulate variable genetic architecture, but keep
#'   the number of SNP markers constant.}
#' }
#' 
#' @examples 
#' # Load some historic data
#' data("s2_cap_haploid")
#' data("s2_snp_info")
#' 
#' # Create a genome
#' map <- lapply(split(s2_snp_info, s2_snp_info$chrom), 
#'               function(chr) structure(chr$cM_pos, names = chr$rs) )
#' 
#' genome <- sim_genome(map = map, type = "hypred")
#' 
#' # Simulate two traits with 30 QTL, 1/3 with pairwise linkage of 0 < x <= 2 cM
#' # and 2/3 with pairwise linkage of 2 < x <= 20. The first set of QTL should have
#' # an LD of at least 0.5, and the second set of QTL should have an LD of at least
#' # 0.1
#' qtl.model <- replicate(2, matrix(nrow = 30, ncol = 4), simplify = FALSE)
#' prob.corr <- cbind(c(2, 20), c(0.333, 0.667), c(0.5, 0.1))
#' 
#' # Simulate the genetic architecture
#' genome <- sim_multi_gen_model(genome = genome, qtl.model = qtl.model, 
#'                               geno = s2_cap_haploid, prob.corr = prob.corr,
#'                               pos.corr = TRUE, add.dist = "geometric", max.qtl = 30)
#' 
#' pop <- create_pop(genome = genome, geno = s2_cap_haploid)
#' 
#' # Genetic correlation prior to adjustment
#' cor(pop$geno_val$trait1, pop$geno_val$trait2)
#' 
#' # Adjust the genome
#' prob.corr <- cbind(2, 1, 0.3)
#' genome <- adj_gen_model(genome = genome, geno = s2_cap_haploid, prob.corr = prob.corr, 
#' pos.corr = T)
#' pop <- create_pop(genome = genome, geno = s2_cap_haploid)
#' 
#' # Genetic correlation after adjustment
#' cor(pop$geno_val$trait1, pop$geno_val$trait2)
#' 
#' @import dplyr
#' @importFrom purrr pmap_chr
#' 
#' @export
#' 
sim_multi_gen_model <- function(genome, qtl.model, geno, prob.corr, pos.corr, ...) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Extract information based on the type
  type <- attr(genome, "type")
  
  # Is there already genetic architecture?
  # If so clear it
  if (!is.null(genome$gen_model)) {
    genome$gen_model <- NULL
  }
  
  # Is the qtl.model a list? If not make it one
  # Also determine the number of traits
  if (!is.list(qtl.model)) {
    qtl.model <- list(qtl.model)
    n_trait <- 1
    
  } else {
    n_trait <- length(qtl.model)
    # Sort the QTL models based on decreasing number of QTL
    n_qtl_per_trait <- sapply(X = qtl.model, FUN = nrow)
    
    qtl.model <- qtl.model[order(n_qtl_per_trait, decreasing = TRUE)]
  }
  
  # Make sure each qtl.model matrix has four columns
  if (any(sapply(qtl.model, ncol) < 4))
    stop("The matrix or matrices in 'qtl.model' must have at least 4 columns.")
  
  # The QTL model must be all NA or no NA
  random_qtl <- sapply(qtl.model, function(mat) all(is.na(mat)))
  
  if (!sum(random_qtl) %in% c(0, length(random_qtl)))
    stop("For one or more traits, the qtl.model must be either totally complete (i.e.
         no NA) or totally missing (i.e. all NA).")
  
  # If the geno input is a list, combine
  if (is.list(geno))
    geno <- do.call("cbind", geno)
  
  # If the geno input is haploid data, convert to geno
  if (is_haploid(geno))
    geno <- haploid_to_geno(geno)
  
  
  # pos.corr must be logicial, if present. If not, set to NULL
  if (missing(pos.corr) | is.null(pos.corr)) {
    pos.corr <- NULL
    
  } else {
    pos.corr <- as.logical(pos.corr)
    
  }
  
  
  ## Error handling of prob.corr
  # Make sure it has three columns
  if (ncol(prob.corr) < 3)
    stop("The 'prob.corr' input must have 3 columns.")
  
  # Are all probabilities between 0 and 1?
  if (!all(prob.corr[,2] >= 0 & prob.corr[,2] <= 1))
    stop("The probabilities in 'prob.corr' are not all between 0 and 1.")
  
  # Are all correlations between 0 and 1?
  if (!all(prob.corr[,3] >= 0 & prob.corr[,3] <= 1))
    stop("The minimum correlations in 'prob.corr' are not all between 0 and 1.")
  
  # Do the probabilities sum to 1
  if (!sum(prob.corr[,2]) == 1)
    stop("The probabilities in 'prob.corr' do not sum to 1.")
  
  # Are any of the levels of p greater than 50 or less than 0
  if (!all(prob.corr[,1] >= 0 & prob.corr[,1] <= 50))
    stop("The distances in 'prob.corr' must be between 0 and 50.")
  
  # Sort the matrix on order of p
  prob.corr <- prob.corr[order(prob.corr[,1]),, drop = FALSE]
  
  # Extract other arguments
  other.args <- list(...)
  
  # Number of chromosomes
  n_chr <- nchr(genome)
  
  # length of chromosomes
  chr_len <- chrlen(genome)
  
  
  # How many traits?
  n_trait <- length(qtl.model)
  
  # If only one trait, return unedited
  if (n_trait < 2)
    stop("The 'qtl.model' must call for two or more traits.")
  
  
  # If any should be randomly generated, continue
  if (any(random_qtl)) {
    
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
      
    } else if (!add.dist %in% c("normal", "geometric")) {
      stop("The argument 'add.dist' must be 'geometric' or 'normal.'")
    }
    
    # If dom.dist is not NULL, it must be "normal"
    if (!is.null(dom.dist)) 
      if(dom.dist != "normal")
        stop("If the argument 'dom.dist' is provided, it must equal 'normal.'")
    
    
    ## Simulate QTL for the first trait
    n_qtl <- nrow(qtl.model[[1]])
    
    # Simulate additive effects
    if (add.dist == "normal") {
      add.eff <- rnorm(n_qtl) 
      
    } else if (add.dist == "geometric") {
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
    
    
    
    # If the max.qtl argument is NULL, set it to the maximum number of QTL
    # for any trait
    max.qtl <- other.args$max.qtl
    
    if (is.null(max.qtl))
      max.qtl <- max(sapply(X = qtl.model, nrow))
    
    
    
    # Sample marker names to become QTL
    marker_sample <- sample(x = markernames(genome, include.qtl = TRUE), size = max.qtl)
    
    # Get the positions of those markers
    marker_sample_pos <- find_markerpos(genome = genome, marker = marker_sample) %>%
      mutate(marker = row.names(.))
    
    # Set the new length of add.eff and dom.eff
    length(add.eff) <- max.qtl
    length(dom.eff) <- max.qtl
    
    # Pad with 0, not NA, and add to the df
    marker_sample_pos1 <- marker_sample_pos %>%
      mutate(add.eff = ifelse(is.na(add.eff), 0, add.eff),
             dom.eff = ifelse(is.na(dom.eff), 0, dom.eff)) %>%
      arrange(chr, pos)
    
    # Extract the chr and pos
    chr <- marker_sample_pos1$chr
    pos <- marker_sample_pos1$pos
    add.eff <- marker_sample_pos1$add.eff
    dom.eff <- marker_sample_pos1$dom.eff
    
    # Assign marker name for these QTL
    qtl_marker_name <- marker_sample_pos1$marker
    
    # Assemble the matrix
    qtl_specs[[1]] <- data.frame(chr = chr, pos = pos, add_eff = add.eff, dom_eff = dom.eff,
                                 qtl_name = qtl_marker_name, qtl1_pair = NA, 
                                 stringsAsFactors = FALSE)
    
    
    for (t in seq(2, n_trait)) {
      
      # Number of QTL
      n_qtl <- nrow(qtl.model[[t]])
      
      # Simulate additive effects
      if (add.dist == "normal") { 
        add.eff <- rnorm(n_qtl)
        
      } else if (add.dist == "geometric") {
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
      
      # Set the new length of add.eff and dom.eff
      length(add.eff) <- max.qtl
      length(dom.eff) <- max.qtl
      
      # If the length of prob.corr is 1, output a vector of that prob.corr
      if (nrow(prob.corr) == 1) {
        qtl_designator <- rep(x = prob.corr[,1], times = n_qtl)
        
      } else {
        # Randomly designated each QTL to share some correlation with QTL of the
        # first trait
        qtl_designator <- sample(prob.corr[,1], size = n_qtl, 
                                 prob = prob.corr[,2], replace = TRUE)
      }
      
      # Pull out the linkage information and desired LD
      corr_level <- prob.corr[,c(1,3), drop = FALSE]
      
      # Create an empty data.frame
      qtl_specs[[t]] <- as.data.frame(
        matrix(data = NA, nrow = max.qtl, ncol = 6, 
               dimnames = list(NULL, names(qtl_specs[[1]])))
      )
      
      ## Create a vector of row numbers from which to draw QTL from trait 1
      ## These must be "real" QTL (i.e. add_eff > 0)
      qtl1_row_vec <- which(qtl_specs[[1]]$add_eff != 0)
      
      # Extract the qtl1 names
      qtl1_names <- qtl_specs[[1]]$qtl_name
      
      
      ## Iterate over the correlation levels
      for (i in seq(nrow(corr_level))) {
        
        # Extract the linkage level p
        p <- corr_level[i,1, drop = TRUE]
        # Extract the LD level r
        r <- corr_level[i,2, drop = TRUE]
        
        # If prob is 0, simulate pleiotropy by drawing chr and pos from the first trait
        if (p == 0) {
          
          qtl_one_sample <- sample(qtl1_row_vec, size = sum(qtl_designator == p))
          
          # Subset chromosomes and positions for these first QTL
          qtl_specs[[t]][which(qtl_designator == p), ] <- 
            qtl_specs[[1]][qtl_one_sample, ]
          qtl_specs[[t]]$qtl1_pair[which(qtl_designator == p)] <- 
            qtl_specs[[1]]$qtl_name[qtl_one_sample]
          
          # Replace the qtl1_row_vec with the the index of non-sampled qtl
          qtl1_row_vec <- setdiff(qtl1_row_vec, qtl_one_sample)
          
          # Sample QTL for linkage
        } else if (p > 0 & p < 50) {
          
          # Sample QTL from trait 1
          qtl_one_sample <- sample(qtl1_row_vec, size = sum(qtl_designator == p))
          
          # If the length of the sample is 0, skip
          if (length(qtl_one_sample) == 0)
            next
          
          # Get those QTL positions
          qtl_one_pos <- qtl_specs[[1]][qtl_one_sample, , drop = FALSE]
          
          # What is the previous level of p? If it is the first, set the minimum 
          # distance to 0.000001
          prev_p <- match(p, prob.corr[,1]) - 1
          min_dist <- ifelse(prev_p == 0, 1e-6, corr_level[prev_p])
          
          # If min_dist is 0, set to 1e-6
          min_dist <- ifelse(min_dist == 0, 1e-6, min_dist)
          
          # Max dist is equal to p
          max_dist <- p
          
          
          # Sample from markers in the range
          prox_mar <- find_proxmarkers(genome = genome, marker = qtl_one_pos$qtl_name,
                                       min.dist = min_dist, max.dist = max_dist)
          
          # Remove those that are from the first trait
          prox_mar_unique <- lapply(prox_mar, setdiff, qtl1_names)
          
          # Iterate over each list name and sample proximal markers
          # After each iteration, remove the sampled marker from the remainder of the list
          prox_mar_sample <- structure(vector("character", length(prox_mar_unique)), 
                                       names = names(prox_mar_unique))
          
          for (j in seq_along(prox_mar_unique)) {
            
            qtl1 <- names(prox_mar_unique)[j]
            prox_qtl <- prox_mar_unique[[j]]
            
            # If prox_qtl has length 0, return NA
            if (length(prox_qtl) == 0)
              prox_mar_sample[j] <- NA
            
            geno_q <- pull_genotype(genome = genome, geno = geno, loci = c(qtl1, prox_qtl))
            # Subset for the qtl1
            qtl1_cor <- suppressWarnings(cor(geno_q)[qtl1, , drop = TRUE])
            
            # Are there markers with correlations above the threshold?
            test_corr <- abs(qtl1_cor[-1]) >= r
            
            if (any(na.omit(test_corr))) {
              # If so, sample one and return
              prox_mar_sample[j] <- sample(names(which(test_corr)), 1)
              
            } else {
              # If not, return NA
              prox_mar_sample[j] <- NA
              
            }
              
              # Remove that marker from sampling in the future
              prox_mar_unique <- sapply(prox_mar_unique, setdiff, prox_mar_sample[j])
            
          }
          
          # Create an empty data.frame
          marker_sample_pos <- as.data.frame(
            matrix(data = NA, nrow = length(prox_mar_sample), ncol = 3, 
                   dimnames = list(NULL, c("chr", "pos", "marker"))))
          
          # Get the positions of those markers
          marker_sample_pos[!is.na(prox_mar_sample),] <- 
            find_markerpos(genome = genome, marker = na.omit(prox_mar_sample)) %>%
            mutate(marker = row.names(.))
          
          
          # Assign this info to the qtl_specs df for trait t
          qtl_specs[[t]][which(qtl_designator == p), c("chr", "pos", "qtl_name")] <- 
            marker_sample_pos
          
          # Extract the QTL names that are paired
          qtl_specs[[t]][which(qtl_designator == p), "qtl1_pair"] <- 
            qtl_specs[[1]]$qtl_name[qtl_one_sample]
          
          # Replace the qtl1_row_vec with the the index of non-sampled qtl
          qtl1_row_vec <- setdiff(qtl1_row_vec, qtl_one_sample)
          
          
          # If p is 50, randomly draw positions and chromosomes, regardless of the first trait
        } else if (p == 50) {
          
          
          # Randomly draw markers to be QTL, excluding those already designated as QTL
          avail_markers <- setdiff(markernames(genome, include.qtl = TRUE), 
                                   unlist(lapply(qtl_specs, "[[", "qtl_name")))
          
          marker_sample <- sample(x = avail_markers, size = sum(qtl_designator == p))
          
          # Get the positions of those markers
          marker_sample_pos <- find_markerpos(genome = genome, marker = marker_sample) %>%
            mutate(marker = row.names(.))
          
          # Assign this info to the qtl_specs df for trait t
          qtl_specs[[t]][which(qtl_designator == p), c("chr", "pos", "qtl_name")] <- 
            marker_sample_pos
          
          
        } # Close the p object if statements 
        
      } # Close the loop per cor level
      
      # Re-add the additive and dominance effects
      qtl_specs[[t]]$add_eff <- ifelse(is.na(add.eff), 0, add.eff)
      qtl_specs[[t]]$dom_eff <- ifelse(is.na(dom.eff), 0, dom.eff)
      
      ## Fill in markers that are NA
      # Randomly draw markers to be QTL, excluding those already designated as QTL
      avail_markers <- setdiff(markernames(genome, include.qtl = TRUE), 
                               unlist(lapply(qtl_specs, "[[", "qtl_name")))
      
      marker_sample <- sample(x = avail_markers, size = sum(is.na(qtl_specs[[t]]$qtl_name)))
      
      # Get the positions of those markers
      marker_sample_pos <- find_markerpos(genome = genome, marker = marker_sample) %>%
        mutate(marker = row.names(.))
      
      # Assign this info to the qtl_specs df for trait t
      qtl_specs[[t]][which(is.na(qtl_specs[[t]]$qtl_name)), c("chr", "pos", "qtl_name")] <- 
        marker_sample_pos
      
      # Extract the df
      qtlmod_t <- qtl_specs[[t]]
      
      
      
      # If poss.corr is not NULL, adjust the effects to reflect that correlation
      # If not, do nothing
      if (!is.null(pos.corr)) {
      
      
        ## Adjust the effects of QTL for trait 2 + to achieve the desired correlation
        # A new vector of additive effects
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
              geno_q <- pull_genotype(genome = genome, geno = geno, loci = c(q$qtl_name, q$qtl1_pair))
              
              # What is the sign of the correlation?
              # If the sign is NA, set to 1
              cor_sign <- suppressWarnings(sign(cor(geno_q)[-1,1]))
              cor_sign <- ifelse(is.na(cor_sign), 1, cor_sign)
              
            }
            
            # What is the sign of the qtl_t add_eff?
            qtl_t_sign <- sign(q$add_eff)
            # What is the sign of the qtl1 pair?
            qtl1_sign <- sign(subset(qtl_specs[[1]], qtl_name == q$qtl1_pair, 
                                     add_eff, drop = TRUE))
            
            # If the intended correlation is positive...
            if (pos.corr) {
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
        
        # Edit the t-th model
        qtl_specs[[t]] <- qtlmod_t %>%
          mutate(chr = factor(chr, levels = chrnames(genome)),
                 add_eff = adj_add_eff)
        
      } # Close the pos.corr if statement
      
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
    
    } # Close if else statement
  
  # Add the genetic model to the genome
  genome[["gen_model"]] <- qtl_specs %>%
    lapply(arrange, chr, pos) %>%
    lapply(mutate, chr = factor(chr, levels = chrnames(genome)))
  
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
  
  # Return the genome
  return(genome)
  
} # Close the function


