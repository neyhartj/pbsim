#' Define the genetic model of a trait
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
#' randomly assigned based on the number of rows in the matrix. See \code{Details}
#' for more information.
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
#' If the \code{qtl.model} imput is completely NA, the following rules apply:
#' \itemize{
#'   \item{QTL positions are randomly drawn, with no regard to uniformity
#'   over chromosomes.}
#'   \item{The distribution of additive effects must be provided. This can
#'   be included as an additional argument: \code{add.dist = "normal"} or
#'   \code{add.dist = "geometric"}. For a distribution of \code{"normal"},
#'   additive effects are generated via the \code{\link[stats]{rnorm}} function.
#'   For a distribution of \code{"geometric"}, additive effects are calculated for
#'   the k-th QTL as \eqn{a^k} where \eqn{a = (1 - L) / (1 + L)} and \eqn{L} is 
#'   the number of QTL (Lande and Thompson, 1990).}
#'   \item{Dominance effects are assumed to be zero, but can be simulated
#'   with the optional argument: \code{dom.dist = "normal"}. Currently, only
#'   normal distributions are supported.}
#'   }
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
#' index <- c(100, 200, 150, 300, 50, 400)
#' a <- c(1, 0.25, 0.5, 0.25, 0.25, 0.5)
#' d <- 0
#' 
#' qtl.model <- cbind(chromosome, index, a, d)
#' 
#' genome <- sim_gen_model(genome, qtl.model)
#' 
#' @export
#' 
sim_gen_model <- function(genome, qtl.model, ...) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Make sure qtl.model has four columns
  if (ncol(qtl.model) < 4)
    stop("The input 'qtl.model' must have at least 4 columns.")
  
  # Extract other arguments
  other.args <- list(...)
  
  # Is the qtl.model a NA model?
  if (all(is.na(qtl.model))) {
    random.qtl <- TRUE
    
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
    
    # Otherwise check the matrix
  } else {
    random.qtl <- FALSE
    
    # The qtl model must not have any NAs
    if (any(is.na(qtl.model)))
      stop("Unless the 'qtl.model' is completely NA, there can be no NA elements.")
    
    # Are the chromosomes in the genome?
    n.chr <- length(genome$len)
    
    if (!all(qtl.model[,1] %in% seq(n.chr)))
      stop("The chromosome numbers in 'qtl.model' are not chromosomes in the 
           genome.")
    
    # Are the marker indices correct?
    qtl.ind <- split(qtl.model[,2], qtl.model[,1])
    
    if (!all(mapply(qtl.ind, genome$n.mar, FUN = function(q, n) q <= n)))
      stop("The QTL indices in 'qtl.model' are not within the number of markers
            on the chromosomes.")
    
    # Are the additive effects positive?
    if (!all(qtl.model[,3] >= 0))
      stop("The additive effects in 'qtl.model' must be positive.")
    
  }
  
  # Number of qtl
  n.qtl <- nrow(qtl.model)
  
  # Number of chromosomes in the genome
  n.chr <- length(genome$len)
  
  # Should QTL be randomized?
  if (random.qtl) {
    
    # Simulate additive effects
    if (add.dist == "normal") 
      add.eff <- rnorm(n.qtl)
    
    if (add.dist == "geometric") {
      a <- (1 - n.qtl) / (1 + n.qtl)
      add.eff <- sample(abs( a ^ (seq(n.qtl)) ))
      # Randomly assign favorable or unfavorable for the first allele
      add.eff <- add.eff * sample(x = c(-1, 1), size = n.qtl, replace = TRUE)
    }
    
    # Simulate dominance effect 
    if (is.null(dom.dist)) {
      dom.eff <- 0
    
    } else {
      dom.eff <- rnorm(n.qtl)
    }
    
    # Randomly draw chromosomes
    chr <- sort(sample(seq_len(n.chr), n.qtl, replace = TRUE))
    
    # Randomly draw SNPs to become QTL.
    qtl.ind <- lapply(X = split(chr, chr), FUN = function(chr.list) { 
      chr.num <- unique(chr.list)
      sort(sample(seq_len(genome$n.mar[chr.num]), length(chr.list))) })

    qtl.ind <- unlist(qtl.ind)
    
    # Assemble the matrix
    qtl.model <- cbind(chr, qtl.ind, add.eff, dom.eff)

  }
  
  # Convert the matrix to a data.frame
  qtl_model_df <- as.data.frame(qtl.model)
    
  # Extract cM positions for the QTL
  qtl.ind <- split(qtl.model[,2], qtl.model[,1])
  
  # Extract the map only for chromosomes in which QTL were designated
  gen_map <- genome$map[as.numeric(unique(names(qtl.ind)))]
  
  qtl.pos <- mapply(qtl.ind, gen_map, FUN = function(i, m) 
    data.frame(snp = names(m[i]), pos = m[i], row.names = NULL), SIMPLIFY = FALSE)
  qtl.pos <- do.call("rbind", qtl.pos)
  
  # Add the position to the qtl model matrix and rearrange
  qtl_model <- cbind(qtl_model_df, qtl.pos)[,c(1, 5:6, 2:4)]
  
  # Change row.names names
  row.names(qtl_model) <- NULL
  
  # Add the genetic model to the genome
  genome[["gen_model"]] <- qtl_model
  
  return(genome)
    
} # Close the function
