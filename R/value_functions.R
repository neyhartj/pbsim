#' Extract values from a population
#' 
#' @description 
#' Functions to extract phenotypic values, genotypic values, or predicted
#' genotypic values.
#' 
#' @param pop A 'pop' object.
#' @param means For phenotypic data, should phenotypic means be returned, or
#' all observations?
#' 
#' @details 
#' 
#' \describe{
#'   \item{pheno}{Extract phenotypic values.}
#'   \item{gv}{Extract genotypic values.}
#'   \item{pgv}{Extract predicted genotypic values.}
#' }
#' 
#' @examples 
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
#' # Simulate the population
#' pop <- sim_pop(genome = genome, n.ind = 100)
#' # Generate phenotypes
#' pop <- sim_phenoval(pop = pop, h2 = 0.5)
#' 
#' 
#' ## Extract genotypic values
#' gv(pop)
#' 
#' ## Extract phenotypic values
#' pheno(pop)
#' 
#' ## Try to get predicted phenotypic values. There are none, so this returns \code{NULL}.
#' pgv(pop)
#' 
#' # Predict genotypic values
#' pop <- pred_geno_val(genome = genome, training.pop = pop, candidate.pop = pop)
#' # Now extract predicted genotypic values
#' pgv(pop)
#' 
#' @export
#' 
pheno <- function(pop, means = TRUE) {
  # Error
  stopifnot(class(pop) == "pop")
  
  if (means) {
    pop$pheno_val$pheno_mean
  } else{
    pop$pheno_val$pheno_obs
  }
  
}

#' @describeIn pheno
#' 
#' @export
#' 
gv <- function(pop) {
  # Error
  stopifnot(class(pop) == "pop")
  pop$geno_val
}

#' @describeIn pheno
#' 
#' @export
#' 
pgv <- function(pop) {
  # Error
  stopifnot(class(pop) == "pop")
  pop$pred_val
}


