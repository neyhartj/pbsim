#' Number of individuals in a population
#' 
#' @param pop An object of class \code{pop}.
#' 
#' @return 
#' Scalar number of individuals.
#' 
#' @export
#' 
nind <- function(pop) {
  
  # Make sure pop inherits the class "pop"
  if (!inherits(pop, "pop"))
    stop("The input 'pop' must be of class 'pop'.")
  
  # Does the pop object have genotypic values
  if (is.null(pop$geno_val))
    stop("The 'pop' object must have the data.frame of genotypic values")
  
  # Number of rows in the genotypic value matrix
  nrow(pop$geno_val)
  
} # Close function


#' Names of individuals
#' 
#' @param pop An object of class \code{pop}.
#' 
#' @return 
#' Character vector of individual names.
#' 
#' @export
#' 
indnames <- function(pop) {
  
  # Make sure pop inherits the class "pop"
  if (!inherits(pop, "pop"))
    stop("The input 'pop' must be of class 'pop'.")
  
  # Does the pop object have genotypic values
  if (is.null(pop$geno_val))
    stop("The 'pop' object must have the data.frame of genotypic values")
  
  # Return the individual names from the genotypic value df
  pop$geno_val$ind
  
} # Close function




#' Genotype a population
#' 
#' @param genome An object of class \code{genome}.
#' @param pop An object of class \code{pop}.
#' @param error.rate The genotyping error rate. This argument is not yet operational.
#' 
#' @return 
#' A matrix of dimensions \code{nind} x \code{nmarkers} in format z {-1, 0, 1}.
#' 
#' @export
#' 
genotype <- function(genome, pop, error.rate = 0) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Make sure there is a genetic model
  if (is.null(genome$gen_model))
    stop("No genetic model has been declared for the genome")
  
  # Make sure pop inherits the class "pop"
  if (!inherits(pop, "pop"))
    stop("The input 'pop' must be of class 'pop'.")
  
  # Get the names of the markers
  marker_names <- markernames(genome)
  
  # Combine the genos per chromosome
  geno <- do.call("cbind", pop$geno)
  
  # Subset only the markers and subtract 1
  subset(x = geno, select = marker_names, drop = FALSE) - 1
  
} # Close the function 



#' Subset a population object for specific individuals
#' 
#' @param pop An object of class \code{pop}.
#' @param individual A character of individuals to subset from the \code{pop}.
#' 
#' @details 
#' If \code{pheno_val} is present in the \code{pop}, the variance components are
#' dropped.
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
#' individual <- c("2ND27380", "06WA-406.6")
#' 
#' pop_subset <- subset_pop(pop = pop, individual = individual)
#' 
#' @import dplyr
#' 
#' @export
#' 
subset_pop <- function(pop, individual) {
  
  # Error handling
  # Make sure pop inherits the class "pop"
  if (!inherits(pop, "pop"))
    stop("The input 'pop' must be of class 'pop'.")

  # Convert individual to character
  individual <- as.character(individual)
    
  # Make sure the individuals specified are in the the pop
  if (!all(individual %in% indnames(pop)))
    stop("Not all of the individuals in 'individual' are in the 'pop' object.")
  
  # Find the element names
  element_names <- names(pop)

  
  # Empty pop object
  new_pop <- structure(vector("list", length(pop)), class = "pop", names = names(pop))
  
  # Subset various components
  new_pop$geno <- lapply(X = pop$geno, FUN = "[", individual, , drop = FALSE)
  new_pop$geno_val <- filter(pop$geno_val, ind %in% individual)

  # Subset phenotypic values, if present
  if ("pheno_val" %in% element_names) {
    # Get rid of the variance component estimate
    new_pop$pheno_val <- pop$pheno_val[-1]
    
    # Subset the phenotypic observations and pheno_mean
    new_pop$pheno_val <- lapply(X = new_pop$pheno_val, filter, ind %in% individual)
    
  }
  
  # Subset haploids, if present
  if ("haploids" %in% element_names) {
    
    new_pop$haploids <- lapply(pop$haploids, "[", ,,individual)
    
  }
  
  # Subset predicted genotypic values, if present
  if ("pred_val" %in% element_names) {
    
    new_pop$pred_val <- filter(pop$pred_val, ind %in% individual)
    
  }
    
  # Return the population
  return(new_pop)
  
} # Close the function



#' Combine a list of populations
#' 
#' @param pop_list A list of objects of class \code{pop}.
#' 
#' @details 
#' If \code{pheno_val} is present in the \code{pop}, the variance components are
#' dropped.
#' 
#' @import dplyr
#' @importFrom purrr pmap
#' @importFrom abind abind
#' 
#' @export
#' 
combine_pop <- function(pop_list) {
  
  # Make sure each element of 'pop_list' is a pop
  if (!all(sapply(X = pop_list, FUN = inherits, "pop")))
    stop("One of more of the elements in 'pop_list' is not a 'pop' object.")
  
  # Combine element names
  element_names <- lapply(pop_list, names) %>% 
    unlist() %>%
    unique()
  
  # Create a new pop object with elements present in any of the pops
  new_pop <- structure(vector("list", length(element_names)), class = "pop", 
                       names = element_names)
  
  # Combine genotypes
  # First extract the 'geno' element from each pop
  geno_list <- lapply(pop_list, "[[", "geno")
  # Combine
  new_pop$geno <- pmap(geno_list, rbind)
  
  # Combine genotypic values
  new_pop$geno_val <- do.call("rbind", lapply(pop_list, "[[", "geno_val"))

  # Combine phenotypic values if present
  if ("pheno_val" %in% element_names) {
    # Get rid of the variance component estimate
    new_pop$pheno_val <- structure(vector("list", 2), names = c("pheno_obs", "pheno_mean"))
    
    # Subset the 'pheno_val' element
    pheno_list <- lapply(X = pop_list, FUN = "[[", "pheno_val")
    
    # Combine
    new_pop$pheno_val$pheno_obs <- do.call("rbind", lapply(X = pheno_list, FUN = "[[", "pheno_obs"))
    new_pop$pheno_val$pheno_mean <- do.call("rbind", lapply(X = pheno_list, FUN = "[[", "pheno_mean"))
    
  }
  
  # Combine haploids if present
  if ("haploids" %in% element_names) {
    
    # Subset the 'haploids' element
    haploid_list <- lapply(X = pop_list, FUN = "[[", "haploids")
    
    # Empty array
    new_pop$haploids <- pmap(haploid_list, abind)
    
  }
  
  # Combine predicted genotypic values, if present
  if ("pred_val" %in% element_names) {
    
    new_pop$pred_val <- do.call("rbind", lapply(pop_list, "[[", "pred_val"))
    
  }
  
  # Return the new pop
  return(new_pop)
  
} # Close the function



#' Make selections from a population
#' 
#' @param pop An object of class \code{pop}.
#' @param intensity Either the prortion of individual in the population to select
#' or the number of individuals in the population to select.
#' @param type The type of selection to perform. If \code{"phenotypic"}, individuals
#' in the population are selected based on phenotypic values, if \code{"genomic"}, 
#' individuals in the population are selected based on predicted genotypic values,
#' and if \code{"random"}, individuals in the population are selected randomly.
#' @param tail The direction in which to select. If \code{"upper"}, the individuals
#' with the highest values are selected, and if \code{"lower"}, the individuals with
#' the lowest values are selected.
#' 
#' @details 
#' If multiple traits are present, "trait1" is used for selection.
#' 
#' @return 
#' An object of class \code{pop} that is a subset of the input \code{pop} for
#' the selections.
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
#' pop <- sim_phenoval(pop, h2 = 0.5)
#' 
#' pop_selected <- selection(pop = pop, intensity = 50)
#' 
#' @import dplyr
#' 
#' @export 
#' 
selection <- function(pop, intensity = 0.1, type = c("phenotypic", "genomic", "random"),
                      tail = c("upper", "tail")) {
  
  # Error handling
  # Make sure pop inherits the class "pop"
  if (!inherits(pop, "pop"))
    stop("The input 'pop' must be of class 'pop'.")
  
  # Match arguments
  type <- match.arg(type)
  tail <- match.arg(tail)
  
  # Number in the pheno.vec
  n_ind <- nind(pop)
  
  # If the sel.intensity is between 0 and 1, find the total number to keep
  if (intensity > 0 & intensity < 1) {
    # Number to select
    intensity_n <- round(nind(pop) * intensity)
  } else {
    intensity_n <- intensity
  }
  
  # If direction is "worst", make intensity_n negative
  intensity_n <- ifelse(tail == "low", intensity_n * -1, intensity_n)
  
  # If phenotypic selection
  if (type == "phenotypic") {
    
    # Check for genotypic values
    if (is.null(pop$pheno_val))
      stop("Phenotypic selection cannot proceed within phenotypes in the population.")
    
    # Select on the means
    selected <- pop$pheno_val$pheno_mean %>% 
      top_n(n = intensity_n, wt = trait1)
    
  } else if (type == "genomic") {
    
    # Check for PGVs
    if (is.null(pop$pred_val))
      stop("Genomic selection cannot proceed within predicted genotypic values in the population.")
  
    # Check for predicted values
    selected <- pop$pred_val %>% 
      top_n(n = intensity_n, wt = trait1)
    
  } else if (type == "random") {
    
    # Random selections
    selected <- pop$geno_val %>% 
      sample_n(size = abs(intensity_n))
    
  }
  
  # Use the individuals in the selection to subset the population and return
  subset_pop(pop = pop, individual = selected$ind)
  
} # Close the function
    






