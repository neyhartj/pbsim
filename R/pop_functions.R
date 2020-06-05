#' Create a population object
#' 
#' @description 
#' Assembles genotype data and into a \code{pop} object.
#' 
#' @param genome An object of class \code{genome}.
#' @param geno Genotype data on a population to phenotype. If the genome type is 
#' \code{"pbsim"}, must be a matrix of dimensions \code{n.ind} x \code{n.loci}, 
#' the elements of which must be z {0, 1, 2}, or a list of such matrices. If the 
#' genome type is \code{"hypred"}, must be an array of dimensions \code{2} x 
#' \code{n.loci} x \code{n.ind}, the elements of which must be z {0, 1}.
#' @param ignore.gen.model Logical - should any gene model be ignored when creating
#' the population? Use this to force a population without a gene model.
#' 
#' 
#' @details 
#' A \code{pop} is similar to a \code{cross} object in \code{\link[qtl]{qtl-package}} 
#' (see \code{\link[qtl]{read.cross}}). The \code{pop} object stores information on the
#' genome, the genotypes at genetic markers, and phenotypes. The \code{pop} object is
#' meant to be a bit more flexible, without the pedigree or family structure required in
#' a \code{cross} object.
#' 
#' @return 
#' An object of class \code{pop} with genotype information for the individuals in
#' that population and the genotypic value of those individuals.
#' 
#' The genotypic value of individuals is calculcated as the sum of the QTL effects
#' carried by each individual. The genetic variance is calculated as the variance
#' of these genotypic values (\eqn{V_G = var(g)}).
#' 
#' @examples 
#' 
#' \dontrun{
#' 
#' # Use data from the PopVar package
#' library(PopVar)
#' data("think_barley")
#' 
#' # Format the map correctly and simulate a genome
#' map_in <- map.in_ex[,-1]
#' row.names(map_in) <- map.in_ex$mkr
#' genome <- sim_genome(map = table_to_map(map_in))
#' 
#' genos <- apply(X = G.in_ex[-1,-1], MARGIN = 2, FUN = as.numeric)
#' dimnames(genos) <- list(as.character(G.in_ex$V1[-1]), as.character(unlist(G.in_ex[1,-1])))
#' 
#' # Impute with the rounded mean
#' genos1 <- apply(X = genos, MARGIN = 2, FUN = function(snp) {
#'   mean <- ifelse(mean(snp, na.rm = T) < 0, -1, 1)
#'   snp[is.na(snp)] <- mean
#'   return(snp) })
#' 
#' ## Create a population without a genetic model
#' pop <- create_pop(genome = genome, geno = genos1 + 1, ignore.gen.model = T)
#' 
#' ## Create a genetic model with 15 QTL
#' qtl.model <- matrix(NA, ncol = 4, nrow = 15)
#' genome <- sim_gen_model(genome = genome, qtl.model = qtl.model, add.dist = "geometric")
#' 
#' pop <- create_pop(genome = genome, geno = genos1 + 1)
#' 
#' }
#' 
#' 
#' @importFrom abind abind
#' 
#' @export
#' 
create_pop <- function(genome, geno, ignore.gen.model = FALSE) {
  
  ## Error handling
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = geno, ignore.gen.model = ignore.gen.model))
    stop("The geno did not pass. See warning for reason.")
  
  
  # Create empty pop list
  pop <- structure(vector("list", 2), class = "pop", names = c("geno", "geno_val"))
  
  
  # Sort the genos on individual names
  geno <- geno[order(row.names(geno)),]
  
  # Split the geno matrix into chromosomes
  geno_split <- split_geno(genome = genome, geno = geno, ignore.gen.model = ignore.gen.model)
  
  # Calculate the genotypic value - ignore if told so
  if (ignore.gen.model) {
    geno_val <- NULL
  } else {
    geno_val <- calc_genoval(genome = genome, geno = geno_split)
  }
  
  
  # Add data to the pop
  pop[["geno"]] <- geno_split
  pop[["geno_val"]] <- geno_val
  
  # Return
  return(pop)
  
} # Close the function




#' Simulate a population
#' 
#' @description
#' Simulates a random population of given population size. SNPs are simulated
#' to be in linkage equilibrium.
#' 
#' @param genome A genome object.
#' @param n.ind The number of individual in the population.
#' @param ignore.gen.model Logical - should the gene model be ignored?
#' 
#' @return 
#' A \code{pop} object.
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
#' 
#' @export
#' 
sim_pop <- function(genome, n.ind, ignore.gen.model = FALSE) {
  
  # Check inputs
  stopifnot(class(genome) == "genome")
  stopifnot(class(n.ind) == "numeric")
  n.ind <- as.integer(n.ind)
  
  # Map
  map <- genome$map
  
  
  # Create individual names
  ind_names <- paste("ind", formatC(x = seq(n.ind), width = nchar(n.ind), flag = 0, format = "d"), sep = "")
  
  # Simulate markers
  geno <- lapply(X = map, FUN = function(m) {
    mar_chr <- length(m)
    M <- replicate(n = mar_chr, ifelse(runif(n = n.ind) < 0.5, 0, 2))
    `dimnames<-`(M, list(ind_names, names(m))) })
  
  geno1 <- do.call("cbind", geno)
  
  # Create the population
  create_pop(genome = genome, geno = geno1, ignore.gen.model = ignore.gen.model)
  
}












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
  
  # Determine the individuals from the genotype data
  nind_count <- sapply(X = pop$geno, nrow)
  # Are these all equal?
  if (length(unique(nind_count)) > 1)
    stop ("The number of individuals in the population is uneven. Please check.")

  # Number of rows in the genotype matrix
  unique(nind_count)
  
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
  
  # Determine the individuals from the genotype data
  nind_count <- sapply(X = pop$geno, nrow)
  # Are these all equal?
  if (length(unique(nind_count)) > 1)
    stop ("The number of individuals in the population is uneven. Please check.")
  
  # Return the row names
  row.names(pop$geno[[1]])
  
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
  
  ## If the pop contains an object called "marker_geno," export those instead of the
  ## actual loci genotypes
  if (!is.null(pop$marker_geno)) {
    geno <- pop$marker_geno
    
  } else {
    # Combine the genos per chromosome
    geno <- do.call("cbind", pop$geno)
    
  }
  

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
#' 
#' # Subset the population
#' subset_pop(pop, c("ind001", "ind100"))
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
  new_pop$geno_val <- subset(pop$geno_val, ind %in% individual)

  # Subset phenotypic values, if present
  if ("pheno_val" %in% element_names) {
    # Get rid of the variance component estimate
    new_pop$pheno_val <- pop$pheno_val[c("pheno_obs", "pheno_mean")]
    
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
#' # Simulate two populations
#' pop1 <- sim_pop(genome = genome, n.ind = 100)
#' pop2 <- sim_pop(genome = genome, n.ind = 100)
#' 
#' # Subset the populations
#' pop1 <- subset_pop(pop1, indnames(pop1)[1:5])
#' pop2 <- subset_pop(pop2, indnames(pop2)[20:25])
#' 
#' # Combine
#' combine_pop(list(pop1, pop2))
#' 
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
  element_names <- unique(unlist(lapply(pop_list, names)))
  
  ## Get the individual names from each population
  pop_list_indiv <- lapply(X = pop_list, FUN = indnames)
  
  # Create a new pop object with elements present in any of the pops
  new_pop <- structure(vector("list", length(element_names)), class = "pop", 
                       names = element_names)
  
  # Combine genotypes
  # First extract the 'geno' element from each pop
  geno_list <- lapply(pop_list, "[[", "geno")
  # Combine
  new_pop$geno <- pmap(geno_list, rbind)
  
  # Combine pedigree information
  new_pop$pedigree <- do.call("rbind", lapply(pop_list, "[[", "pedigree"))
  
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
    
    # Add factor levels to the phenotype means; missing levels will become NA; remove factors
    # Subset for easier use
    pheno_mean_new <- new_pop$pheno_val$pheno_mean
    pheno_mean_new$ind <- factor(pheno_mean_new$ind, levels = unlist(pop_list_indiv))
    
    # Make missing cases explicit
    pheno_mean_new1 <- as.data.frame(complete(data = pheno_mean_new, ind))
    pheno_mean_new1$ind <- as.character(pheno_mean_new1$ind)
    
    # Replace
    new_pop$pheno_val$pheno_mean <- pheno_mean_new1
    
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
  
  # Combine marker genotypess, if present
  if ("marker_geno" %in% element_names) {
    
    new_pop$marker_geno <- do.call("rbind", lapply(pop_list, "[[", "marker_geno"))
    
  }
  
  # Return the new pop
  return(new_pop)
  
} # Close the function



#' Make selections from a population
#' 
#' @param pop An object of class \code{pop}.
#' @param intensity Either the prortion of individual in the population to select
#' or the number of individuals in the population to select.
#' @param index The coefficients for the selection index. Positive coefficients 
#' equate to selection on higher trait values, and negative coefficients equate 
#' to selection on lower trait value. Must be a vector if length \code{n_trait}.
#' If one trait is present, the coefficient is 1 or -1.
#' @param type The type of selection to perform. If \code{"phenotypic"}, individuals
#' in the population are selected based on phenotypic values, if \code{"genomic"}, 
#' individuals in the population are selected based on predicted genotypic values,
#' and if \code{"random"}, individuals in the population are selected randomly.
#' 
#' @details 
#' If one trait is present, selection is performed on that one trait.
#' 
#' If two traits are present, an index is calculated using the index.
#' 
#' If there is a tie in the phenotypic or predicted genotypic values, individuals
#' are randomly chosen.
#' 
#' @return 
#' An object of class \code{pop} that is a subset of the input \code{pop} for
#' the selections.
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
#' # Simulate the population
#' pop <- sim_pop(genome = genome, n.ind = 100)
#' # Select randomly
#' select_pop(pop = pop, intensity = 0.1, type = "random")
#' 
#' # Phenotype the population
#' pop <- sim_phenoval(pop = pop, h2 = 0.5)
#' # Select on phenotypes
#' select_pop(pop = pop, intensity = 0.1, type = "phenotypic")
#' 
#' # Predict genotypic values
#' pop <- pred_geno_val(genome = genome, training.pop = pop, candidate.pop = pop)
#' # Select on pgvs
#' select_pop(pop = pop, intensity = 0.1, type = "genomic")
#' 
#' 
#' @import dplyr
#' 
#' @export 
#' 
select_pop <- function(pop, intensity = 0.1, index = 1, type = c("phenotypic", "genomic", "random")) {
  
  # Error handling
  # Make sure pop inherits the class "pop"
  if (!inherits(pop, "pop"))
    stop("The input 'pop' must be of class 'pop'.")
  
  # Match arguments
  type <- match.arg(type)

  # Number in the pheno.vec
  n_ind <- nind(pop)
  
  # If the sel.intensity is between 0 and 1, find the total number to keep
  if (intensity > 0 & intensity < 1) {
    # Number to select
    intensity_n <- round(nind(pop) * intensity)
  } else {
    intensity_n <- intensity
  }
  
  # How many traits
  n_trait <- ncol(pop$geno_val) - 1
  
  # Make sure that type is of appropriate length
  if (length(index) != n_trait)
    stop("The number of elements in the input 'type' must be the same as the
         number of trait.")
  
  # Rescale the index
  index <- scale(index, scale = abs(sum(index)), center = FALSE)
  
  # If phenotypic selection
  if (type == "phenotypic") {
    
    # Check for genotypic values
    if (is.null(pop$pheno_val))
      stop("Phenotypic selection cannot proceed without phenotypic information on the population.")
    
    # Recode the value
    selected <- pop$pheno_val$pheno_mean
    selected[,-1] <- selected[,-1, drop = FALSE] * matrix(index, nrow = n_ind, ncol = n_trait, byrow = T)
    
    # Calculate an index and take the top n
    selected <- selected %>% 
      mutate(index = rowSums(select(., -1))) %>% 
      top_n(n = intensity_n, wt = index)
    
    # Is the df greater than the number of intended selections?
    if (nrow(selected) > intensity_n) {
      
      # Separate those selections with the greatest value
      top_selected <- selected %>% 
        filter(index != min(index))
      
      # How many?
      n_top <- nrow(top_selected)
      
      # Separate those selections with the lowest value
      bot_selected <- selected %>%
        filter(index == min(index))
      
      # Sample among the bottom randomly to bring the number of selections up to
      # the intended number
      bot_selected_sample <- bot_selected %>% 
        sample_n(size = intensity_n - n_top)
      
      # Bind rows and sort
      selected <- bind_rows(top_selected, bot_selected_sample) %>% 
        arrange(ind)
      
    }
    
  } else if (type == "genomic") {
    
    # Check for PGVs
    if (is.null(pop$pred_val))
      stop("Genomic selection cannot proceed withouit predicted genotypic values for the population.")
      
    # Recode the value
    selected <- pop$pred_val
    selected[,-1] <- selected[,-1, drop = FALSE] * matrix(index, nrow = n_ind, ncol = n_trait, byrow = T)
    
    # Calculate an index, select the top_n, then sort on the index
    selected <- selected %>% 
      mutate(index = rowSums(select(., -1))) %>% 
      top_n(n = intensity_n, wt = index) %>%
      arrange(desc(index))
    
    # Is the df greater than the number of intended selections?
    if (nrow(selected) > intensity_n) {
      
      # Separate those selections with the greatest value
      top_selected <- selected %>% 
        filter(index != min(index))
      
      # How many?
      n_top <- nrow(top_selected)
      
      # Separate those selections with the lowest value
      bot_selected <- selected %>%
        filter(index == min(index))
      
      # Sample among the bottom randomly to bring the number of selections up to
      # the intended number
      bot_selected_sample <- bot_selected %>% 
        sample_n(size = intensity_n - n_top)
      
      # Bind rows and sort
      selected <- bind_rows(top_selected, bot_selected_sample) %>% 
        arrange(ind)
      
    }
    
  } else if (type == "random") {
    
    # Random selections
    selected <- pop$geno_val %>% 
      sample_n(size = intensity_n)
    
  }
  
  # Use the individuals in the selection to subset the population and return
  subset_pop(pop = pop, individual = selected$ind)
  
} # Close the function
    






