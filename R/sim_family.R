#' Simulate a family from a pedigree
#' 
#' @param genome An object of class \code{genome}.
#' @param pedigree A \code{pedigree} detailing the scheme to develop the family.
#' Use \code{\link{sim_pedigree}} to generate.
#' @param founder.pop An object of class \code{pop} with the geno information for
#' the founders. Use the \code{\link{subset.pop}} function to subset a \code{pop}
#' object.
#' @param map.function The mapping function used to convert genetic distance to
#' recombination fractions. Can be one of "haldane", "kosambi", "c-f", or "morgan".
#' @param ignore.gen.model Logical - should the gene model be ignored?
#' 
#' @details 
#' Other arguments can be passed to internal functions in \code{sim_family}:
#' \describe{
#'   \item{\code{dh}}{A logical indicating if double-haploid lines should be created
#'   at the end of selfing. Doubled-haploids are generated at the last generation
#'   of selfing (e.g. if \emph{F_3} individuals are specified in the pedigree, DH
#'   lines are induced after the \emph{F_2}). \strong{Currently ignored.}}
#'   \item{\code{marker.gen}}{Numeric indicating the selfing generation at which to genotype
#'   individuals in the family. For instance, if \code{marker.gen = 2}, individuals
#'   are genotyped at the F_3 stage. \strong{Currently ignored.}}
#'   \item{\code{m}}{The crossover interference parameter. See 
#'   \code{\link[simcross]{sim_from_pedigree}} for more information. }
#'   \item{\code{p}}{The proportion of crosses from non-interference process.
#'   See \code{\link[simcross]{sim_from_pedigree}} for more information.}
#'   \item{\code{family.num}}{A integer designator of the family.}
#'   \item{\code{cycle.num}}{A integer designator of the breeding cycle.}
#' }
#' 
#' @return 
#' An object of class \code{pop}.
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
#' ## 2-way population
#' 
#' # Simulate the founder genotypes
#' founder.pop <- sim_founders(genome, pat.freq = c(0, 1))
#' 
#' # Create a pedigree with 100 individuals selfed to the F_3 generation
#' ped <- sim_pedigree(n.par = 2, n.ind = 100, n.selfgen = 2)
#' fam <- sim_family(genome = genome, pedigree = ped, founder.pop = founder.pop)
#' 
#' # Create an F_6 family, genotyped at the F_3
#' ped <- sim_pedigree(n.par = 2, n.ind = 100, n.selfgen = 5)
#' fam <- sim_family(genome = genome, pedigree = ped, founder.pop = founder.pop,
#'                   marker.gen = 2)
#' 
#' # Create a pedigree with 100 RIL individuals
#' ped <- sim_pedigree(n.par = 2, n.ind = 100, n.selfgen = Inf)
#' fam <- sim_family(genome = genome, pedigree = ped, founder.pop = founder.pop)
#' 
#' # Create a pedigree with 100 doubled-haploid individuals induced after the F_2
#' # generation.
#' ped <- sim_pedigree(n.par = 2, n.ind = 100, n.selfgen = 2)
#' fam <- sim_family(genome = genome, pedigree = ped, founder.pop = founder.pop, dh = TRUE)
#' 
#' 
#' ## 4-way population
#' 
#' # Simulate founders
#' founder.pop <- sim_founders(genome, n.str = 4, pat.freq = c(0, 0, 1))
#' 
#' #' # Create a pedigree with 100 individuals selfed to the F_3 generation
#' ped <- sim_pedigree(n.par = 4, n.ind = 100, n.selfgen = 2)
#' fam <- sim_family(genome = genome, pedigree = ped, founder.pop = founder.pop)
#' 
#' # Create a pedigree with 100 RIL individuals
#' ped <- sim_pedigree(n.par = 4, n.ind = 100, n.selfgen = Inf)
#' fam <- sim_family(genome = genome, pedigree = ped, founder.pop = founder.pop)
#' 
#' # Create a pedigree with 100 doubled-haploid individuals induced after the F_2
#' # generation.
#' ped <- sim_pedigree(n.par = 4, n.ind = 100, n.selfgen = 2)
#' fam <- sim_family(genome = genome, pedigree = ped, founder.pop = founder.pop, dh = TRUE)
#' 
#' 
#'                   
#' @importFrom simcross check_pedigree
#' @importFrom mpMap2 detailedPedigree simulateMPCross
#' @importFrom qtl sim.cross
#' 
#' @export
#' 
sim_family <- function(genome, pedigree, founder.pop, map.function = c("haldane","kosambi","c-f","morgan"),
                       ignore.gen.model = FALSE, ...) {
  
  # Error
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # founder.pop needs to be a pop object
  if (!inherits(founder.pop, "pop"))
    stop("The input 'founder.pop' must be of class 'pop'")
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = founder.pop$geno, ignore.gen.model = ignore.gen.model))
    stop("The geno did not pass. See warning for reason.")
  
  # Check the pedigree
  if (!check_pedigree(pedigree, ignore_sex = TRUE))
    stop("The pedigree is not formatted correctly.")
  
  # How many founders?
  n_founders <- nind(founder.pop)
  
  # Are the number of founders correct vis a vis the pedigree?
  if (n_founders != sum(pedigree$gen == 0))
    stop("The number of founders in the input 'founders' is not equal to the
         number of founders in the 'pedigree.'")
  
  # Parse other arguments
  other.args <- list(...)
  map.function <- match.arg(map.function)
  
  ####
  # For doubled-haploids
  dh <- ifelse(is.null(other.args$dh), FALSE, other.args$dh)
  ## DH capability is not working; add warning
  if (dh) warning("Doubled-haploid generation is currently not supported. The 'dh' argument is ignored.")
  ####
  
  
  # Marker genotyping generations
  # Set to the max pedigree generation
  marker_gen <- other.args$marker.gen
  # Marker gen 
  if (!is.null(marker_gen)) {
    if (marker_gen <= 0) stop ("marker.gen cannot be less than or equal to 0.")
  } else {
    marker_gen <- max(pedigree$gen)
  }
  
  # Extract the individual ids of the finals
  final_id <- subset(pedigree, gen == max(gen))$id
  
  # Generate new progeny names
  n_ind <- length(final_id) # Number of individuals
  # Find the maximum generation in the pedigree
  max_gen <- max(pedigree$gen)
  
  # For naming
  family_num <- ifelse(is.null(other.args$family.num), 1, other.args$family.num)
  cycle_num <- ifelse(is.null(other.args$cycle.num), 1, other.args$cycle.num)
  
  # New names
  new_names <- paste0("C", cycle_num, "_", max_gen, formatC(family_num, width = 3, flag = "0"), 
                      "-", formatC(seq_len(n_ind), width = 3, flag = "0"))
  
  
  # Determine the selfing level
  selfing <- attr(pedigree, "selfing")
  
  # Extract the map
  map <- genome$map
  
  # For xodata
  m <- ifelse(is.null(other.args$m), 0, other.args$m)
  p <- ifelse(is.null(other.args$p), 1, other.args$p)
  
  # Extract founder genotypes
  founder_geno <- do.call("cbind", founder.pop$geno)
  
  
  if (selfing == "finite") {
    
    ## Control flow for marker_gen
    if ( marker_gen == max(pedigree$gen) ) {
    
      ## Option 1
      # Turn the pedigree into a detailedPedigree object
      det_ped <- detailedPedigree(lineNames = as.character(pedigree$id), mother = pedigree$mom, 
                                  father = pedigree$dad, initial = seq_len(n_founders), 
                                  observed = pedigree$observed, selfing = selfing)
      
      # Simulate using mpMap2
      cross_sim <- simulateMPCross(map = map, pedigree = det_ped, mapFunction = map.function)
      
      # Extract the progeny genotypes
      prog_multipoint <- cross_sim@geneticData@.Data[[1]]@finals
      prog_multipoint[prog_multipoint > n_founders] <- NA
      
      ## Convert the progeny genotypes to the parental states
      prog_geno <- t(apply(X = prog_multipoint, MARGIN = 1, FUN = function(prog) {
        founder_geno[cbind(prog, seq(ncol(founder_geno)))] }))
      # Replace NA with het
      prog_geno[is.na(prog_geno)] <- 1
      # Filler
      prog_marker_geno <- NULL
      
    } else if ( marker_gen < max(pedigree$gen) ) {
      
      stop("Marker generations less than the inbreeding generations are currently not supported.")
    
      # # Create a founder pop from the parents of the family
      # founderPop <- newMapPop(genMap = structure(genome$map, class = "numeric"), haplotypes = founder.pop$geno,
      #                         inbred = TRUE, ploidy = 2)
      # # Simulation parameter object
      # simParamTemp <- SimParam$new(founderPop = founderPop)
      # parentPop <- newPop(rawPop = founderPop, simParam = simParamTemp)
      # 
      # 
      # ## Split by number of founders
      # ## 2
      # if (n_founders == 2) {
      # 
      #   # Create the cross - F1
      #   pop1 <- makeCross(pop = parentPop, crossPlan = cbind(1,2), nProgeny = n_ind, simParam = simParamTemp)  
      #   
      # ## 4
      # } else {
      #   
      #   # Create the first cross, then the second
      #   cross1 <- makeCross(pop = parentPop, crossPlan = cbind(1,2), nProgeny = 1, simParam = simParamTemp)
      #   cross2 <- makeCross(pop = parentPop, crossPlan = cbind(3,4), nProgeny = 1, simParam = simParamTemp)
      #   # F1
      #   pop1 <- makeCross2(females = cross1, males = cross2, nProgeny = n_ind, crossPlan = cbind(1,1),
      #                      simParam = simParamTemp)
      #   
      # }
      # 
      # # Inbreed to marker generation
      # for (s in seq(2, marker_gen + 1)) pop1 <- self(pop = pop1, nProgeny = 1, simParam = simParamTemp)
      # # Genotype
      # prog_marker_geno <- pullSegSiteGeno(pop = pop1, simParam = simParamTemp)
      # # continue to inbreed
      # for (s in seq(s+1, max(pedigree$gen))) pop1 <- self(pop = pop1, nProgeny = 1, simParam = simParamTemp)
      # 
      # # Get segregating sites
      # prog_geno <- pullSegSiteGeno(pop = pop1, simParam = simParamTemp)
      # dimnames(prog_marker_geno) <- list(new_names, markernames(genome, include.qtl = TRUE))
      
      
      
    } else {
      stop ("marker.gen cannot be greater than the inbreeding generation.")
      
    }
    
    
    # Else infinite inbreeding
  } else {
    
    # Otherwise use sim.cross from qtl
    if (n_founders == 4) {
      
      # Create a list of founder genotypes - all alleles are unique (1-4)
      founderGenos <- lapply(map, function(chr) matrix(data = 1:4, nrow = length(chr), ncol = 4, byrow = TRUE))
      
      cross_sim <- sim.cross(map = genome$map, n.ind = length(final_id), type = "ri4self", 
                             map.function = map.function, random.cross = FALSE,
                             founderGeno = founderGenos)
      
      # Extract progeny genos
      prog_multipoint <- do.call("cbind", lapply(X = cross_sim$geno, FUN = "[[", "data"))
      prog_multipoint[prog_multipoint == 4] <- 3
      prog_multipoint[prog_multipoint == 8] <- 4

      
    } else if (n_founders == 2) {
      
      cross_sim <- sim.cross(map = genome$map, n.ind = length(final_id), type = "riself", 
                             map.function = map.function)
      
      prog_multipoint <- do.call("cbind", lapply(X = cross_sim$geno, FUN = "[[", "data"))
      
    }
    
    ## Convert the progeny genotypes to the parental states
    prog_geno <- t(apply(X = prog_multipoint, MARGIN = 1, FUN = function(prog) {
      founder_geno[cbind(prog, seq(ncol(founder_geno)))] }))
    # Replace NA with het
    prog_geno[is.na(prog_geno)] <- 1
    # Filler
    prog_marker_geno <- NULL
    
  }
  
  
  dimnames(prog_geno) <- list(new_names, markernames(genome, include.qtl = TRUE))
    
  # Create the pop
  pop_out <- create_pop(genome = genome, geno = prog_geno, ignore.gen.model = ignore.gen.model)
  # Add markers, if available
  pop_out$marker_geno <- prog_marker_geno
  return(pop_out)
  
} # Close the function




#' Simulate multiple families from a crossing block
#' 
#' @param genome An object of class \code{genome}.
#' @param pedigree A \code{pedigree} detailing the scheme to develop the family.
#' Use \code{\link{sim_pedigree}} to generate. May be a single pedigree object 
#' or (in the case of different family sizes) a list of length \code{nrow(crossing.block)}.
#' @param founder.pop An object of class \code{pop} with the geno information for
#' the parents. Additional individuals can be present, but they
#' will be filtered according to the parents in the \code{crossing.block}.
#' @param crossing.block A crossing block detailing the crosses to make. Must be a
#' \code{data.frame} with 2 columns (for 2-way crosses) or 4 columns (for 4-way crosses).
#' See \code{\link{sim_crossing.block}}.
#' @param ... Additional arguments passed to \code{\link{sim_family}}.
#' 
#' @return 
#' An object of class \code{pop} with the information for all individuals in 
#' the families specified by the crossing block.
#' 
#' @examples 
#' 
#' # Simulate a genome
#' n.mar  <- rep(500, 7)
#' len <- runif(n = 7, min = 140, 170)
#' 
#' genome <- sim_genome(len, n.mar)
#' 
#' # Simulate a quantitative trait influenced by 200 QTL
#' L <- 200
#' qtl.model <- matrix(NA, L, 4)
#' genome <- sim_gen_model(genome = genome, qtl.model = qtl.model, 
#'                         add.dist = "geometric", max.qtl = L)
#' 
#' # Simulate the genotypes for 8 founders
#' founder.pop <- sim_founders(genome = genome, n.str = 8)
#' 
#' ## 2-way population
#' 
#' # Generate a crossing block with 5 crosses
#' cb <- sim_crossing_block(parents = indnames(founder.pop), n.crosses = 25)
#' 
#' # Create a bi-parental pedigree with 100 individuals selfed to the F_3 generation
#' ped <- sim_pedigree(n.par = 2, n.ind = 100, n.selfgen = 2)
#' 
#' # Simulate a group of families from the crossing block
#' fam_cb <- sim_family_cb(genome = genome, pedigree = ped, founder.pop = founder.pop, 
#'                         crossing.block = cb)
#'                         
#' # Create a bi-parental pedigree with 100 individuals selfed to the F_3 generation
#' ped <- sim_pedigree(n.par = 2, n.ind = 100, n.selfgen = 5)
#'                         
#' # Create many F_6 families, genotyped at the F_3
#' fam_cb <- sim_family_cb(genome = genome, pedigree = ped, founder.pop = founder.pop, 
#'                         crossing.block = cb, marker.gen = 2)
#' 
#'                         
#' ## 4-way population
#' 
#' # Generate a crossing block with 5 crosses
#' cb <- sim_crossing_block(parents = indnames(founder.pop), n.crosses = 5, type = "4way")
#' 
#' # Create a pedigree with 100 individuals selfed to the F_3 generation
#' ped <- sim_pedigree(n.par = 4, n.ind = 100, n.selfgen = 2)
#' 
#' # Simulate a group of families from the crossing block
#' fam_cb <- sim_family_cb(genome = genome, pedigree = ped, founder.pop = founder.pop, 
#'                         crossing.block = cb)
#'                         
#' ## Simulate RIL familes of different sizes
#' 
#' # Generate a crossing block with 5 crosses
#' cb <- sim_crossing_block(parents = indnames(founder.pop), n.crosses = 5, type = "2way")
#' # Number of individuals per cross
#' nIndCross <- c(20, 30, 40, 50, 60)
#' 
#' # Create a list of pedigrees
#' pedList <- sapply(nIndCross, sim_pedigree, n.par = 2, simplify = FALSE)
#' 
#' # Simulate a group of families from the crossing block
#' fam_cb <- sim_family_cb(genome = genome, pedigree = pedList, founder.pop = founder.pop, 
#'                         crossing.block = cb)
#'
#'
#'
#' @importFrom simcross check_pedigree
#' @import dplyr
#' 
#' @export
#' 
sim_family_cb <- function(genome, pedigree, founder.pop, crossing.block, ...) {
  
  # Error handling
  if (!inherits(genome, "genome")) stop("The input 'genome' must be of class 'genome.'")
  
  # founder.pop needs to be a pop object
  if (!inherits(founder.pop, "pop")) stop("The input 'founder.pop' must be of class 'pop'")
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = founder.pop$geno)) stop("The geno did not pass. See warning for reason.")
  
  # If pedigree is a list of pedigrees, make sure it is
  # the same length as nrows of the crossing.block
  if (inherits(pedigree, "list")) {
    stopifnot(length(pedigree) == nrow(crossing.block))
    # check each pedigree
    if (any(!sapply(pedigree, simcross::check_pedigree, ignore_sex = TRUE))) stop("One or more pedigrees are not formatted correctly.")
    
  } else {
    # Check the pedigree
    if (!simcross::check_pedigree(pedigree, ignore_sex = TRUE)) stop("The pedigree is not formatted correctly.")
    
  }
  
  # Are all of the parents in the crossing block in the founder.pop?
  if (!all(unique(unlist(crossing.block)) %in% indnames(founder.pop))) {
    stop("Not all of the parents in the crossing block are in the 'founder.pop'.")
  }
  
  fam_cb <- vector("list", nrow(crossing.block))
  
  # Separate flow by whether pedigree is a list
  if (inherits(pedigree, "list")) {
    
    # Seq along the crossing block
    for (i in seq_along(fam_cb)) {
      
      parents <- unlist(crossing.block[i,])
      founder_geno <- subset_pop(pop = founder.pop, individual = parents)
      fam_cb[[i]] <- sim_family(genome = genome, pedigree = pedigree[[i]], founder.pop = founder_geno, 
                                family.num = i, ... = ...)
      
    }
    
    
  } else {
  
    # Seq along the crossing block
    for (i in seq_along(fam_cb)) {
      
      parents <- unlist(crossing.block[i,])
      founder_geno <- subset_pop(pop = founder.pop, individual = parents)
      fam_cb[[i]] <- sim_family(genome = genome, pedigree = pedigree, founder.pop = founder_geno, 
                                family.num = i, ... = ...)
      
    }
    
  }
    

  # Combine the populations and return
  combine_pop(pop_list = fam_cb)

} # Close the function
  
  
  
  
  


