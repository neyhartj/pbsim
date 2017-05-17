#' Simulate a family from a pedigree
#' 
#' @param genome An object of class \code{genome}.
#' @param pedigree A \code{pedigree} detailing the scheme to develop the family.
#' Use \code{\link{sim_pedigree}} to generate.
#' @param founder.pop An object of class \code{pop} with the geno information for
#' the founders. Use the \code{\link{subset.pop}} function to subset a \code{pop}
#' object.
#' 
#' @details 
#' Other arguments can be passed to internal functions in \code{sim_family}:
#' \describe{
#'   \item{\code{dh}}{A logical indicating if double-haploid lines should be created
#'   at the end of selfing. Doubled-haploids are generated at the last generation
#'   of selfing (e.g. if \emph{F_3} individuals are specified in the pedigree, DH
#'   lines are induced after the \emph{F_2}).}
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
#' # Simulate the founder genotypes
#' founder.pop <- sim_founders(genome)
#' 
#' # Create a pedigree with 100 individuals selfed to the F_3 generation
#' ped <- sim_pedigree(n.ind = 100, n.selfgen = 2)
#' 
#' fam <- sim_family(genome = genome, pedigree = ped, founder.pop = founder.pop)
#' 
#' # Create a pedigree with 100 RIL individuals
#' ped <- sim_pedigree(n.ind = 100, n.selfgen = Inf)
#' 
#' fam <- sim_family(genome = genome, pedigree = ped, founder.pop = founder.pop)
#' 
#' # Create a pedigree with 100 doubled-haploid individuals induced after the F_2
#' # generation.
#' ped <- sim_pedigree(n.ind = 100, n.selfgen = 2)
#' 
#' fam <- sim_family(genome = genome, pedigree = ped, founder.pop = founder.pop, dh = TRUE)
#'                   
#' ## The above commands can be run using a hypred genome
#' genome <- sim_genome(len, n.mar, type = "hypred")
#' 
#' # Simulate a quantitative trait influenced by 50 QTL
#' qtl.model <- matrix(NA, 50, 4)
#' genome <- sim_gen_model(genome = genome, qtl.model = qtl.model, 
#'                         add.dist = "geometric", max.qtl = 50)
#' 
#' # Simulate the founder genotypes
#' founder.pop <- sim_founders(genome)
#' 
#' ped <- sim_pedigree(n.ind = 100, n.selfgen = 2)
#' 
#' fam <- sim_family(genome = genome, pedigree = ped, founder.pop = founder.pop)
#' 
#' fam <- sim_family(genome = genome, pedigree = ped, founder.pop = founder.pop, dh = TRUE)
#' 
#'  
#' @import simcross
#' @importFrom purrr map
#' @importFrom purrr pmap
#' @importFrom stringr str_pad
#' @importFrom stringr str_c
#' @importFrom qtl sim.cross
#' 
#' @export
#' 
sim_family <- function(genome, pedigree, founder.pop, ...) {
  
  # Error
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # founder.pop needs to be a pop object
  if (!inherits(founder.pop, "pop"))
    stop("The input 'founder.pop' must be of class 'pop'")
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = founder.pop$geno))
    stop("The geno did not pass. See warning for reason.")
  
  # Check the pedigree
  if (!simcross::check_pedigree(pedigree, ignore_sex = TRUE))
    stop("The pedigree is not formatted correctly.")
  
  # Get the genome type
  type <- attr(genome, "type")
  
  
  # How many founders?
  n_founders <- nind(founder.pop)
  
  # Are the number of founders correct vis a vis the pedigree?
  if (n_founders != sum(pedigree$gen == 0))
    stop("The number of founders in the inpute 'fouders' is not equal to the
         number of founders in the 'pedigree.'")
  
  # Parse other arguments
  other.args <- list(...)
  
  # For doubled-haploids
  dh <- ifelse(is.null(other.args$dh), FALSE, other.args$dh)
  
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
  new_names <- str_c("C", cycle_num, "_", max_gen, str_pad(family_num, width = 3, pad = 0), 
                     "-", str_pad(seq_len(n_ind), width = 3, pad = 0)) 
  
  
  # Determine the selfing level
  selfing <- attr(pedigree, "selfing")
  
  
  # Split the stream by genome type
  if (type == "pbsim") {
    
    # Extract the map
    map <- genome$map
    
    # For xodata
    m <- ifelse(is.null(other.args$m), 0, other.args$m)
    p <- ifelse(is.null(other.args$p), 1, other.args$p)
    
    if (selfing == "partial") {
    
      # Generate cross-over data
      xo_data <- sim_from_pedigree_allchr(pedigree = pedigree, map = map, m = m, p = p)
      
      # Simulate DH if called for
      if (dh) {
        xo_data <- induce_dh(xodat = xo_data, pedigree = pedigree)
      }
      
      # Simulate genotypic data
      prog_genos <- convert2geno_allchr(xodat = xo_data, map = map, id = final_id)
      
    } else {
      
      # Otherwise use sim.cross
      cross_sim <- qtl::sim.cross(map = genome$map, n.ind = length(final_id), type = "riself", 
                                  m = m, p = p)
      
      # Extract progeny genos
      prog_genos <- lapply(X = cross_sim$geno, FUN = function(geno_chr) 
        ifelse(geno_chr$data == 1, 1, ifelse(geno_chr$data == 2, 3, NA)) )
      
      # Bind
      prog_genos <- do.call("cbind", prog_genos)
      
      # Add row.names
      row.names(prog_genos) <- final_id
      
    }
    
    # Extract founder genotypes
    founder_geno <- do.call("cbind", founder.pop$geno)
    
    # Create a matrix with the founder1, het, and founder2 genotypes
    founder_multipoint <- rbind(founder_geno[1,], colMeans(founder_geno), founder_geno[2,])
    
    # Convert the progeny genotypes to the parental states
    prog_geno <- t(apply(X = prog_genos, MARGIN = 1, FUN = function(prog) {
      founder_multipoint[cbind(prog, seq(nrow(founder_multipoint)))] }))
    
    dimnames(prog_geno) <- list(new_names, markernames(genome, include.qtl = TRUE))
    
    
    
  } else if (type == "hypred") {
    
    # If selfing is complete, error out
    if (selfing == "complete")
      stop("Complete selfing is not supported in 'hypred'")
      
    # Extract the haploids as generation 0
    gen0 <- founder.pop$haploids
    
    # Iterate from generation 0 to max_gen
    for (g in seq(0, max_gen)) {
      
      if (g == 0) {
        
        pedigree_append <- subset(pedigree, gen == g) %>%
          group_by_(.dots = names(.)) %>%
          do(hap = {map(gen0, function(chr_array) chr_array[,,.$id])} )
        
      } else {
        
        pedigree_append <- subset(pedigree, gen == g) %>%
          group_by_(.dots = names(.)) %>%
          do(hap = {
            
            # Gamete1
            mom_id <- match(.$mom, pedigree_append$id)
            gamete1 <- pbsim:::recombine_hypred(genome = genome, haploids = pedigree_append$hap[[mom_id]])
            
            # Gamete2
            dad_id <- match(.$dad, pedigree_append$id)
            gamete2 <- pbsim:::recombine_hypred(genome = genome, haploids = pedigree_append$hap[[dad_id]])
            
            # Combine and return
            pmap(list(gamete1, gamete2), rbind) })
        
      }}

    # If doubled haploids should be induced, sample one haploid from each
    # chromosome from each individual and double it
    if (dh) {
      resample_haps <- map(pedigree_append$hap, function(ind)
        map(ind, function(chr) chr[rep(sample(c(1,2), 1), 2),] ))
      
    } else {
      resample_haps <- pedigree_append$hap
      
    }
    
    
    # Convert to list of haploids per chromosome
    haploids <- pmap(resample_haps, list) %>%
      map(simplify2array)
    
    # Add new line names
    prog_geno <- map(haploids, function(chr_haps) 
      structure(chr_haps, dimnames = list(NULL, colnames(chr_haps), new_names)) )
    
  }
    
  # Create the pop
  create_pop(genome = genome, geno = prog_geno)
  
} # Close the function




#' Simulate multiple families from a crossing block
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
#' @param ... Additional arguments passed to \code{\link{sim_family}}.
#' 
#' @return 
#' An object of class \code{pop} with the information for all individuals in 
#' the families specified by the crossing block.
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
#' # Simulate a group of families from the crossing block
#' fam_cb <- sim_family_cb(genome = genome, pedigree = ped, founder.pop = founder.pop, crossing.block = cb)
#' 
#' 
#' ## The above commands can be run using a hypred genome
#' genome <- sim_genome(len, n.mar, type = "hypred")
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
#' # Simulate a group of families from the crossing block
#' fam_cb <- sim_family_cb(genome = genome, pedigree = ped, founder.pop = founder.pop, crossing.block = cb)
#' 
#' @importFrom simcross check_pedigree
#' @import dplyr
#' 
#' @export
#' 
sim_family_cb <- function(genome, pedigree, founder.pop, crossing.block, ...) {
  
  # Error habdling
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # founder.pop needs to be a pop object
  if (!inherits(founder.pop, "pop"))
    stop("The input 'founder.pop' must be of class 'pop'")
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = founder.pop$geno))
    stop("The geno did not pass. See warning for reason.")
  
  # Check the pedigree
  if (!simcross::check_pedigree(pedigree, ignore_sex = TRUE))
    stop("The pedigree is not formatted correctly.")
  
  # Check the crossing block
  if (ncol(crossing.block) != 2) {
    stop("The crossing block should have two columns.")
  } else {
    crossing.block <- as.data.frame(crossing.block)
  }
  
  # Are all of the parents in the crossing block in the founder.pop?
  if (!all(unique(unlist(crossing.block)) %in% indnames(founder.pop)))
    stop("Not all of the parents in the crossing block are in the 'founder.pop'.")
  
  # Simulate families for each cross
  fam_cb <- crossing.block %>% 
    # Add family number
    mutate(fam_num = row_number()) %>%
    rowwise() %>%
    do(fam = {
      founder_geno <- subset_pop(pop = founder.pop, individual = c(.$parent1, .$parent2))
      sim_family(genome = genome, pedigree = pedigree, founder.pop = founder_geno, 
                 family.num = .$fam_num, ... = ...) })
  
  fam_cb <- crossing.block %>% 
    # Add family number
    mutate(fam_num = row_number()) %>%
    rowwise() %>%
    do(fam = {
      founder_geno <- subset_pop(pop = founder.pop, individual = c(.$parent1, .$parent2))
      sim_family(genome = genome, pedigree = pedigree, founder.pop = founder_geno, 
                 family.num = .$fam_num, other.args) })

  
  # Combine the populations and return
  combine_pop(pop_list = fam_cb$fam)

} # Close the function
  
  
  
  
  


