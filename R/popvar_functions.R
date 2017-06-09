#' Convert data for use in PopVar
#' 
#' @param genome An object of class \code{genome}.
#' @param geno Genotype data on a population to phenotype. Must be a matrix of dimensions
#' \code{n.ind} x \code{n.loci}, the elements of which must be z {0, 1, 2}.
#' 
#' @return 
#' A \code{data.frame} of marker genotypes ready for use in \code{\link[PopVar]{pop.predict}}.
#' 
#' @export
#' 
geno_to_popvar <- function(genome, geno) {
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = geno))
    stop("The geno did not pass. See warning for reason.")
  
  # If the geno input is a list, recombine
  if (is.list(geno))
    geno <- do.call("cbind", geno)
  
  # Get the individual names from the geno
  ind_names <- row.names(geno)
  # Marker names
  marker_names <- markernames(genome, include.qtl = FALSE)
  
  # Is the length of the marker names the same as the columns in 'geno'?
  # If not, subset the geno for markers, not QTL
  if (length(marker_names) != ncol(geno)) {
    geno_to_convert <- pull_genotype(genome = genome, geno = geno, loci = marker_names)
  } else {
    geno_to_convert <- geno
  }

  # Recode
  geno_recode <- geno_to_convert - 1
  
  # Create data.frame for output
  as.data.frame(cbind( c("", ind_names), rbind(marker_names, geno_recode)) )
  
} # Close the function


#' Convert data for use in PopVar
#' 
#' @describeIn geno_to_popvar 
#' 
#' @import dplyr
#' @importFrom qtl sim.cross
#' 
#' @export
#' 
map_to_popvar <- function(genome) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Convert the map to a table
  map_table <- map_to_table(genome)
  
  map_table %>% 
    mutate(marker = row.names(.)) %>% 
    select(marker, chr, pos)
  
}



#' Quicker procedures for PopVar
#' 
#' @param genome An object of class \code{genome}.
#' @param map An object of class map to simulate the virtual bi-parental populations.
#' @param training.pop An object of class \code{pop} with the elements \code{geno} and 
#' \code{pheno_val}. This is used as the training population.
#' @param founder.pop An object of class \code{pop} with the geno information for
#' the parents. Additional individuals can be present in \code{parent_pop}. They
#' will be filtered according to the parents in the \code{crossing_block}.
#' @param crossing.block A crossing block detailing the crosses to make. Must be a
#' \code{data.frame} with 2 columns: the first gives the name of parent 1, and the 
#' second gives the name of parent 2. See \code{\link{sim_crossing_block}}.
#' @param k.sp The standardized selection coefficient (e.g. a selection intensity of 0.1 means k_sp ~ 1.745).
#' @param map.function The map function for converting genetic distance to recomination fraction.
#' 
#' @examples 
#' 
#' # Simulate a genome
#' n.mar  <- c(505, 505, 505)
#' len <- c(120, 130, 140)
#' 
#' genome <- sim_genome(len, n.mar)
#' 
#' # Simulate another genome for use in creating the bi-parental populations
#' genome_use <- sim_genome(len, n.mar)
#' 
#' # Simulate a quantitative trait influenced by 50 QTL
#' qtl.model <- matrix(NA, 50, 4)
#' genome <- sim_gen_model(genome = genome, qtl.model = qtl.model, 
#'                         add.dist = "geometric", max.qtl = 50)
#'                         
#' # Remove QTL from the second genome map
#' map_use <- lapply(genome_use$map, function(chr) chr[setdiff(names(chr), qtlnames(genome))])
#' 
#' # Simulate the genotypes of eight founders
#' founder_pop <- sim_founders(genome, n.str = 8)
#' founder_pop <- sim_phenoval(pop = founder_pop, h2 = 0.5)
#' 
#' ped <- sim_pedigree(n.ind = 100, n.selfgen = 2)
#' 
#' # Extract the founder names
#' parents <- indnames(founder_pop)
#' 
#' # Generate a crossing block with 5 crosses
#' cb <- sim_crossing_block(parents = parents, n.crosses = 5)
#' 
#' # Use the map in the second genome to predict the genetic variance in the 
#' # bi-parental populations
#' pop_predict_quick(genome = genome, map = map_use, training.pop = founder_pop, 
#'                   founder.pop = founder_pop, crossing.block = cb, k.sp = 1.745, 
#'                   map.function = "haldane")
#' 
#' 
#' 
#' @import dplyr
#' @importFrom qtl sim.cross mf.h mf.k mf.m mf.cf
#' @importFrom purrr map pmap
#' 
#' @export
#' 
pop_predict_quick <- function(genome, map, training.pop, founder.pop, crossing.block,
                              k.sp, map.function = c("haldane","kosambi","cf","morgan")) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Set the map function
  map.function <- match.arg(map.function)
  
  # First predict marker effects
  mar_eff <- pred_mar_eff(genome = genome, training.pop = training.pop)$mar_eff
  
  # Predict genotypic values
  pgv <- pred_geno_val(genome = genome, training.pop = training.pop, candidate.pop = founder.pop,
                       method = "RRBLUP")
  
  # Convert the usable map to a df
  map_df <- qtl::map2table(map) %>%
    data.frame(marker = row.names(.), ., stringsAsFactors = FALSE)
  
  # Combine marker name, position, and effect
  mar_specs <- full_join(x = map_df, mar_eff, by = "marker")
  # Add marker names to row.names
  row.names(mar_specs) <- mar_specs$marker
  
  # Calculate the pairwise distance between markers
  mar_specs_pair <- mar_specs %>%
    split(.$chr) %>%
    map(function(mar_chr) {
      
      # Pairwise combinations of markers
      mar_pairs <- combn(x = mar_chr$marker, m = 2) %>% 
        t() %>% 
        data.frame(stringsAsFactors = FALSE) %>% 
        structure(names = c("mar1", "mar2"))
      
      # Gather the marker effects and positions for those markers
      mar_pairs %>%
        mutate(mar1_eff = mar_chr[mar1, "trait1"], 
               mar2_eff = mar_chr[mar2, "trait1"], 
               mar1_pos = mar_chr[mar1, "pos"], 
               mar2_pos = mar_chr[mar2, "pos"],
               d = abs(mar1_pos - mar2_pos),
               c = switch(map.function,
                          haldane = mf.h(d),
                          kosambi = mf.k(d),
                          cf = mf.cf(d),
                          morgan = mf.m(d)) ) %>%
        select(mar1:mar2_eff, c) })

  
  # Iterate over the parents in the crossing block
  crossing.block %>%
    by_row(function(pars) {
      
      # Extract parental genotypes and recode
      parents <- subset_pop(pop = founder.pop, individual = pars)$geno %>%
        lapply(function(geno) geno - 1)
      
      # Iterate over the chromosomes of the parents and map
      # Which markers are segregating?
      mar_seg <- list(map, parents) %>%
        pmap(function(map_chr, geno_chr)
          geno_chr[,names(map_chr)] %>% 
            apply(MARGIN = 2, FUN = function(snp) n_distinct(snp) > 1) %>%
            which() %>%
            names() )
      
      # Subset markers that are segregating and add the genotype of the first parent
      mar_specs_pair_seg <- list(mar_specs_pair, parents, mar_seg) %>%
        pmap(function(specs_chr, par_chr, seg_chr) {
          specs_chr %>% 
            filter(mar1 %in% seg_chr & mar2 %in% seg_chr) %>% 
            mutate(mar1_par1 = par_chr[1,mar1], mar2_par1 = par_chr[1,mar2]) })
      
      
      # Calculate the variance of each segregating marker
      mar_var <- mar_specs_pair_seg %>%
        bind_rows() %>%
        distinct(mar1, .keep_all = TRUE) %>%
        mutate(var_a = (mar1_eff^2)) %>%
        summarize(var_a = sum(var_a))
      
      # Calculate the covariance between marker pairs
      mar_covar <- mar_specs_pair_seg %>%
        bind_rows() %>%
        # Calculate covariance
        mutate(frac = (1 - (2 * c)) / (1 + (2 * c)),
               eff_prod = (mar1_eff * mar1_par1) * (mar2_eff * mar2_par1),
               cov = frac * eff_prod ) %>%
        summarize(covar = sum(cov))
      
      # Calculate the mean PGV based on the markers
      pred_mu <- pgv$pred_val %>% 
        filter(ind %in% pars) %>% 
        summarize(mean_pgv = mean(trait1)) %>%
        as.numeric()
      
      # Calculate genetic variance from the expected equation
      pred_varG <- as.numeric(mar_var + (2 * mar_covar)) 
      
      # Calculate the mu_sp (high and low)
      # From Zhong and Jannink 2007
      pred_mu_sp_high <- pred_mu + (k.sp * pred_varG)
      pred_mu_sp_low <- pred_mu - (k.sp * pred_varG)
      
      # Return vector
      data.frame(pred_mu, pred_varG, pred_mu_sp_high, pred_mu_sp_low) }, .collate = "cols") %>%
    
    # Rename
    rename(pred_mu = pred_mu1, pred_varG = pred_varG1, pred_mu_sp_high = pred_mu_sp_high1, 
           pred_mu_sp_low = pred_mu_sp_low1)
  
} # Close the function
        
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

