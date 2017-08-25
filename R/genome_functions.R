#' Summarize a genome object
#' 
#' @param x An object of class \code{genome}.
#' 
#' @export
#' 
summary.genome <- function(x) {
  
  # Summarize by hypred or pbsim
  type <- attr(x, "type")
  
  n_chr <- nchr(x)
  len <- chrlen(x)
  n_mar <- nmar(x, by.chr = TRUE)
  
  # Extract the genetic model
  gen_model <- x$gen_model
  
  # Show
  cat("\nGenome summary \n\n")
  cat("Genome type: ", type, "\n\n")
  cat("Number of chromosomes: ", n_chr, "\n")
  cat("Length of chromosomes (in cM): ", len, "\n")
  cat("Markers per chromosome: ", n_mar, "\n")
  
  # Report the genetic model
  if (is.null(gen_model)) {
    cat("\nGenetic model summary: No genetic model")
    
  } else {
    
    # Number of traits
    n_trait <- length(gen_model)
    cat("\nNumber of traits with genetic models: ", n_trait, "\n")
    
    # Iterate over traits
    for (t in seq(n_trait)) {
      
      cat("\n\nSummary for trait ", t, ":\n")
      
      n_qtl <- nrow(gen_model[[t]])
      eff_qtl <- subset(gen_model[[t]], add_eff != 0 | dom_eff != 0)
      n_eff_qtl <- nrow(eff_qtl)
      qtl_chr <- table(eff_qtl$chr)
      
      
      cat("Total QTL: ", n_qtl, "\n")
      cat("Number of effective QTL: ", n_eff_qtl, "\n")
      cat("Effective QTL per chromosome: ", qtl_chr, "\n\n") 
      cat("Distribution of additive effects of effective QTL:\n")
      print(summary(eff_qtl$add_eff))
      cat("Distribution of dominance effects of effective QTL: \n")
      print(summary(eff_qtl$dom_eff))
      
    }
    
  }
  
}

#' Printing genomes
#' 
#' @rdname summary.genome
#' 
#' @export
#' 
print.genome <- function(x) summary(x)



#' Extract the number of markers in the genome
#' 
#' @param genome An object of class \code{genome}.
#' @param by.chr Logical. Should the number of markers per chromosome be returned?
#' 
#' @return 
#' If \code{by.chr = TRUE}, a vector of markers per chromosomes. If \code{by.chr = FALSE},
#' a scalar of the total number of markers.
#' 
#' @export
#' 
nmar <- function(genome, by.chr = FALSE) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Extract information based on the type
  type <- attr(genome, "type")
  
  if (type == "pbsim") {
    
    # All loci names
    loci_names <- lapply(genome$map, names)
    
  } else if (type == "hypred") {
    # Get all loci names
    loci_names <- lapply(genome$hypredGenomes, slot, "pos.snp") %>% lapply(names)
    
  }
  
  # Get the QTL names if present
  if (!is.null(genome$gen_model)) {
    # Get QTL names
    qtl_names <- qtlnames(genome)
    # Find the loci that are not QTL
    marker_names <- lapply(X = loci_names, FUN = setdiff, y = qtl_names)
    
  } else {
    # Otherwise all of the loci are markers
    marker_names <- loci_names
    
  }
  
  # Number of markers
  n_marker <- sapply(X = marker_names, FUN = length)

  # If by.chr is true, give results by chromosome. If not, sum
  if (by.chr) {
    return(n_marker)
  } else{
    return(sum(n_marker))
  }
    
  
} # Close the function


#' Extract marker names
#' 
#' @param genome An object of class \code{genome}.
#' @param chr See \code{\link[qtl]{markernames}}.
#' @param include.qtl Logical. Should the names of QTL be returned?
#' 
#' @return 
#' A vector of character strings (the marker names).
#' 
#' @import dplyr
#' 
#' @export
#' 
markernames <- function(genome, chr, include.qtl = FALSE) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # If chr is missing, assume all chromosomes
  if (missing(chr)) {
    chr <- chrnames(genome)
  }
  
  # Extract information based on the type
  type <- attr(genome, "type")
  
  if (type == "pbsim") {
    
    # All loci names
    loci_names <- lapply(genome$map, names)
    
  } else if (type == "hypred") {
    # Get all loci names
    loci_names <- lapply(genome$hypredGenomes, slot, "pos.snp") %>% lapply(names)
    
  }
  
  # Get the QTL names if present
  if (!include.qtl) {
    # Get QTL names
    qtl_names <- qtlnames(genome)
    # Find the loci that are not QTL
    marker_names <- lapply(X = loci_names, FUN = setdiff, y = qtl_names)
    
  } else {
    # Otherwise all of the loci are markers
    marker_names <- loci_names
    
  }

  # What chromosomes to return?
  structure(unlist(marker_names[chr]), names = NULL)
  
} # Close function
  

#' Extract the number of chromosomes in the genome
#' 
#' @param genome An object of class \code{genome}.
#' 
#' @return 
#' A scalar of the total number of chromosomes.
#' 
#' @export
#' 
nchr <- function(genome) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Extract information based on the type
  type <- attr(genome, "type")
  
  if (type == "pbsim") {
    n_chr <- length(genome$map)
    
  } else if (type == "hypred") {
    n_chr <- length(genome$hypredGenomes)
  
  }
  
  # Extract and return chromosomes
  return(n_chr)
  
} # Close the function



#' Extract chromosome names
#' 
#' @param genome An object of class \code{genome}.
#' 
#' @return 
#' A vector of character strings (the chromosome names).
#' 
#' @export
#' 
chrnames <- function(genome, chr) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Extract information based on the type
  type <- attr(genome, "type")
  
  if (type == "pbsim") {
    chr_names <- names(genome$map)
    
  } else if (type == "hypred") {
    chr_names <- names(genome$hypredGenomes)
    
  }
  
  return(chr_names)
  
} # Close function


#' Extract the length of chromosomes
#' 
#' @param genome An object of class \code{genome}.
#' 
#' @return 
#' A vector of chromosome lengths.
#' 
#' @importFrom methods slot
#' 
#' @export
#' 
#' 
chrlen <- function(genome) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Extract information based on the type
  type <- attr(genome, "type")
  
  if (type == "pbsim") {
    chr_len <- genome$len
    
  } else if (type == "hypred") {
    chr_len <- sapply(genome$hypredGenomes, slot, "len.chr") * 100
    
  }
  
  return(chr_len)
  
} # Close the function



#' Extract the QTL from a genome
#' 
#' @param genome An object of class \code{genome}.
#' @param unique Logical. Should information of only the unique QTL be returned?
#' 
#' @import dplyr
#' 
#' @export
#' 
pull_qtl <- function(genome, unique = TRUE) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Make sure there is a genetic model
  if (is.null(genome$gen_model))
    stop("No genetic model has been declared for the genome")
  
  # Unique qtl
  all_qtl <- genome$gen_model %>% 
    mapply(seq_along(.), FUN = function(qtlmod, traitn) 
      mutate(qtlmod, trait = traitn), SIMPLIFY = FALSE) %>%
    bind_rows()
  
  if (unique) {
    return(distinct(all_qtl, chr, pos, .keep_all = TRUE))
  } else {
    return(all_qtl)
  }
  
} # Close the function

#' Extract the number of QTL in the genome
#' 
#' @param genome An object of class \code{genome}.
#' @param by.chr Logical. Should the number of QTL per chromosome be returned?
#' 
#' @details
#' Only unique QTL are counted. That is, pleiotropic QTL only count once.
#' 
#' @return 
#' If \code{by.chr = TRUE}, a vector of the QTL per chromosomes. If 
#' \code{by.chr = FALSE}, a scalar of the total number of QTL.
#' 
#' @import dplyr
#' 
#' @export
#' 
nqtl <- function(genome, by.chr = FALSE) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # If there is no genetic model, there are no QTL
  if (is.null(genome$gen_model)) {
    n_qtl <- structure(rep(0, nchr(genome)), names = chrnames(genome))
    
  } else {
    # Number of unique qtl per chromosome
    qtl_unique <- pull_qtl(genome = genome)
    
    n_qtl <- structure(table(qtl_unique$chr), names = chrnames(genome))
    
  }

  # If by.chr is true, give results by chromosome. If not, sum
  if (by.chr) {
    return(n_qtl)
  } else{
    return(sum(n_qtl))
  }
  
} # Close the function



#' Extract the names of QTL in the genome
#' 
#' @param genome An object of class \code{genome}.
#' 
#' @details
#' Only unique QTL are counted. That is, pleiotropic QTL only count once.
#' 
#' @return 
#' If \code{by.chr = TRUE}, a vector of the QTL names.
#' 
#' @export
#' 
qtlnames <- function(genome, chr) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Make sure there is a genetic model
  if (is.null(genome$gen_model))
    stop("No genetic model has been declared for the genome")
  
  # If chr is missing, assume all chromosomes
  if (missing(chr)) {
    chr <- chrnames(genome)
  }
  
  # Get the QTL names
  all_qtl <- pull_qtl(genome, unique = FALSE)
  # Return
  subset(x = all_qtl, chr %in% chr, qtl_name, drop = TRUE)

} # Close the function


#' Extract the pleiotropic QTL from the genome
#' 
#' @param genome An object of class \code{genome}.
#' 
#' @return 
#' A \code{data.frame} of the pleiotropic QTL per trait.
#'  
#' @import dplyr
#' 
#' @export
#' 
pull_plei_qtl <- function(genome) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Make sure there is a genetic model
  if (is.null(genome$gen_model))
    stop("No genetic model has been declared for the genome")
  
  # Combind genetic models for traits
  all_mods <- lapply(seq_along(genome$gen_model), FUN = function(i) 
    genome$gen_model[[i]] %>% mutate(trait = i) ) %>%
    bind_rows() %>%
    select(trait, chr:dom_eff)

  # Find duplications by finding the pairwise union
  dup_mod <- duplicated(subset(all_mods, select = c(chr, pos))) | 
    duplicated(subset(all_mods, select = c(chr, pos)), fromLast = TRUE)
  
  all_mods %>% 
    slice(which(dup_mod))
  
} # Close the function


#' Extract the perfect markers from the genome
#' 
#' @description 
#' Extract SNP markers that are also QTL
#' 
#' @param genome An object of class \code{genome}.
#' @param by.chr Logical. Should the number of perfect markers per chromosome be returned?
#' 
#' @return 
#' A \code{list} of the perfect QTL.
#' 
#' @examples
#' n.mar  <- c(3, 3, 3)
#' len <- c(120, 130, 140)
#' 
#' map <- lapply(X = seq(3), FUN = function(i) 
#'   structure(c(20, 50, 80), names = paste("C", i, "M", seq(3), sep = "")))
#' 
#' genome <- sim_genome(len, n.mar, map)
#' 
#' # Simulate genetic architecture with one perfect marker per chromosome
#' chromosome <- seq(3)
#' pos <- rep(50, 3)
#' a <- rep(1, 3)
#' d <- 0
#' 
#' qtl.model <- cbind(chromosome, pos, a, d)
#' 
#' genome <- sim_gen_model(genome, qtl.model)
#'  
#' @import dplyr
#' 
#' @export
#' 
pull_perf_mar <- function(genome, by.chr = FALSE) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Make sure there is a genetic model
  if (is.null(genome$gen_model))
    stop("No genetic model has been declared for the genome")
  
  # Pull out the unique QTL
  unique_qtl <- pull_qtl(genome = genome, unique = TRUE)
  
  qtl_ann <- unique_qtl %>% 
    rowwise() %>% 
    mutate(perf_mar = list(genome$map[[chr]] %in% pos)) %>%
    mutate(perf_mar = ifelse(any(perf_mar), names(genome$map[[chr]])[which(perf_mar)], NA)) %>%
    as.data.frame()
  
  # Reduce
  perf_count <- qtl_ann %>% 
    filter(!is.na(perf_mar))
    
  # Split by chromsome?
  if (by.chr) {
    split(perf_count, perf_count$chr)
  } else {
    return(perf_count)
  }
  
} # Close the function


#' Extract the total number of loci in the genome
#' 
#' @description 
#' Calculate the total number of loci (markers + QTL)
#' 
#' @param genome An object of class \code{genome}.
#' @param by.chr Logical. Should the number of loci per chromosome be returned?
#' 
#' @export
#' 
nloci <- function(genome, by.chr = FALSE) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Extract the type
  type <- attr(genome, "type")
  
  # Pipe by type
  if (type == "pbsim") {
  
    # Use the map to determine the number of loci
    n_loci <- sapply(genome$map, length)
    
  } else if (type == "hypred") {
    
    # Pull out the slots
    n_loci <- sapply(genome$hypredGenomes, slot, "num.snp.chr")
  
  }
  
  # By chromosome?
  if (by.chr) {
    return(n_loci)
  } else {
    return(sum(n_loci))
  }
  
} # Close the function


#' Convert the genetic map in the genome to a table
#' 
#' @param genome An object of class \code{genome}.
#' 
#' @return 
#' A \code{data.frame} with two columns: chromosome name and position, with row names
#' equal to marker names
#' 
#' @import dplyr
#' 
#' @export
#' 
map_to_table <- function(genome) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Extract the type
  type <- attr(genome, "type")
  
  # Pipe by type
  if (type == "pbsim") {
    
    pos_df <- lapply(X = genome$map, structure, class = NULL) %>% 
      lapply(FUN = as.data.frame) %>% 
      lapply(structure, names = "pos") %>% 
      mapply(chrnames(genome), FUN = function(pos, chrname)
        data.frame(chr = chrname, pos, stringsAsFactors = FALSE), SIMPLIFY = FALSE) %>%
      do.call("rbind", .)
    
    row.names(pos_df) <- markernames(genome, include.qtl = TRUE)
    
  } else if (type == "hypred") {
    
    # Create a data.frame of chr and pos
    pos_df <- lapply(genome$hypredGenomes, slot, "pos.snp") %>% 
      lapply(FUN = as.data.frame) %>% 
      lapply(structure, names = "pos") %>%
      mapply(chrnames(genome), FUN = function(pos, chrname)
        data.frame(chr = chrname, pos * 100, stringsAsFactors = FALSE), SIMPLIFY = FALSE) %>%
      do.call("rbind", .)
    
    row.names(pos_df) <- markernames(genome, include.qtl = TRUE)
    
  }
  
  # Return the data.frame
  return(pos_df)
  
}

#' Convert a table of genetic positions to a map
#' 
#' @description This function is a simple wrapper around \code{\link[qtl]{table2map}}.
#' 
#' @param df A \code{data.frame} of genetic positions. The first column is the 
#' chromosome and the second column is the genetic position (in cM). Row names
#' should be the marker names.
#' 
#' @return 
#' A \code{list} with length equal to the number of chromosomes. Elements in the
#' list are named numeric vectors of genetic positions and the names are the
#' marker names. 
#' 
#' @importFrom qtl table2map
#' 
#' @export
#' 
table_to_map <- function(df) table2map(df)




#' Find the position of markers in the genome
#' 
#' @param genome An object of class \code{genome}.
#' @param marker See \code{\link[qtl]{find.markerpos}}.
#' 
#' @import dplyr
#' 
#' @export
#' 
find_markerpos <- function(genome, marker) {
  
  # Convert the map to a table
  pos_df <- map_to_table(genome)

  # Subset the data.frame for the marker
  marker_pos <- pos_df[marker,, drop = FALSE]
  
  return(marker_pos)
  
} # Close the function



#' Find the position of markers near a locus
#' 
#' @description 
#' Finds the position of markers near a particular locus. By default, the function
#' returns the flanking markers, however markers within a certain range can also
#' be returned.
#' 
#' @param genome An object of class \code{genome}.
#' @param marker See \code{\link[qtl]{find.markerpos}}.
#' @param ... Additional arguments. Markers within a position range of \code{marker}
#' can be specified using the \code{min.dist} (returns markers at \code{min.dist}
#' from \code{marker}) and \code{max.dist} (returns markers no more than \code{min.dist}
#' from \code{marker}).
#' 
#' @examples 
#' 
#' # Simulate the genome
#' n.mar  <- c(505, 505, 505)
#' len <- c(120, 130, 140)
#' 
#' genome <- sim_genome(len, n.mar)
#' 
#' # Sample marker names to lookup
#' sample_markers <- sample(markernames(genome), size = 3)
#' 
#' @importFrom qtl map2table
#' @importFrom qtl sim.cross
#' @importFrom qtl chrlen
#' @importFrom qtl find.flanking
#' 
#' @export
#' 
find_proxmarkers <- function(genome, marker, ...) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Are 'marker' in the total marker names?
  if (!all(marker %in% markernames(genome = genome, include.qtl = TRUE)))
    stop("The markers in 'marker' are not in the genome.")
  
  # Extract the type
  type <- attr(genome, "type")
  
  # Extract the other arguments
  other.args <- list(...)
  
  min.dist <- other.args$min.dist
  max.dist <- other.args$max.dist
  
  # Find the position of the marker
  marker_pos <- find_markerpos(genome = genome, marker = marker)
  
  ## Are min.dist and max.dist both NULL?
  # If so, return flanking markers
  if (all(is.null(min.dist), is.null(max.dist))) {
    
    # Iterate over the markers
    flanking_pos <- marker_pos %>% 
      mutate(marker = row.names(.)) %>% 
      group_by(marker) %>% 
      do({
        
        # Markers in the chr
        chr_markers <- markernames(genome, chr = .$chr)
        # Find the position of the marker on the map      
        map_ind <- match(.$marker, chr_markers)
        # Find the position of the marker to the left and right
        left_marker <- ifelse(map_ind == 1, NA, chr_markers[map_ind - 1])
        right_marker <- chr_markers[map_ind + 1]
        # Return data.frame
        data.frame(left = left_marker, right = right_marker, stringsAsFactors = FALSE)
        
      }) %>%
      as.data.frame()
    
    return(flanking_pos)
    
  } else if (is.null(max.dist)) {
    
    # Find all markers on the chromosome at least min.dist away
    
    # First declare the range of positions in which the markers can be
    lower_range <- data.frame(chr = marker_pos$chr, min = 0, 
                              max = marker_pos$pos - min.dist)
    
    upper_range <- data.frame(chr = marker_pos$chr, min = marker_pos$pos + min.dist,
                              max = chr_len[marker_pos$chr])
    
  } else if (is.null(min.dist)) {
    # Find markers within the maximum distance
    
    lower_range <- data.frame(chr = marker_pos$chr, min = marker_pos$pos - max.dist, 
                              max = marker_pos$pos)
    
    upper_range <- data.frame(chr = marker_pos$chr, min = marker_pos$pos,
                              max = marker_pos$pos + max.dist)
    
  } else {
    
    # Find markers in the min/max range 
    
    lower_range <- data.frame(chr = marker_pos$chr, min = marker_pos$pos - max.dist, 
                              max = marker_pos$pos - min.dist)
    
    upper_range <- data.frame(chr = marker_pos$chr, min = marker_pos$pos + min.dist,
                              max = marker_pos$pos + max.dist)
    
  }
    
  # List of marker names
  lower_mar <- upper_mar <- vector("list", nrow(marker_pos))
  
  # Convert the map to a table
  snp_table <- map_to_table(genome)
  
  # Iterate over markers
  for (i in 1:nrow(marker_pos)) {
    range_mar <- subset(snp_table, chr == lower_range$chr[i] & pos >= lower_range$min[i] & pos <= lower_range$max[i])
    
    # Add the row.names
    lower_mar[[i]] <- row.names(range_mar)
  }
    
  # Iterate over markers
  for (i in 1:nrow(marker_pos)) {
    
    range_mar <- subset(snp_table, chr == upper_range$chr[i] & pos >= upper_range$min[i] & pos <= upper_range$max[i])
    
    # Add the row.names
    upper_mar[[i]] <- row.names(range_mar)
  }
  
  # Combine
  prox_mar <- structure(mapply(lower_mar, upper_mar, FUN = c, SIMPLIFY = FALSE),
                        names = marker)
  
  # Return
  return(prox_mar)
    
} # Close the function
    
    
  
  
  
  
  
  
  
  
  
  
  
  