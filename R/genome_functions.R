#' Create a simulated genome
#' 
#' @description
#' Creates a list containing information on a simulated genome
#' 
#' @param len A vector specifying the chromosome lengths (in cM).
#' @param n.mar A vector specifying the umber of markers per chromosome.
#' @param map A list of marker positions (in cM), where each list element is a named vector of marker positions that
#' are located on a single chromosome. The names of these elements in the list are chromosome \emph{numbers}. If 
#' \code{NULL} (default), marker positions are drawn from a uniform distribution. See \code{\link[qtl]{sim.map}} for more
#' information.
#' @param eq.spacing If TRUE, markers will be equally spaced. See \code{\link[qtl]{sim.map}}.
#' @param type Deprecated. 'hyped' genomes are no longer supported.
#' 
#' @return 
#' Object of class "genome" with the length of each chromosome, the number of
#' markers per chromosome, and the genetic map.
#' 
#' @examples
#' 
#' \dontrun{
#' n.mar  <- c(505, 505, 505)
#' len <- c(120, 130, 140)
#' 
#' genome <- sim_genome(len, n.mar)
#' 
#' exp_map <- replicate(1, runif(10), simplify = FALSE)
#' # This will fail because the vector elements in the list are not mapped
#' sim_genome(map = exp_map)
#' 
#' # Name the markers
#' names(exp_map[[1]]) <- paste0("M", seq_along(exp_map[[1]]))
#' # This will fail because the list is not named.
#' sim_genome(map = exp_map)
#' 
#' # Give the chromosomes names
#' names(exp_map) <- "chr1"
#' # This will fail because the chromosome names are not coercible to numbers
#' sim_genome(map = exp_map)
#' 
#' # This will work
#' names(exp_map) <- 1
#' sim_genome(map = exp_map)
#' }
#' 
#' 
#' @importFrom qtl sim.map
#' 
#' @export 
#' 
#' 
sim_genome <- function(len, n.mar, map, eq.spacing = FALSE) {
  
  # Error handling
  
  # len and n.mar must be present if map is not
  if (missing(map)) {
    len <- as.numeric(len)
    n.mar <- as.numeric(n.mar)
    
  } else {
    if (missing(len))
      len <- sapply(map, max)
    
    if (missing(n.mar))
      n.mar <- sapply(map, length)
    
  }
  
  # The length of the len vector must be the same as that of n.mar
  if (length(len) != length(n.mar)) {
    stop("The length of the input 'len' vector does not equal the length of
         'n.mar.'")
  }
  
  # Create an empty list
  genome <- structure(vector("list"), class = "genome")
  
  # If genetic map is NULL, sample from uniform distribution
  if (missing(map)) {
    # Use R/qtl
    map <- sim.map(len = len, n.mar = n.mar, include.x = FALSE, eq.spacing = eq.spacing)
    
    # If map is not null, check compatability
  } else {
    
    # First check that the markers are named
    if (any(sapply(map, function(m) is.null(names(m)))))
      stop("Marker positions in the input 'map' must be named.")
    
    # Make sure the map has the chromosome class attribute
    map <- lapply(map, structure, class = "A")
    # Reclass the map
    class(map) <- "map"
    
  }
  
  # Assemble the genome
  genome[["len"]] <- len
  genome[["n.mar"]] <- n.mar
  genome[["map"]] <- map
  
  
  
  
  
  # Return the genome
  return(genome)
  
} # Close the function


#' Summarize a genome object
#' 
#' @param x An object of class \code{genome}.
#' 
#' @export
#' 
summary.genome <- function(x) {
  
  
  n_chr <- nchr(x)
  len <- chrlen(x)
  n_mar <- nmar(x, by.chr = TRUE)
  
  # Extract the genetic model
  gen_model <- x$gen_model
  
  # Show
  cat("\nGenome summary \n\n")
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
  
  # All loci names
  loci_names <- lapply(genome$map, names)
  
  
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
  
  # All loci names
  loci_names <- lapply(genome$map, names)

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
  
  # Extract and return chromosomes
  length(genome$map)

  
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
chrnames <- function(genome) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  names(genome$map)
  
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
  
  genome$len
  
  
} # Close the function



#' Extract the QTL from a genome
#' 
#' @param genome An object of class \code{genome}.
#' @param unique Logical. Should information of only the unique QTL be returned?
#' 
#' @importFrom dplyr distinct
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
  all_qtl <- mapply(genome$gen_model, seq_along(genome$gen_model), FUN = function(qtlmod, traitn) {
    qtlmod$trait <- traitn
    return(qtlmod)
  }, SIMPLIFY = FALSE)
  
  all_qtl <- do.call("rbind", all_qtl)
  
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
  all_mods <- lapply(seq_along(genome$gen_model), FUN = function(i) {
    genome$gen_model[[i]]$trait <- i
    return(genome$gen_model[[i]]) })
  
  all_mods <- do.call("rbind", all_mods)[,c("trait", "chr", "pos", "add_eff", "dom_eff")]

  # Find duplications by finding the pairwise union
  dup_mod <- duplicated(subset(all_mods, select = c(chr, pos))) | 
    duplicated(subset(all_mods, select = c(chr, pos)), fromLast = TRUE)
  
  all_mods[which(dup_mod),,drop = FALSE]
  
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
  
  # Use the map to determine the number of loci
  n_loci <- sapply(genome$map, length)
  
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
#' @export
#' 
map_to_table <- function(genome) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Extract the genome map
  genome_map <- genome$map
  
  ## apply a function over the chromosomes in the map
  pos_df <- mapply(seq_along(genome_map), genome_map, FUN = function(chr, map) 
    cbind(chr = chr, as.data.frame(structure(map, class = NULL), nm = "pos")), SIMPLIFY = FALSE)
  
  pos_df <- do.call("rbind", pos_df)
  
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
#' @param ignore.gen.model Logical - should the gene model be ignored?
#' @param ... Additional arguments. Markers within a position range of \code{marker}
#' can be specified using the \code{min.dist} (returns markers at \code{min.dist}
#' from \code{marker}) and \code{max.dist} (returns markers no more than \code{min.dist}
#' from \code{marker}).
#' @param include.qtl Logical. Should QTL be included in the results? If FALSE, the genome
#' should have a genetic model.
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
#' sample_markers <- sample(markernames(genome, include.qtl = TRUE), size = 3)
#' 
#' find_proxmarkers(genome = genome, marker = sample_markers, include.qtl = TRUE)
#' 
#' @importFrom qtl map2table
#' @importFrom qtl sim.cross
#' @importFrom qtl chrlen
#' @importFrom qtl find.flanking
#' 
#' @export
#' 
find_proxmarkers <- function(genome, marker, ..., include.qtl = FALSE) {
  
  # Make sure genome inherits the class "genome."
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Are 'marker' in the total marker names?
  if (!all(marker %in% markernames(genome = genome, include.qtl = TRUE)))
    stop("The markers in 'marker' are not in the genome.")
  
  # Extract the other arguments
  other.args <- list(...)
  
  min.dist <- other.args$min.dist
  max.dist <- other.args$max.dist
  
  # Find the position of the marker
  marker_pos <- find_markerpos(genome = genome, marker = marker)
  
  # Chromosome lengths
  chr_len <- chrlen(genome)
  
  ## Are min.dist and max.dist both NULL?
  # If so, return flanking markers
  if (all(is.null(min.dist), is.null(max.dist))) {
    
    flanking_pos <- list()
    
    for (i in seq(nrow(marker_pos))) {
      
      mkr <- marker_pos[i,,drop = FALSE]
      
      # Markers in the chr
      chr_markers <- markernames(genome, chr = mkr$chr, include.qtl = include.qtl)
      # Find the position of the marker on the map      
      map_ind <- match(row.names(mkr), chr_markers)
      # Find the position of the marker to the left and right
      left_marker <- ifelse(map_ind == 1, NA, chr_markers[map_ind - 1])
      right_marker <- chr_markers[map_ind + 1]
      # Return data.frame
      flanking_pos[[i]] <- data.frame(left = left_marker, right = right_marker, stringsAsFactors = FALSE)
      
    }
    
    flanking_pos <- cbind(marker_pos, do.call("rbind", flanking_pos))
    
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
    
    
  
  
  
  
  
  
  
  
  
  
  
  