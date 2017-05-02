#' Summarize a genome object
#' 
#' @param x An object of class \code{genome}.
#' 
#' @export
#' 
summary.genome <- function(x) {
  
  # Extract information
  n_chr <- nchr(x)
  len <- x$len
  n_mar <- nmar(x, by.chr = TRUE)
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
      
      cat("\nSummary for trait: ", t, "\n")
      
      n_qtl <- nrow(gen_model[[t]])
      qtl_chr <- table(gen_model[[t]][,1])

      cat("Number of QTL: ", n_qtl, "\n")
      cat("QTL per chromosome: ", qtl_chr, "\n") 
      cat("Distribution of additive effects:\n")
      print(summary(gen_model[[t]][,3]))
      cat("Distribution of dominance effects: \n")
      print(summary(gen_model[[t]][,4]))
      
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
