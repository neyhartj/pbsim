#' Summary
#' @export
#' 
summary.genome <- function(x) {
  
  # Extract information
  n.chr <- length(x$len)
  len <- x$len
  n.mar <- x$n.mar
  gen_model <- x$gen_model
  
  # Show
  cat("\nGenome summary \n\n")
  cat("Number of chromosomes: ", n.chr, "\n")
  cat("Length of chromosomes (in cM): ", len, "\n")
  cat("Markers per chromosome: ", n.mar, "\n")
  
  # Report the genetic model
  if (is.null(gen_model)) {
    cat("\nGenetic model summary: No genetic model")
    
  } else {
    
    n.qtl <- nrow(gen_model)
    qtl.chr <- table(gen_model[,1])
    
    cat("\nGenetic model summary\n")
    cat("Number of QTL: ", n.qtl, "\n")
    cat("QTL per chromosome: ", qtl.chr, "\n") 
    cat("Distribution of additive effects:\n")
    print(summary(gen_model[,3]))
    cat("Distribution of dominance effects: \n")
    print(summary(gen_model[,4]))
    
  }
  
}

#' Print
#' @export
#' 
print.genome <- function(x) summary(x)
