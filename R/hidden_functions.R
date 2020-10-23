#' Hidden functions; generally not to be called by the user
#' 
#' 
calc_genoval <- function(genome, geno) {
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = geno)) stop("The geno did not pass. See warning for reason.")
  
  # If geno is a list, recombine
  if (is.list(geno)) geno <- do.call("cbind", geno)
  
  ## Iterate over traits in the genetic model
  geno_val <- lapply(X = genome$gen_model, FUN = function(qtlmod) {
    qtl_geno <- pull_genotype(genome = genome, geno = geno, loci = qtlmod$qtl_name)
    # Subtract 1
    qtl_geno1 <- qtl_geno - 1
    
    # Additive effect
    bv <- qtl_geno1 %*% as.matrix(qtlmod$add_eff)
    # Dominance - first convert matrix to dominance matrix, then calculate deviations
    qtl_dom_geno <- ifelse(qtl_geno1 != 0, 0, 1)
    dd <- qtl_dom_geno %*% as.matrix(qtlmod$dom_eff)
    
    # Sum
    gv <- bv + dd 
    
    # If qxe effects are present, calculate the genotype-specific slope
    if (!is.null(qtlmod$gxe_slope)) {
      g_beta <- qtl_geno1 %*% as.matrix(qtlmod$gxe_slope)
    
    } else {
      g_beta <- NULL
      
    }
    
    # return a list
    list(gv = gv, g_beta = g_beta)
    
  })
  
  # Bind columns
  geno_val1 <- do.call("cbind", lapply(X = geno_val, "[[", "gv"))
  g_beta1 <- do.call("cbind", lapply(X = geno_val, "[[", "g_beta"))
  
  # Add trait names
  colnames(geno_val1) <- paste("trait", seq(ncol(geno_val1)), sep = "")
  
  # Convert the row.names to a column and return the data.frame
  gv <- data.frame(ind = row.names(geno_val1), geno_val1, row.names = NULL, stringsAsFactors = FALSE)
  gv <- gv[order(gv$ind),]
  
  # Do the same for the betas, if they exist
  if (!is.null(g_beta1)) {
    colnames(g_beta1) <- paste("trait", seq(ncol(g_beta1)), sep = "")
    
    gbeta <- data.frame(ind = row.names(g_beta1), g_beta1, row.names = NULL, stringsAsFactors = FALSE)
    gbeta <- gbeta[order(gbeta$ind),]
    
  } else{
    gbeta <- NULL
    
  }
  
  list(gv = gv, gbeta = gbeta)

  
} # Close the function






