#' Determine what happens when you call on a "genome" or "chromosome" object
#' 
#' 

# First set genome
setMethod("show", signature(object = "genome"), function(object) {
  
  genome.summary <- list(n.chr = object@n.chr,
                         genome.len = object@genome.len,
                         n.snps = object@n.snps,
                         n.qtl = object@n.qtl,
                         qtl.add.eff = unlist(lapply(X = object@chromosomes, FUN = function(chr) chr@qtl.effects$a)),
                         qtl.dom.eff = unlist(lapply(X = object@chromosomes, FUN = function(chr) chr@qtl.effects$d)),
                         chr.len = sapply(X = object@chromosomes, FUN = function(chr) chr@chr.len),
                         n.perfect.snps = object@n.perfect.snps,
                         chromosomes = object@chromosomes
  )
  
  ##
  cat("\nThe genome object has the following attributes:\n")
  cat(paste("\nTotal genome length of ", genome.summary$genome.len, " M.\n", sep = ""))
  cat(paste(genome.summary$n.chr, "chromosomes of lengths"), genome.summary$chr.len, " M.\n\n")
  cat(genome.summary$n.snps," SNP loci plus ", genome.summary$n.perfect.snps," perfect SNP markers.\n", sep = "")
  cat(genome.summary$n.qtl," QTL across the genome.\n\n", sep = "")
  
  if(!all(genome.summary$qtl.add.eff == 0)) {
    cat("Summary of distribution of additive effects:\n\n")
    print(summary(genome.summary$qtl.add.eff))
    cat("\n")
  }
  
  if(!all(genome.summary$qtl.dom.eff == 0)) {
    cat("Summary of distribution of dominance effects:\n\n")
    print(summary(genome.summary$qtl.dom.eff))
    cat("\n")
  }
}
)

# Now set a chromosome
setMethod("show", signature(object = "chromosome"), function(object) {
  
  chr.summary <- list(chr.len = object@chr.len,
                      n.chr.snps = object@n.snps,
                      n.chr.markers = object@n.markers,
                      n.chr.qtl = object@n.add.qtl,
                      qtl.add.eff = object@qtl.effects$a,
                      qtl.dom.eff = object@qtl.effects$d)
  
  ##
  cat("\nThe chromosome object has the following attributes:\n")
  cat(paste("\nChromosome length of ", chr.summary$chr.len, " M.\n", sep = ""))
  cat(chr.summary$n.chr.snps, "SNP loci on the chromosome.\n")
  cat(chr.summary$n.chr.qtl, "QTL on the chromosome.\n\n")
  
  if(!all(chr.summary$qtl.add.eff == 0)) {
    cat("Summary of distribution of additive effects:\n\n")
    print(summary(chr.summary$qtl.add.eff))
    cat("\n")
  }
  
  if(!all(chr.summary$qtl.dom.eff == 0)) {
    cat("Summary of distribution of dominance effects:\n\n")
    print(summary(chr.summary$qtl.dom.eff))
    cat("\n")
  }
}
)
