#' Class "genome"
#' 
#' @description A class of the genome. Used to store information for each chromosome
#' 
#' @slot n.chr An object of class \code{integer} giving the number of chromosomes in the genome.
#' @slot n.loci An object of class \code{integer} giving the total number of SNP loci in the genome.
#' @slot genome.len An object of class \code{numeric} giving the total Morgan length of the genome.
#' @slot n.qtl An object of class \code{integer} giving the total number of QTL in the genome.
#' @slot chromosomes An object of class \code{list} giving all the \code{chromosome} objects in the genome.
#' 
genome <- setClass("genome", 
                   representation(n.chr = "integer",
                                  n.snps = "integer",
                                  genome.len = "numeric",
                                  n.qtl = "integer",
                                  n.perfect.snps = "integer",
                                  chromosomes = "list"))