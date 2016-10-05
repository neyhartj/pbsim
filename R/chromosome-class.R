#' Class "chromosome"
#' 
#' @description A class of one or many chromosomes found within the class \code{genome}. Each chromosome has information as to the length, number of loci, and number and effect of qtl.
#' 
#' @slot chr.len An object of class \code{numeric} giving the length of the chromosome in Morgans.
#' @slot n.snps An object of class \code{integer} giving the number of SNP loci found on the chromosome.
#' @slot pos.snps An object of class \code{numeric} giving the Morgan position of each SNP locus found on the chromosome.
#' @slot n.markers An object of class \code{integer} giving the number of SNP loci that are designated as markers (i.e. not QTL).
#' @slot pos.markers An object of class \code{numeric} giving the Morgan position of each SNP marker locus found on the chromosome.
#' @slot n.add.qtl An object of class \code{integer} giving the number of QTL that exhibit an additive effect. Since all QTL exhibit an additive effect, this is a proxy for the number of QTL.
#' @slot n.dom.qtl An object of class \code{integer} giving the number of QTL that exhibit some dominance effect.
#' @slot pos.add.qtl An object of class \code{numeric} giving the Morgan position of each QTL that exhibits an additive effect. Since all QTL exhibit an additive effect, this is a proxy for the location of all QTL.
#' @slot pos.dom.qtl An object of class \code{numeric} giving the Morgan position of each QTL that exhibits a dominance effect.
#' @slot qtl.effects An object of class \code{list} giving the additive and domiance effects of all QTL.
#' 
chromosome <- setClass("chromosome", 
                   representation(
                     chr.len = "numeric",
                     n.snps = "integer",
                     pos.snps = "numeric",
                     n.markers = "integer",
                     pos.markers = "list",
                     n.add.qtl = "integer",
                     n.dom.qtl = "integer",
                     pos.add.qtl = "list",
                     pos.dom.qtl = "list",
                     qtl.effects = "list") )

