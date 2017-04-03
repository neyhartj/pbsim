#' Create a crossing block from a vector of parents
#' 
#' @description 
#' Create a crossing block based on parent names, with different options for
#' mating scheme.
#' 
#' @param parents A \code{character} of line names to use as parents. If
#' \code{second.parents} is not provided, crosses are assigned from randomly sampling
#' the entries in \code{parents}.
#' @param second.parents A \code{character} of line names to use as parents. If
#' not \code{NULL}, must be the same length as \code{parents}. Crosses are formed
#' by randomly pairing the entries in \code{parents} with those in 
#' \code{second.parents}, unless \code{scheme = "pass"} is specified.
#' @param n.crosses The number of crosses to generate. Cannot be more than the
#' possible number of crosses.
#' @param scheme The mating scheme. Can be one of \code{"random"}, 
#' \code{"all.pairwise"}, \code{"chain"}, or \code{"pass"}. See \code{Details} 
#' for more information on the rules of these schemes.
#' @param use.parents.once \code{Logical} - should parents be used only once?
#' 
#' @details 
#' Several options are available to generate crossing blocks from a list of
#' parents. Here are the rules used for generating different crossing blocks.
#' 
#' If \code{second.parents = NULL}:
#' 
#' \itemize{
#'   \item{If \code{scheme = "random"} (default), crosses are randomly created
#'   using the parents in \code{parents}. Reciprocal crosses are excluded. If
#'   \code{n.crosses} is less than the total number of possible crosses
#'   (\eqn{total crosses = (n * (n - 1)) / 2}), then \code{n.crosses} crosses are randomly
#'   sampled. If \code{use.parents.once = TRUE}, then the list of crosses is
#'   trimmed such that parents that were already used once are not used again. This
#'   is ignored if the final number of crosses is less than \code{n.crosses}.}
#'   \item{If \code{scheme = "random"}}
#' 
#' }
#' 
#' 
#' @import dplyr
#' 
#' @export
#' 
sim_crossing_block <- function(parents, second.parents = NULL, n.crosses,
                               scheme = c("random", "all.pairwise", "chain", "pass"), 
                               use.parents.once = FALSE) {
  
  # Error
  if (!is.character(first.parents)) 
    stop("The input 'first.parents' must be a character.")
  
  if (!is.null(second.parents)) 
    if (!is.character(second.parents))
      stop("The input 'second.parents' must be a character.")
  
  # Is method a correct choice?
  if (!method %in% c("random", all.pairwise))
    stop("The argument 'method' must be 'random' or 'all.pairwise.'")
  
  # Options if the second.parents are NULL
  if (is.null(second.parents)) {
    
    # Generate all pairwise crosses minus reciprocals
    sample.crosses <- t(combn(x = first.parents, m = 2))
    
    # If the method is all pairwise, return all pairwise
    if (method == "all.pairwise") {
      chosen.crosses <- sample.crosses
      
    } else {
      # Randomly sample n.crosses
      crosses.ind <- sort(sample(seq_len(nrow(sample.crosses)), n.crosses))
      chosen.crosses <- sample.crosses[crosses.ind,]
      
    }
    
  }
    
  
} # Close the function
    
#   
#   
#   
#   
#   # If using parents once, make sure that there are enough crosses
#   if (use.parents.once) {
#     
#     # Find duplicated parents
#     
#     
#     
#     
#   }
#     
# 
#     
#     
#     
#   
#   # Generate all pairwise crossses
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   parent1.lines, parent2.lines, n.crosses, 
#                                 method = "random", use.parents.once = FALSE) {
#   
#   # Find the intersection of parent1.lines and parent2.lines
#   par.intersect <- intersect(parent1.lines, parent2.lines)
#   n.intersect <- length(par.intersect)
#   # Find the number of parent1 and parent2 lines
#   n.par1 <- length(parent1.lines)
#   n.par2 <- length(parent2.lines)
#   
#   if ( all(n.intersect == n.par1, n.intersect == n.par2) ) {
#     n.possible.crosses <- (n.par1 * (n.par2 - 1)) / 2
#     same.lines = T
#   } else {
#     n.possible.crosses <- n.par1 * n.par2
#     same.lines = F
#   }
#   
#   # Create all pairwise crosses
#   sample.crosses <- expand.grid(parent1.lines, parent2.lines) %>%
#     select(Var2, Var1)
#   
#   # Remove selfs
#   sample.crosses1 <- sample.crosses %>%
#     filter(apply(X = sample.crosses, MARGIN = 1, FUN = function(cross) length(unique(cross)) > 1))
#   
#   # If statements for methods
#   if (method == "all.pairwise") {
#     
#     return(sample.crosses1)
#     
#   } else {
#     
#     if (n.crosses > n.possible.crosses) stop("n.crosses is more than possible.")
#     
#     # If parents should only be used once, sample without replacement all lines 
#     # and put them into a matrix
#     
#     if (use.parents.once) {
#       
#       # Sample into a matrix
#       if (same.lines) {
#         
#         # Create random crosses from the parent lines
#         random.crosses <- sample(parent1.lines) %>%
#           matrix(ncol = 2, byrow = T) %>%
#           as.data.frame()
#         
#         # Sample n.crosses from the random crosses
#         selected.crosses <- random.crosses %>%
#           sample_n(size = n.crosses)
#         
#       } else {
#         
#         # Are there any overlapping lines in the two parent vectors?
#         if (n.intersect == 0) {
#           
#           # If not just sample each vector into a matrix
#           random.crosses <- cbind( sample(parent1.lines), sample(parent2.lines) ) %>%
#             as.data.frame()
#           
#           # Sample n.crosses from the random crosses
#           selected.crosses <- random.crosses %>%
#             sample_n(size = n.crosses)
#           
#         } else {
#           stop("Function not available.")
#         }
#       }
#       
#     } else {
#       
#       # Remove reciprocal
#       sample.crosses2 <- sample.crosses1 %>%
#         filter(!apply(X = sample.crosses1, MARGIN = 1, FUN = sort) %>% 
#                  t() %>% 
#                  duplicated() )
#       
#       # Sample n.crosses from the sample crosses
#       selected.crosses <- sample.crosses2 %>%
#         sample_n(size = n.crosses)
#       
#     }
#     
#     # Rename
#     colnames(selected.crosses) <- c("Parent1", "Parent2")
#     row.names(selected.crosses) <- NULL
#     return(selected.crosses)
#   }
#   
# } # Close the function