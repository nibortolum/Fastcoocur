#' Compute Adjacency Matrix from Cooccur object
#'
#' @param cooccur A cooccur object retruned by fast_cooccur
#' @param p.adj a correction method for the p.value
#' @param verbose TRUE or FALSE. If TRUE, display meaningful messages about the progression
#'
#' @return A signed adjacency matrix
#' @export
#'
#' @examples
generate_adjmat <- function(cooccur, p.adj = "none", verbose = TRUE){
  if(class(cooccur)!= "cooccur"){
    stop("coccur must be a cooccur object")
  }
  if(p.adj != "none"){

    if(verbose) {cat("correction pvalues \n")}
    cooccur$results$p_gt <- p.adjust(cooccur$results$p_gt, method = p.adj)
    cooccur$results$p_lt <- p.adjust(cooccur$results$p_lt, method = p.adj)
  }
  if(verbose){cat("Rendering adjacency matrix\n")}
  positive <- as.matrix(cooccur$results[which(cooccur$results$p_gt < 0.05) ,c("spp", "spp_next")])
  negative <- as.matrix(cooccur$results[which(cooccur$results$p_lt < 0.05) ,c("spp", "spp_next")])

  #edgelist (a=origin, b=destination)
  adjmat <- matrix(0, cooccur$species, cooccur$species)
  adjmat[positive] <- 1
  adjmat[negative] <- -1

  adjmat[lower.tri(adjmat)] <- t(adjmat)[lower.tri(adjmat)]
  colnames(adjmat) <- row.names(adjmat) <- cooccur$spp.names

  return(adjmat)
}
