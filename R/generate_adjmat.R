#' Compute Adjacency Matrix from Cooccur object
#'
#' @param cooccur A cooccur object retruned by fast_cooccur
#'
#' @return A signed adjacency matrix
#' @export
#'
#' @examples
generate_adjmat <- function(cooccur){
  if(class(cooccur)!= "cooccur"){
    stop("coccur must be a cooccur object")
  }
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
