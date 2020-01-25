#' Fast Cooccur - compute exact cooccurence betwenn species pairs
#'
#' @details Fast Cooccur relies on \code{\link[future]{future}} and \code{\link[furrr]{future_map_dfr}} to handle parallel computing.
#' This way fast_cooccur is backend - agnostic. You can specify the bakend with the function \code{\link[future]{plan}}.
#' By default, it runs with \code{plan(sequential)}, which run as a single R process. To go parallel on a desktop machine, call \code{plan(multiprocess)}
#' before calling \code{fast_cooccur}. On windows machines, it will fork the current R Session. On Unix machines (OSX nd Linux), it will go multicore.
#'
#' If you need to process very large matrices in a short amount of time, it is possible to run the function on a HPC configured to run the package future.
#'
#' fast_coocur() will perform a strict sanity check on spp_site_mat before launching computation, and will stop execution if any of the conditions are not met.
#' Please ensure that :
#' \enumerate{
#'   \item spp_site_mat is a matrix (not a dataframe)
#'   \item all species are present at least once in the dataset
#'   \item you can give row.names to yourmatrix. It might be useful for downstream analysis.
#' }
#'
#' @param spp_site_mat A apecies (rows) x sites (column) matrix.
#' @param chunks The number of chunks into which the data should be split.
#' If yorking in parallel mode, it is wise to keep it a multiple of the computing units (cores or workers).
#' @param verbose TRUE or FALSE. If TRUE, display some info about the computation
#' @param progress TRUE or FALSE. If TRUE, display a progress bar in parallel mode. Only meaningful if chunks is greather
#' than the number of computing units
#'
#' @return
#' A cooccur object
#' @export
#'
#' @examples
fast_cooccur <- function(spp_site_mat, chunks = 2, verbose = TRUE, progress = TRUE) {

  if(verbose) cat("Sanity check ...\n")
  sanity_check(spp_site_mat, TRUE)

  if(verbose) cat("Preparing for analysis \n")
  spp_site_mat[spp_site_mat>0] <- 1
  # HANDLE ARGUEMENTS
  true_rand_classifier = 0.1
  spp_key <- data.frame(num=1:nrow(spp_site_mat),spp=row.names(spp_site_mat))

  # ORGANIZE & INITIALIZE FOR ANALYSIS
  spp_site_mat[spp_site_mat>0] <- 1
  tsites <- ncol(spp_site_mat)
  nspp <- nrow(spp_site_mat)
  spp_pairs <- choose(nspp,2)

  if(verbose) cat("Generating", spp_pairs, "species pairs\n")

  sp.df <- t(utils::combn(nspp,2, simplify = TRUE))
  sp.df <- data.frame(spp = sp.df[,1], spp_next = sp.df[,2])
  ncores <- future::availableCores() -1
  myl <- spp_pairs
  chunksize <- floor(myl/chunks)

  # chunks <- 1:(ncores - 1)
  split.group <- c(rep(1:(chunks -1), each = chunksize), rep(chunks, myl-(chunks-1)*chunksize))

  if(verbose) cat("Splitting in", chunks, "chunks of roughly", chunksize, "pairs\n\n")

  #define fonction to compute probas to be passed to purrr::map
  get.proba <- function(df, mat){
    out <- as.data.frame(get_probaC(df, mat))
    colnames(out) <- c("spp", "spp_next", "sp1_inc", "sp2_inc", "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt")

    return(out)
  }

  #launching the computation
  if(verbose) cat("Computing Cooccurences\n")

  output <- sp.df %>%
    split(split.group) %>%
    furrr::future_map_dfr(get.proba, mat = spp_site_mat, .progress = progress) %>%
    dplyr::filter(.data$exp_cooccur >= 1)

  n_omitted <- spp_pairs - nrow(output)

  # clasifying true randoms
  true_rand <- (nrow(output[(output$p_gt >= 0.05 & output$p_lt >= 0.05) & (abs(output$obs_cooccur - output$exp_cooccur) <= (tsites * true_rand_classifier)),]))

  # PREPARE AND DELIVER OUTPUT AS CLASS "cooccur"
  output_list <- list(call = match.call(),
                      results = output,
                      positive = nrow(output[output$p_gt < 0.05,]),
                      negative = nrow(output[output$p_lt < 0.05,]),
                      co_occurrences = (nrow(output[output$p_gt < 0.05 | output$p_lt < 0.05,])), #nrow(output[output$p_exactly_obs <= 0.05,]),
                      pairs = nrow(output),
                      random = true_rand,
                      unclassifiable = nrow(output) - (true_rand + nrow(output[output$p_gt < 0.05,]) + nrow(output[output$p_lt < 0.05,])),
                      sites = tsites,
                      species = nspp,
                      percent_sig = (((nrow(output[output$p_gt < 0.05 | output$p_lt < 0.05,]))) / (nrow(output))) * 100,
                      true_rand_classifier = true_rand_classifier
  )


  output_list$spp_key <- spp_key
  output_list$spp.names = row.names(spp_site_mat)

  output_list$omitted <- n_omitted
  output_list$pot_pairs <- spp_pairs

  class(output_list) <- "cooccur"
  return(output_list)

}
