#' Fast Cooccur - compute exact cooccurence betwenn species pairs
#'
#' @param mat A apecies (rows) x sites (column) matrix.
#' @param chunks The number of chunks into which the data should be split.
#'
#' @return
#' A cooccur object
#' @export
#'
#' @examples
fast_cooccur <- function(mat, chunks = 2) {

  cat("Preparing for analysis \n")
  spp_site_mat <- mat
  # HANDLE ARGUEMENTS
  true_rand_classifier = 0.1
  spp_key <- data.frame(num=1:nrow(spp_site_mat),spp=row.names(spp_site_mat))

  # ORGANIZE & INITIALIZE FOR ANALYSIS
  spp_site_mat[spp_site_mat>0] <- 1
  tsites <- ncol(spp_site_mat)
  nspp <- nrow(spp_site_mat)
  spp_pairs <- choose(nspp,2)

  cat("Generating", spp_pairs, "species pairs\n")

  sp.df <- t(combn(nspp,2, simplify = TRUE))
  sp.df <- data.frame(spp = sp.df[,1], spp_next = sp.df[,2])
  ncores <- availableCores() -1
  myl <- spp_pairs
  chunksize <- floor(myl/chunks)

  # chunks <- 1:(ncores - 1)
  split.group <- c(rep(1:(chunks -1), each = chunksize), rep(chunks, myl-(chunks-1)*chunksize))

  cat("Splitting in", chunks, "chunks of roughly", chunksize, "pairs\n\n")

  #define fonction to compute probas to be passed to purrr::map
  get.proba <- function(df, mat){
    out <- as.data.frame(get_probaC(df, mat))
    colnames(out) <- c("spp", "spp_next", "sp1_inc", "sp2_inc", "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt")

    return(out)
  }

  #launching the computation
  cat("Computing probas\n")

  output <- sp.df %>%
    split(split.group) %>%
    furrr::future_map_dfr(get.proba, mat = mat, .progress = TRUE) %>%
    filter(exp_cooccur >= 1)

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
  output_list$spp.names = row.names(mat)

  output_list$omitted <- n_omitted
  output_list$pot_pairs <- spp_pairs

  class(output_list) <- "cooccur"
  return(output_list)

}
