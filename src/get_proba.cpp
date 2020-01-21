#include <Rcpp.h>
using namespace Rcpp;
//' Compute cooccurences probabilities over slices
//'
//' @param df a DataFrame of species pairs
//' @param mat a species matrix
//'
//' @details Internal function. Not to be used by itself
//' @useDynLib FastCooccur
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericMatrix get_probaC(DataFrame df, NumericMatrix mat){

  NumericVector spp = df[0];
  NumericVector spp_next = df[1];
  int chunkSize = spp.length();
  NumericMatrix results(chunkSize, 9);
  int tsites = mat.ncol();


  //declare variables
  int Mspp;
  int Mspp_next;
  double sp1_inc;
  double sp2_inc;
  double prob_cooccur;
  int obs_cooccur;
  double exp_cooccur;
  int max_inc;
  int min_inc;
  int psite;
  NumericVector prob_share_site;
  NumericVector allProb;
  double p_lt;
  double p_gt;

  for(int i = 0; i < chunkSize; i++){

    Mspp = spp(i);
    Mspp_next = spp_next(i);
    sp1_inc = sum(mat(Mspp -1, _)); // ORGIGINAL incidence[incidence[,1]==sp1,2]
    sp2_inc = sum(mat(Mspp_next -1 , _)); // ORGIGINAL incidence[incidence[,1]==sp2,2]
    prob_cooccur = (sp1_inc/tsites) * (sp2_inc/tsites);

    obs_cooccur = sum((mat(Mspp - 1, _ ) + mat(Mspp_next - 1, _ )) > 1);
    exp_cooccur = prob_cooccur * tsites;
    max_inc = max(NumericVector::create(sp1_inc,sp2_inc));
    min_inc = min(NumericVector::create(sp1_inc,sp2_inc));
    psite = tsites+1;
    prob_share_site = rep(0, psite);
    allProb = Rcpp::phyper(seq(0, min_inc), min_inc, tsites - min_inc, max_inc, true, false);
    prob_share_site[0] = allProb[0];

    for(int j = 1; j < allProb.size(); j++){
      prob_share_site[j] = allProb[j]-allProb[j-1];
    }

    p_lt = 0.0;
    p_gt = 0.0;

    for (int j = 0 ; j < tsites; j++){

      if(j <= obs_cooccur){
        p_lt += prob_share_site[(j)];
      }
      if (j >= obs_cooccur){
        p_gt += prob_share_site[(j)];
      }
      if (j == obs_cooccur){
        double p_exactly_obs = prob_share_site[(j)];

      }
    }

    //p_lt <- round(p_lt,5)
    //p_gt <- round(p_gt,5)
    //p_exactly_obs <- round(p_exactly_obs,5)

    results(i,0) = Mspp;
    results(i,1) = Mspp_next;
    results(i,2) = sp1_inc;
    results(i,3) = sp2_inc;
    results(i,4) = obs_cooccur;
    results(i,5) = prob_cooccur;
    results(i,6) = exp_cooccur;
    results(i,7) = p_lt;
    results(i,8) = p_gt;
  }

  return  results;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
ncol = 20
nrow = 100


mat <- matrix(sample(c(0, 1), ncol*nrow, replace = TRUE), ncol = ncol, nrow = nrow)
row.names(mat) <- 1:nrow
sp.df <- t(combn(nrow,2, simplify = TRUE))
df <- data.frame(spp = sp.df[,1], spp_next = sp.df[,2])
get_probaC(df, mat)
*/
