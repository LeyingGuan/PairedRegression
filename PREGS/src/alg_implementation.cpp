#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;

arma::mat dcSVD(arma::mat & X) {
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, X, "dc");
  
  int K = 0;
  for(int k = 0; k < S.n_elem; k++){
    if(S(k) > 1E-10){
      K +=1;
    }
  }
  arma::mat U0(U.n_rows, K, arma::fill::zeros);
  for(int k=0; k < K; k++){
    for(int i=0; i < U.n_rows; i++){
      U0(i,k) += U(i,k);
    }
    
  }
  return U0;
}

// [[Rcpp::export]]
List permutation_PREGSseparate(const arma::vec& x, const arma::mat& Y, arma::mat& Z,
                               const int B){
  arma::mat U = dcSVD(Z);
  arma::vec xresid = x - U*(U.t()*x);
  int N = Y.n_rows;
  int M = Y.n_cols;
  arma::mat yresid_for_perm = Y - U * (U.t() * Y);
  arma::uvec seq0 =arma::conv_to<arma::uvec>::from(arma::linspace<arma::vec>(0, N - 1, N));
  arma::uvec seq = arma::shuffle(seq0);
  arma::vec pvals(M, arma::fill::zeros);
  arma::vec pvals_pos(M, arma::fill::zeros);
  arma::vec pvals_neg(M, arma::fill::zeros);
  arma::vec minusB(M, arma::fill::zeros);
  arma::vec minusB1(M, arma::fill::zeros);
  for(int b = 0; b < B; b++){
    arma::mat Uperm = U.rows(seq) * 1.0;
    arma::vec xresid_perm = xresid(seq) * 1.0;
    arma::mat yresid_for_original = Y - Uperm * (Uperm.t() * Y);
    arma::vec prodb0 = yresid_for_original.t() * xresid;
    arma::vec prodb =  yresid_for_perm.t() * xresid_perm;
    yresid_for_original =yresid_for_original - U * (U.t() * yresid_for_original);
    yresid_for_perm = yresid_for_perm - Uperm * (Uperm.t() * yresid_for_perm);
    //calculate F-stat
    arma::mat yresid_for_original_x(N, M, arma::fill::zeros);
    arma::mat yresid_for_perm_x(N, M, arma::fill::zeros);
    double s2 = arma::accu(arma::square(xresid));
    arma::vec fstat_orginal(M, arma::fill::zeros);
    arma::vec fstat_perm(M, arma::fill::zeros);
    for(int m = 0; m < M; m++){
      if(s2 > 0.0){
        yresid_for_original_x.col(m) = yresid_for_original.col(m) - prodb0(m) * xresid/s2;
        yresid_for_perm_x.col(m) = yresid_for_perm.col(m) - prodb(m) * xresid_perm/s2;
      }else{
        yresid_for_original_x.col(m) = yresid_for_original.col(m);
        yresid_for_perm_x.col(m) = yresid_for_perm.col(m);
      }
      double tmp01 = arma::accu(arma::square(yresid_for_original.col(m)));
      double tmp02 = arma::accu(arma::square(yresid_for_perm.col(m)));
      double tmp1 = arma::accu(arma::square(yresid_for_original_x.col(m)));
      double tmp2 = arma::accu(arma::square(yresid_for_perm_x.col(m)));
      fstat_orginal(m) = tmp01 - tmp1;
      fstat_perm(m) = tmp02 - tmp2;
      if(tmp1 > 0.0 && tmp2 > 0.0){
        fstat_orginal(m) /= tmp1;
        fstat_perm(m) /= tmp2;
      }else if(tmp01 > 0.0 && tmp02 > 0.0){
        fstat_orginal(m) /= tmp01;
        fstat_perm(m) /= tmp02;
      }
      if(fstat_perm(m) > fstat_orginal(m)){
        pvals(m) += 1;
      }else if(fstat_perm(m) == fstat_orginal(m)){
        pvals(m) += .5;
        minusB(m) -= .5;
      }
    }
    seq = arma::shuffle(seq0);
  }
  pvals = (pvals + 1) / (B +  minusB + 1);
  pvals_neg = (pvals_neg + 1) / (B + minusB1 + 1);
  pvals_pos = (pvals_pos + 1) / (B + minusB1 + 1);
  
  return List::create(Named("unsigned") = pvals,
                      Named("pos") = pvals_pos,
                      Named("neg") = pvals_neg);
}

// [[Rcpp::export]]
List permutation_PREGSjoint(const arma::vec& x, const arma::mat& Y, const arma::mat& Z,
                               const int B){
  int N = Y.n_rows;
  int M = Y.n_cols;
  int p = Z.n_cols;
  arma::uvec seq0 =arma::conv_to<arma::uvec>::from(arma::linspace<arma::vec>(0, N - 1, N));
  arma::uvec seq = arma::shuffle(seq0);
  arma::vec pvals(M, arma::fill::zeros);
  arma::vec pvals_pos(M, arma::fill::zeros);
  arma::vec pvals_neg(M, arma::fill::zeros);
  arma::vec minusB(M, arma::fill::zeros);
  arma::vec minusB1(M, arma::fill::zeros);
  for(int b = 0; b < B; b++){
    arma::mat Zperm = Z.rows(seq) * 1.0;
    arma::mat Zconcate(N, p*2);
    for(int j = 0; j < p; j++){
      Zconcate.col(j) = Z.col(j) * 1.0;
      Zconcate.col(j+p) = Zperm.col(j) * 1.0;
    }
    arma::mat U = dcSVD(Zconcate);
    arma::mat yresid = Y - U * (U.t() * Y);
    arma::vec xresid = x - U * (U.t() * x);
    arma::vec xresid_perm = x(seq) -  U * (U.t() * x(seq));
    //*** mutually orthogonal
    arma::vec prodb0 = yresid.t() * xresid;
    arma::vec prodb =  yresid.t() * xresid_perm;
    //calculate F-stat
    arma::mat yresid_for_original_x(N, M, arma::fill::zeros);
    arma::mat yresid_for_perm_x(N, M, arma::fill::zeros);
    double s2_original = arma::accu(arma::square(xresid));
    double s2_perm = arma::accu(arma::square(xresid_perm));
    arma::vec fstat_orginal(M, arma::fill::zeros);
    arma::vec fstat_perm(M, arma::fill::zeros);
    for(int m = 0; m < M; m++){
      if(s2_original > 0.0){
        yresid_for_original_x.col(m) = yresid.col(m) - prodb0(m) * xresid/s2_original;
      }else{
        yresid_for_original_x.col(m) = yresid.col(m) * 1.0;
      }
      if(s2_perm > 0.0){
        yresid_for_perm_x.col(m) = yresid.col(m) - prodb(m) * xresid_perm/s2_perm;
      }else{
        yresid_for_perm_x.col(m) = yresid.col(m) * 1.0;
      }
      fstat_orginal(m) = arma::accu(arma::square(yresid.col(m)));
      fstat_perm(m)= fstat_orginal(m) * 1.0;
      fstat_orginal(m) -= arma::accu(arma::square(yresid_for_original_x.col(m)));
      fstat_perm(m) -= arma::accu(arma::square(yresid_for_perm_x.col(m)));
      if(fstat_perm(m) > fstat_orginal(m)){
        pvals(m) += 1;
      }else if(fstat_perm(m) == fstat_orginal(m)){
        pvals(m) += .5;
        minusB(m) -= .5;
      }
    }
    seq = arma::shuffle(seq0);
  }
  pvals = (pvals + 1) / (B +  minusB + 1);
  pvals_neg = (pvals_neg + 1) / (B + minusB1 + 1);
  pvals_pos = (pvals_pos + 1) / (B + minusB1 + 1);
  
  return List::create(Named("unsigned") = pvals,
                      Named("pos") = pvals_pos,
                      Named("neg") = pvals_neg);
}

// [[Rcpp::export]]
List permutation_FL(const arma::vec& x, const arma::mat& Y, arma::mat& Z,
                    const int B){
  arma::mat U = dcSVD(Z);
  arma::vec xresid = x - U*(U.t()*x);
  int N = Y.n_rows;
  int M = Y.n_cols;
  arma::mat yresid= Y - U * (U.t() * Y);
  arma::uvec seq0 =arma::conv_to<arma::uvec>::from(arma::linspace<arma::vec>(0, N - 1, N));
  arma::uvec seq = arma::shuffle(seq0);
  arma::vec pvals(M, arma::fill::zeros);
  arma::vec pvals_pos(M, arma::fill::zeros);
  arma::vec pvals_neg(M, arma::fill::zeros);
  arma::vec minusB(M, arma::fill::zeros);
  arma::vec minusB1(M, arma::fill::zeros);
  for(int b = 0; b < B; b++){
    arma::mat Uperm = U.rows(seq) * 1.0;
    arma::vec xresid_perm = xresid(seq) * 1.0;
    arma::mat yresid = Y - U * (U.t() * Y);
    arma::vec prodb0 = yresid.t() * xresid;
    arma::vec prodb =  yresid.t() * xresid_perm;
    arma::mat yresid_for_original =yresid * 1.0;
    arma::mat yresid_for_perm = yresid - Uperm * (Uperm.t() * yresid);
    //calculate F-stat
    arma::mat yresid_for_original_x(N, M, arma::fill::zeros);
    arma::mat yresid_for_perm_x(N, M, arma::fill::zeros);
    double s2 = arma::accu(arma::square(xresid));
    arma::vec fstat_orginal(M, arma::fill::zeros);
    arma::vec fstat_perm(M, arma::fill::zeros);
    for(int m = 0; m < M; m++){
      if(s2 > 0.0){
        yresid_for_original_x.col(m) = yresid_for_original.col(m) - prodb0(m) * xresid/s2;
        yresid_for_perm_x.col(m) = yresid_for_perm.col(m) - prodb(m) * xresid_perm/s2;
      }else{
        yresid_for_original_x.col(m) = yresid_for_original.col(m);
        yresid_for_perm_x.col(m) = yresid_for_perm.col(m);
      }
      double tmp01 = arma::accu(arma::square(yresid_for_original.col(m)));
      double tmp02 = arma::accu(arma::square(yresid_for_perm.col(m)));
      double tmp1 = arma::accu(arma::square(yresid_for_original_x.col(m)));
      double tmp2 = arma::accu(arma::square(yresid_for_perm_x.col(m)));
      fstat_orginal(m) = tmp01-tmp1;
      fstat_perm(m) = tmp02 - tmp2;
      if(tmp1 > 0.0 && tmp2 > 0.0){
        fstat_orginal(m) /= tmp1;
        fstat_perm(m) /= tmp2;
      }else if(tmp01 > 0.0 && tmp02 > 0.0){
        fstat_orginal(m) /= tmp01;
        fstat_perm(m) /= tmp02;
      }
      if(fstat_perm(m) >= fstat_orginal(m)){
        pvals(m) += 1;
      }
    }
    seq = arma::shuffle(seq0);
  }
  pvals = (pvals + 1) / (B +  minusB + 1);
  pvals_neg = (pvals_neg + 1) / (B + minusB1 + 1);
  pvals_pos = (pvals_pos + 1) / (B + minusB1 + 1);
  
  return List::create(Named("unsigned") = pvals,
                      Named("pos") = pvals_pos,
                      Named("neg") = pvals_neg);
}


// [[Rcpp::export]]
List permutation_vanilla(const arma::vec& x, const arma::mat& Y, arma::mat& Z,
                    const int B){
  arma::mat U = dcSVD(Z);
  arma::vec xresid = x - U*(U.t()*x);
  int N = Y.n_rows;
  int M = Y.n_cols;
  arma::mat yresid= Y - U * (U.t() * Y);
  arma::uvec seq0 =arma::conv_to<arma::uvec>::from(arma::linspace<arma::vec>(0, N - 1, N));
  arma::uvec seq = arma::shuffle(seq0);
  arma::vec pvals(M, arma::fill::zeros);
  arma::vec pvals_pos(M, arma::fill::zeros);
  arma::vec pvals_neg(M, arma::fill::zeros);
  arma::vec minusB(M, arma::fill::zeros);
  arma::vec minusB1(M, arma::fill::zeros);
  for(int b = 0; b < B; b++){
    arma::mat Uperm = U.rows(seq) * 1.0;
    arma::vec xresid_perm = x(seq) -U * (U.t() * x(seq));
    arma::mat yresid = Y - U * (U.t() * Y);
    arma::vec prodb0 = yresid.t() * xresid;
    arma::vec prodb =  yresid.t() * xresid_perm;
    arma::mat yresid_for_original =yresid * 1.0;
    arma::mat yresid_for_perm = yresid * 1.0;
    //calculate F-stat
    arma::mat yresid_for_original_x(N, M, arma::fill::zeros);
    arma::mat yresid_for_perm_x(N, M, arma::fill::zeros);
    double s2_original = arma::accu(arma::square(xresid));
    double s2_perm = arma::accu(arma::square(xresid_perm));
    arma::vec fstat_orginal(M, arma::fill::zeros);
    arma::vec fstat_perm(M, arma::fill::zeros);
    for(int m = 0; m < M; m++){
      if(s2_original > 0.0){
        yresid_for_original_x.col(m) = yresid.col(m) - prodb0(m) * xresid/s2_original;
      }else{
        yresid_for_original_x.col(m) = yresid.col(m) * 1.0;
      }
      if(s2_perm > 0.0){
        yresid_for_perm_x.col(m) = yresid.col(m) - prodb(m) * xresid_perm/s2_perm;
      }else{
        yresid_for_perm_x.col(m) = yresid.col(m) * 1.0;
      }
      double tmp01 = arma::accu(arma::square(yresid_for_original.col(m)));
      double tmp02 = arma::accu(arma::square(yresid_for_perm.col(m)));
      double tmp1 = arma::accu(arma::square(yresid_for_original_x.col(m)));
      double tmp2 = arma::accu(arma::square(yresid_for_perm_x.col(m)));
      fstat_orginal(m) = tmp01-tmp1;
      fstat_perm(m) = tmp02 - tmp2;
      if(tmp1 > 0.0 && tmp2 > 0.0){
        fstat_orginal(m) /= tmp1;
        fstat_perm(m) /= tmp2;
      }else if(tmp01 > 0.0 && tmp02 > 0.0){
        fstat_orginal(m) /= tmp01;
        fstat_perm(m) /= tmp02;
      }
      if(fstat_perm(m) > fstat_orginal(m)){
        pvals(m) += 1;
      }
    }
    seq = arma::shuffle(seq0);
  }
  pvals = (pvals + 1) / (B +  minusB + 1);
  pvals_neg = (pvals_neg + 1) / (B + minusB1 + 1);
  pvals_pos = (pvals_pos + 1) / (B + minusB1 + 1);
  
  return List::create(Named("unsigned") = pvals,
                      Named("pos") = pvals_pos,
                      Named("neg") = pvals_neg);
}