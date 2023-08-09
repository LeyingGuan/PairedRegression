#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;

// [[Rcpp::export]]
List permutation_conformal_C(arma::mat Xresid, arma::vec y,  arma::mat U, arma::umat perm_idx, std::string type = "coef") {
  
  arma::vec sum2 = sum(square(Xresid), 0).t();
  int p = Xresid.n_cols;
  int B = perm_idx.n_cols;
  int n = Xresid.n_rows;
  
  arma::vec pvals(p, arma::fill::zeros);
  arma::vec pvals_pos(p, arma::fill::zeros);
  arma::vec pvals_neg(p, arma::fill::zeros);
  arma::vec minusB(p, arma::fill::zeros);
  arma::vec minusB1(p, arma::fill::zeros);
  arma::vec yresid = y - U * (U.t() * y);
  for(int b = 0; b < B; b++){
    arma::uvec current_perm_idx = perm_idx.col(b);
    arma::mat Uperm = U.rows(current_perm_idx);
    
    arma::vec yresid_perm = y - Uperm * (Uperm.t() * y);
    arma::mat X2 = Xresid.rows(current_perm_idx);
    arma::vec prodb0 = Xresid.t() * yresid_perm;
    arma::vec prodb = X2.t() * yresid;
    
    arma::vec yresid_ = yresid - Uperm * (Uperm.t() * yresid);
    arma::vec yresid_perm_ = yresid_perm - U * (U.t() * yresid_perm);
    arma::vec tval_b(p, arma::fill::zeros), tval_b0(p, arma::fill::zeros);
    if(type == "coef"){
      tval_b += prodb;
      tval_b0 += prodb0;
    }else{
      arma::mat aa0(n, p);
      arma::mat aa1(n, p);
      
      for(int j = 0; j < p; j++){
        if(sum2(j) > 0){
          aa0.col(j) = yresid_perm_ - prodb0[j] / sum2(j) * Xresid.col(j);
          aa1.col(j) = yresid_ - prodb[j] / sum2(j) * X2.col(j);
        }else{
          aa0.col(j) = yresid_perm_;
          aa1.col(j) = yresid_;
        }
      }
      
      arma::vec aa0_vec = sum(square(aa0), 0).t();
      arma::vec aa1_vec = sum(square(aa1), 0).t();
      
      tval_b += prodb;
      tval_b0 += prodb0;
      
      for(int j = 0; j < p; j++){
        if(std::abs(prodb(j)) > 0 && std::abs(prodb0(j)) > 0 && (type != "coef")){
          tval_b(j) = prodb(j) / std::sqrt(aa1_vec(j));
          tval_b0(j) = prodb0(j)/ std::sqrt(aa0_vec(j));
        }
      } 
    }
    arma::uvec ii = find(abs(tval_b0) < abs(tval_b));
    arma::uvec ii3 = find(abs(tval_b0) == abs(tval_b));
    arma::uvec ii1 = find(tval_b0 < tval_b);
    arma::uvec ii2 = find(tval_b0 > tval_b);
    arma::uvec ii4 = find(tval_b0 == tval_b);
    if(!ii.is_empty()){
      pvals(ii) += 1;
    }
    if(!ii3.is_empty()){
      pvals(ii3) += 0.5;
      minusB(ii3) -= 0.5;
    }
    
    if(!ii1.is_empty()){
      pvals_pos(ii1) += 1;
    }
    
    if(!ii2.is_empty()){
      pvals_neg(ii2) += 1;
    }
    if(!ii4.is_empty()){
      pvals_pos(ii4) += 0.5;
      pvals_neg(ii4) += 0.5;
      minusB1(ii4) -= 0.5;
    }
  }
  pvals = (pvals + 1) / (B +  minusB + 1);
  pvals_neg = (pvals_neg + 1) / (B + minusB1 + 1);
  pvals_pos = (pvals_pos + 1) / (B + minusB1 + 1);
  
  return List::create(Named("unsigned") = pvals,
                      Named("pos") = pvals_pos,
                      Named("neg") = pvals_neg);
}

// // [[Rcpp::export]]
// List permutation_conformal_C(arma::mat Xresid, arma::vec yfitted, arma::vec yresid, arma::mat U, arma::umat perm_idx, std::string type = "coef") {
//   
//   arma::vec sum2 = sum(square(Xresid), 0).t();
//   int p = Xresid.n_cols;
//   int B = perm_idx.n_cols;
//   int n = Xresid.n_rows;
//   
//   arma::vec pvals(p, arma::fill::zeros);
//   arma::vec pvals_pos(p, arma::fill::zeros);
//   arma::vec pvals_neg(p, arma::fill::zeros);
//   arma::vec minusB(p, arma::fill::zeros);
//   arma::vec minusB1(p, arma::fill::zeros);
//   for(int b = 0; b < B; b++){
//     arma::uvec current_perm_idx = perm_idx.col(b);
//     arma::mat Uperm = U.rows(current_perm_idx);
//     
//     arma::vec ybresid = yresid - Uperm * (Uperm.t() * yresid);
//     arma::mat X2 = Xresid.rows(current_perm_idx);
//     arma::vec prodb0 = Xresid.t() * ybresid;
//     arma::vec prodb = X2.t() * ybresid;
//     
//     arma::vec tval_b(p, arma::fill::zeros), tval_b0(p, arma::fill::zeros);
//     
//     if(type == "coef"){
//       tval_b += prodb;
//       tval_b0 += prodb0;
//     }else{
//       arma::mat aa0(n, p);
//       arma::mat aa1(n, p);
//       
//       for(int j = 0; j < p; j++){
//         if(sum2(j) > 0){
//           aa0.col(j) = ybresid - prodb0[j] / sum2(j) * Xresid.col(j);
//           aa1.col(j) = ybresid - prodb[j] / sum2(j) * X2.col(j);
//         }else{
//           aa0.col(j) = ybresid;
//           aa1.col(j) = ybresid;
//         }
//       }
//       
//       arma::vec aa0_vec = sum(square(aa0), 0).t();
//       arma::vec aa1_vec = sum(square(aa1), 0).t();
//       
//       tval_b += prodb;
//       tval_b0 += prodb0;
//       
//       for(int j = 0; j < p; j++){
//         if(std::abs(prodb(j)) > 0 && std::abs(prodb0(j)) > 0 && (type != "coef")){
//           tval_b(j) = prodb(j) / std::sqrt(aa1_vec(j));
//           tval_b0(j) = prodb0(j)/ std::sqrt(aa0_vec(j));
//         }
//       } 
//     }
//     
//     arma::uvec ii = find(abs(tval_b0) < abs(tval_b));
//     arma::uvec ii3 = find(abs(tval_b0) == abs(tval_b));
//     arma::uvec ii1 = find(tval_b0 <= tval_b);
//     arma::uvec ii2 = find(tval_b0 >= tval_b);
//     arma::uvec ii4 = find(tval_b0 == tval_b);
//     if(!ii.is_empty()){
//       pvals(ii) += 1;
//     }
//     if(!ii3.is_empty()){
//       pvals(ii3) += 0.5;
//       minusB(ii3) -= 0.5;
//     }
//     
//     if(!ii1.is_empty()){
//       pvals_pos(ii1) += 1;
//     }
//     
//     if(!ii2.is_empty()){
//       pvals_neg(ii2) += 1;
//     }
//     // if(!ii4.is_empty()){
//     //   pvals_pos(ii4) += 1.0;//0.5;
//     //   pvals_neg(ii4) += 1.0;//0.5;
//     //   minusB1(ii4) -= 0.0;//0.5;
//     // }
//   }
//   
//   pvals = (pvals + 1) / (B +  minusB + 1);
//   pvals_neg = (pvals_neg + 1) / (B + minusB1 + 1);
//   pvals_pos = (pvals_pos + 1) / (B + minusB1 + 1);
//   
//   return List::create(Named("unsigned") = pvals,
//                       Named("pos") = pvals_pos,
//                       Named("neg") = pvals_neg);
// }



// [[Rcpp::export]]
List permutation_FL_C(arma::mat Xresid, arma::vec yfitted, arma::vec yresid, arma::mat U, arma::umat perm_idx, std::string type = "coef") {
  
  arma::vec sum2 = sum(square(Xresid), 0).t();
  int p = Xresid.n_cols;
  int B = perm_idx.n_cols;
  int n = Xresid.n_rows;
  
  arma::vec pvals(p, arma::fill::zeros);
  arma::vec pvals_pos(p, arma::fill::zeros);
  arma::vec pvals_neg(p, arma::fill::zeros);
  arma::vec prodb0 = Xresid.t() * yresid;
  
  arma::mat aa0(n, p);
  for(int j = 0; j < p; j++){
    if(sum2(j) > 0){
      aa0.col(j) = yresid - prodb0[j] / sum2(j) * Xresid.col(j);
    }else{
      aa0.col(j) = yresid;
    }
  }
  arma::vec aa0_vec = sum(square(aa0), 0).t();

  for(int b = 0; b < B; b++){
    arma::uvec current_perm_idx = perm_idx.col(b);
    arma::mat Uperm = U.rows(current_perm_idx);
    arma::mat X2 = Xresid.rows(current_perm_idx);
    
    arma::vec ybresid = yresid - Uperm * (Uperm.t() * yresid);
    arma::vec prodb = X2.t() * ybresid;
    
    arma::vec tval_b(p), tval_b0(p);
    
    if(type == "coef"){
      tval_b = prodb;
      tval_b0 = prodb0 * 1.0;
    }else{
      arma::mat aa1(n, p);
      for(int j = 0; j < p; j++){
        if(sum2(j) > 0){
          aa1.col(j) = ybresid - prodb[j] / sum2(j) * X2.col(j);
        }else{
          aa1.col(j) = ybresid;
        }
      }
      arma::vec aa1_vec = sum(square(aa1), 0).t();
      tval_b = prodb;
      tval_b0 = prodb0 * 1.0;
      for(int j = 0; j < p; j++){
        if(std::abs(prodb(j)) > 0 && std::abs(prodb0(j)) > 0){
          tval_b(j) = prodb(j) / std::sqrt(aa1_vec(j));
          tval_b0(j) = prodb0(j)/ std::sqrt(aa0_vec(j));
        }
      } 
    }
    arma::uvec ii = find(abs(tval_b0) <= abs(tval_b));
    if(!ii.is_empty()){
      pvals(ii) += 1;
    }
    
    arma::uvec ii1 = find(tval_b0 <= tval_b);
    arma::uvec ii2 = find(tval_b0 >= tval_b);
    
    if(!ii1.is_empty()){
      pvals_pos(ii1) += 1;
    }
    
    if(!ii2.is_empty()){
      pvals_neg(ii2) += 1;
    }
  }
  
  pvals = (pvals + 1) / (B + 1);
  pvals_neg = (pvals_neg + 1) / (B + 1);
  pvals_pos = (pvals_pos + 1) / (B + 1);
  
  return List::create(Named("unsigned") = pvals,
                      Named("pos") = pvals_pos,
                      Named("neg") = pvals_neg);
}


// [[Rcpp::export]]
List permutation_simple_C(arma::mat X, arma::vec yresid, arma::mat U, arma::umat perm_idx) {
  int p = X.n_cols;
  int B = perm_idx.n_cols;
  int n = X.n_rows;
  
  arma::mat Xresid = X -  U * (U.t() * X);
  arma::vec prod = Xresid.t() *  yresid;
  arma::vec sum2 = sqrt(sum(square(Xresid), 0).t());
  
  arma::vec pvals(p, arma::fill::zeros);
  arma::vec pvals_pos(p, arma::fill::zeros);
  arma::vec pvals_neg(p, arma::fill::zeros);
  
  for(int b = 0; b < B; b++){
    arma::uvec current_perm_idx = perm_idx.col(b);
    arma::mat Xperm = X.rows(current_perm_idx);
    arma::mat Xpermresid = Xperm -  U * (U.t() * Xperm);
    arma::vec prodb = Xpermresid.t() *  yresid;
    arma::vec sum2_perm = sqrt(sum(square(Xpermresid), 0).t());
    arma::vec tval = prod * 1.0;
    arma::vec tval_b = prodb * 1.0;
    for(int j = 0; j < p; j++){
      if(std::abs(prodb(j)) > 0 && std::abs(prod(j)) > 0){
        tval_b(j) = prodb(j) / sum2_perm(j);
        tval(j) = prod(j)/ sum2(j);
      }
    }
    arma::uvec ii = find(abs(tval) <= abs(tval_b));
    if(!ii.is_empty()){
      pvals(ii) += 1;
    }
    
    arma::uvec ii1 = find(tval < tval_b);
    arma::uvec ii2 = find(tval > tval_b);
    if(!ii1.is_empty()){
      pvals_pos(ii1) += 1;
    }
    
    if(!ii2.is_empty()){
      pvals_neg(ii2) += 1;
    }
  }
  
  pvals = (pvals + 1) / (B + 1);
  pvals_neg = (pvals_neg + 1) / (B + 1);
  pvals_pos = (pvals_pos + 1) / (B + 1);
  
  return List::create(Named("unsigned") = pvals,
                      Named("pos") = pvals_pos,
                      Named("neg") = pvals_neg);
}


