// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// permutation_conformal_C
List permutation_conformal_C(arma::mat Xresid, arma::vec yfitted, arma::vec yresid, arma::mat U, arma::umat perm_idx, std::string type);
RcppExport SEXP _PREGS_permutation_conformal_C(SEXP XresidSEXP, SEXP yfittedSEXP, SEXP yresidSEXP, SEXP USEXP, SEXP perm_idxSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xresid(XresidSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yfitted(yfittedSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yresid(yresidSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::umat >::type perm_idx(perm_idxSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(permutation_conformal_C(Xresid, yfitted, yresid, U, perm_idx, type));
    return rcpp_result_gen;
END_RCPP
}
// permutation_FL_C
List permutation_FL_C(arma::mat Xresid, arma::vec yfitted, arma::vec yresid, arma::mat U, arma::umat perm_idx, std::string type);
RcppExport SEXP _PREGS_permutation_FL_C(SEXP XresidSEXP, SEXP yfittedSEXP, SEXP yresidSEXP, SEXP USEXP, SEXP perm_idxSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xresid(XresidSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yfitted(yfittedSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yresid(yresidSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::umat >::type perm_idx(perm_idxSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(permutation_FL_C(Xresid, yfitted, yresid, U, perm_idx, type));
    return rcpp_result_gen;
END_RCPP
}
// permutation_simple_C
List permutation_simple_C(arma::mat X, arma::vec yresid, arma::mat U, arma::umat perm_idx);
RcppExport SEXP _PREGS_permutation_simple_C(SEXP XSEXP, SEXP yresidSEXP, SEXP USEXP, SEXP perm_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yresid(yresidSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::umat >::type perm_idx(perm_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(permutation_simple_C(X, yresid, U, perm_idx));
    return rcpp_result_gen;
END_RCPP
}
