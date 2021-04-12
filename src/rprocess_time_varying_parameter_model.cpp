#include <RcppDist.h>
#include <string>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// [[Rcpp::export]]
Rcpp::NumericVector rprocess_time_varying_parameter_model(
    arma::mat cov_mat, arma::vec m_betas,
    double cr, double mys, double fr, double ms, double gdp, double sigma_u
) {
    arma::colvec z = {1, cr, mys, fr, ms, gdp};
    arma::rowvec betas = rmvnorm(1, m_betas, cov_mat).row(0);

    double e_lpd = (betas * z).eval()(0, 0) + R::rnorm(0, exp(sigma_u));

    return Rcpp::NumericVector::create(
        Rcpp::Named("beta_0", betas(0)),
        Rcpp::Named("beta_1", betas(1)),
        Rcpp::Named("beta_2", betas(2)),
        Rcpp::Named("beta_3", betas(3)),
        Rcpp::Named("beta_4", betas(4)),
        Rcpp::Named("beta_5", betas(5)),
        Rcpp::Named("e_lpd", e_lpd)
    );
}
