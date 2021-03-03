#include <RcppDist.h>
#include <string>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// [[Rcpp::export]]
double dm_cpp2(
        double e_lpd,
        double g,
        double g_p,
        double gamma_0,
        double gamma_1,
        double sigma_e,
        double sigma_d,
        double rho_dp,
        double lpd,
        double d_e,
        bool as_log
) {
    arma::mat corr_mat = {
        {1.0, rho_dp},
        {rho_dp, 1}
    };
    arma::mat std_matrix = {
        {sigma_e, 0.0},
        {0.0, sigma_d}
    };
    arma::mat cov_matrix = std_matrix * corr_mat * std_matrix;

    double rho_t = 1 / (1 + exp(-e_lpd));
    double k_t = -log(rho_t) - ((1.0 - rho_t) * log((1.0 / rho_t) - 1.0));

    double m_1 = (
        (k_t / (1.0 - rho_t))
        + (gamma_0 / (1.0 - rho_t))
        + ((g - gamma_0) / (1.0 - gamma_1 * rho_t))
    );
    double m_2 = g_p;
    arma::vec mu = {m_1, m_2};
    arma::mat vals = {{lpd, d_e}};
    return(dmvnorm(vals, mu, cov_matrix, as_log).eval()(0, 0));
}
