# Rcpp function for rw_laten_lpd_model
Rcpp::cppFunction(
    '
    double rw_laten_lpd_dmeasure_cpp(
        double e_lpd, double sigma_e,
        Rcpp::List other_par,
        bool as_log
    ) {
        arma::mat A(5, 5);
        arma::rowvec h = {0.0, 1.0, -1.0, 0.0, 0.0};
        arma::colvec y_t = {
            other_par["lpd"], other_par["d_e"], other_par["r_e"],
            other_par["inflation"], 1
        };
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 5; j++) {
                A(i, j) = other_par[
                    "param_" + std::to_string((i*5) + (j + 1))
                ];
            }
        }

        for (int i = 0; i < 5; i++) {
            A(4, i) = i == 4 ? 1 : 0;
        }

        double rho_t = 1 / (1 + exp(-e_lpd));
        double k_t = -log(rho_t) - ((1.0 - rho_t) * log((1.0 / rho_t) - 1.0));
        double m_t = (
            (k_t / (1.0 - rho_t))
            + (
                h * A
                * arma::inv(arma::eye(5, 5) - rho_t * A)
                * y_t
            ).eval()(0, 0)
        );
        return(R::dnorm4(y_t(0), m_t, exp(sigma_e), as_log));
    }
    ',
    depends = c("RcppArmadillo, RcppDist")
)

# pomp model dmeasure for latent lpd model
rw_latent_lpd_dmeasure <- function(e_lpd, sigma_e, ..., log) {
    values <- list(...)
    return(rw_laten_lpd_dmeasure_cpp(e_lpd, sigma_e, values, log))
}
