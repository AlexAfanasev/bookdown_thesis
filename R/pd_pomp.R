# Rcpp function for rw_laten_lpd_model
Rcpp::sourceCpp(here::here("src", "pd_pomp_dmeasure.cpp"))

# pomp model dmeasure for latent lpd model
rw_latent_lpd_dmeasure <- function(e_lpd, sigma_e, ..., log) {
    values <- list(...)
    return(rw_laten_lpd_dmeasure_cpp(e_lpd, sigma_e, values, log))
}
