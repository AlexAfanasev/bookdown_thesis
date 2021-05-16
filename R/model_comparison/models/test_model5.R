# SETUP
source(here::here("R", "pd_pomp.R"))
y <- read.csv(here::here("data", "final_dataset.csv"))
source(here::here("R", "covars.R"))
Rcpp::sourceCpp(here::here("src", "rprocess_time_varying_parameter_model.cpp"))

# QUESTION: WHICH MODEL TO CHOSE???
# DO NOT USE THIS
# MODEL 5: MODEL WITH Time varying COVARIATES (stationary)
model_6 <- pomp::pomp(
    data = y[, c(1, 2)], times = "time", t0 = 0,
    rinit = function(beta_0_m, beta_1_m, beta_2_m, ...) {
        return(c(e_lpd = 3.5, beta_0 = beta_0_m, beta_1 = beta_1_m,
                 beta_2 = beta_2_m))
    },
    rprocess = pomp::discrete_time(
        function(
            e_lpd, beta_0, beta_1, beta_2,
            sigma_0, sigma_1, sigma_2,
            cr, mys, fr, ms, gdp, sigma_u, ...
        ){
            std_mat <- diag(
                exp(c(sigma_0, sigma_1, sigma_2))
            )^2

            m_betas <- c(beta_0, beta_1, beta_2)

            return(rprocess_time_varying_parameter_model(
                cov_mat = std_mat, cr = cr, mys = mys,
                fr = fr, ms = ms, gdp = gdp, m_betas = m_betas,
                sigma_u = sigma_u
            ))
        },
        delta.t = 1
    ),
    dmeasure = rw_latent_lpd_dmeasure,
    statenames = c("e_lpd", "beta_0", "beta_1", "beta_2"),
    paramnames = c("sigma_e", "e_lpd_0", "beta_0_m", "beta_1_m", "sigma_u",
                   "beta_2_m", "sigma_0", "sigma_1", "sigma_2"),
    covar = pomp::covariate_table(covars, times = "time"),
    covarnames = colnames(covars[, -1])
)

theta <- c(
    sigma_u = log(0.05),
    sigma_e = log(0.025),
    beta_0_m = 2.2,
    beta_1_m = 1.7,
    beta_2_m = 0.2,
    sigma_0 = log(0.01),
    sigma_1 = log(0.01),
    sigma_2 = log(0.01)
)
res <- pomp::pmcmc(
    model_6, Nmcmc = 50, Np = 250,
    proposal = pomp::mvn.diag.rw(
        c(sigma_e = 0.015, sigma_u = 0.015,
          beta_0_m = 0.015, beta_1_m = 0.015, beta_2_m = 0.015,
          sigma_0 = 0.015, sigma_1 = 0.015, sigma_2 = 0.015)
    ),
    params = theta,
    verbose = TRUE
)

pf <- pomp::pfilter(model_6, Np = 500, params = res@params, filter.mean = TRUE)
