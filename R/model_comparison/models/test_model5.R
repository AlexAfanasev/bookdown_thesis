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
    rinit = function(
        e_lpd_0, beta_0_m, beta_1_m, beta_2_m, beta_3_m, beta_4_m,
        beta_5_m, ...
    ) {
        return(c(e_lpd = e_lpd_0, beta_0 = beta_0_m, beta_1 = beta_1_m,
                 beta_2 = beta_2_m, beta_3 = beta_3_m, beta_4 = beta_4_m,
                 beta_5 = beta_5_m))
    },
    rprocess = pomp::discrete_time(
        function(
            e_lpd, beta_0_m, beta_1_m, beta_2_m, beta_3_m, beta_4_m, beta_5_m,
            sigma_0, sigma_1, sigma_2, sigma_3, sigma_4, sigma_5,
            cr, mys, fr, ms, gdp, sigma_u, ...
        ){
            std_mat <- diag(
                exp(c(sigma_0, sigma_1, sigma_2, sigma_3, sigma_4, sigma_5))
            )^2

            m_betas <- c(beta_0_m, beta_1_m, beta_2_m, beta_3_m, beta_4_m,
                         beta_5_m)

            return(rprocess_time_varying_parameter_model(
                cov_mat = std_mat, cr = cr, mys = mys,
                fr = fr, ms = ms, gdp = gdp, m_betas = m_betas,
                sigma_u = sigma_u
            ))
        },
        delta.t = 1
    ),
    dmeasure = rw_latent_lpd_dmeasure,
    statenames = c("e_lpd", "beta_0", "beta_1", "beta_2", "beta_3", "beta_4",
                   "beta_5"),
    paramnames = c("sigma_e", "e_lpd_0", "beta_0_m", "beta_1_m", "sigma_u",
                   "beta_2_m", "beta_3_m", "beta_4_m", "beta_5_m", "sigma_0",
                   "sigma_1", "sigma_2", "sigma_3", "sigma_4", "sigma_5"),
    covar = pomp::covariate_table(covars, times = "time"),
    covarnames = colnames(covars[, -1])
)

i_theta <- as.vector(result[[2]]@params)
theta <- c(
    e_lpd_0 = i_theta[9],
    sigma_u = i_theta[7],
    sigma_e = i_theta[8],
    beta_0_m = i_theta[1],
    beta_1_m = i_theta[2],
    beta_2_m = i_theta[3],
    beta_3_m = i_theta[4],
    beta_4_m = i_theta[5],
    beta_5_m = i_theta[6],
    sigma_0 = log(0.01),
    sigma_1 = log(0.01),
    sigma_2 = log(0.01),
    sigma_3 = log(0.01),
    sigma_4 = log(0.01),
    sigma_5 = log(0.01)
)
res <- pomp::pmcmc(
    model_6, Nmcmc = 50, Np = 250,
    proposal = pomp::mvn.diag.rw(
        c(e_lpd_0 = 0.03, sigma_e = 0.03, sigma_u = 0.03,
          beta_0_m = 0.03, beta_1_m = 0.03, beta_2_m = 0.03,
          beta_3_m = 0.03, beta_4_m = 0.03, beta_5_m = 0.03,
          sigma_0 = 0.03, sigma_1 = 0.03, sigma_2 = 0.03, sigma_3 = 0.03,
          sigma_4 = 0.03, sigma_5 = 0.03)
    ),
    params = theta,
    verbose = TRUE
)

pf <- pomp::pfilter(model_6, Np = 500, params = res@params, filter.mean = TRUE)
