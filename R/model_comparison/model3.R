# SETUP
source(here::here("R", "pd_pomp.R"))
y <- read.csv(here::here("data", "final_dataset.csv"))
source(here::here("R", "covars.R"))

# QUESTION: WHICH MODEL TO CHOSE???
# MODEL 3: MODEL WITH AR AND LAGS & COVARIATES
model_3 <- pomp::pomp(
    data = y[, c(1, 2)], times = "time", t0 = 0,
    rinit = function(e_lpd_0, ...) {
        return(c(e_lpd = e_lpd_0))
    },
    rprocess = pomp::discrete_time(
        pomp::Csnippet(
            "
            e_lpd = (
                a0
                + tanh(phi)*e_lpd
                + (tanh(phi)-1)*(
                    -beta_0
                    - beta_1*l_cr
                    - beta_2*l_mys
                    - beta_3*l_fr
                    - beta_4*l_ms
                    - beta_5*l_gdp
                )
                + a1*(cr-l_cr)
                + a2*(mys-l_mys)
                + a3*(fr-l_fr)
                + a4*(ms-l_ms)
                + a5*(gdp-l_gdp)
                + rnorm(0, exp(sigma_u))
            );
            "
        ),
        delta.t = 1
    ),
    dmeasure = rw_latent_lpd_dmeasure,
    statenames = c("e_lpd"),
    paramnames = c("sigma_u", "sigma_e", "e_lpd_0", "beta_0", "phi",
                   "beta_1", "beta_2", "beta_3", "beta_4", "beta_5",
                   "a0", "a1", "a2", "a3", "a4", "a5"),
    covar = pomp::covariate_table(covars, times = "time"),
    covarnames = colnames(covars[, -1])
)

theta <- c(
    e_lpd_0 = 3.5, sigma_e = log(0.05), sigma_u = log(0.05), phi = atanh(0.95),
    beta_0 = 0.0, beta_1 = 0, beta_2 = 0, beta_3 = 0, beta_4 = 0, beta_5 = 0,
    a0 = 0.175, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0
)
res <- pomp::pmcmc(
    model_3, Nmcmc = 10000, Np = 1000,
    proposal = pomp::mvn.diag.rw(
        c(e_lpd_0 = 0.01, sigma_e = 0.01, sigma_u = 0.01, phi = 0.01,
          beta_0 = 0.01, beta_1 = 0.01, beta_2 = 0.01,
          beta_3 = 0.01, beta_4 = 0.01, beta_5 = 0.01,
          a0 = 0.01, a1 = 0.01, a2 = 0.01, a3 = 0.01, a4 = 0.01, a5 = 0.01)
    ),
    params = theta
)
