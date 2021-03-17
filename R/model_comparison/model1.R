# SETUP
source(here::here("R", "pd_pomp.R"))
y <- read.csv(here::here("data", "final_dataset.csv"))
source(here::here("R", "covars.R"))

# QUESTION: WHICH MODEL TO CHOSE???
# MODEL 2: only with covariates
model_1 <- pomp::pomp(
    data = y[, c(1, 2)], times = "time", t0 = 0,
    rinit = function(e_lpd_0, ...) {
        return(c(e_lpd = e_lpd_0))
    },
    rprocess = pomp::discrete_time(
        pomp::Csnippet(
            "
            e_lpd = (
              beta_0
              + beta_1*cr
              + beta_2*mys
              + beta_3*fr
              + beta_4*ms
              + beta_5*gdp
              + rnorm(0, exp(sigma_u))
            );
            "
        ),
        delta.t = 1
    ),
    dmeasure = rw_latent_lpd_dmeasure,
    statenames = c("e_lpd"),
    paramnames = c("sigma_u", "sigma_e", "e_lpd_0", "beta_0",
                   "beta_1", "beta_2", "beta_3", "beta_4", "beta_5"),
    covar = pomp::covariate_table(covars, times = "time"),
    covarnames = colnames(covars[, -1])
)

theta <- c(
    beta_0 = 2.8, beta_1 = 0.0, beta_2 = 0.0, beta_3 = 0.0,
    beta_4 = 0.0, beta_5 = 0.0,
    sigma_u = log(0.05), sigma_e = log(0.05), e_lpd_0 = 3.4
)
res <- pomp::pmcmc(
    model_1, Nmcmc = 5000, Np = 350,
    proposal = pomp::mvn.diag.rw(
        c(e_lpd_0 = 0.02, sigma_e = 0.02, sigma_u = 0.02,
          beta_0 = 0.02, beta_1 = 0.02, beta_2 = 0.02,
          beta_3 = 0.02, beta_4 = 0.02, beta_5 = 0.02)
    ),
    params = theta
)
