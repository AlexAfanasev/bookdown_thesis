# SETUP
source(here::here("R", "pd_pomp.R"))
y <- read.csv(here::here("data", "final_dataset.csv"))
source(here::here("R", "covars.R"))

# QUESTION: WHICH MODEL TO CHOSE???
# MODEL 2: MODEL WITH AR & COVARIATES
model_2 <- pomp::pomp(
    data = y[, c(1, 2)], times = "time", t0 = 0,
    rinit = function(e_lpd_0, ...) {
        return(c(e_lpd = e_lpd_0))
    },
    rprocess = pomp::discrete_time(
        pomp::Csnippet(
            "
            e_lpd = (
              beta_0
              + tanh(phi)*e_lpd
              + beta_1*mys
              + beta_2*ms
              + rnorm(0, exp(sigma_u))
            );
            "
        ),
        delta.t = 1
    ),
    dmeasure = rw_latent_lpd_dmeasure,
    statenames = c("e_lpd"),
    paramnames = c("sigma_u", "sigma_e", "e_lpd_0", "beta_0", "phi",
                   "beta_1", "beta_2"),
    covar = pomp::covariate_table(covars, times = "time"),
    covarnames = colnames(covars[, -1])
)
rm(covars, y)

theta <- c(
    e_lpd_0 = 3.5, sigma_e = log(0.02), sigma_u = log(0.05), phi = atanh(0.95),
    beta_0 = 0.12, beta_1 = 0.08, beta_2 = 0.02
)
res <- pomp::pmcmc(
    model_2, Nmcmc = 500, Np = 2000,
    proposal = pomp::mvn.diag.rw(
        c(e_lpd_0 = 0.001, sigma_e = 0.015, sigma_u = 0.015, phi = 0.015,
          beta_0 = 0.001, beta_1 = 0.001, beta_2 = 0.001,
          beta_3 = 0.001, beta_4 = 0.001, beta_5 = 0.001)
    ),
    params = theta,
    dprior = function(sigma_u,
                      sigma_e,
                      phi,
                      e_lpd_0,
                      beta_0,
                      beta_1,
                      beta_2,
                      beta_3,
                      beta_4,
                      beta_5,
                      ...,
                      log) {
        p_sigma_u <- exp(sigma_e)
        p_sigma_e <- exp(sigma_u)
        p_phi <- (1 - tanh(phi) ^ 2)
        lik <- p_sigma_u * p_sigma_e * p_phi
        return(ifelse(log, log(lik), lik))
    }
)

a <- replicate(
    10, pomp::logLik(pomp::pfilter(model_2, Np = 2000, params = theta))
)
