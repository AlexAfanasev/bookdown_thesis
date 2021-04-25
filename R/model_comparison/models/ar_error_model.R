# SETUP
source(here::here("R", "pd_pomp.R"))
y <- read.csv(here::here("data", "final_dataset.csv"))
source(here::here("R", "covars.R"))

# QUESTION: WHICH MODEL TO CHOSE???
# MODEL 4: MODEL WITH AR Errors & COVARIATES
model_2 <- pomp::pomp(
    data = y[, c(1, 2)], times = "time", t0 = 0,
    rinit = function(e_lpd_0, ...) {
        return(c(e_lpd = e_lpd_0, u = 0))
    },
    rprocess = pomp::discrete_time(
        pomp::Csnippet(
            "
            u = tanh(phi) * u + rnorm(0, exp(sigma_u));
            e_lpd = (
              beta_0
              + beta_1*cr
              + beta_2*mys
              + beta_3*ms
              + beta_4*fr
              + beta_5*gdp
              + u
            );
            "
        ),
        delta.t = 1
    ),
    dmeasure = rw_latent_lpd_dmeasure,
    statenames = c("e_lpd", "u"),
    paramnames = c("sigma_u", "sigma_e", "e_lpd_0", "beta_0", "phi",
                   "beta_1", "beta_2", "beta_3", "beta_4", "beta_5"),
    covar = pomp::covariate_table(covars, times = "time"),
    covarnames = colnames(covars[, -1])
)
rm(covars, y)

theta <- c(
    e_lpd_0 = 3.5, sigma_e = log(0.05), sigma_u = log(0.05), phi = atanh(0.95),
    beta_0 = 2.75, beta_1 = 0.18, beta_2 = 1.86, beta_3 = 0.75, beta_4 = 0.69,
    beta_5 = -0.15
)
res <- pomp::pmcmc(
    model_2, Nmcmc = 500, Np = 1000,
    proposal = pomp::mvn.diag.rw(
        c(sigma_e = 0.02, sigma_u = 0.02, phi = 0.02,
          beta_0 = 0.02, beta_1 = 0.02, beta_2 = 0.02,
          beta_3 = 0.02, beta_4 = 0.02, u_0 = 0.02, beta_5 = 0.02)
    ),
    params = theta,
    dprior = function(sigma_u,
                      sigma_e,
                      phi,
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
