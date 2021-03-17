# SETUP
source(here::here("R", "pd_pomp.R"))
y <- read.csv(here::here("data", "final_dataset.csv"))

# QUESTION: WHICH MODEL TO CHOSE???
# BASIC rw model:
rw_latent_lpd_pomp <- pomp::pomp(
    data = y[, c(1, 2)], times = "time", t0 = 0,
    rinit = function(e_lpd_0, ...) {
        return(c(e_lpd = e_lpd_0))
    },
    rprocess = pomp::discrete_time(
        pomp::Csnippet(
            "e_lpd = e_lpd + rnorm(0, exp(sigma_u));"
        ),
        delta.t = 1
    ),
    dmeasure = rw_latent_lpd_dmeasure,
    statenames = "e_lpd",
    paramnames = c("e_lpd_0", "sigma_u", "sigma_e"),
    covar = pomp::covariate_table(rbind(0, y[, -c(2, 6, 7, 8, 9)]),
                                  times = "time"),
    covarnames = colnames(rbind(0, y[, -c(1, 2, 6, 7, 8, 9)]))
)

theta <- c(e_lpd_0 = 3.5, sigma_u = log(0.05), sigma_e = log(0.05))
res <- pomp::pmcmc(
    rw_latent_lpd_pomp,
    Np = 350,
    Nmcmc = 5000,
    params = theta,
    proposal = pomp::mvn.diag.rw(
        c(e_lpd_0 = 0.02, sigma_u = 0.02, sigma_e = 0.02)
    )
)
