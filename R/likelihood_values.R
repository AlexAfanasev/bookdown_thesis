source(here::here("R", "pd_pomp.R"))
y <- read.csv(here::here("data", "final_dataset.csv"))
source(here::here("R", "covars.R"))

r <- readRDS(
    here::here("data", "results", "model_study_result_basic_model.rds")
)
N <- ncol(r[[1]]@data)
rw_model <- lapply(r, function(x){x@traces[10000:50000, ]})
rw_model <- do.call(rbind, rw_model)

rw_model_pomp <- pomp::pomp(
    data = y[, c(1, 2)],
    times = "time",
    t0 = 0,
    rinit = function(e_lpd_0, ...) {
        return(c(e_lpd = e_lpd_0))
    },
    rprocess = pomp::discrete_time(
        pomp::Csnippet("e_lpd = e_lpd + rnorm(0, exp(sigma_u));"),
        delta.t = 1
    ),
    dmeasure = rw_latent_lpd_dmeasure,
    statenames = c("e_lpd"),
    paramnames = c("e_lpd_0", "sigma_u", "sigma_e"),
    covar = pomp::covariate_table(rbind(0, y[, -c(2, 6, 7, 8, 9)]),
                                  times = "time"),
    covarnames = colnames(rbind(0, y[, -c(1, 2, 6, 7, 8, 9)]))
)

theta <- c(e_lpd_0 = mean(rw_model[, 3]),
           sigma_u = log(mean(exp(rw_model[, 4]))),
           sigma_e = log(mean(exp(rw_model[, 5]))))

l_rw_model <- mean(
    replicate(10, pomp::logLik(pomp::pfilter(rw_model_pomp, Np = 1000, params = theta)))
)

model_1_pomp <- pomp::pomp(
    data = y[, c(1, 2)], times = "time", t0 = 0,
    rinit = function(e_lpd_0, ...) {
        return(c(e_lpd = e_lpd_0))
    },
    rprocess = pomp::discrete_time(
        pomp::Csnippet(
            "
            e_lpd = (
              beta_0
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
    paramnames = c("sigma_u", "sigma_e", "e_lpd_0", "beta_0",
                   "beta_1", "beta_2"),
    covar = pomp::covariate_table(covars, times = "time"),
    covarnames = colnames(covars[, -1])
)

model_1 <- lapply(
    readRDS(
        here::here("data", "results", "model_study_result_static_model.rds")
    ), function(x){x@traces[10000:50000, ]}
)
model_1 <- do.call(rbind, model_1)

theta <- c(
    beta_0 = mean(model_1[, 3]), beta_1 = mean(model_1[, 4]), beta_2 = mean(model_1[, 5]), sigma_u = log(mean(exp(model_1[, 6]))), sigma_e = log(mean(exp(model_1[, 7]))), e_lpd_0 = 3.5
)

l_cov_model <- mean(
    replicate(10, pomp::logLik(pomp::pfilter(model_1_pomp, Np = 1000, params = theta)))
)

saveRDS(l_rw_model, file = here::here("data", "results", "rw_lik.rds"))
saveRDS(l_cov_model, file = here::here("data", "results", "cov_lik.rds"))
