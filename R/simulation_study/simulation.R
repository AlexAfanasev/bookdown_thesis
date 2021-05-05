source(here::here("R", "simulation_study", "simulation_setup.R"))

nonlinear_simulation <- function(N_obs, true_params) {
    data_simulation <- generate_simulation_data(N_obs, true_params)
    theta <- generate_start_params()

    # transform to state x and process y
    data <- data.frame(cbind(1:N, data_simulation$observed))
    colnames(data) <- c("time", "y", "d", "x1", "x2")
    y <- data[, c("time", "y")]
    covar <- data[, c("time", "d", "x1", "x2")]
    covar_2 <- covar
    covar_2[1:(N - 1), 3:4] <- covar_2[2:N, 3:4]

    # create pomp object
    pomp_simulation <- pomp::pomp(
        # init data
        data = y,
        times = "time",
        t0 = 1,
        # init process / state
        rinit = function(x_0, phi, ...) {
            c(x = x_0)
        },
        rprocess = pomp::discrete_time(
            pomp::Csnippet(
                "
            x = (
                beta_0 + tanh(phi) * x + beta_1 * x1 + beta_2 * x2
                + rnorm(0, exp(sigma_x))
            );
            "
            ),
            delta.t = 1
        ),
        # init meas / obs
        dmeasure = pomp::Csnippet(
            "
        double rho_t = 1.0 / (1.0 + exp(-x));
        double k_t = -log(rho_t) - ((1.0 - rho_t) * log((1.0 / rho_t) - 1.0));
        double m = (k_t / (1 - rho_t)) + d;
        lik = dnorm(y, m, exp(sigma_y), give_log);
        "
        ),
        # state names
        statenames = "x",
        # obs names
        obsnames = "y",
        # param names
        paramnames = c(
            "beta_1",
            "beta_2",
            "beta_0",
            "sigma_x",
            "sigma_y",
            "phi",
            "x_0"
        ),
        # covar
        covar = pomp::covariate_table(covar_2, times = "time"),
        covarnames = c("d", "x1", "x2")
    )

    # can use optim only for finding MLE --> not for CI
    optim_result <- replicate(10,
                              {
                                  op_res <- optim(
                                      par = theta,
                                      fn = function(par) {
                                          # print(par)
                                          pomp::logLik(pomp::pfilter(
                                              pomp_simulation,
                                              params = par,
                                              Np = 1500
                                          ))
                                      },
                                      control = list("fnscale" = -1)
                                  )
                                  return(c(op_res$par, value = op_res$value))
                              })

    k <- which.max(optim_result["value",])
    est_params <- optim_result[1:length(theta), k]
    pf <- pomp::pfilter(
        pomp_simulation,
        params = est_params,
        Np = 1500,
        filter.mean = TRUE
    )

    # pmcmc test - flat prior
    pmmh_noninformative <- pomp::pmcmc(
        pomp_simulation,
        Nmcmc = 20000,
        Np = 500,
        proposal = pomp::mvn.diag.rw(
            c(
                beta_0 = 0.025,
                sigma_x = 0.025,
                sigma_y = 0.025,
                beta_1 = 0.025,
                beta_2 = 0.025,
                phi = 0.025,
                x_0 = 0.025
            )
        ),
        params = theta,
        verbose = FALSE,
        dprior = function(beta_0,
                          sigma_x,
                          sigma_y,
                          beta_1,
                          beta_2,
                          phi,
                          x_0,
                          ...,
                          log) {
            p_x_0 <- 1
            p_beta_0 <- 1
            p_sigma_x <- exp(sigma_x)
            p_sigma_y <- exp(sigma_y)
            p_beta_1 <- 1
            p_beta_2 <- 1
            p_phi <- (1 - tanh(phi) ^ 2)
            lik <- (
                p_x_0 * p_beta_0 * p_sigma_x * p_sigma_y * p_beta_1 *
                    p_beta_2 * p_phi
            )
            return(ifelse(log, log(lik), lik))
        }
    )

    # pmcmc test - true informative prior
    pmmh_true_informative <- pomp::pmcmc(
        pomp_simulation,
        Nmcmc = 20000,
        Np = 500,
        proposal = pomp::mvn.diag.rw(
            c(
                beta_0 = 0.025,
                sigma_x = 0.025,
                sigma_y = 0.025,
                beta_1 = 0.025,
                beta_2 = 0.025,
                phi = 0.025,
                x_0 = 0.025
            )
        ),
        params = theta,
        verbose = FALSE,
        dprior = function(beta_0,
                          sigma_x,
                          sigma_y,
                          beta_1,
                          beta_2,
                          phi,
                          x_0,
                          ...,
                          log) {
            p_beta_0 <- dnorm(beta_0, 0.4, 0.1)
            p_beta_1 <- dnorm(beta_1, 0.5, 0.2)
            p_beta_2 <- dnorm(beta_2, -0.3, 0.2)
            p_sigma_x <- (
                dgamma(exp(sigma_x), shape = 12, scale = 0.005) * exp(sigma_x)
            )
            p_sigma_y <- (
                dgamma(exp(sigma_y), shape = 6, scale = 0.005) * exp(sigma_y)
            )
            p_phi <- (
                truncnorm::dtruncnorm(tanh(phi), -1, 1, 0.8, 0.1) *
                    (1 - tanh(phi) ^ 2)
            )
            p_x_0 <- dnorm(x_0, 3.5, 0.1)
            lik <- (
                p_beta_0 * p_sigma_x * p_sigma_y * p_beta_1 * p_beta_2 *
                    p_phi * p_x_0
            )
            return(ifelse(log, log(lik), lik))
        }
    )

    # pmcmc test - misspecified informative prior
    pmmh_false_informative <- pomp::pmcmc(
        pomp_simulation,
        Nmcmc = 20000,
        Np = 500,
        proposal = pomp::mvn.diag.rw(
            c(
                beta_0 = 0.025,
                sigma_x = 0.025,
                sigma_y = 0.025,
                beta_1 = 0.025,
                beta_2 = 0.025,
                phi = 0.025,
                x_0 = 0.025
            )
        ),
        params = theta,
        verbose = FALSE,
        dprior = function(beta_0,
                          sigma_x,
                          sigma_y,
                          beta_1,
                          beta_2,
                          phi,
                          x_0,
                          ...,
                          log) {
            p_beta_0 <- dnorm(beta_0, 0.8, 0.1)
            p_beta_1 <- dnorm(beta_1, 0.0, 0.1)
            p_beta_2 <- dnorm(beta_2, 0.0, 0.1)
            p_sigma_x <- (
                dgamma(exp(sigma_x), shape = 10, scale = 0.01) * exp(sigma_x)
            )
            p_sigma_y <- (
                dgamma(exp(sigma_y), shape = 10, scale = 0.005) * exp(sigma_y)
            )
            p_phi <- (
                truncnorm::dtruncnorm(tanh(phi), -1, 1, 0.3, 0.1) *
                    (1 - tanh(phi) ^ 2)
            )
            p_x_0 <- dnorm(x_0, 3.0, 0.1)
            lik <- (
                p_beta_0 * p_sigma_x * p_sigma_y * p_beta_1 * p_beta_2 *
                    p_phi * p_x_0
            )
            return(ifelse(log, log(lik), lik))
        }
    )

    return(
        list(
            data_observed = data_simulation$observed,
            data_latent = data_simulation$latent,
            theta = theta,
            smc_result = list(
                par = est_params,
                lik = pf@loglik,
                pf = pf@filter.mean[1, ],
                optim_result = optim_result
            ),
            pmmh_noninformative = list(
                state_mean = colMeans(
                    pmmh_noninformative@filter.traj[1, 10000:20000, ]
                ),
                state_lower = apply(
                    pmmh_noninformative@filter.traj[1, 10000:20000, ],
                    MARGIN = 2,
                    quantile, probs = 0.025
                ),
                state_upper = apply(
                    pmmh_noninformative@filter.traj[1, 10000:20000, ],
                    MARGIN = 2,
                    quantile, probs = 0.975
                ),
                traces = pmmh_noninformative@traces
            ),
            pmmh_true_informative = list(
                state_mean = colMeans(
                    pmmh_true_informative@filter.traj[1, 10000:20000, ]
                ),
                state_lower = apply(
                    pmmh_true_informative@filter.traj[1, 10000:20000, ],
                    MARGIN = 2,
                    quantile, probs = 0.025
                ),
                state_upper = apply(
                    pmmh_true_informative@filter.traj[1, 10000:20000, ],
                    MARGIN = 2,
                    quantile, probs = 0.975
                ),
                traces = pmmh_true_informative@traces
            ),
            pmmh_false_informative = list(
                state_mean = colMeans(
                    pmmh_false_informative@filter.traj[1, 10000:20000, ]
                ),
                state_lower = apply(
                    pmmh_false_informative@filter.traj[1, 10000:20000, ],
                    MARGIN = 2,
                    quantile, probs = 0.025
                ),
                state_upper = apply(
                    pmmh_false_informative@filter.traj[1, 10000:20000, ],
                    MARGIN = 2,
                    quantile, probs = 0.975
                ),
                traces = pmmh_false_informative@traces
            )
        )
    )
}

# bake this stuff
system.time(
    simulation_result <- pomp::bake(
        here::here("data", "results", "simulation_study_result.rds"),
        {
            library(parallel)
            library(doParallel)
            library(foreach)
            library(doRNG)
            cl <- makeCluster(3)
            registerDoParallel(cl)
            registerDoRNG(1234)
            results <- foreach(i = 1:N_sim) %dopar% {
                nonlinear_simulation(
                    N, true_params_simulation
                )
            }
            stopCluster(cl)
            results
        }
    )
)

