# set true parameters
N <- 50
true_params_simulation <- c(
    beta_1 = 0.4,
    beta_2 = -0.2,
    sigma_y = log(0.015),
    sigma_x = log(0.05),
    phi = 0.85
)

# generate data
generated_data <- generate_simulation_data(N, true_params_simulation)

# transform to state x and process y
data <- data.frame(cbind(1:N, generated_data$observed))
colnames(data) <- c("time", "y", "d", "x1", "x2")
y <- data[, c("time", "y")]
covar <- data[, c("time", "d", "x1", "x2")]
covar_2 <- covar
covar_2[1:(N - 1), 3:4] <- covar_2[2:N, 3:4]
x <- generated_data$latent
head(y)
tail(y)

# visualize data
plot(y[, 2], type = "l")
lines(x, col = "red")
legend(
    "topleft",
    legend = c("obs", "latent"),
    col = c("black", "red"),
    lty = c(1, 1)
)

# TODO: SHOULD I USE COVAR OR COVAR_2??? --> use covar_2
# TODO: Justify process for d_t ???? --> assume it is given as an covariate (simplify simulation)

# create pomp object
pomp_simulation <- pomp::pomp(
    # init data
    data = y, times = "time", t0 = 1,
    # init process / state
    rinit = function(mu, ...){
        c(x = mu)
    },
    rprocess = pomp::discrete_time(
        pomp::Csnippet(
            "
            x = (
                mu + tanh(phi) * (x - mu) + beta_1 * x1 + beta_2 * x2
                + rnorm(0, exp(sigma_x))
            );
            "
        ), delta.t = 1
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
    paramnames = c("beta_1", "beta_2", "mu", "sigma_x",
                   "sigma_y", "phi"),
    # covar
    covar = pomp::covariate_table(covar_2, times = "time"),
    covarnames = c("d", "x1", "x2")
)

# test optim parameter estimation
theta <- c(
    beta_1 = 0.0,
    beta_2 = 0.0,
    mu = 3.0,
    sigma_x = log(0.1),
    sigma_y = log(0.05),
    phi = atanh(0.8)
)
# can use optim only for finding MLE --> not for CI
system.time(
    optim_result <- replicate(
        10,
        {
            op_res <- optim(
                par = theta,
                fn = function(par){
                    # print(par)
                    pomp::logLik(pomp::pfilter(
                        pomp_simulation, params = par, Np = 500
                    ))
                },
                control = list(
                    "fnscale" = -1
                )
            )
            return(c(op_res$par, value = op_res$value))
        }
    )
)

k <- which.max(optim_result["value", ])
(est_params <- optim_result[1:length(theta), k])
true_params_simulation

# compare pfilter result
pf1 <- pomp::pfilter(
    pomp_simulation,
    params = c(
        beta_1 = as.vector(true_params_simulation["beta_1"]),
        beta_2 = as.vector(true_params_simulation["beta_2"]),
        sigma_y = as.vector(true_params_simulation["sigma_y"]),
        sigma_x = as.vector(true_params_simulation["sigma_x"]),
        phi = atanh(as.vector(true_params_simulation["phi"])),
        mu = 3.5
    ),
    Np = 500,
    filter.mean = TRUE
)
pomp::logLik(pf1)

pf2 <- pomp::pfilter(
    pomp_simulation,
    params = est_params,
    Np = 500,
    filter.mean = TRUE
)
pomp::logLik(pf2)

# pmcmc test
system.time(
    a <- pomp::pmcmc(
        pomp_simulation,
        Nmcmc = 6000,
        Np = 500,
        proposal = pomp::mvn.diag.rw(
            c(mu = 0.03, sigma_x = 0.03, sigma_y = 0.03, beta_1 = 0.03,
              beta_2 = 0.03, phi = 0.03)
        ),
        params = theta,
        verbose = FALSE
    )
)

# pomp::plot(a)
