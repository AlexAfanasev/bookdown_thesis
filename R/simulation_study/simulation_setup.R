# Simulation functions
generate_simulation_data <- function(n, params) {
    # simulate d_t
    d_t <- arima.sim(
        model = list(
            order = c(1, 0, 1),
            ar = 0.85,
            ma = 0.15
        ),
        n = n,
        sd = sqrt(0.003)
    )
    d_t <- d_t - 1.0

    # simulate latent mean
    latent_lpd_sim <- vector("numeric", length = n)
    x_1 <- vector("numeric", length = n)
    x_2 <- vector("numeric", length = n)
    for (i in 1:n) {
        if (i == 1) {
            x_1[i] <- 0
            x_2[i] <- 0
            latent_lpd_sim[i] <- params["x_0"]
        } else {
            while (TRUE) {
                val_2 <- x_1[i - 1] + rnorm(1, 0, 0.05)
                val_3 <- x_2[i - 1] + rnorm(1, 0, 0.05)
                val_1 <- (
                    params["beta_0"]
                    + tanh(params["phi"]) * latent_lpd_sim[i - 1]
                    + params["beta_1"] * val_2
                    + params["beta_2"] * val_3
                    + rnorm(1, 0, exp(params["sigma_x"]))
                )
                if (val_1 > 2.775 & val_1 < 4.487) {
                    x_1[i] <- val_2
                    x_2[i] <- val_3
                    latent_lpd_sim[i] <- val_1
                    break
                }
            }
        }
    }

    # compute observed lpd
    observed_lpd <- vector("numeric", length = n)
    for (i in 1:n) {
        rho <- 1 / (1 + exp(-latent_lpd_sim[i]))
        k <- -log(rho) - ((1 - rho) * log((1 / rho) - 1))
        observed_lpd[i] <- ((k / (1 - rho)) + d_t[i]
                            + rnorm(1, 0, exp(params["sigma_y"])))
    }
    return(list(
        observed = cbind(observed_lpd, d_t, x_1, x_2),
        latent = latent_lpd_sim
    ))
}

generate_start_params <- function() {
    return(
        c(
            beta_0 = runif(1, 0.0, 1.0),
            beta_1 = runif(1, -0.5, 0.5),
            beta_2 = runif(1, -0.5, -0.5),
            sigma_x = log(runif(1, 0.05, 0.1)),
            sigma_y = log(runif(1, 0.02, 0.05)),
            phi = atanh(runif(1, 0.8, 0.99)),
            x_0 = runif(1, 3.0, 3.2)
        )
    )
}

# Setup simulation
N <- 76
N_sim <- 1000
true_params_simulation <- c(
    beta_0 = 0.525,
    beta_1 = 0.4,
    beta_2 = -0.2,
    sigma_x = log(0.05),
    sigma_y = log(0.02),
    phi = atanh(0.85),
    x_0 = 3.5
)
