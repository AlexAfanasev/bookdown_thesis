# function for computing probability of observing drawn state
prob_state_simulation_nonlinear <- function(
    y, state, start_point, params = params
) {
    last_obs <- y[nrow(y), ]
    latent_mean <- state[1, dim(state)[2], ]
    rho_t <- 1 / (1 + exp(-latent_mean))
    k_t <- -log(rho_t) - ((1 - rho_t) * log((1 / rho_t) - 1))
    means <- (k_t / (1 - rho_t)) + last_obs[2]

    probs <- dnorm(
        last_obs[1], means, exp(params[["sigma_e"]])
    )
    return(probs)
}

# function for sampling states for simulation1
sample_state_simulation_nonlinear <- function(
    previous_state, params, time_point, y
) {
    dim_previous_state <- dim(previous_state)
    if (time_point == 0) {
        return(rep(params[["eta_0"]], dim_previous_state[3]))
    } else {
        previous_sample <- previous_state[1, dim_previous_state[2], ]

        sampled_data <- vector("numeric", length = dim_previous_state[3])
        for (i in 1:dim_previous_state[3]) {
            while (TRUE) {
                val <- (
                    previous_sample[i]
                    + rnorm(1, 0, exp(params[["sigma_u"]]))
                )
                if (val > 2.9 & val < 4.35) {
                    sampled_data[i] <- val
                    break
                }
            }
        }

        return(sampled_data)
    }
}
