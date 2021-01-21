# function for computing probability of observing drawn state
prob_state_simulation_linear <- function(y, state, start_point, params = params) {
    means <- (
        state[1, dim(state)[2], ]
        + tanh(params[["theta"]]) * state[2, dim(state)[2], ]
    )

    probs <- dnorm(
        as.vector(tail(y, 1)), means, exp(params[["sigma_e"]])
    )
    return(probs)
}

# function for sampling states for simulation1
sample_state_simulation_linear <- function(previous_state, params, time_point, y) {
    dim_previous_state <- dim(previous_state)
    if (time_point == 0) {
        sampled_data <- matrix(
            c(
                rep(params[["eta_0"]], dim_previous_state[3]),
                rep(0, dim_previous_state[3])
            ), nrow = 2, ncol = dim_previous_state[3], byrow = TRUE
        )
        return(sampled_data)
    } else {
        previous_sample <- previous_state[, dim_previous_state[2], ]

        sampled_data <- matrix(nrow = 2, ncol = dim_previous_state[3])
        for (i in 1:dim_previous_state[3]) {
            while (TRUE) {
                val <- (
                    previous_sample[1, i]
                    + rnorm(1, 0, exp(params[["sigma_u"]]))
                )
                if (val > 2.9 & val < 4.35) {
                    sampled_data[1, i] <- val
                    if (dim_previous_state[2] == 1) {
                        sampled_data[2, i] <- 0.0
                    } else {
                        sampled_data[2, i] <- (
                            y[time_point - 1, 1]
                            - previous_sample[1, i]
                            - tanh(params[["theta"]]) * previous_sample[2, i]
                        )
                    }
                    break
                }
            }
        }

        return(sampled_data)
    }
}
