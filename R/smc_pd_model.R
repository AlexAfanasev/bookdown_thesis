# function for computing mean of lpd_t using campbell-shiller PV model
p_v <- function(rho, a_mat, last_obs) {
    return(
        tryCatch(
            {
                k <- -log(rho) - ((1 - rho) * log((1 / rho) - 1))
                (
                    (k / (1 - rho))
                    + t(c(0, 1, -1, 0, 0))
                    %*% a_mat
                    %*% solve(diag(5) - rho * a_mat)
                    %*% c(last_obs, 1)
                )
            },
            error = function(cond){
                return(Inf)
            }
        )
    )
}

# function for computing probability of observing drawn state
prob_state <- function(y, state, start_point, params = params) {
    previous_state <- tail(state[1, , ], 1)

    rho_t <- 1 / (1 + exp(-previous_state))

    var_model <- vars::VAR(
        y[(nrow(y) - start_point + 1):nrow(y), 1:4], p = 1
    )
    coef_matrix <- do.call("rbind", lapply(var_model$varresult, coef))
    coef_matrix <- rbind(coef_matrix, c(rep(0, ncol(coef_matrix) - 1), 1))

    last_obs <- as.vector(tail(y[, 1:4], 1))

    means <- sapply(rho_t, p_v, a_mat = coef_matrix, last_obs = last_obs)
    probs <- dnorm(last_obs[1], means, exp(params[["sigma_e"]]))

    return(probs)
}

# function for sampling state with covariates
sample_state_with_covariates <- function(
    previous_state, params, time_point, y
) {
    dim_previous_state <- dim(previous_state)
    if (time_point == 0) {
        return(rep(params[["eta_0"]], dim_previous_state[3]))
    } else if (dim_previous_state[2] == 1) {
        return(rep(params[["eta_1"]], dim_previous_state[3]))
    } else {
        eta_t_1 <- previous_state[, dim_previous_state[2], ]
        eta_t_2 <- ifelse(
            dim_previous_state[2] <= 1 , 0,
            previous_state[, dim_previous_state[2] - 1, ]
        )
        y_t_1 <- as.vector(y[time_point - 1, ])

        new_sample <- (
            params[["beta_0"]]
            + params[["beta_eta_1"]] * eta_t_1
            + params[["beta_eta_2"]] * eta_t_2
            + params[["beta_cr"]] * y_t_1[5]
            + params[["beta_my"]] * y_t_1[6]
            + params[["beta_fr"]] * y_t_1[7]
            + params[["beta_ms"]] * y_t_1[8]
            + rnorm(dim_previous_state[3], 0, exp(params[["sigma_u"]]))
        )

        return(new_sample)
    }
}

# function for sampling state without covariates
sample_state_without_covariates <- function(
    previous_state, params, time_point, y
) {
    dim_previous_state <- dim(previous_state)
    if (time_point == 0) {
        return(rep(params[["eta_0"]], dim_previous_state[3]))
    } else {
        new_sample <- (
            previous_state[1, dim_previous_state[2], ]
            + rnorm(dim_previous_state[3], 0, exp(params[["sigma_u"]]))
        )

        return(new_sample)
    }
}
