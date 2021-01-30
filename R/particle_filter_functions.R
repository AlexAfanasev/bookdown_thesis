#' Sequential Monte Carlo // Particle Filter
#'
#' This function runs the bootstrap particle filter for given parameters.
#'
#' @param params A list of parameter values.
#' @param y A matrix of observations and covariates.
#' @param start_point A scalar defining the start point
#' @param prob_state Function for evaluating probability of drawing particular
#'  state: prob_state(y, state, params)
#' @param sample_state Function for sampling states:
#'  sample_state(previous_state, params, time_point, y)
#' @param n_particles A scalar defining the number of particles.
#' @param n_states Number of state variables.
#' @param more_information Boolean if more information should be printed.
#'
#' @return Returns a list containing the estimated negative log likelihood and
#' an array containing the estimated state and 95% CI.
smc <- function(params,
                y,
                start_point,
                prob_state,
                sample_state,
                n_particles,
                n_states,
                more_information = FALSE) {
    if (isTRUE(more_information)) {
        print("Parameters:")
        print(paste(paste(
            names(params), round(unlist(params), 4), sep = ": "
        ),
        collapse = " | "))
    }

    n_t <- nrow(y[start_point:nrow(y), , drop = FALSE])
    sampled_states <-
        array(0, dim = c(n_states, n_t + 1, n_particles))
    state_weights <- matrix(0, nrow = n_t, ncol = n_particles)
    estimated_states <- array(0, dim = c(n_states, n_t, 3))

    # Initialize drawing states
    sampled_states[, 1,] <- sample_state(
        previous_state = sampled_states[, 1, , drop = F],
        params = params,
        time_point = 0,
        y = y[1:(start_point - 1), , drop = FALSE]
    )

    prob_obs <- 0

    for (time_point in 1:n_t) {
        # Sample states
        sampled_states[, time_point + 1,] <- sample_state(
            previous_state = sampled_states[, 1:time_point, , drop = F],
            params = params,
            time_point = time_point + start_point - 1,
            y = y[1:(start_point + time_point - 1), , drop = FALSE]
        )

        # Compute Probability of drawn states
        state_weights[time_point,] <- prob_state(
            y = y[1:(start_point + time_point - 1), , drop = FALSE],
            state = sampled_states[, 1:(time_point + 1), , drop = F],
            start_point = start_point,
            params = params
        )

        if (isTRUE(any(is.na(state_weights[time_point,])))) {
            print(paste("NaN weight! Iteration:", time_point))
            return(list(prob = -Inf, state = NULL))
        } else if (all(state_weights[time_point,] == 0.0)) {
            print(paste("All weights are equal to 0 at iteration:", time_point))
            return(list(prob = -Inf, state = NULL))
        }

        # Resample particles (complete paths)
        sampled_indices <- sample(n_particles,
                                  n_particles,
                                  replace = TRUE,
                                  prob = state_weights[time_point,])
        sampled_states <-
            sampled_states[, , sampled_indices, drop = F]

        # Estimate state
        estimated_states[, time_point,] <- rowMeans(
            sampled_states[, time_point + 1, , drop = FALSE]
        )

        # Estimated state 95% CI
        estimated_states[, time_point, 2:3] <- apply(
            sampled_states[, time_point + 1, , drop = FALSE], 1,
            quantile, probs = c(0.025, 0.975)
        )

        # Compute Log Likelihood
        prob_obs <-
            prob_obs + log(mean(state_weights[time_point,]))
        state_weights <- state_weights[, sampled_indices, drop = F]
    }

    pmmh_particle <- sample(n_particles, 1)
    pmmh_path <- sampled_states[, -1, pmmh_particle]
    pmmh_weights <- state_weights[, pmmh_particle]

    return(
        list(
            prob = -1 * prob_obs,
            state = estimated_states,
            pmmh_prob = -1 * sum(log(pmmh_weights)),
            pmmh_path = pmmh_path
        )
    )
}


#' Particle Marginal Metropolis Hastings
#'
#' This function runs the pmmh algorithm and returns samples from posterior.
#'
#' @param params A list of parameter values.
#' @param y A matrix of observations and covariates.
#' @param start_point A scalar defining the start point
#' @param prob_state Function for evaluating probability of drawing particular
#'  state: prob_state(y, state, params)
#' @param sample_state Function for sampling states:
#'  sample_state(previous_state, params, time_point, y)
#' @param likelihood_prior Function that returns the likelihood of drawing the
#'   sampled values from the prior distributions:
#'   likelihood_prior(values, log = FALSE)
#' @param proposal Function that returns new proposal values given the previous
#'   proposals: proposal(previous_proposal)
#' @param likelihood_proposal Function that returns likelihood of drawing new
#'   proposal given the previous proposal:
#'   likelihood_proposal(new_proposal, previous_proposal, log = FALSE)
#' @param n_states Number of state variables.
#' @param n_met Scalar as number of samples using metropolis hasting algorithm
#' @param n_particles Scalar as number of particles for smc algorithm
#' @param burn_in Scalar as number of burn in samples
#' @param thining Scalar to use each nth sample
#' @param more_information Boolean if more information should be printed.
#'
#' @return Returns matrix containing the samples from the posterior. The burn in
#'   iterations are removed and the rest is thinned. The matrix contains
#'   number of params columns.
pmmh <- function(params,
                 y,
                 start_point,
                 prob_state,
                 sample_state,
                 likelihood_prior,
                 proposal,
                 likelihood_proposal,
                 n_states = 1,
                 n_met = 6500,
                 n_particles = 1000,
                 burn_in = 499,
                 thining = 1,
                 more_information = FALSE,
                 adaptive_n = NULL) {
    thetas <- matrix(nrow = n_met, ncol = length(params))
    paths <-
        matrix(nrow = n_met, ncol = nrow(y[start_point:nrow(y),]))
    theta_start <- params
    acceptance_rate <- rep(0, n_met)

    if (isTRUE(more_information)) {
        print("Iteration 0.")
    }

    # compute probability of y1:T using smc
    smc_result_start <- smc(
        params = theta_start,
        y = y,
        start_point = start_point,
        prob_state = prob_state,
        sample_state = sample_state,
        n_particles = n_particles,
        n_states = n_states,
        more_information = more_information
    )

    # breaks if start params result in -Inf probability
    if (smc_result_start$prob == -Inf) {
        stop("Start params resulted in -Inf probability.")
    }

    for (met_iteration in 1:n_met) {
        if (isTRUE(more_information)) {
            print(paste("Iteration:", met_iteration))
        }

        if (!is.null(adaptive_n)) {
            if (met_iteration %% adaptive_n == 0) {
                prop_cov <- var(var(thetas[1:(met_iteration - 1), ]))
            } else {
                prop_cov <- NULL
            }
        } else {
            prop_cov <- NULL
        }

        # sample new proposals and compute prob y1:T as long as not -Inf
        while (TRUE) {
            theta_proposal <- proposal(theta_start, prop_cov = prop_cov)

            smc_result_proposal <- smc(
                params = theta_proposal,
                y = y,
                start_point = start_point,
                prob_state = prob_state,
                sample_state = sample_state,
                n_particles = n_particles,
                n_states = n_states,
                more_information = more_information
            )

            if (smc_result_proposal$prob != -Inf) {
                break()
            } else if (isTRUE(more_information)) {
                print("New proposal params resulted in -Inf probability.")
            }
        }

        # computed acceptance probability
        prob <- (
            (smc_result_start$prob - smc_result_proposal$prob)
            + (sum(
                likelihood_prior(theta_proposal, log = TRUE)
            )
            - sum(
                likelihood_prior(theta_start, log = TRUE)
            ))
            + (
                sum(
                    likelihood_proposal(theta_start, theta_proposal,
                                        log = TRUE, prop_cov = prop_cov)
                )
                - sum(
                    likelihood_proposal(theta_proposal, theta_start,
                                        log = TRUE, prop_cov = prop_cov)
                )
            )
        )

        # overwrite initial proposal with new if accepted
        if (prob > log(runif(1))) {
            acceptance_rate[met_iteration] <- 1
            theta_start <- theta_proposal
            smc_result_start <- smc_result_proposal
        }

        # store sampled parameters
        thetas[met_iteration,] <- unlist(theta_start)
        # store specific path
        paths[met_iteration,] <- smc_result_start$pmmh_path
    }

    return(list(
        theta = thetas[seq(burn_in + 1, n_met, by = thining), ,
                       drop = FALSE],
        paths = paths[seq(burn_in + 1, n_met, by = thining), ,
                      drop = FALSE],
        acceptance_rate = acceptance_rate
    ))
}
