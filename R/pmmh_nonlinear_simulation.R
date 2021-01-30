# function for computing likelihood of values from prior distribution
likli_prior_nl_sim_ui <- function(values, log = FALSE) {
    probs <- list(
        beta_1 = 0.0,
        beta_2 = 0.0,
        sigma_e = 0.0,
        sigma_u = 0.0
    )
    # beta_1
    probs$beta_1 <- dnorm(values$beta_1, 0.0, 10.0, log = log)
    # beta_2
    probs$beta_2 <- dnorm(values$beta_2, 0.0, 10.0, log = log)
    # log(sigma_e)
    probs$sigma_e <-
        dnorm(values$sigma_e, log(0.1), sd = 10.0, log = log)
    # log(sigma_u)
    probs$sigma_u <-
        dnorm(values$sigma_u, log(0.1), sd = 10.0, log = log)
    return(unlist(probs))
}

# function for computing new proposal given previous proposals
proposal_nl_sim_ui <- function(previous_proposal, prop_cov = NULL) {
    while (TRUE) {
        new_proposal <- list(
            beta_1 = 0.0,
            beta_2 = 0.0,
            sigma_e = 0.0,
            sigma_u = 0.0
        )
        prop_mean <- c(
            previous_proposal$beta_1,
            previous_proposal$beta_2,
            previous_proposal$sigma_e,
            previous_proposal$sigma_u
        )
        if (is.null(prop_cov)) {
            prop_cov <- diag(c(0.025, 0.025, 0.05, 0.05))^2
        } else {
            prop_cov <- (
                0.95 * ((2.38^2) / 4) * prop_cov
                + 0.05 * diag(c(0.025, 0.025, 0.05, 0.05))^2
            )
        }

        prop <- mvtnorm::rmvnorm(1, prop_mean, prop_cov)

        # beta_1
        new_proposal$beta_1 <- prop[1]
        # beta_2
        new_proposal$beta_2 <- prop[2]
        # log(sigma_e)
        new_proposal$sigma_e <- prop[3]
        # log(sigma_u)
        new_proposal$sigma_u <- prop[4]
        if (all(likli_prior_nl_sim_ui(new_proposal) > 0.0)) {
            return(new_proposal)
        }
    }
}

# function for computing likelihood of new proposal given previous proposals
likli_proposal_nl_sim_ui <-
    function(new_proposal, previous_proposal, log = FALSE, prop_cov = NULL) {
        prop_mean <- c(
            previous_proposal$beta_1,
            previous_proposal$beta_2,
            previous_proposal$sigma_e,
            previous_proposal$sigma_u
        )
        if (is.null(prop_cov)) {
            prop_cov <- diag(c(0.025, 0.025, 0.05, 0.05))^2
        } else {
            prop_cov <- (
                0.95 * ((2.38^2) / 4) * prop_cov
                + 0.05 * diag(c(0.025, 0.025, 0.05, 0.05))^2
            )
        }
        new_obs <- c(
            new_proposal$beta_1,
            new_proposal$beta_2,
            new_proposal$sigma_e,
            new_proposal$sigma_u
        )
        return(as.vector(mvtnorm::dmvnorm(new_obs, prop_mean, prop_cov,
                                          log = log)))
    }
