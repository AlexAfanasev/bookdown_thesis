# function for computing likelihood of values from prior distribution
likelihood_prior <- function(values, log = FALSE) {
    probs <- list(
        eta_0 = 0.0, beta_0 = 0.0, beta_eta = 0.0,
        beta_cr = 0.0, beta_my = 0.0, beta_fr = 0.0,  beta_ms = 0.0,
        sigma_e = 0.0, sigma_u = 0.0
    )
    # eta_0
    probs[["eta_0"]] <- dnorm(values[["eta_0"]], 3.5, 0.1, log = log)
    # beta_0
    probs[["beta_0"]] <- dnorm(values[["beta_0"]], 0.0, 0.1, log = log)
    # beta_eta
    probs[["beta_eta"]] <- dnorm(values[["beta_eta"]], 1.7, 0.1, log = log)
    # beta_cr
    probs[["beta_cr"]] <- dnorm(values[["beta_cr"]], 0.0, 0.05, log = log)
    # beta_my
    probs[["beta_my"]] <- dnorm(values[["beta_my"]], 0.15, 0.1, log = log)
    # beta_fr
    probs[["beta_fr"]] <- dnorm(values[["beta_fr"]], 0.0, 0.1, log = log)
    # beta_ms
    probs[["beta_ms"]] <- dnorm(values[["beta_ms"]], 0.0, 0.05, log = log)
    # log(sigma_e)
    probs[["sigma_e"]] <- dnorm(
        values[["sigma_e"]], log(0.01), sd = 0.25, log = log
    )
    # log(sigma_u)
    probs[["sigma_u"]] <- dnorm(
        values[["sigma_u"]], log(0.1), sd = 0.25, log = log
    )
    return(unlist(probs))
}

# function for computing new proposal given previous proposals
proposal <- function(previous_proposal) {
    while (TRUE) {
        new_proposal <- list(
            eta_0 = 0.0, beta_0 = 0.0, beta_eta = 0.0,
            beta_cr = 0.0, beta_my = 0.0, beta_fr = 0.0,  beta_ms = 0.0,
            sigma_e = log(0.1), sigma_u = log(0.1)
        )

        # eta_0
        new_proposal[["eta_0"]] <- rnorm(1, previous_proposal[["eta_0"]], 0.001)
        # beta_0
        new_proposal[["beta_0"]] <- rnorm(
            1, previous_proposal[["beta_0"]], 0.0001
        )
        # beta_eta_1
        new_proposal[["beta_eta"]] <- rnorm(
            1, previous_proposal[["beta_eta"]], 0.0001
        )
        # beta_cr
        new_proposal[["beta_cr"]] <- rnorm(
            1, previous_proposal[["beta_cr"]], 0.0001
        )
        # beta_my
        new_proposal[["beta_my"]] <- rnorm(
            1, previous_proposal[["beta_my"]], 0.0001
        )
        # beta_fr
        new_proposal[["beta_fr"]] <- rnorm(
            1, previous_proposal[["beta_fr"]], 0.0001
        )
        # beta_ms
        new_proposal[["beta_ms"]] <- rnorm(
            1, previous_proposal[["beta_ms"]], 0.0001
        )
        # log(sigma_e)
        new_proposal[["sigma_e"]] <- rnorm(
            1, previous_proposal[["sigma_e"]], 0.0001
        )
        # log(sigma_u)
        new_proposal["sigma_u"] <- rnorm(
            1, previous_proposal[["sigma_u"]], 0.0001
        )
        if (all(likelihood_prior(new_proposal) > 0.0)) {
            return(new_proposal)
        }
    }
}

# function for computing likelihood of new proposal given previous proposals
likelihood_proposal <- function(new_proposal, previous_proposal, log = FALSE) {
    return(
        c(
            dnorm(
                new_proposal[["eta_0"]],
                previous_proposal[["eta_0"]], 0.0001, log = log
            ),
            dnorm(
                new_proposal[["beta_0"]],
                previous_proposal[["beta_0"]], 0.0001, log = log
            ),
            dnorm(
                new_proposal[["beta_eta"]],
                previous_proposal[["beta_eta"]], 0.0001, log = log
            ),
            dnorm(
                new_proposal[["beta_cr"]],
                previous_proposal[["beta_cr"]], 0.0001, log = log
            ),
            dnorm(
                new_proposal[["beta_my"]],
                previous_proposal[["beta_my"]], 0.0001, log = log
            ),
            dnorm(
                new_proposal[["beta_fr"]],
                previous_proposal[["beta_fr"]], 0.0001, log = log
            ),
            dnorm(
                new_proposal[["beta_ms"]],
                previous_proposal[["beta_ms"]], 0.0001, log = log
            ),
            dnorm(
                new_proposal[["sigma_e"]],
                previous_proposal[["sigma_e"]], 0.0001, log = log
            ),
            dnorm(
                new_proposal[["sigma_u"]],
                previous_proposal[["sigma_u"]], 0.0001, log = log
            )
        )
    )
}
