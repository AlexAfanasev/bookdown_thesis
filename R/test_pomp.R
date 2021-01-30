dm <- function(e_lpd, sigma_e, ..., log) {
    values <- list(...)
    A_mat <- matrix(
        c(
            values$param_1,
            values$param_2,
            values$param_3,
            values$param_4,
            values$param_5,
            values$param_6,
            values$param_7,
            values$param_8,
            values$param_9,
            values$param_10,
            values$param_11,
            values$param_12,
            values$param_13,
            values$param_14,
            values$param_15,
            values$param_16,
            values$param_17,
            values$param_18,
            values$param_19,
            values$param_20
        ),
        nrow = 4,
        byrow = TRUE
    )
    A_mat <- rbind(A_mat, c(rep(0, ncol(A_mat) - 1), 1))
    rho_t <- 1 / (1 + exp(-e_lpd))
    k_t <- -log(rho_t) - ((1 - rho_t) * log((1 / rho_t) - 1))
    m_t <- (
        (k_t / (1 - rho_t))
        + t(c(0, 1, -1, 0, 0))
        %*% A_mat
        %*% solve(diag(5) - rho_t * A_mat)
        %*% c(values$lpd, values$d_e, values$r_e, values$inflation, 1)
    )
    return(dnorm(values$lpd, m_t, sigma_e, log = log))
}


dm2 <- function(e_lpd, sigma_e, ..., log) {
    values <- list(...)

}
