source(here::here("R", "simulation_study", "simulation.R"))


sim_4 <- simulation_result[[4]]

res_smc <- matrix(nrow = N_sim, ncol = length(true_params_simulation))
mse_smc <- vector("numeric", length = N_sim)
res_pmmh_noninformative <- matrix(
    nrow = N_sim, ncol = length(true_params_simulation)
)
mse_pmmh_noninformative <- vector("numeric", length = N_sim)
res_pmmh_true_informative <- matrix(
    nrow = N_sim, ncol = length(true_params_simulation)
)
mse_pmmh_true_informative <- vector("numeric", length = N_sim)
res_pmmh_false_informative <- matrix(
    nrow = N_sim, ncol = length(true_params_simulation)
)
mse_pmmh_false_informative <- vector("numeric", length = N_sim)

for (i in 1:N_sim) {
    res_smc[i, ] <- simulation_result[[i]]$smc_result$par
    mse_smc[i] <- sqrt(mean((
        simulation_result[[i]]$smc_result$pf -
            simulation_result[[i]]$data_latent
    )[2:N]^2))
    for (j in 1:length(true_params_simulation)) {
        res_pmmh_noninformative[i, j] <- mean(
            simulation_result[[i]]$pmmh_noninformative$traces[25000:50000, j+2]
        )
        mse_pmmh_noninformative[i] <- sqrt(mean((
            simulation_result[[i]]$pmmh_noninformative$state_mean -
                simulation_result[[i]]$data_latent
        )[2:N]^2))
        res_pmmh_true_informative[i, j] <- mean(
            simulation_result[[i]]$pmmh_true_informative$traces[25000:50000, j+2]
        )
        mse_pmmh_true_informative[i] <- sqrt(mean((
            simulation_result[[i]]$pmmh_true_informative$state_mean -
                simulation_result[[i]]$data_latent
        )[2:N]^2))
        res_pmmh_false_informative[i, j] <- mean(
            simulation_result[[i]]$pmmh_false_informative$traces[25000:50000, j+2]
        )
        mse_pmmh_false_informative[i] <- sqrt(mean((
            simulation_result[[i]]$pmmh_false_informative$state_mean -
                simulation_result[[i]]$data_latent
        )[2:N]^2))
    }
}

rm(simulation_result)
save(file = here::here("data", "results", "simulation_results.RData"))
