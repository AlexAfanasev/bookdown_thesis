# get covars
read.csv(here::here("data", "final_dataset.csv"))
covars <- rbind(0, y[, -2])
covars <- cbind(
    covars[, 1:9], covars[, 5:9], covars[, 10:29]
)
covars[1:(nrow(covars) - 1), colnames(covars)[5:9]] <- covars[
    2:nrow(covars), colnames(covars)[5:9]
]
colnames(covars)[10:14] <- c("l_cr", "l_mys", "l_fr", "l_ms", "l_gdp")
