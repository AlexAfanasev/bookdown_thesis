# SETUP
source(here::here("R", "pd_pomp.R"))
y <- read.csv(here::here("data", "final_dataset.csv"))
source(here::here("R", "covars.R"))
Rcpp::sourceCpp(here::here("src", "rprocess_time_varying_parameter_model.cpp"))

# QUESTION: WHICH MODEL TO CHOSE???
# DO NOT USE THIS
# MODEL 4: MODEL WITH Time varying COVARIATES
model_5 <- pomp::pomp(
    data = y[, c(1, 2)], times = "time", t0 = 0,
    rinit = function(
        e_lpd_0, beta_0_0, beta_1_0, beta_2_0, beta_3_0, beta_4_0,
        beta_5_0, ...
    ) {
        return(c(e_lpd = e_lpd_0, beta_0 = beta_0_0, beta_1 = beta_1_0,
                 beta_2 = beta_2_0, beta_3 = beta_3_0, beta_4 = beta_4_0,
                 beta_5 = beta_5_0))
    },
    rprocess = pomp::discrete_time(
        function(
            e_lpd, beta_0, beta_1, beta_2, beta_3, beta_4, beta_5,
            sigma_0, sigma_1, sigma_2, sigma_3, sigma_4, sigma_5,
            cor_01, cor_02, cor_03, cor_04, cor_05, cor_12, cor_13, cor_14,
            cor_15, cor_23, cor_24, cor_25, cor_34, cor_35, cor_45,
            cr, mys, fr, ms, gdp, ...
        ){
            cor_vector <- tanh(c(
                cor_01, cor_02, cor_03, cor_04, cor_05,
                cor_12, cor_13, cor_14, cor_15,
                cor_23, cor_24, cor_25,
                cor_34, cor_35,
                cor_45
            ))
            cor_mat <- matrix(1, nrow = 6, ncol = 6)
            cor_index <- matrix(1:36, nrow = 6, ncol = 6)
            cor_mat[cor_index[upper.tri(cor_index)]] <- cor_vector
            cor_mat[lower.tri(cor_mat)] <- t(cor_mat)[lower.tri(cor_mat)]
            std_mat <- diag(
                exp(c(sigma_0, sigma_1, sigma_2, sigma_3, sigma_4, sigma_5))
            )

            m_betas <- c(beta_0, beta_1, beta_2, beta_3, beta_4, beta_5)

            return(rprocess_time_varying_parameter_model(
                cor_mat = cor_mat, std_mat = std_mat, cr = cr, mys = mys,
                fr = fr, ms = ms, gdp = gdp, m_betas = m_betas
            ))
        },
        delta.t = 1
    ),
    dmeasure = rw_latent_lpd_dmeasure,
    statenames = c("e_lpd", "beta_0", "beta_1", "beta_2", "beta_3", "beta_4",
                   "beta_5"),
    paramnames = c("sigma_e", "e_lpd_0", "beta_0_0", "beta_1_0",
                   "beta_2_0", "beta_3_0", "beta_4_0", "beta_5_0", "sigma_0",
                   "sigma_1", "sigma_2", "sigma_3", "sigma_4", "sigma_5",
                   "cor_01", "cor_02", "cor_03", "cor_04", "cor_05",
                   "cor_12", "cor_13", "cor_14", "cor_15",
                   "cor_23", "cor_24", "cor_25",
                   "cor_34", "cor_35",
                   "cor_45"),
    covar = pomp::covariate_table(covars, times = "time"),
    covarnames = colnames(covars[, -1])
)

i_theta <- as.vector(result[[2]]@params)
theta <- c(
    e_lpd_0 = i_theta[9],
    sigma_e = i_theta[8],
    beta_0_0 = i_theta[1],
    beta_1_0 = i_theta[2],
    beta_2_0 = i_theta[3],
    beta_3_0 = i_theta[4],
    beta_4_0 = i_theta[5],
    beta_5_0 = i_theta[6],
    sigma_0 = log(0.01),
    sigma_1 = log(0.01),
    sigma_2 = log(0.01),
    sigma_3 = log(0.01),
    sigma_4 = log(0.01),
    sigma_5 = log(0.01),
    cor_01 = 0, cor_02 = 0, cor_03 = 0, cor_04 = 0, cor_05 = 0,
    cor_12 = 0, cor_13 = 0, cor_14 = 0, cor_15 = 0,
    cor_23 = 0, cor_24 = 0, cor_25 = 0,
    cor_34 = 0, cor_35 = 0,
    cor_45 = 0
)
res <- pomp::pmcmc(
    model_5, Nmcmc = 50, Np = 1000,
    proposal = pomp::mvn.diag.rw(
        c(e_lpd_0 = 0.03, sigma_e = 0.03,
          beta_0_0 = 0.03, beta_1_0 = 0.03, beta_2_0 = 0.03,
          beta_3_0 = 0.03, beta_4_0 = 0.03, beta_5_0 = 0.03,
          sigma_0 = 0.03, sigma_1 = 0.03, sigma_2 = 0.03, sigma_3 = 0.03,
          sigma_4 = 0.03, sigma_5 = 0.03,
          cor_01 = 0.03, cor_02 = 0.03, cor_03 = 0.03, cor_04 = 0.03,
          cor_05 = 0.03,
          cor_12 = 0.03, cor_13 = 0.03, cor_14 = 0.03, cor_15 = 0.03,
          cor_23 = 0.03, cor_24 = 0.03, cor_25 = 0.03,
          cor_34 = 0.03, cor_35 = 0.03,
          cor_45 = 0.03)
    ),
    params = res@params,
    verbose = TRUE
)
