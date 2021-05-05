# QUESTION: WHICH MODEL TO CHOSE???
system.time(
    result <- pomp::bake(
        here::here("data", "results", "model_study_result_ar_error_model.rds"),
        {
            models <- rep("ar_error_model.R", 3)

            library(parallel)
            library(doParallel)
            library(foreach)
            library(doRNG)
            cl <- makeCluster(3)
            registerDoParallel(cl)
            registerDoRNG(42)
            out <- foreach(i = 1:3) %dopar% {
                source(here::here("R", "model_comparison", "models", models[i]))
                res
            }
            stopCluster(cl)
            out
        }
    )
)
