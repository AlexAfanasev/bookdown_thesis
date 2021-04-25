# QUESTION: WHICH MODEL TO CHOSE???
system.time(
    result <- pomp::bake(
        here::here("data", "results", "model_study_result_ardl_model.rds"),
        {
            models <- rep("ardl_model.R", 3)

            doRNG::registerDoRNG(1598260027L)
            cl <- parallel::makeCluster(3)
            doParallel::registerDoParallel(cl)
            out <- foreach::`%dopar%`(
                foreach::foreach(i = 1:3),
                {
                    source(here::here("R", "model_comparison", "models",
                                      models[i]))
                    res
                }
            )
            parallel::stopCluster(cl)
            out
        }
    )
)
