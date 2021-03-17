# QUESTION: WHICH MODEL TO CHOSE???
system.time(
    result <- pomp::bake(
        "model_study_result.rds",
        {
            models <- rep(
                c("basic_model.R", "model1.R", "model2.R", "model3.R"), 2
            )

            doRNG::registerDoRNG(1598260027L)
            cl <- parallel::makeCluster(4)
            doParallel::registerDoParallel(cl)
            out <- foreach::`%dopar%`(
                foreach::foreach(i = 1:8),
                {
                    source(here::here("R", models[i]))
                    res
                }
            )
            parallel::stopCluster(cl)
            out
        }
    )
)
