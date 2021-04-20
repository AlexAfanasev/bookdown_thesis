# prepares raw dataset for analysis
dataset <- readODS::read_ods(here::here("data", "final_dataset.ods"))
relevant_columns <- c(
    "Date",
    "Real Consumption (per Capita)",
    "Consumption (per Capita)",
    "USPOPTO",
    "USPOP24Y",
    "USPOP29Y",
    "USPOP44Y",
    "USPOP49Y",
    "S&P",
    "Dividend",
    "Earnings",
    "CPI",
    "rf",
    "fund_rate",
    "real_m1",
    "real_gdp"
)
# filter relevant columns
dataset <- dataset[, relevant_columns]

# rename columns
renamed_columns <- c(
    "date",
    "real_consumption",
    "consumption",
    "us_pop_total",
    "us_pop_20_24",
    "us_pop_25_29",
    "us_pop_40_44",
    "us_pop_45_49",
    "price",
    "dividend",
    "earnings",
    "cpi",
    "rf",
    "fund_rate",
    "real_m1",
    "real_gdp"
)
colnames(dataset) <- renamed_columns
row.names(dataset) <- dataset$date
dataset <- dataset[, colnames(dataset)[2:length(colnames(dataset))]]

# transform to xts object
dataset <- xts::xts(
    dataset, order.by = as.Date(row.names(dataset), "%Y-%m-%d")
)
dataset <- dataset["1949-12-31/"]
# dataset <- dataset["1949-12-31/2019-12-31"]

# linear interpolation for population data
pop_columns <- c(
    "us_pop_total",
    "us_pop_20_24",
    "us_pop_25_29",
    "us_pop_40_44",
    "us_pop_45_49"
)
dataset[
    -dataset[xts:::endof(dataset, "years"), which.i = TRUE],
    pop_columns
] <- NA
storage.mode(dataset) <- "numeric"
dataset[, pop_columns] <- zoo::na.approx(dataset[, pop_columns])
dataset <- dataset[complete.cases(dataset), ]
# dataset <- dataset[xts:::endof(dataset, "years")]

# consumption volatility
con <- dataset$real_consumption
quar_c_growth <- diff(con) / stats::lag(con, k = 1)
quar_c_growth <- quar_c_growth[!is.na(quar_c_growth)]
ar_mod <- arima(quar_c_growth, order = c(1, 0, 0))
con_vol <- log(
    zoo::rollapply(
        zoo::zoo(abs(residuals(ar_mod)), order.by = zoo::index(con)),
        width = 12, FUN = sum, align = "right"
    )
)

# middle aged young ratio
middle_aged <- dataset$us_pop_40_44 + dataset$us_pop_45_49
young <- dataset$us_pop_20_24 + dataset$us_pop_25_29
middle_aged_young_ratio <- middle_aged / young

# rolling VAR parameters + other covariates
r_e <- (
    log((dataset$price + dataset$dividend) /
            stats::lag(dataset$price, k = 1))
    - log(1 + (dataset$rf / 100))
)
d_e <- (
    log(dataset$dividend / stats::lag(dataset$dividend, k = 1))
    - log(1 + (dataset$rf / 100))
)
con_vol <- con_vol
my <- middle_aged_young_ratio
fr <- log(1 + (dataset$fund_rate/100))
ms <- log(dataset$real_m1/1000)
gdp <- log(dataset$real_gdp/10000)
inflation <- log(dataset$cpi / stats::lag(dataset$cpi, k = 1))
lpd <- log((dataset$price / dataset$dividend))

var_data <- cbind(lpd, d_e, r_e, inflation, con_vol, my, fr, ms, gdp)
var_data <- var_data[complete.cases(var_data),]["1963-09-30/"]

## dynamic var
starting_point <- 30
start_w <- 10
params <- matrix(nrow = nrow(var_data) - starting_point + 1, ncol = 20)
for (i in starting_point:nrow(var_data)) {
    residuals_model <- vector("numeric")
    for (w in start_w:(starting_point - 1)) {
        var_model <- vars::VAR(var_data[(i - w + 1):i, 1:4], p = 1,
                               type = "const")
        # test root
        e <- abs(eigen(vars::Acoef(var_model)[[1]])$values)
        if (any(e >= 1)) {
            residuals_model[w - start_w + 1] <- Inf
        } else {
            residuals_j <- tail(residuals(var_model), 5)[, c(2, 3)]
            residuals_model[w - start_w + 1] <- sqrt(
                mean((residuals_j[, 1] - residuals_j[, 2])^2)
            )
        }
    }
    best_model <- ifelse(all(residuals_model == Inf),
                         start_w,
                         which.min(residuals_model) + start_w - 1)
    model <- vars::VAR(var_data[(i - best_model + 1):i, 1:4], type = "const")
    e <- abs(eigen(vars::Acoef(model)[[1]])$values)
    params[i - starting_point + 1, ] <- as.vector(t(vars::Bcoef(model)))
}

var_data <- cbind(var_data[starting_point:nrow(var_data), ], params)

## static var
# starting_point <- 21
# params <- matrix(nrow = nrow(var_data) - starting_point + 1, ncol = 20)
# for (i in starting_point:nrow(var_data)) {
#   model <- vars::VAR(var_data[(i - starting_point + 1):i, 1:4], type = "const")
#   params[i - starting_point + 1, ] <- as.vector(t(vars::Bcoef(model)))
# }
# var_data <- cbind(var_data[starting_point:nrow(var_data), ], params)

## complete var
# starting_point <- 21
# params <- matrix(nrow = nrow(var_data) - starting_point + 1, ncol = 20)
# for (i in starting_point:nrow(var_data)) {
#   model <- vars::VAR(var_data[1:i, 1:4], type = "const")
#   params[i - starting_point + 1, ] <- as.vector(t(vars::Bcoef(model)))
# }
# var_data <- cbind(var_data[starting_point:nrow(var_data), ], params)

y <- var_data[complete.cases(var_data), ]
y <- cbind(1:nrow(y), y)
colnames(y) <- c("time", "lpd", "d_e", "r_e", "inflation", "cr", "mys", "fr",
                 "ms", "gdp", sapply(1:20, function(x){sprintf("param_%i", x)}))
dates <- zoo::index(y)
y <- as.data.frame(y, row.names = NULL)

head(y, 10)

write.csv(y, here::here("data", "final_dataset.csv"), row.names = FALSE)
