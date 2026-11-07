# ---------------------------------------------------------------------------
# Doctor visit counts: Poisson regression and bootstrap inference
# Original report:
# A statistical analysis of doctor visit counts using Poisson Regression and Bootstrap methods-1.pdf
#
# Data: data/dvisus.csv, or pass a CSV path as the first argument.
# Required columns: docvis, age, lincome, female, educ, married, white,
#                   hispanic, private, chronic
# Packages: boot, lmtest, sandwich
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else file.path("data", "dvisus.csv")
if (!file.exists(data_file)) {
  stop(
    "Data file not found: ", data_file,
    "\nPlace dvisus.csv in data/ or pass its path as the first argument."
  )
}

required_packages <- c("boot", "lmtest", "sandwich")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install required package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

visits <- read.csv(data_file, stringsAsFactors = FALSE)
required_columns <- c(
  "docvis", "age", "lincome", "female", "educ",
  "married", "white", "hispanic", "private", "chronic"
)
missing_columns <- setdiff(required_columns, names(visits))
if (length(missing_columns) > 0) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "))
}
if (nrow(visits) == 0) {
  stop("The doctor-visits data file contains no observations.")
}
if (!all(vapply(visits[required_columns], is.numeric, logical(1)))) {
  stop("All analysis columns must be numeric.")
}
if (any(!vapply(visits[required_columns], function(x) all(is.finite(x)), logical(1)))) {
  stop("Analysis columns must contain only finite values; remove or resolve missing values.")
}
binary_columns <- c("female", "married", "white", "hispanic", "private", "chronic")
invalid_binary <- binary_columns[
  !vapply(visits[binary_columns], function(x) all(x %in% c(0, 1)), logical(1))
]
if (length(invalid_binary) > 0) {
  stop("These columns must be coded as 0/1: ", paste(invalid_binary, collapse = ", "))
}
if (any(visits$docvis < 0)) {
  stop("docvis must be a non-negative count.")
}
if (any(visits$docvis != floor(visits$docvis))) {
  stop("docvis must contain integer counts.")
}

cat("Summary statistics:\n")
print(summary(visits))
cat(
  "Mean doctor visits:", mean(visits$docvis),
  "\nVariance of doctor visits:", var(visits$docvis), "\n"
)

par(mfrow = c(2, 2))
hist(visits$docvis, breaks = 30, col = "lightblue",
     main = "Doctor Visits", xlab = "Number of visits")
hist(visits$age, breaks = 30, col = "lightgreen",
     main = "Age", xlab = "Age (decades)")
hist(visits$lincome, breaks = 30, col = "lightcoral",
     main = "Log Income", xlab = "Log income (/1000)")
hist(visits$educ, breaks = 30, col = "wheat",
     main = "Education", xlab = "Years of schooling")
par(mfrow = c(1, 1))

poisson_formula <- docvis ~
  age + lincome + female + educ + married + white + hispanic + private + chronic
poisson_model <- glm(
  poisson_formula,
  family = poisson(link = "log"),
  data = visits
)
print(summary(poisson_model))

robust_vcov <- sandwich::vcovHC(poisson_model, type = "HC0")
cat("\nCoefficient tests with HC0 sandwich standard errors:\n")
print(lmtest::coeftest(poisson_model, vcov. = robust_vcov))

standard_errors <- sqrt(diag(vcov(poisson_model)))
robust_standard_errors <- sqrt(diag(robust_vcov))
standard_error_comparison <- data.frame(
  Estimate = coef(poisson_model),
  Standard_SE = standard_errors,
  Robust_SE = robust_standard_errors
)

# Row bootstrap for coefficient standard errors.
bootstrap_coefficient <- function(data, indices) {
  sample_data <- data[indices, , drop = FALSE]
  fit <- glm(poisson_formula, family = poisson(link = "log"), data = sample_data)
  coef(fit)
}

bootstrap_replicates <- 1000L
set.seed(42069)
coefficient_bootstrap <- boot::boot(
  data = visits,
  statistic = bootstrap_coefficient,
  R = bootstrap_replicates
)
coefficient_draws <- coefficient_bootstrap$t
colnames(coefficient_draws) <- names(coef(poisson_model))
bootstrap_standard_errors <- apply(coefficient_draws, 2, sd)
standard_error_comparison$Bootstrap_SE <- bootstrap_standard_errors
cat("\nConventional, robust, and bootstrap standard errors:\n")
print(round(standard_error_comparison, 5))

# Joint Wald test of income, marital status, and Hispanic ethnicity.
tested_parameters <- c("lincome", "married", "hispanic")
beta_tested <- coef(poisson_model)[tested_parameters]
standard_vcov_tested <- vcov(poisson_model)[tested_parameters, tested_parameters]
robust_vcov_tested <- robust_vcov[tested_parameters, tested_parameters]
bootstrap_vcov <- cov(coefficient_draws)
bootstrap_vcov_tested <- bootstrap_vcov[tested_parameters, tested_parameters]

wald_test <- function(beta, covariance) {
  statistic <- drop(crossprod(beta, solve(covariance, beta)))
  c(
    statistic = statistic,
    p_value = pchisq(statistic, df = length(beta), lower.tail = FALSE)
  )
}

wald_standard <- wald_test(beta_tested, standard_vcov_tested)
wald_robust <- wald_test(beta_tested, robust_vcov_tested)
wald_bootstrap <- wald_test(beta_tested, bootstrap_vcov_tested)
wald_table <- data.frame(
  Method = c("Standard", "Sandwich (HC0)", "Bootstrap covariance"),
  Statistic = c(wald_standard["statistic"], wald_robust["statistic"],
                wald_bootstrap["statistic"]),
  P_Value = c(wald_standard["p_value"], wald_robust["p_value"],
              wald_bootstrap["p_value"])
)
cat("\nJoint Wald test: lincome = married = hispanic = 0\n")
print(wald_table)

# Empirical bootstrap p-value following the report's resampling procedure.
# This resamples the observed rows; it does not impose the null hypothesis.
bootstrap_wald_statistic <- function(data, indices) {
  sample_data <- data[indices, , drop = FALSE]
  fit <- glm(poisson_formula, family = poisson(link = "log"), data = sample_data)
  beta <- coef(fit)[tested_parameters]
  covariance <- sandwich::vcovHC(fit, type = "HC0")[
    tested_parameters, tested_parameters
  ]
  drop(crossprod(beta, solve(covariance, beta)))
}

set.seed(69420)
wald_bootstrap_distribution <- boot::boot(
  data = visits,
  statistic = bootstrap_wald_statistic,
  R = bootstrap_replicates
)
observed_robust_wald <- unname(wald_robust["statistic"])
empirical_bootstrap_p_value <- mean(
  wald_bootstrap_distribution$t[, 1] >= observed_robust_wald
)
cat("Empirical bootstrap p-value:", empirical_bootstrap_p_value, "\n")

# Studentized bootstrap intervals for the report's restricted model.
restricted_formula <- docvis ~ age + female + educ + white + private + chronic
restricted_model <- glm(
  restricted_formula,
  family = poisson(link = "log"),
  data = visits
)
original_estimates <- coef(restricted_model)
original_standard_errors <- sqrt(diag(vcov(restricted_model)))

studentized_statistic <- function(data, indices) {
  sample_data <- data[indices, , drop = FALSE]
  fit <- glm(restricted_formula, family = poisson(link = "log"), data = sample_data)
  estimates <- coef(fit)
  standard_errors <- sqrt(diag(vcov(fit)))
  (estimates - original_estimates) / standard_errors
}

set.seed(21)
studentized_bootstrap <- boot::boot(
  data = visits,
  statistic = studentized_statistic,
  R = bootstrap_replicates
)
studentized_draws <- studentized_bootstrap$t
alpha <- 0.05
coefficient_intervals <- lapply(seq_along(original_estimates), function(i) {
  t_star <- studentized_draws[, i]
  t_star <- t_star[is.finite(t_star)]
  quantiles <- quantile(t_star, probs = c(alpha / 2, 1 - alpha / 2))
  c(
    lower = original_estimates[i] - quantiles[2] * original_standard_errors[i],
    upper = original_estimates[i] - quantiles[1] * original_standard_errors[i]
  )
})
coefficient_intervals <- do.call(rbind, coefficient_intervals)
rownames(coefficient_intervals) <- names(original_estimates)
interval_table <- data.frame(
  Coefficient = names(original_estimates),
  Estimate = unname(original_estimates),
  Multiplicative_Effect = exp(unname(original_estimates)),
  CI_Lower = exp(coefficient_intervals[, "lower"]),
  CI_Upper = exp(coefficient_intervals[, "upper"]),
  row.names = NULL
)
cat("\n95% studentized bootstrap intervals (multiplicative scale):\n")
print(round(interval_table, 4))
