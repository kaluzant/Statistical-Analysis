# ---------------------------------------------------------------------------
# How smoking status is affected by other variables
# Original report: How smoking status is affected by other variables.pdf
#
# Data: data/smoking.csv, or pass a CSV path as the first argument.
# Required columns: smoker, ban, age, age2, education, aahisp, female
# Packages: ggplot2
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else file.path("data", "smoking.csv")
if (!file.exists(data_file)) {
  stop(
    "Data file not found: ", data_file,
    "\nPlace smoking.csv in data/ or pass its path as the first argument."
  )
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
}

smoke <- read.csv(data_file, stringsAsFactors = FALSE)
required_columns <- c(
  "smoker", "ban", "age", "age2", "education", "aahisp", "female"
)
missing_columns <- setdiff(required_columns, names(smoke))
if (length(missing_columns) > 0) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "))
}

complete_rows <- complete.cases(smoke[, required_columns])
if (any(!complete_rows)) {
  cat("Excluding", sum(!complete_rows), "rows with missing analysis values.\n")
}
smoke <- smoke[complete_rows, required_columns, drop = FALSE]
rownames(smoke) <- NULL
if (nrow(smoke) == 0) {
  stop("No complete observations remain in the smoking data.")
}

numeric_columns <- c("smoker", "ban", "age", "age2", "aahisp", "female")
not_numeric <- numeric_columns[
  !vapply(smoke[numeric_columns], is.numeric, logical(1))
]
if (length(not_numeric) > 0) {
  stop("These columns must be numeric: ", paste(not_numeric, collapse = ", "))
}
if (!all(vapply(smoke[numeric_columns], function(x) all(is.finite(x)), logical(1)))) {
  stop("Numeric analysis columns must contain only finite values.")
}

binary_columns <- c("smoker", "ban", "aahisp", "female")
invalid_binary <- binary_columns[
  !vapply(smoke[binary_columns], function(x) all(x %in% c(0, 1)), logical(1))
]
if (length(invalid_binary) > 0) {
  stop("These columns must be coded as 0/1: ", paste(invalid_binary, collapse = ", "))
}

education_raw <- trimws(as.character(smoke$education))
education_codes <- suppressWarnings(as.numeric(education_raw))
education_levels <- c(
  "HS Dropout", "HS Graduate", "Some College", "College Graduate", "Masters +"
)
if (!anyNA(education_codes) && all(education_codes %in% seq_along(education_levels))) {
  smoke$education <- factor(
    education_codes,
    levels = seq_along(education_levels),
    labels = education_levels
  )
} else {
  education <- rep(NA_character_, length(education_raw))
  education[grepl("drop", education_raw, ignore.case = TRUE)] <- "HS Dropout"
  education[grepl("some", education_raw, ignore.case = TRUE)] <- "Some College"
  education[grepl("master|post.?grad", education_raw, ignore.case = TRUE)] <- "Masters +"
  education[
    grepl("college|bachelor", education_raw, ignore.case = TRUE) &
      is.na(education)
  ] <- "College Graduate"
  education[
    grepl("grad", education_raw, ignore.case = TRUE) & is.na(education)
  ] <- "HS Graduate"
  if (anyNA(education)) {
    stop(
      "Education must use codes 1-5 or labels identifying the five report categories: ",
      paste(education_levels, collapse = ", "),
      "."
    )
  }
  smoke$education <- factor(education, levels = education_levels)
}

if (any(!is.finite(smoke$age)) || any(!is.finite(smoke$age2))) {
  stop("Age and age2 must contain only finite values.")
}

cat("Observations:", nrow(smoke), "\n")
print(summary(smoke))

# Exploratory summaries and plots
print(
  ggplot2::ggplot(smoke, ggplot2::aes(x = age)) +
    ggplot2::geom_histogram(bins = 30, fill = "skyblue", colour = "white") +
    ggplot2::labs(title = "Age Distribution", x = "Age (years)", y = "Count") +
    ggplot2::theme_minimal()
)
print(
  ggplot2::ggplot(smoke, ggplot2::aes(x = education)) +
    ggplot2::geom_bar(fill = "skyblue") +
    ggplot2::labs(title = "Education Levels", x = "Education", y = "Count") +
    ggplot2::theme_minimal()
)

plot_smoking_proportion <- function(data, group, title, x_label) {
  print(
    ggplot2::ggplot(
      data,
      ggplot2::aes(x = factor(.data[[group]]), fill = factor(smoker))
    ) +
      ggplot2::geom_bar(position = "fill") +
      ggplot2::labs(
        title = title,
        x = x_label,
        y = "Proportion",
        fill = "Smoking status"
      ) +
      ggplot2::theme_minimal()
  )
}

plot_smoking_proportion(
  smoke, "education", "Smoking Status by Education", "Education"
)
plot_smoking_proportion(
  smoke, "ban", "Smoking Status by Work Area Ban", "Work area ban"
)
plot_smoking_proportion(
  smoke, "aahisp", "Smoking Status by Race/Ethnicity", "African-American or Hispanic"
)
plot_smoking_proportion(smoke, "female", "Smoking Status by Gender", "Female")

# Full model from the report and its written-up selected model.
full_model <- glm(
  smoker ~ aahisp * (ban + age + age2 + education + female),
  family = binomial(link = "logit"),
  data = smoke
)
print(summary(full_model))

selected_formula <- smoker ~
  aahisp + ban + age + age2 + education + female +
  aahisp:age2 + aahisp:education + aahisp:female
selected_model <- glm(
  selected_formula,
  family = binomial(link = "logit"),
  data = smoke
)

cat("\nLikelihood-ratio test: selected model versus full interaction model\n")
print(anova(selected_model, full_model, test = "Chisq"))
cat("\nSelected model summary\n")
print(summary(selected_model))
cat("AIC of full model:", AIC(full_model), "\n")
cat("AIC of selected model:", AIC(selected_model), "\n")

# Odds ratios and profile-likelihood confidence intervals
coefficient_summary <- summary(selected_model)$coefficients
profile_intervals <- suppressMessages(confint(selected_model))
odds_ratio_table <- data.frame(
  Term = rownames(coefficient_summary),
  Odds_Ratio = exp(coefficient_summary[, "Estimate"]),
  CI_Lower = exp(profile_intervals[, 1]),
  CI_Upper = exp(profile_intervals[, 2]),
  P_Value = coefficient_summary[, "Pr(>|z|)"],
  row.names = NULL
)
print(odds_ratio_table)

# Model diagnostics and an influential-observation sensitivity fit
par(mfrow = c(2, 2))
plot(
  fitted(selected_model), rstandard(selected_model, type = "deviance"),
  xlab = "Fitted probability", ylab = "Standardized deviance residual",
  main = "Residuals vs fitted"
)
abline(h = 0, col = "red")
qqnorm(rstandard(selected_model, type = "deviance"))
qqline(rstandard(selected_model, type = "deviance"), col = "red")
plot(
  hatvalues(selected_model), ylab = "Leverage",
  main = "Leverage"
)
cook_distance <- cooks.distance(selected_model)
cook_threshold <- 4 / nobs(selected_model)
plot(cook_distance, ylab = "Cook's distance", main = "Influential observations")
abline(h = cook_threshold, col = "red", lty = 2)
par(mfrow = c(1, 1))

influential_rows <- which(cook_distance > cook_threshold)
cat("Cook's distance threshold:", cook_threshold, "\n")
cat("Influential observations:", length(influential_rows), "\n")
print(influential_rows)
if (length(influential_rows) > 0) {
  sensitivity_data <- smoke[-influential_rows, , drop = FALSE]
  sensitivity_model <- glm(
    selected_formula,
    family = binomial(link = "logit"),
    data = sensitivity_data
  )
  cat("\nSensitivity model after excluding flagged observations\n")
  print(summary(sensitivity_model))
}

# Compare categorical education with the report's ordinal numeric alternative.
smoke$education_numeric <- as.numeric(smoke$education)
continuous_education_model <- glm(
  smoker ~
    aahisp + ban + age + age2 + education_numeric + female +
    aahisp:age2 + aahisp:education_numeric + aahisp:female,
  family = binomial(link = "logit"),
  data = smoke
)
cat("\nEducation coding comparison\n")
cat("AIC (categorical education):", AIC(selected_model), "\n")
cat("AIC (numeric education):", AIC(continuous_education_model), "\n")
print(anova(continuous_education_model, selected_model, test = "Chisq"))
