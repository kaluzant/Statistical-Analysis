# ---------------------------------------------------------------------------
# Sydney Airport weekly maximum temperatures: 1973 vs 2003
# Original report: Is there any difference in temperatures in Sydney 1973 vs 2003.pdf
#
# Data: data/sydney_temperatures.csv, or pass a CSV path as the first argument.
# Required columns: Week, Temp1973, Temp2003
# Packages: base R (stats)
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) {
  args[1]
} else {
  file.path("data", "sydney_temperatures.csv")
}
if (!file.exists(data_file)) {
  stop(
    "Data file not found: ", data_file,
    "\nPlace sydney_temperatures.csv in data/ or pass its path as the first argument."
  )
}

sydney <- read.csv(data_file, stringsAsFactors = FALSE)
required_columns <- c("Week", "Temp1973", "Temp2003")
missing_columns <- setdiff(required_columns, names(sydney))
if (length(missing_columns) > 0) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "))
}
if (!all(vapply(sydney[required_columns], is.numeric, logical(1)))) {
  stop("Week, Temp1973, and Temp2003 must be numeric.")
}
sydney <- sydney[complete.cases(sydney[, required_columns]), required_columns]
if (nrow(sydney) == 0) {
  stop("No complete paired temperature observations remain.")
}
if (!all(vapply(sydney, function(x) all(is.finite(x)), logical(1)))) {
  stop("Week and temperature columns must contain only finite values.")
}
sydney <- sydney[order(sydney$Week), ]
if (anyDuplicated(sydney$Week)) {
  stop("Week must identify one paired observation per week.")
}

cat("Temperature summary:\n")
print(summary(sydney[, c("Temp1973", "Temp2003")]))

plot(
  sydney$Week, sydney$Temp1973,
  type = "l", col = "blue",
  ylim = range(sydney[, c("Temp1973", "Temp2003")]),
  xlab = "Week", ylab = "Maximum weekly temperature (C)",
  main = "Sydney Airport Maximum Weekly Temperatures"
)
lines(sydney$Week, sydney$Temp2003, col = "red")
legend(
  "topright", legend = c("1973", "2003"),
  col = c("blue", "red"), lty = 1
)

boxplot(
  sydney$Temp1973, sydney$Temp2003,
  names = c("1973", "2003"),
  main = "Maximum Weekly Temperatures",
  ylab = "Temperature (C)"
)

temperature_difference <- sydney$Temp2003 - sydney$Temp1973
hist(
  temperature_difference, breaks = 20,
  main = "Weekly Temperature Differences (2003 - 1973)",
  xlab = "Temperature difference (C)"
)
abline(v = mean(temperature_difference), col = "blue", lwd = 2)

# Paired, two-sided Wilcoxon signed-rank test with continuity correction.
wilcoxon_result <- wilcox.test(
  sydney$Temp2003,
  sydney$Temp1973,
  paired = TRUE,
  alternative = "two.sided",
  exact = FALSE,
  correct = TRUE
)
print(wilcoxon_result)

# Leave-one-out cross-validation error for the normal-kernel smoother.
loocv_error <- function(x, y, bandwidth) {
  squared_errors <- vapply(seq_along(x), function(i) {
    prediction <- stats::ksmooth(
      x[-i], y[-i],
      kernel = "normal",
      bandwidth = bandwidth,
      x.points = x[i]
    )$y
    (y[i] - prediction)^2
  }, numeric(1))
  mean(squared_errors)
}

week <- sydney$Week
temperature_1973 <- sydney$Temp1973
temperature_2003 <- sydney$Temp2003
bandwidth_candidates <- seq(1, 20, by = 0.5)
cv_error_1973 <- vapply(
  bandwidth_candidates,
  function(h) loocv_error(week, temperature_1973, h),
  numeric(1)
)
cv_error_2003 <- vapply(
  bandwidth_candidates,
  function(h) loocv_error(week, temperature_2003, h),
  numeric(1)
)

cat(
  "LOOCV-minimizing bandwidths:",
  "1973 =", bandwidth_candidates[which.min(cv_error_1973)],
  "; 2003 =", bandwidth_candidates[which.min(cv_error_2003)], "\n"
)
par(mfrow = c(1, 2))
plot(
  bandwidth_candidates, cv_error_1973, type = "b", col = "blue",
  xlab = "Bandwidth", ylab = "LOOCV mean squared error",
  main = "Bandwidth Selection: 1973"
)
abline(v = bandwidth_candidates[which.min(cv_error_1973)], col = "red", lty = 2)
plot(
  bandwidth_candidates, cv_error_2003, type = "b", col = "red",
  xlab = "Bandwidth", ylab = "LOOCV mean squared error",
  main = "Bandwidth Selection: 2003"
)
abline(v = bandwidth_candidates[which.min(cv_error_2003)], col = "blue", lty = 2)
par(mfrow = c(1, 1))

# The report chose h = 8 for smoother seasonal curves rather than the
# bandwidths that minimize LOOCV error.
chosen_bandwidth <- 8
kernel_1973 <- stats::ksmooth(
  week, temperature_1973, kernel = "normal",
  bandwidth = chosen_bandwidth, x.points = week
)
kernel_2003 <- stats::ksmooth(
  week, temperature_2003, kernel = "normal",
  bandwidth = chosen_bandwidth, x.points = week
)

plot(
  week, temperature_1973, type = "l", col = "blue",
  ylim = range(c(temperature_1973, temperature_2003, kernel_1973$y, kernel_2003$y)),
  xlab = "Week", ylab = "Maximum temperature (C)",
  main = "Kernel-Smoothed Temperatures"
)
lines(week, kernel_1973$y, col = "blue", lwd = 2, lty = 2)
lines(week, temperature_2003, col = "red")
lines(week, kernel_2003$y, col = "red", lwd = 2, lty = 2)
legend(
  "topright",
  legend = c("1973 data", "1973 kernel", "2003 data", "2003 kernel"),
  col = c("blue", "blue", "red", "red"),
  lty = c(1, 2, 1, 2),
  lwd = c(1, 2, 1, 2)
)

# Local polynomial alternative from the report.
loess_span <- 0.3
loess_1973 <- stats::loess(temperature_1973 ~ week, span = loess_span)
loess_2003 <- stats::loess(temperature_2003 ~ week, span = loess_span)
loess_prediction_1973 <- predict(loess_1973, newdata = data.frame(week = week))
loess_prediction_2003 <- predict(loess_2003, newdata = data.frame(week = week))

par(mfrow = c(1, 2))
plot(week, temperature_1973, col = "blue", xlab = "Week", ylab = "Temperature (C)")
lines(week, kernel_1973$y, col = "blue", lty = 2, lwd = 2)
lines(week, loess_prediction_1973, col = "darkgreen", lwd = 2)
legend("topright", c("Data", "Kernel", "LOESS"), col = c("blue", "blue", "darkgreen"),
       lty = c(1, 2, 1), lwd = c(1, 2, 2))
plot(week, temperature_2003, col = "red", xlab = "Week", ylab = "Temperature (C)")
lines(week, kernel_2003$y, col = "red", lty = 2, lwd = 2)
lines(week, loess_prediction_2003, col = "purple", lwd = 2)
legend("topright", c("Data", "Kernel", "LOESS"), col = c("red", "red", "purple"),
       lty = c(1, 2, 1), lwd = c(1, 2, 2))
par(mfrow = c(1, 1))
