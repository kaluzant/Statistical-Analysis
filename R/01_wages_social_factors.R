# ---------------------------------------------------------------------------
# Social factors affecting weekly wages
# Original report: Social factors affecting wages.pdf
#
# Data: data/wages.csv
#   Region, Location, Age, EducCateg, EducCode, Race, Wage
#
# Packages: ggplot2, dplyr, MASS
# ---------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(MASS)  # boxcox()

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else file.path("data", "wages.csv")
if (!file.exists(data_file)) {
  stop(
    "Data file not found: ", data_file,
    "\nPlace wages.csv in data/ or pass its path as the first argument."
  )
}

wages <- read.csv(data_file, stringsAsFactors = TRUE)

required_columns <- c(
  "Region", "Location", "Age", "EducCateg", "EducCode", "Race", "Wage"
)
missing_columns <- setdiff(required_columns, names(wages))
if (length(missing_columns) > 0) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "))
}
if (any(!is.finite(wages$Wage[!is.na(wages$Wage)]))) {
  stop("Wage must contain only finite values or NA.")
}
if (any(wages$Wage <= 0, na.rm = TRUE)) {
  stop("Wage must be positive wherever it is observed because the analysis uses log(Wage).")
}

# ---------------------------------------------------------------------------
# 1. Exploratory data analysis
# ---------------------------------------------------------------------------

str(wages)
summary(wages)

mean_wage <- mean(wages$Wage, na.rm = TRUE)
cat("Overall mean weekly wage:", round(mean_wage, 2), "\n")

# Boxplot of wages by a grouping variable, with group means and the overall
# mean shown for reference.
create_boxplot <- function(group_var, title) {
  means <- wages %>%
    group_by(.data[[group_var]]) %>%
    summarise(MeanWage = mean(Wage, na.rm = TRUE), .groups = "drop")

  ggplot(wages, aes(x = .data[[group_var]], y = Wage,
                    fill = .data[[group_var]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                 fill = "red") +
    geom_text(
      data = means,
      aes(x = .data[[group_var]], y = MeanWage, label = round(MeanWage, 1)),
      vjust = -0.5, colour = "blue", size = 3, inherit.aes = FALSE
    ) +
    geom_hline(yintercept = mean_wage, colour = "red", linetype = "dashed",
               linewidth = 1) +
    labs(title = title, x = group_var, y = "Weekly Earnings") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 10),
      legend.position = "none"
    )
}

plot_race <- create_boxplot("Race", "Weekly Earnings by Race")
plot_location <- create_boxplot("Location", "Weekly Earnings by Location")
plot_region <- create_boxplot("Region", "Weekly Earnings by Region")
plot_educ <- create_boxplot("EducCateg", "Weekly Earnings by Education Category") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

print(plot_race)
print(plot_location)
print(plot_region)
print(plot_educ)

# Education code (ordered) vs wage
educ_means <- wages %>%
  group_by(EducCode) %>%
  summarise(MeanWage = mean(Wage, na.rm = TRUE), .groups = "drop")

print(
  ggplot(wages, aes(x = factor(EducCode), y = Wage)) +
    geom_boxplot(outlier.shape = NA, fill = "skyblue", alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                 fill = "red") +
    geom_text(
      data = educ_means,
      aes(x = factor(EducCode), y = MeanWage, label = round(MeanWage, 1)),
      vjust = -0.5, colour = "blue", size = 3, inherit.aes = FALSE
    ) +
    geom_hline(yintercept = mean_wage, colour = "red", linetype = "dashed",
               linewidth = 1) +
    labs(title = "Weekly Earnings vs Education Codes",
         x = "Education Code", y = "Weekly Earnings") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
)

# Average wage by age
average_wage_by_age <- wages %>%
  group_by(Age) %>%
  summarise(Wage = mean(Wage, na.rm = TRUE), .groups = "drop")

print(
  ggplot(average_wage_by_age, aes(x = Age, y = Wage)) +
    geom_line(colour = "blue") +
    geom_point(colour = "red") +
    geom_hline(yintercept = mean_wage, colour = "red", linetype = "dashed",
               linewidth = 1) +
    labs(title = "Average Wage by Age", x = "Age", y = "Average Wage") +
    theme_minimal()
)

# Quartiles and IQR for each grouping variable
calculate_quartiles_iqr <- function(group_var) {
  wages %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      Q1 = quantile(Wage, 0.25, na.rm = TRUE),
      Median = median(Wage, na.rm = TRUE),
      Q3 = quantile(Wage, 0.75, na.rm = TRUE),
      IQR = IQR(Wage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(.data[[group_var]])
}

print(calculate_quartiles_iqr("Race"))
print(calculate_quartiles_iqr("Location"))
print(calculate_quartiles_iqr("EducCateg"))
print(calculate_quartiles_iqr("Region"))

# ---------------------------------------------------------------------------
# 2. Initial normal linear model and diagnostics
# ---------------------------------------------------------------------------

lm1 <- lm(Wage ~ Age + EducCode + Location + Region + Race, data = wages)
summary(lm1)

par(mfrow = c(2, 2))
plot(lm1)
par(mfrow = c(1, 1))

studentised_residuals <- rstudent(lm1)

par(mfrow = c(1, 2))
plot(fitted(lm1), studentised_residuals,
     xlab = "Fitted Values", ylab = "Studentised Residuals",
     main = "Studentised Residuals vs Fitted Values")
abline(h = 0, col = "red", lty = 2)
qqnorm(studentised_residuals,
       main = "Normal Q-Q Plot of Studentised Residuals")
qqline(studentised_residuals, col = "red")
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------
# 3. Box-Cox transformation and log-wage model
# ---------------------------------------------------------------------------

boxcox_result <- boxcox(lm1, lambda = seq(-2, 2, by = 0.1))
best_lambda <- boxcox_result$x[which.max(boxcox_result$y)]
cat("Best lambda for Box-Cox transformation:", round(best_lambda, 4), "\n")

# The report chooses lambda = 0 (log wage) for interpretability.
wages$LogWage <- log(wages$Wage)
lm2 <- lm(LogWage ~ Age + EducCode + Location + Race + Region, data = wages)
summary(lm2)

studentised_residuals <- rstudent(lm2)

par(mfrow = c(1, 2))
plot(fitted(lm2), studentised_residuals,
     xlab = "Fitted Values", ylab = "Studentised Residuals",
     main = "Studentised Residuals vs Fitted Values (log model)")
abline(h = 0, col = "red", lty = 2)
abline(h = c(-3, 3), col = "red", lty = 2)
qqnorm(studentised_residuals,
       main = "Normal Q-Q Plot of Studentised Residuals (log model)")
qqline(studentised_residuals, col = "red")
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------
# 4. Influential observations (Cook's distance)
# ---------------------------------------------------------------------------

n <- nrow(wages)
p <- length(coef(lm2))
cook_dist <- cooks.distance(lm2)
cook_threshold <- 8 / (n - 2 * p)
influential_obs <- which(cook_dist > cook_threshold)

plot(cook_dist, type = "h",
     main = "Cook's Distance with Threshold",
     xlab = "Observation", ylab = "Cook's Distance",
     ylim = c(0, max(c(cook_dist, cook_threshold))))
abline(h = cook_threshold, col = "red", lty = 2, lwd = 2)

cat("Removing", length(influential_obs), "influential observations\n")
print(influential_obs)

wages_clean <- wages[-influential_obs, ]
lm2_clean <- lm(LogWage ~ Age + EducCode + Location + Race + Region,
                data = wages_clean)
summary(lm2_clean)

studentised_residuals <- rstudent(lm2_clean)

par(mfrow = c(1, 2))
plot(fitted(lm2_clean), studentised_residuals,
     xlab = "Fitted Values", ylab = "Studentised Residuals",
     main = "Studentised Residuals vs Fitted Values (cleaned)")
abline(h = 0, col = "red", lty = 2)
abline(h = c(-3, 3), col = "red", lty = 2)
qqnorm(studentised_residuals,
       main = "Normal Q-Q Plot of Studentised Residuals (cleaned)")
qqline(studentised_residuals, col = "red")
par(mfrow = c(1, 1))

# Coefficient table for the cleaned and transformed model
coef_table <- summary(lm2_clean)$coefficients
print(data.frame(
  Variable = rownames(coef_table),
  Coefficient = coef_table[, "Estimate"],
  Std.Error = coef_table[, "Std. Error"],
  P.Value = coef_table[, "Pr(>|t|)"]
))

# ---------------------------------------------------------------------------
# 5. Education: numerical code vs categorical factor
# ---------------------------------------------------------------------------

# The report compares the two education codings on the full sample.
lm_edu_num <- lm(LogWage ~ Age + EducCode + Location + Race + Region,
                 data = wages)
lm_edu_cat <- lm(LogWage ~ Age + EducCateg + Location + Race + Region,
                 data = wages)

summary(lm_edu_num)
summary(lm_edu_cat)

cat("BIC (numerical education):", BIC(lm_edu_num), "\n")
cat("BIC (categorical education):", BIC(lm_edu_cat), "\n")
