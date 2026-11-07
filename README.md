# Statistical Analysis

A collection of statistical reports and independent R scripts. Each analysis
explores a different method or model, with the original PDF write-up retained
alongside its runnable script.

## Repository contents

| Analysis | Standalone script | Original report |
|---|---|---|
| Social factors and weekly wages | [`R/01_wages_social_factors.R`](./R/01_wages_social_factors.R) | [Social factors affecting wages](./Social%20factors%20affecting%20wages.pdf) |
| Smoking status | [`R/02_smoking_status.R`](./R/02_smoking_status.R) | [How smoking status is affected by other variables](./How%20smoking%20status%20is%20affected%20by%20other%20variables.pdf) |
| Sydney temperatures, 1973 vs 2003 | [`R/03_sydney_temperatures.R`](./R/03_sydney_temperatures.R) | [Is there any difference in temperatures in Sydney?](./Is%20there%20any%20difference%20in%20temperatures%20in%20Sydney%201973%20vs%202003.pdf) |
| Doctor visit counts | [`R/04_doctor_visits.R`](./R/04_doctor_visits.R) | [Poisson regression and bootstrap methods](./A%20statistical%20analysis%20of%20doctor%20visit%20counts%20using%20Poisson%20Regression%20and%20Bootstrap%20methods-1.pdf) |

The scripts are independent: run one analysis at a time with its corresponding
input data. The reports are the source for the questions, interpretations, and
reported results; the R files organize the analyses into executable scripts,
add input checks, and correct incomplete or inconsistent code from the
appendices.

## Requirements and running the analyses

Install R, then install the packages used by the analyses:

```r
install.packages(c("ggplot2", "dplyr", "MASS", "sandwich", "lmtest", "boot"))
```

Run commands from the repository root. Each script looks for its documented
default file under `data/`; alternatively, pass the path to a CSV as the first
argument:

```sh
Rscript R/01_wages_social_factors.R data/wages.csv
Rscript R/02_smoking_status.R data/smoking.csv
Rscript R/03_sydney_temperatures.R data/sydney_temperatures.csv
Rscript R/04_doctor_visits.R data/dvisus.csv
```

The scripts print summaries and model results to the console and draw plots on
the active R graphics device. They report an error if the input file or required
columns are missing. The data files are not included in this repository, so
download or prepare them separately before running the analyses.

## Input data

All inputs are CSV files. Column names are case-sensitive.

| Script | Default input | Required columns and notes |
|---|---|---|
| Wages | `data/wages.csv` | `Region`, `Location`, `Age`, `EducCateg`, `EducCode`, `Race`, `Wage`. `Wage` is weekly earnings and must be positive for the log-wage model. |
| Smoking | `data/smoking.csv` | `smoker`, `ban`, `age`, `age2`, `education`, `aahisp`, `female`. The four indicator columns use `0`/`1`; `age2` is age squared. Education can be coded 1-5 or use labels for high-school dropout, high-school graduate, some college, college graduate, and master's degree or higher. |
| Sydney temperatures | `data/sydney_temperatures.csv` | `Week`, `Temp1973`, `Temp2003`. Each week should occur once and have a paired temperature in both years. |
| Doctor visits | `data/dvisus.csv` | `docvis`, `age`, `lincome`, `female`, `educ`, `married`, `white`, `hispanic`, `private`, `chronic`. `docvis` is a non-negative count; the indicator columns use `0`/`1`. |

The wage report describes a 2011 U.S. sample and the doctor-visits report
describes U.S. private-sector workers from 2001 who were not publicly insured
or self-employed. The smoking report describes a U.S. sample but does not give
a downloadable data source or filename. The Sydney report cites the Australian
Bureau of Meteorology and uses Sydney Airport observations. The reports do not
provide complete download links or license information for every dataset;
check the original source and its terms before redistributing any data.

## Analyses and reported findings

### 1. Social factors affecting wages

The script starts with descriptive summaries and wage distributions by
education, race, location, and region. It fits a multiple linear regression
for weekly wages, checks residual diagnostics, examines a Box-Cox
transformation, and fits a log-wage model. Cook's distance is used to flag
influential observations for a refit, and the script compares numerical and
categorical encodings of education with BIC.

The report notes an overall mean weekly wage of about `$1,056.40` in its sample.
Its cleaned log-wage model finds associations with age, education, location,
race, and some regions. These are sample associations, not causal effects; the
report also notes that its sample contains substantially more white than Black
respondents.

### 2. Smoking status

The script explores smoking prevalence by education, workplace smoking ban,
race/ethnicity, and gender. It fits a binomial logistic regression with
interactions between `aahisp` (African-American or Hispanic indicator) and the
other predictors, compares the full model with the selected model using a
likelihood-ratio test and AIC, and reports odds ratios with profile-likelihood
confidence intervals. It also checks model diagnostics, refits after flagging
high-influence observations as a sensitivity analysis, and compares categorical
with numeric education.

The report's selected model associates workplace smoking bans and higher
education with lower smoking odds. It reports interactions between
race/ethnicity and education, age squared, and gender, so main-effect odds
ratios should not be read as applying uniformly to every subgroup.

### 3. Sydney temperatures

The script compares paired weekly maximum temperatures for Sydney Airport in
1973 and 2003. It includes descriptive plots, a two-sided paired Wilcoxon
signed-rank test with continuity correction, leave-one-out cross-validation
(LOOCV) for normal-kernel smoothing, and LOESS as an alternative smoother.

The report gives a Wilcoxon p-value of `0.6132`; at the 5% level it does not
reject the no-difference null hypothesis. It uses a kernel bandwidth of `8`
for smoother seasonal curves and a LOESS span of `0.3`. The report also notes
that weekly temperatures are serially related and that symmetry of paired
differences is an assumption of the signed-rank test. Failure to reject the
null is not proof that the years' temperature distributions are identical.

### 4. Doctor visit counts

The script fits a Poisson generalized linear model with a log link, calculates
HC0 sandwich standard errors, and uses row-resampling bootstraps for coefficient
standard errors, a joint Wald test, and studentized confidence intervals for a
restricted model. It also reports the mean and variance of the count outcome
and plots the main variables.

The report's sample has 1,084 observations, a mean of `3.39` visits, and a
variance of `24.50`, indicating overdispersion relative to a standard Poisson
model. Its joint test of the income, marital-status, and Hispanic coefficients
has a reported p-value of `0.015` with conventional model-based covariance,
`0.662` with the HC0 sandwich covariance, and `0.677` with bootstrap
covariance. The report also gives an empirical bootstrap p-value of `0.798`.
The code follows the report's row-resampling procedure for that empirical
p-value; it does not impose the null hypothesis when generating the bootstrap
distribution. Interpret that result with this limitation in mind.

## Reproducibility and interpretation

The reports contain the original explanations and output, while the repository
does not contain the underlying datasets. Exact results therefore depend on
using the same source data, data cleaning, and category coding described in
each report. Bootstrap analyses use 1,000 resamples and fixed seeds in the
doctor-visits script to make runs reproducible for a given input.

These are exploratory analyses of observational data. Regression coefficients,
odds ratios, and group differences describe associations under each model's
assumptions; they do not establish causal effects. Review the assumptions and
diagnostics in the linked reports before using the results for decisions.
