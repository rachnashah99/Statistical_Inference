## ----install_pckgs--------------------------------------------------------------------------------------------------------------------
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(gridExtra)) install.packages("gridExtra")
if(!require(pheatmap)) install.packages("pheatmap")
if(!require(carat)) install.packages("caret")
if(!require(pscl)) install.packages("pscl")
if(!require(Metrics)) install.packages("Metrics")
if(!require(ggstats)) install.packages("ggstats")
if(!require(DescTools)) install.packages("DescTools")
# Packages required to install DMwR
if(!require(zoo) | !require(xts) | !require(quantmod)) install.packages(c("zoo","xts","quantmod"))
if(!require(abind) | !require(ROCR)) install.packages(c("abind", "ROCR"))
if(!require(DMwR)) install.packages( "https://cran.r-project.org/src/contrib/Archive/DMwR/DMwR_0.4.1.tar.gz", repos=NULL, type="source" )


## ----imp_pckgs------------------------------------------------------------------------------------------------------------------------
#library import
library(tidyverse)
library(gridExtra)
library(pheatmap)
library(zoo)
library(knitr)
library(caret)
library(DMwR)
library(pscl)
library(Metrics)
library(ggstats)
library(DescTools)


## ----import_data----------------------------------------------------------------------------------------------------------------------
# Load the patient data from GitHub
patient.data <- read.csv("https://raw.github.com/Neatherblok/Pregnancy_Caused_Diabetes_Predictor/main/data/patients.csv")

# Read first few rows of dataset
head(patient.data)


## ----rm_duplicates--------------------------------------------------------------------------------------------------------------------
# Data Dimension before deleting duplicates
dim(patient.data)

# Select only all unique observations
unique.data<-unique(patient.data)

# Data Dimension after deleting duplicates
dim(unique.data)


## ----target_pie_chart_1---------------------------------------------------------------------------------------------------------------
# Target distribution pie chart
ggplot(patient.data, aes(x = "", fill = factor(Diagnosis))) +
  geom_bar(width = 1) +
  coord_polar(theta = "y") +
  labs(fill = "Diagnosis") +
  ggtitle("Diagnosis Distribution") +
  theme_void()


## ----balancing_target-----------------------------------------------------------------------------------------------------------------
# Create subset of Diagnosis = 0 observations
no.diabetes <- patient.data[patient.data$Diagnosis == 0, ]

#  Finding the odd rows: Check for 0 values in specified columns and remove the rows
filtered.data <- no.diabetes[!(no.diabetes$Glucose == 0 | no.diabetes$BloodPressure == 0 | no.diabetes$SkinThickness == 0 | no.diabetes$Insulin == 0 | no.diabetes$BMI == 0), ]

# Combine the filtered rows with the remaining rows that have Diagnosis value not equal to 0
balanced.patient.data <- rbind(filtered.data, patient.data[patient.data$Diagnosis != 0, ])

# row index are reset
row.names(balanced.patient.data) <- NULL

# Show first few rows of data set
head(balanced.patient.data)

# Show last few rows of data set
tail(balanced.patient.data)


## ----target_pie_chart_2---------------------------------------------------------------------------------------------------------------
# Distribution of target variable post downscaling
ggplot(balanced.patient.data, aes(x = "", fill = factor(Diagnosis))) +
  geom_bar(width = 1) +
  coord_polar(theta = "y") +
  labs(fill = "Diagnosis") +
  ggtitle("Diagnosis_Distribution") +
  theme_void()

# Diagnosis distribution in numbers which follows approx 50%:50%
diagnosis.counts <- table(balanced.patient.data$Diagnosis)
print(diagnosis.counts)


## ----descriptive_stat-----------------------------------------------------------------------------------------------------------------
# Summarize descriptive stats
summary(balanced.patient.data[, 1:8])


## ----box_plot-------------------------------------------------------------------------------------------------------------------------
# Box Plot
independent.features <- balanced.patient.data[, 1:8]
long_data <- pivot_longer(independent.features, everything(), names_to = "variable", values_to = "value")
ggplot(long_data, aes(x = variable, y = value)) +
  geom_boxplot() +
  xlab("") +
  ylab("Value") +
  theme_minimal()


## ----histogram------------------------------------------------------------------------------------------------------------------------
# Histogram
Pregnancies <- ggplot(independent.features, aes(x = Pregnancies)) + geom_histogram() + ggtitle("Pregnancies")
Glucose <- ggplot(independent.features, aes(x = Glucose)) + geom_histogram() + ggtitle("Glucose")
BloodPressure <- ggplot(independent.features, aes(x = BloodPressure)) + geom_histogram() + ggtitle("Blood Pressure")
SkinThickness <- ggplot(independent.features, aes(x = SkinThickness)) + geom_histogram() + ggtitle("Skin Thickness")
Insulin <- ggplot(independent.features, aes(x = Insulin)) + geom_histogram() + ggtitle("Insulin")
BMI <- ggplot(independent.features, aes(x = BMI)) + geom_histogram() + ggtitle("BMI")
Pedigree <- ggplot(independent.features, aes(x = Pedigree)) + geom_histogram() + ggtitle("Pedigree")
Age <- ggplot(independent.features, aes(x = Age)) + geom_histogram() + ggtitle("Age")
grid.arrange(
  Pregnancies, Glucose, BloodPressure, SkinThickness,
  Insulin, BMI, Pedigree, Age,
  nrow = 3
)


## ----miss_val_handling----------------------------------------------------------------------------------------------------------------
# All zeros are treated as missing values and replaced with median
for (i in 2:8) {
  balanced.patient.data[, i] <- ifelse(balanced.patient.data[, i] == 0, median(balanced.patient.data[balanced.patient.data[, i] != 0, i], na.rm = TRUE), balanced.patient.data[, i])
}

# Change name of data set
cleaned.patient.data <- balanced.patient.data

# Read First Few rows to see if 0's were replaced by median of the column
head(cleaned.patient.data)


## ----outlier_handling-----------------------------------------------------------------------------------------------------------------
# Removing outliers using capping method
treating.outliers <- function(x) {
  qntls <- quantile(x, probs = c(0.25, 0.75))
  iqr <- qntls[2] - qntls[1]
  lb <- qntls[1] - 1.5 * iqr
  ub <- qntls[2] + 1.5 * iqr
  x[x < lb] <- lb
  x[x > ub] <- ub
  return(x)
}

# Execute function for all independent variables
cleaned.patient.data[, 1:8] <- apply(cleaned.patient.data[, 1:8], 2, treating.outliers)

# By the revised box plot you can see that the outliers are treated and not present in the data that we are going to further analyse
cleaned.independent.variables <- cleaned.patient.data[, 1:8]
long_data <- pivot_longer(cleaned.independent.variables, everything(), names_to = "variable", values_to = "value")
ggplot(long_data, aes(x = variable, y = value)) +
  geom_boxplot() +
  xlab("") +
  ylab("Value") +
  theme_minimal()


## -------------------------------------------------------------------------------------------------------------------------------------
# check summary post data cleaning
summary(cleaned.patient.data)


## ----correlogram----------------------------------------------------------------------------------------------------------------------
# Check for correlation post data cleaning
Correlation.Mtrx <- cor(cleaned.patient.data)
Correlation.Percentages <- round(Correlation.Mtrx * 100, 2)
pheatmap(Correlation.Percentages,
         main = "Correlation Matrix (%)",
         fontsize = 8,
         fontsize_row = 6,
         fontsize_col = 6,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         show_rownames = TRUE,
         show_colnames = TRUE,
         display_numbers = TRUE,
         number_color = "black",
         number_format = "%.2f%%")



## ----rand_split_1---------------------------------------------------------------------------------------------------------------------
# pre-set properties for random split 
train.size = floor(0.8*nrow(cleaned.patient.data))
set.seed(130)

# randomly split data in r according to generator seed 130
picked = sample(seq_len(nrow(cleaned.patient.data)),size = train.size)

# Divide the data set up in two data set according to picked sample
train.data.not.normal = cleaned.patient.data[picked,]
test.data.not.normal = cleaned.patient.data[-picked,]

# Reset row index
rownames(train.data.not.normal) <- NULL
rownames(test.data.not.normal) <- NULL

# Show the two split data sets
head(train.data.not.normal)
head(test.data.not.normal)


## ----normalize------------------------------------------------------------------------------------------------------------------------
# Set up normalization function
process <- preProcess(cleaned.patient.data, method=c("range"))

# Execute Min-Max normalization on the cleaned data set
norm.data <- predict(process, cleaned.patient.data)

# Show normalized data set
head(norm.data)


## ----rand_split_2---------------------------------------------------------------------------------------------------------------------
# Divide the data set up in two data set according to picked sample
train.data.normal = norm.data[picked,]
test.data.normal = norm.data[-picked,]

# Reset row index
rownames(train.data.normal) <- NULL
rownames(test.data.normal) <- NULL

# Show the two split data sets
head(train.data.normal)
head(test.data.normal)


## ----train_lrm_1----------------------------------------------------------------------------------------------------------------------
# Train Logistic Regression model 2
model.1 <- glm(Diagnosis ~ Pregnancies + Age + Glucose + BMI, data = train.data.normal, family = binomial(link='logit'))


## ----summ_lrm_1-----------------------------------------------------------------------------------------------------------------------
# Summarize the findings of the first Logistic Regression model
summary(model.1)


## ----pred_lrm_1-----------------------------------------------------------------------------------------------------------------------
# Predict the outcome using Logistic Regression Model 1
mdl.1.fitted.results <- predict(model.1, test.data.normal[,c("Pregnancies", "Age", "Glucose", "BMI")], type='response')

# Transform predicted probability values with a probability higher than 0.5 to 1
# Transform values with probability lower than 0.5 to 0
predictions.mdl.1 <- ifelse(mdl.1.fitted.results >= 0.5, 1, 0)


## ----conf_matrix_lrm_1----------------------------------------------------------------------------------------------------------------
# Create a classification report for Logistic Regression Model 1
# This includes a confusion matrix, accuracy score, ACI score, confidence interval and other scores
classification.report.mdl.1 <- caret::confusionMatrix(factor(predictions.mdl.1), factor(test.data.normal$Diagnosis))
print(classification.report.mdl.1)


## ----train_lrm_2----------------------------------------------------------------------------------------------------------------------
# Train Logistic Regression model 2
model.2 <- glm(Diagnosis ~.,family=binomial(link='logit'),data=train.data.normal)


## ----summ_lrm_2-----------------------------------------------------------------------------------------------------------------------
# Summarize the findings of the second Logistic Regression model
summary(model.2)


## ----pred_lrm_2-----------------------------------------------------------------------------------------------------------------------
# Predict the outcome using Logistic Regression Model 2
mdl.2.fitted.results <- predict(model.2, test.data.normal, type='response')

# Transform predicted probability values with a probability higher than 0.5 to 1
# Transform values with probability lower than 0.5 to 0
predictions.mdl.2 <- ifelse(mdl.2.fitted.results >= 0.5, 1, 0)


## ----conf_matrix_lrm_2----------------------------------------------------------------------------------------------------------------
# Create a classification report
# This includes a confusion matrix, accuracy score, ACI score, confidence interval and other scores
classification.report.mdl.2 <- caret::confusionMatrix(factor(predictions.mdl.2), factor(test.data.normal$Diagnosis))
print(classification.report.mdl.2)


## ----train_lrm_3----------------------------------------------------------------------------------------------------------------------
# Train Logistic Regression model 3
model.3 <- glm(Diagnosis ~Pregnancies + Age + Glucose + BMI, data = train.data.not.normal, family = binomial(link='logit'))


## ----summ_lrm_3-----------------------------------------------------------------------------------------------------------------------
# Summarize the findings of the third Logistic Regression model
summary(model.3)


## ----pred_lrm_3-----------------------------------------------------------------------------------------------------------------------
# Predict the outcome using Logistic Regression Model 3
mdl.3.fitted.results <- predict(model.3, test.data.not.normal[,c("Pregnancies", "Age", "Glucose", "BMI")], type='response')

# Transform predicted probability values with a probability higher than 0.5 to 1
# Transform values with probability lower than 0.5 to 0
predictions.mdl.3 <- ifelse(mdl.3.fitted.results >= 0.5, 1, 0)


## ----conf_matrix_lrm_3----------------------------------------------------------------------------------------------------------------
# Create a classification report for Logistic Regression Model 3
# This includes a confusion matrix, accuracy score, ACI score, confidence interval and other scores
classification.report.mdl.3 <- caret::confusionMatrix(factor(predictions.mdl.3), factor(test.data.not.normal$Diagnosis))
print(classification.report.mdl.3)


## ----train_lrm_4----------------------------------------------------------------------------------------------------------------------
# Train Logistic Regression model 4
model.4 <- glm(Diagnosis ~ ., data = train.data.not.normal, family = binomial(link='logit'))


## ----summ_lrm_4-----------------------------------------------------------------------------------------------------------------------
# Summarize the findings of the fourth Logistic Regression model
summary(model.4)


## ----pred_lrm_4-----------------------------------------------------------------------------------------------------------------------
# Predict the outcome using Logistic Regression Model 4
mdl.4.fitted.results <- predict(model.4, test.data.not.normal,type='response')

# Transform predicted probability values with a probability higher than 0.5 to 1
# Transform values with probability lower than 0.5 to 0
predictions.mdl.4 <- ifelse(mdl.4.fitted.results >= 0.5, 1, 0)


## ----conf_matrix_lrm_4----------------------------------------------------------------------------------------------------------------

# Create a classification report for Logistic Regression Model 4
# This includes a confusion matrix, accuracy score, ACI score, confidence interval and other scores
classification.report.mdl.4 <- caret::confusionMatrix(factor(predictions.mdl.4), factor(test.data.not.normal$Diagnosis))
print(classification.report.mdl.4)


## ----train_lrm_5----------------------------------------------------------------------------------------------------------------------
# Train Logistic Regression model 5
model.5 <- glm(Diagnosis ~ Pregnancies + Glucose + BMI, data = train.data.normal, family = binomial(link='logit'))


## ----summ_lrm_5-----------------------------------------------------------------------------------------------------------------------
# Summarize the findings of the fifth Logistic Regression model
summary(model.5)


## ----pred_lrm_5-----------------------------------------------------------------------------------------------------------------------
# Predict the outcome using Logistic Regression Model 5
mdl.5.fitted.results <- predict(model.5, test.data.normal,type='response')

# Transform predicted probability values with a probability higher than 0.5 to 1
# Transform values with probability lower than 0.5 to 0
predictions.mdl.5 <- ifelse(mdl.5.fitted.results >= 0.5, 1, 0)


## ----conf_matrix_lrm_5----------------------------------------------------------------------------------------------------------------
# Create a classification report for Logistic Regression Model 5
# This includes a confusion matrix, accuracy score, ACI score, confidence interval and other scores
classification.report.mdl.5 <- caret::confusionMatrix(factor(predictions.mdl.5), factor(test.data.normal$Diagnosis))
print(classification.report.mdl.5)


## ----comparison-----------------------------------------------------------------------------------------------------------------------
df <- data.frame(Metrics=c("Accuracy", "AIC", "Confidence Interval Lower", "Confidence Interval Higher"),
                  select_feature_normal=c(classification.report.mdl.1$overall[1], AIC(model.1),
                              classification.report.mdl.1$overall[3], classification.report.mdl.1$overall[4]),
                 
                  all_normal=c(classification.report.mdl.2$overall[1], AIC(model.2), 
                               classification.report.mdl.2$overall[3], classification.report.mdl.2$overall[4]),
                 
                  select_feature_not_normal=c(classification.report.mdl.3$overall[1], AIC(model.3),
                                classification.report.mdl.3$overall[3], classification.report.mdl.3$overall[4]),
                 
                  all_not_normal=c(classification.report.mdl.4$overall[1], AIC(model.4), 
                                classification.report.mdl.4$overall[3], classification.report.mdl.4$overall[4]),
                 
                  without_age_normal=c(classification.report.mdl.5$overall[1], AIC(model.5), 
                                classification.report.mdl.5$overall[3], classification.report.mdl.5$overall[4]))

# Set metrics column as index
rownames(df) <- df$Metrics

# Remove original metrics column from data frame
df$Metrics <- NULL

# Show results
df


## ----coef_plot------------------------------------------------------------------------------------------------------------------------
# Visualize logistic regression coefficients
ggcoef_model(model.1)
ggcoef_model(model.2)
ggcoef_model(model.3)
ggcoef_model(model.4)
ggcoef_model(model.5)


## ----anova_runs-----------------------------------------------------------------------------------------------------------------------
# Number of times one anova comparison has been run
anova.runs <- list(
  anova.1 = 0,
  anova.2 = 0,
  anova.3 = 0,
  anova.4 = 0,
  anova.5 = 0
)


## ----hypothesis-----------------------------------------------------------------------------------------------------------------------
# Testing significance of model 1 compared to model 3
# H0: Model 1 (min-max scaled;best features selected) is more significant than model 3 (original scale; best features selected)
execute.anova.1 <- function() {
  eval.parent(substitute(anova.runs$anova.1 <- anova.runs$anova.1 + 1))
  return(anova(model.1,model.3, test='LRT'))
}

# Testing significance of model 2 compared to model 4
# H0: Model 2 (min-max scaled; all features) is more significant than model 4 (original scale; all features)
execute.anova.2 <- function() {
  eval.parent(substitute(anova.runs$anova.2 <- anova.runs$anova.2 + 1))
  return(anova(model.2, model.4, test='LRT'))
}

# Testing significance of model 1 compared to model 2
# H0: Model 1 (min-max scaled; best features selected) is more significant than model 2 (min-max scaled; all features)
execute.anova.3 <- function() {
  eval.parent(substitute(anova.runs$anova.3 <- anova.runs$anova.3 + 1))
  return(anova(model.1, model.2, test='LRT'))
}

# Testing significance of model 1 compared to model 5
# H0: Model 1 (min-max scaled; best features selected) is more significant than model 5 (min-max scaled; best features without Age)
execute.anova.4 <- function() {
  eval.parent(substitute(anova.runs$anova.4 <- anova.runs$anova.4 + 1))
  return(anova(model.1, model.5, test='LRT'))
}

# Testing significance of model 3 compared to model 4
# H0: Model 3 (original scale; best features selected) is more significant than model 4 (original scale; all features)
execute.anova.5 <- function() {
  eval.parent(substitute(anova.runs$anova.5 <- anova.runs$anova.5 + 1))
  return(anova(model.3, model.4, test='LRT'))
}

# Execute the four Likeliness Ratio tests
# Print Likeliness Ratio test results
anova.1 <- execute.anova.1()
print(anova.1)
anova.2 <- execute.anova.2()
print(anova.2)
anova.3 <- execute.anova.3()
print(anova.3)
anova.4 <- execute.anova.4()
print(anova.4)
anova.5 <- execute.anova.5()
print(anova.5)



## ----run_more-------------------------------------------------------------------------------------------------------------------------
# Run Likeliness Ratio tests 3 more times, to use bonferroni correction
for (i in 1:3) {
  anova.1 <- execute.anova.1()
  anova.2 <- execute.anova.2()
  anova.3 <- execute.anova.3()
  anova.4 <- execute.anova.4()
  anova.5 <- execute.anova.5()
  
}
print("All four Likeliness Ratio tests have run 3 more times.")


## ----bonferroni_correction------------------------------------------------------------------------------------------------------------

bonferroni_corr <- function(p_value, num_run_test) {
  
  # Apply Bonferroni correction
  adjusted.alpha <- 0.05 / num_run_test
  
  # Adjust the p-value using Bonferroni correction
  adjusted.p.value <- p.adjust(p_value, method = "bonferroni")
  
  # Compare the adjusted p-value to the adjusted alpha
  if(is.na(adjusted.p.value)) {
    # The two models are the same and have no difference
    cat("There is no significant difference between the two models, as the two models are the same.", "\nThe adjusted significance level is:", adjusted.alpha, "\nThe adjusted p-value is:", adjusted.p.value, "\n")
  }
  else if (adjusted.p.value >= adjusted.alpha) {
    # Fail to reject the null hypothesis
    cat("There is no significant difference between the two models (after Bonferroni correction).", "\nThe adjusted significance level is:", adjusted.alpha, "\nThe adjusted p-value is:", adjusted.p.value, "\n")
  } else {
    # Reject the null hypothesis and conclude significant difference
    cat("There is a significant difference between the two models (after Bonferroni correction).", "\nThe adjusted significance level is:", adjusted.alpha, "\nThe adjusted p-value is:", adjusted.p.value, "\n")

  }
}

# Execute Bonferroni correction for all four Likeliness Ratio tests
bonferroni_corr(anova.1$Pr[2], anova.runs$anova.1)
bonferroni_corr(anova.2$Pr[2], anova.runs$anova.2)
bonferroni_corr(anova.3$Pr[2], anova.runs$anova.3)
bonferroni_corr(anova.4$Pr[2], anova.runs$anova.4)
bonferroni_corr(anova.5$Pr[2], anova.runs$anova.5)

