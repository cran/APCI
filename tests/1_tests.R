# install the package and use this script to test the package
rm(list=ls())
library("APCI")
test_data <- APCI::women9017
test_data$acc <- as.factor(test_data$acc)
test_data$pcc <- as.factor(test_data$pcc)
test_data$educc <- as.factor(test_data$educc)
test_data$educr <- as.factor(test_data$educr)

APC_I <- APCI::apci(outcome = "inlfc",
                    age = "acc",
                    period = "pcc",
                    cohort = "ccc",
                    weight = "wt",
                    data = test_data,dev.test=FALSE,
                    print = F,
                    family = "gaussian")

# check results
summary(APC_I)

APC_I$model
summary(APC_I$model)
APC_I$dev_global
APC_I$dev_local
APC_I$intercept
APC_I$age_effect
APC_I$period_effect
APC_I$cohort_average
APC_I$cohort_slope

apci.bar(model = APC_I, age = "acc",period = "pcc")
