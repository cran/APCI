\name{temp_model}
\alias{temp_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate APC-I model
}
\description{
Estimate the APCI original model. This is a generalized linear regression model.
}
\usage{
temp_model(data, outcome = "inlfc", age = "acc",
period = "pcc", cohort = NULL, weight = NULL,
covariate = NULL, family = "quasibinomial", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing the outcome variable, age group indicator, period group indicator, and covariates to be used in the model. If the variable(s) are not found in data, there will be an error message reminding the users to check the input data again.
}
  \item{outcome}{
An object of class character containing the name of the outcome variable. The outcome variable can be continuous, categorical, or count.
}
  \item{age}{
An object of class character representing the age group index taking on a small number of distinct values in the data. Usually, the vector should be converted to a factor (or the terms of "category" and "enumerated type").
}
  \item{period}{
An object of class character, similar to the argument of age, representing the time period index in the data.
}
  \item{cohort}{
An optional object of class character representing cohort membership index in the data. Usually, the cohort index can be generated from the age group index and time period index in the data because of the intrinsic relationship among these three time-related indices.
}
  \item{weight}{
An optional vector of sample weights to be used in the model fitting process. If non-NULL, the weights will be used in the first step to estimate the model. Observations with negative weights will be automatically dropped in modeling.
}
  \item{covariate}{
An optional vector of characters, representing the name(s) of the user-specified covariate(s) to be used in the model. If the variable(s) are not found in data, there will be an error message reminding the users to check the data again.
}
  \item{family}{
Used to specify the statistical distribution of the error term and link function to be used in the model. Usually, it is a character string naming a family function. For example, family can be "binomial", "multinomial"", or "gaussian". Users could also check R package glm for more details of family functions.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
