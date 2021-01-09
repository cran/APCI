\name{cohortdeviation}
\alias{cohortdeviation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate cohort deviation
}
\description{
calculate cohort deviation
}
\usage{
cohortdeviation(A, P, C, model = temp6,
weight = "wt", covariate, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A, P, C}{
The numbers of age groups, period groups, and cohort groups separately.
}
  \item{model}{
A generalized linear regression model generated from the internal function temp_model
}
  \item{weight}{
An optional vector of sample weights to be used in the model fitting process. If non-NULL, the weights will be used in the first step to estimate the model. Observations with negative weights will be automatically dropped in modeling.
}
  \item{covariate}{
An optional vector of characters, representing the name(s) of the user-specified covariate(s) to be used in the model. If the variable(s) are not found in data, there will be an error message reminding the users to check the data again.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
