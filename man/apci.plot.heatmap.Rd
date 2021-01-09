\name{apci.plot.heatmap}
\alias{apci.plot.heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot the heatmap for APC-I model
}
\description{
plot the heatmap to show cohort effects
}
\usage{
apci.plot.heatmap(model, age, period, color_map = NULL,
color_scale = NULL, quantile = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A list, inheriting the corresponding results generated by function apci.
}
  \item{age}{
A vector, representing the age group index taking on a small number of distinct values in the data. Usually, the vector should be converted to a factor (or the terms of "category" and "enumerated type").
}
  \item{period}{
An object of class character, similar to the argument of age, representing the time period index in the data.
}
  \item{color_map}{
A vector, representing the color palettes to be used in the figure. The default setting is greys if color_map is NULL. Alternations, for example, can be c("blue", "yellow"), blues, etc.
}
  \item{color_scale}{
A vector including two numbers indicating the limit of the values to be plotted. The first number is the minimum value to be visualized and the second is the maximum value to be visualized. If NULL, the algorithm will automatically select the limits from the data (estimation results) to set up the scale.
}
  \item{quantile}{
A number between 0 and 1, representing the percentiles to be used in visualizing the data or model. If NULL, the original magnitude will be used.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}

\examples{
library("APCI")
## load data
test_data <- APCI::women9017
test_data$acc <- as.factor(test_data$acc)
test_data$pcc <- as.factor(test_data$pcc)

## run APCI model
APC_I <- apci(outcome = "inlfc",
              age = "acc",
              period = "pcc",
              cohort = "ccc",
              weight = "wt",
              data = test_data,dev.test=FALSE,
              family = "gaussian")

## plot heatmap
apci.plot.heatmap(model = APC_I, age = "acc",period = 'pcc',
                  color_map = c('blue','yellow'))

}

