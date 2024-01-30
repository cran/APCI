#' Estimate APC-I model
#'
#' Estimate the APCI original model. This is a generalized linear regression model.
#'
#' @param data A data frame containing the outcome variable, age group
#' indicator, period group indicator, and covariates to be used in the model.
#' If the variable(s) are not found in data, there will be an error message
#' reminding the users to check the input data again.
#' @param outcome An object of class character containing the name of the
#' outcome variable. The outcome variable can be continuous, categorical,
#' or count.
#' @param age An object of class character representing the age group index
#' taking on a small number of distinct values in the data. Usually, the vector
#' should be converted to a factor (or the terms of "category" and "enumerated
#' type").
#' @param period An object of class character, similar to the argument of age,
#' representing the time period index in the data.
#' @param cohort An optional object of class character representing cohort
#' membership index in the data. Usually, the cohort index can be generated
#' from the age group index and time period index in the data because of the
#' intrinsic relationship among these three time-related indices.
#' @param weight An optional vector of sample weights to be used in the model
#' fitting process. If non-NULL, the weights will be used in the first step to
#' estimate the model. Observations with negative weights will be automatically
#' dropped in modeling.
#' @param covariate An optional vector of characters, representing the name(s)
#' of the user-specified covariate(s) to be used in the model. If the
#' variable(s) are not found in data, there will be an error message reminding
#' the users to check the data again.
#' @param family Used to specify the statistical distribution of the error
#' term and link function to be used in the model. Usually, it is a character
#' string naming a family function. For example, family can be "binomial",
#' "multinomial"", or "gaussian". Users could also check R package glm for
#' more details of family functions.
#' @param gee Logical, indicating if the data is cross-sectional data or
#' longitudinal/panel data. If \code{TRUE}, the generalized estimating equation
#' will be used to correct the standard error estimates. The default is
#' \code{FALSE}, indicating that the data are cross-sectional.
#' @param id A vector of character, specifying the cluster index in longitudinal
#' data. It is required when \code{gee} is \code{TRUE}. The length of the vector
#' should be the same as the number of observations.
#' @param corstr A character string, specifying a possible correlation
#' structure in the error terms when \code{gee} is \code{TRUE}. The following
#' are allowed: \code{independence}, \code{fixed}, \code{stat\_M\_dep},
#' \code{non\_stat\_M\_dep}, \code{exchangeable}, \code{AR-M} and
#' \code{unstructured}. The default value is \code{exchangeable}.
#' @param \dots Additional arguments to be passed to the function.
#'
#' @return A list containing:
#' \item{A}{Age group index.}
#' \item{P}{Period group index.}
#' \item{C}{Cohort group index.}
#' \item{model}{Fitted APCI models of outcome on predictors.}
#'
#' @export


# load packages: if packages is not successfully loaded, install the corresponding
# packages
# pkgs <- c("haven","survey","tidyverse","magrittr",'data.table')
# installpackages <- lapply(pkgs,function(x){
#   if(x %in% rownames(installed.packages()) == FALSE) {install.packages(`x`)}
# })
# loadpackages <- lapply(pkgs,function(x){
#   library(`x`,character.only = T)
# })
# rm(list = c("installpackages","loadpackages","pkgs"))



# change variable names ()

# get temp6
temp_model <- function(data,
                       outcome = "inlfc",
                       age = "acc",
                       period = "pcc",
                       cohort = NULL,
                       weight = NULL,
                       covariate = NULL,
                       family = "quasibinomial",
                       gee = FALSE,
                       id = NULL,
                       corstr = "exchangeable",
                       ...){
  ###############
  #no missing data
  ###############
  data2 <- as.data.frame(data)

  if(family=="binomial"&gee==FALSE){
    family <- "quasibinomial"
  }

  # data2 = na.omit( data[ , unique(c(outcome,age,period,cohort,weight,covariate,id)[!is.null(c(outcome,age,period,cohort,weight,covariate,id))]) ] )

  if(is.null(weight)){
    weight <- 1
  }

  # A = nlevels(data2[,age])
  # P = nlevels(data2[,period])
  A = length(unique(data2[,age]))
  P = length(unique(data2[,period]))
  C = A + P - 1

  data2$acc <- data2[,age]
  data2$pcc <- data2[,period]

  data2$acc <- as.factor(data2$acc)
  data2$pcc <- as.factor(data2$pcc)

  # if(age!='acc'&length(str_which(covariate,age))==0 ){
  # data2[,age] <- NULL
  # }
  # if(period!='pcc'&length(str_which(covariate,period))==0){
  # data2[,period] <- NULL
  # }
  #
  # if(length(str_which(covariate,age))>0){
  # covariate <- covariate[-str_which(covariate,age)]
  # }
  # if(length(str_which(covariate,period))>0){
  # covariate <- covariate[-str_which(covariate,period)]
  # }

  options(contrasts=c("contr.sum","contr.poly"), na.action = na.omit)
  wtdata2 = survey::svydesign(id=~1, strata=NULL, data=data2,
                              weights=as.formula(paste0("~",weight)))

  if(!is.null(covariate)&length(covariate)>0){
    temp6_formula <- as.formula(paste0(outcome,"~",paste(c(paste(covariate,collapse = "+"),
                                                           paste0(age,"*",period)),collapse = "+")))
  }else{
    temp6_formula <- as.formula(paste0(outcome,"~",paste(c(paste0(age,"*",period)),collapse = "+")))
  }

  print(temp6_formula)

  if(!is.null(covariate)&length(covariate)>0){
    temp6_formula <- as.formula(paste0(outcome,"~",paste(c(paste(covariate,collapse = "+"),
                                                           paste0('acc',"*",'pcc')),collapse = "+")))
  }else{
    temp6_formula <- as.formula(paste0(outcome,"~",paste(c(paste0('acc',"*",'pcc')),collapse = "+")))
  }

  if(gee==TRUE){
    temp6 = gee::gee(temp6_formula,
                     id = id,
                     # id = get(id),
                     # id = eval(parse(text = id)),
                     data = data2,
                     family = get(family),
                     corstr = corstr)
    temp6$model <- data2

  }else{
  temp6 = survey::svyglm(temp6_formula,
                         wtdata2,
                         family = get(family))
  if(temp6$df.residual == 0){
    temp6 = glm(temp6_formula,
                 data2,
                 family = get(family))
  }
  }

  # temp6
  # output
  list(A=A,P=P,C=C,model=temp6)

}

