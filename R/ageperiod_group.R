#' Get the cohort index matrix for any age and period groups
#'
#' This function returns the cohort index matrix for any age and period
#' groups. The cohort index matrix will then be used to extract the cohort
#' effects.
#'
#' @param age_range,period_range Numeric vector indicating the actual
#' age and period range (e.g., 10 to 59 years old from 2000 to 2019).
#' @param age_interval,period_interval,age_group,period_group Numeric
#' values or character vectors indicating how age and period are
#' grouped. \code{age_interval} and \code{period_interval} are numbers
#' indicating the width of age and period groups respectively.
#' \code{age_group} and \code{period_group} are character vectors
#' explicitly listing all potential age and period groups. Either
#' \code{age_interval}(\code{period_interval}) or \code{age_group}
#' (\code{period_group}) have to be defined when \code{unequal_interval}
#' is \code{TRUE}.
#'
#' @return a matrix respresenting the relationship among age, period,
#' and cohort groups under the current setting.
#'
#' @examples
#' ## age and period groups have equal width
#' ageperiod_group(age_range = 10:59, period_range = 2000:2019,
#'                 age_interval = 5, period_interval = 5)
#' ageperiod_group(age_range = 10:59, period_range = 2000:2019,
#'                 age_group = c("10-14","15-19","20-24","25-29",
#'                               "30-34","35-39","40-44","45-49",
#'                               "50-54","55-59"),
#'                 period_group = c("2000-2004","2005-2009","2010-2014",
#'                                  "2015-2019"))
#'
#' ## age and period groups have unequal width
#' ageperiod_group(age_range = 10:59, period_range = 2000:2019,
#'                 age_interval = 10, period_interval = 5)
#' ageperiod_group(age_range = 10:59, period_range = 2000:2019,
#'                 age_group = c("10-19","20-29","30-39","40-49","50-59"),
#'                 period_group = c("2000-2004","2005-2009",
#'                                  "2010-2014","2015-2019"))
#'
#' @export




# Define age and period groups
# one drawback of this function: it does not allow for seeting up
# the initial grouping point.
# also, if unequal, the cohort group numbers will not be A + P -1
# therefore, in calculating the cohort average, this can be a problem
# since argument C should be changed

# if the users define some groups not always the equal, this function should
# allow for the customized grouping strategy to make everything more
# flexible

ageperiod_group <- function(age_range,
                            period_range,
                            age_interval=NULL,
                            period_interval=NULL,
                            age_group = NULL,
                            period_group = NULL){



  # define expanded data frame for specific age and period ranges
  df <- expand.grid(A=age_range,P=period_range)
  if(!is.null(age_group)){
    # age_group <- c("10-19","20-29","30-39","40-49","50-59")
    age_group_split <- stringr::str_extract_all(age_group, "[0-9]+",simplify = T)
    df$acc <- cut(df$A,
        breaks = as.numeric(c(age_group_split[1],age_group_split[,2])),
        right = T, include.lowest = T,
        labels = 1:length(age_group))
  }else{
  df$acc <- floor((df$A-min(age_range))/age_interval) + 1
  }

  if(!is.null(period_group)){
    # period_group <- c("2000-2004","2005-2009","2010-2014","2015-2019")
    period_group_split <- stringr::str_extract_all(period_group, "[0-9]+",simplify = T)
    df$pcc <- cut(df$P,
                  breaks = as.numeric(c(period_group_split[1],period_group_split[,2])),
                  right = T, include.lowest = T,
                  labels = 1:length(period_group))
  }else{
  df$pcc <- floor((df$P-min(period_range))/period_interval) + 1
  }
  # define cohort
  df$cohort <- df$P - df$A

  df <- dplyr::mutate(dplyr::group_by(df,acc,pcc),crange = paste(min(cohort),max(cohort)) )
  df <- dplyr::ungroup(df)
  df$ccc <- dplyr::group_indices(dplyr::group_by(df,crange))

  matrix(unique(dplyr::select(.data = df, acc,pcc,ccc))$ccc,
         nrow = length(unique(df$acc)),
         ncol = length(unique(df$pcc)))
}



