#' @rawNamespace useDynLib(tsarma, .registration=TRUE); useDynLib(tsarma_TMBExports)
#' @keywords internal
#' @importFrom TMB MakeADFun sdreport
#' @import tsmethods
#' @import data.table
#' @import methods
#' @importFrom stats na.omit printCoefmat ar coef dnorm lm logLik residuals fitted pnorm var vcov confint qnorm sigma AIC BIC nobs simulate arima ar.ols sd ARMAtoMA embed approxfun predict runif
#' @importFrom zoo na.fill coredata index is.zoo as.zoo
#' @importFrom xts xts as.xts is.xts merge.xts
#' @importFrom sandwich estfun bwNeweyWest vcovHAC vcovOPG bread
#' @importFrom numDeriv grad jacobian
#' @importFrom nloptr nloptr
#' @importFrom flextable flextable as_flextable set_caption add_footer_row add_footer_lines append_chunks as_chunk as_equation as_paragraph compose colformat_double set_header_labels padding bold align autofit hline width
#' @importFrom future.apply future_lapply
#' @importFrom future %<-%
#' @importFrom progressr handlers progressor
#' @importFrom lubridate %m+% tz days weeks years
#' @importFrom tsdistributions ddist rdist qdist pdist
#' @importFrom graphics abline box curve grid hist layout legend lines par points title
#' @importFrom Rdpack reprompt
#' @importFrom utils tail head data
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
