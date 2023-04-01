#' ARMA Model Specification
#'
#' @description Specifies an ARMA model prior to estimation.
#' @details The specification object holds the information and data which is
#' then passed to the maximum likelihood estimation routines.
#' @param y an xts vector.
#' @param constant whether to estimate a constant (mean) for y,
#' @param order the (p,q) ARMA order.
#' @param xreg an optional xts matrix of regressors in the conditional mean
#' equation.
#' @param distribution a valid distribution from the available
#' re-parameterized distributions of the package.
#' @param ... not used.
#' @return An object of class \dQuote{tsarma.spec}.
#' @details
#' The arma specification is meant for stationary data, not integrated. The constant
#' represents the mean of y so that the equation becomes:
#' \deqn{y[t] = mu[t] + phi[1](y[t-1] - mu[t-1]) + \dots + phi[p](y[t-p]-mu[t-p]) +
#' theta[1]e[t-1] + \dots + theta[q]e[t-q] + e[t]}
#' where \deqn{mu[t]} is the estimated mean of y plus any additional pre-lagged
#' regressors passed via the xreg argument. If the absence of regressors,
#' then \deqn{mu[t]} is constant across all time indices.
#' Initialization of the parameters is based on Hannan and Rissanen (1982).
#' @aliases arma_modelspec
#' @rdname arma_modelspec
#' @author Alexios Galanos
#' @export
#'
#'
arma_modelspec <- function(y, order = c(0,0), constant = TRUE, xreg = NULL, distribution = "norm", ...)
{
    if  (!is.xts(y)) {
        stop("y must be an xts object")
    }
    if (any(is.na(y))) stop("\ntsarma does not currently support missing (NA) values. Try tsissm instead.")
    spec <- initialize_data(y)
    order <- c(abs(as.integer(order[1])),abs(as.integer(order[2])))
    distribution <- match.arg(distribution[1], choices = valid_distributions())
    parmatrix <- .parameters_arma(y, constant = constant, order = order, xreg = xreg, distribution = distribution)
    cmodel <- c(max(order), order[1], order[2], distribution_class(distribution))
    spec$model$model <- "arma"
    spec$model$order <- order
    if (is.null(xreg)) {
        spec$xreg$xreg <- matrix(0, ncol = 1, nrow = NROW(y))
        spec$xreg$include_xreg <- FALSE
    } else {
        spec$xreg$xreg <- coredata(xreg)
        spec$xreg$include_xreg <- TRUE
    }
    spec$distribution <- distribution
    spec$parmatrix <- parmatrix
    spec$model_options <- cmodel
    spec$model$constant <- constant
    # initialize parameters
    class(spec) <- "tsarma.spec"
    return(spec)
}
