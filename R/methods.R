#' Estimates an ARMA model given a specification object using maximum likelihood and autodiff
#'
#' @param object an object of class tsarma.spec.
#' @param control solver control parameters.
#' @param solver only \dQuote{nloptr} is currently supported (see \code{\link[nloptr]{nloptr}}).
#' @param ... not currently used.
#' @return An object of class \dQuote{tsarma.estimate}.
#' @details The underlying code is written using the TMB framework which uses
#' automatic differentiation and hence allows the generation of analytic
#' derivatives.
#' Stationarity of the process is controlled through a set of inequality constraints
#' on the characteristic roots of the AR/MA
#' The estimation makes 2 passes to the solver. The first pass uses no parameter
#' scaling, whilst in the second pass the parameters (as well as bounds) are scaled
#' making use of the estimated hessian from the first pass in order to generate
#' a hopefully more robust solution.
#' @export estimate.tsarma.spec
#' @aliases estimate
#' @rdname estimate
#' @author Alexios Galanos
#' @export
#'
estimate.tsarma.spec <- function(object, solver = "nloptr", control = NULL, ...)
{
    if (is.null(control)) control <- nloptr_fast_options(trace = FALSE)
    estimate_pars <- sum(object$parmatrix$estimate)
    if (estimate_pars == 0) {
        warnings("\nall parameters are fixed...returning filtered object instead.")
        out <- tsfilter(object)
    } else {
        start_timer <- Sys.time()
        out <- .estimate_arma(object, solver, control, ...)
        end_timer <- Sys.time()
        elapsed <- end_timer - start_timer
        out$elapsed <- elapsed
    }
    return(out)
}

#' Extract Model Coefficients
#'
#' @description Extract the estimated coefficients of a model.
#' @param object an object of class \dQuote{tsarma.estimate}.
#' @param ... not currently used.
#' @return A numeric named vector of estimated coefficients.
#' @aliases coef
#' @method coef tsarma.estimate
#' @rdname coef
#' @export
#'
#'
coef.tsarma.estimate <- function(object, ...)
{
    out <- object$parmatrix[estimate == 1]$value
    names(out) <- object$parmatrix[estimate == 1]$parameter
    return(out)
}

#' Extract Model Standard Deviation
#'
#' @description Extract the standard deviation from a ARMA model.
#' @param object an object of class \dQuote{tsarma.estimate}.
#' @param ... not currently used.
#' @return An numeric vector of length 1.
#' @aliases sigma
#' @method sigma tsarma.estimate
#' @rdname sigma
#' @export
#'
#'
sigma.tsarma.estimate <- function(object, ...)
{
    parameter <- NULL
    return(object$parmatrix[parameter == "sigma"]$value)
}

#' Extract Model Fitted Values
#'
#' @description Extract the fitted values of the estimated model.
#' @param object an object of class \dQuote{tsarma.estimate}.
#' @param ... not currently used.
#' @return An xts vector of the fitted values. Since only a constant is supported
#' in the conditional mean equation this is either a vector with a constant else
#' a vector with zeros.
#' @aliases fitted
#' @method fitted tsarma.estimate
#' @rdname fitted
#' @export
#'
#'
fitted.tsarma.estimate <- function(object, ...)
{
    idx <- object$spec$target$index
    f <- xts(object$fitted, idx)
    return(f)
}

#' Extract Model Residuals
#'
#' @description Extract the residuals of the estimated model.
#' @param object an object of class \dQuote{tsarma.estimate}.
#' @param standardize logical. Whether to standardize the residuals by the
#' conditional volatility.
#' @param ... not currently used.
#' @return An xts vector of the residuals. If the model had no constant in
#' the conditional mean equation then this just returns the original data (which
#' is assumed to be zero mean noise).
#' @aliases residuals
#' @method residuals tsarma.estimate
#' @rdname residuals
#' @export
#'
#'
residuals.tsarma.estimate <- function(object, standardize = FALSE, ...)
{
    res <- xts(object$residuals, object$spec$target$index)
    if (standardize) {
        res <- res/sigma(object)
    }
    return(res)
}

#' The Covariance Matrix of the Estimated Parameters
#'
#' @param object an object of class tsarma.estimate.
#' @param adjust logical. Should a finite sample adjustment be made? This amounts
#' to multiplication with n/(n-k) where n is the number of observations and k
#' the number of estimated parameters.
#' @param type valid choices are \dQuote{H} for using the analytic hessian
#' for the bread, \dQuote{OP} for the outer product of gradients, \dQuote{QMLE}
#' for the Quasi-ML sandwich estimator (Huber-White), and \dQuote{NW} for the Newey-West
#' adjusted sandwich estimator (a HAC estimator).
#' @param ... additional parameters passed to the Newey-West bandwidth function to
#' determine the optimal lags.
#' @return The variance-covariance matrix of the estimated parameters.
#' @method vcov tsarma.estimate
#' @aliases vcov
#' @rdname vcov
#' @export
#'
vcov.tsarma.estimate <- function(object, adjust = FALSE, type = c("H","OP","QMLE","NW"), ...)
{
    type <- match.arg(type[1],c("H","OP","QMLE","NW"))
    N <- nrow(estfun(object))
    if (type == "H") {
        V <- solve(bread(object))
    } else if (type == "QMLE") {
        bread. <- bread(object)
        meat. <- meat_tsarma(object, adjust = adjust)
        bread. <- solve(bread.)
        V <- bread. %*% meat. %*% bread.
    } else if (type == "OP") {
        V <- vcovOPG(object, adjust = adjust)
    } else if (type == "NW") {
        bread. <- bread(object)
        meat. <- meatHAC_tsarma(object, adjust = adjust, ...)
        bread. <- solve(bread.)
        V <- bread. %*% meat. %*% bread.
    }
    par_names <- object$parmatrix[estimate == 1]$parameter
    colnames(V) <- rownames(V) <- par_names
    return(V)
}

#' Confidence Intervals for Model Parameters
#'
#' @param object an object of class tsarma.estimate.
#' @param parm a specification of which parameters are to be given confidence intervals,
#' either a vector of numbers or a vector of names. If missing, all parameters
#' are considered.
#' @param level the confidence level required.
#' @param vcov_type valid choices are \dQuote{H} for using the analytic hessian
#' for the bread, \dQuote{OP} for the outer product of gradients, \dQuote{QMLE}
#' for the Quasi-ML sandwich estimator (Huber-White), and \dQuote{NW} for the Newey-West
#' adjusted sandwich estimator (a HAC estimator).
#' @param ... additional parameters passed to the Newey-West bandwidth function to
#' determine the optimal lags.
#' @return A matrix (or vector) with columns giving lower and upper confidence
#' limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2
#' in % (by default 2.5% and 97.5%).
#' @method confint tsarma.estimate
#' @aliases confint
#' @rdname confint
#' @export
#'
confint.tsarma.estimate <- function(object, parm, level = 0.95, vcov_type = "H", ...)
{
    # extract the names of the estimated parameters
    par_names <- object$parmatrix[estimate == 1]$parameter
    if (length(par_names) == 0) return(NULL)
    coefficients <- coef(object)
    if (missing(parm)) {
        parm <- par_names
    } else if (is.numeric(parm)) {
        parm <- par_names[parm]
    } else if (is.character(parm)) {
        parm <- par_names[which(par_names %in% parm)]
    }
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
    vc <- vcov(object, type = vcov_type)
    ses <- sqrt(diag(vcov(object, type = vcov_type)))
    ses <- ses[parm]
    ci[] <- coefficients[parm] + ses %o% fac
    return(ci)
}


#' Extract Log-Likelihood
#'
#' @description Extract the log likelihood of the model at the estimated optimal
#' parameter values.
#' @param object an object of class \dQuote{tsarma.estimate}.
#' @param ... not currently used.
#' @return An object of class \dQuote{logLik} with attributes for \dQuote{nobs} and
#' \dQuote{df}. The latter is equal to the number of estimated parameters
#' plus sum(ARMA order) (the initialization values).
#' @aliases logLik
#' @method logLik tsarma.estimate
#' @rdname logLik
#' @export
#'
#'
logLik.tsarma.estimate <- function(object, ...)
{
    out <- -1.0 * object$loglik
    attr(out,"nobs") <- object$nobs
    attr(out,"df") <- object$npars
    class(out) <- "logLik"
    return(out)
}


#' ARMA Model Estimation Summary
#'
#' @description Summary method for class \dQuote{tsarma.estimate}
#' @param object an object of class \dQuote{tsarma.estimate}.
#' @param vcov_type the type of standard errors based on the vcov estimate (see \code{\link{vcov}}).
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param ... not currently used.
#' @return A list with summary information of class \dQuote{tsarma.summary}.
#' @aliases summary
#' @method summary tsarma.estimate
#' @rdname summary
#' @export
#'
#'
summary.tsarma.estimate <- function(object, digits = 4, vcov_type = "H", ...)
{
    estimate_pars <- sum(object$parmatrix$estimate)
    if (estimate_pars == 0) return(NULL)
    V <- vcov(object, type = vcov_type)
    est <- object$parmatrix[estimate == 1]$value
    par_names <- object$parmatrix[estimate == 1]$parameters
    se <- sqrt(diag(V))
    tval <- est/se
    coefficients <- cbind(Estimate = est, `Std. Error` = se,`t value` = tval, `Pr(>|t|)` = 2*(1 - pnorm(abs(tval))))
    n_obs <- nobs(object)
    n_parameters <- length(coef(object))
    llh <- -object$loglik
    elapsed <- object$elapsed
    conditions <- object$conditions[c("kkt1","kkt2","evratio")]
    distribution <- object$spec$distribution
    equation <- tsequation(object)
    coefficients <- as.data.table(coefficients, keep.rownames = TRUE)
    setnames(coefficients, "rn","term")
    syms <- object$parmatrix[estimate == 1]$symbol
    out <- list(coefficients = coefficients, distribution = distribution,
                loglikelihood = llh, n_obs = n_obs, n_parameters = n_parameters,
                df = object$df,
                AIC = AIC(object),
                BIC = BIC(object),
                elapsed = elapsed, conditions = conditions, equation = equation,
                model = object$spec$model$model, symbol = syms,
                order = object$spec$model$order,
                equation = object$parmatrix[estimate == 1]$equation)
    class(out) <- "summary.tsarma"
    return(out)
}

#' Model Estimation Summary Print method
#'
#' @description Print method for class \dQuote{summary.tsarma}
#' @param x an object of class \dQuote{summary.tsarma}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param signif.stars logical. If TRUE, ‘significance stars’ are printed for each coefficient.
#' @param include.symbols logical. If TRUE, replaces parameter names with their symbols (if they exist).
#' @param include.equation logical. If TRUE, adds a section with the symbolic model equation.
#' @param include.statistics logical. If TRUE, adds a section with summary statistics on the model.
#' @param table.caption an optional string for the table caption.
#' @param format either prints to \dQuote{console} or prints and returns a \dQuote{flextable} object.
#' @param ... additional arguments passed to flextable print method.
#' @return Invisibly returns a flextable object if output was set to \dQuote{flextable} else
#' the original summary object if output was set to \dQuote{console}.
#' @aliases print.summary.tsarma
#' @method print summary.tsarma
#' @rdname print
#' @export
#'
#'
print.summary.tsarma <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"),
                                  include.symbols = TRUE, include.equation = TRUE,
                                  include.statistics = TRUE, table.caption = paste0(toupper(x$model)," Model Summary"),
                                  format = c("console","flextable"), ...)
{
    format <- match.arg(format[1], c("console","flextable"))
    if (format == "console") {
        .print_screen(x, digits = digits, signif.stars = signif.stars, table.caption = table.caption, ...)
    } else {
        out <- .print_flextable(x, digits = digits, signif.stars = signif.stars,
                                include.symbols = include.symbols, include.equation = include.equation,
                                include.statistics = include.statistics,
                                table.caption = table.caption, ...)
        print(out)
        return(invisible(out))
    }
}


#' Model Equation (LaTeX)
#'
#' @description Generates a list of model equations in LaTeX.
#' @param object an object of class \dQuote{tsarma.estimate}.
#' @param ... not currently used.
#' @return A list of equations in LaTeX which can be used in documents.
#' @details This method is called in the summary when the format output option
#' chosen is \dQuote{flextable}.
#' @aliases tsequation
#' @method tsequation tsarma.estimate
#' @rdname tsequation
#' @export
#'
tsequation.tsarma.estimate <- function(object, ...)
{
    if (object$spec$xreg$include_xreg) {
        x <- extract_model_values(object, object_type = "estimate", value_name = "xreg")
    } else {
        x <- NULL
    }
    out <- .equation_arma(object$spec$model$order, xreg = x, distribution = object$spec$distribution)
    return(out)
}


#' Unconditional Value
#'
#' @description General method the unconditional value of a model.
#' @param object an object.
#' @param ... additional parameters passed to the method.
#' @return A scalar of the unconditional value. For ARMA models this represents
#' the unconditional mean.
#' @aliases unconditional
#' @method unconditional tsarma.estimate
#' @rdname unconditional
#' @export
#'
#
unconditional.tsarma.estimate <- function(object, ...)
{
    x <- extract_model_values(object, object_type = "estimate", value_name = "xreg")
    xi <- extract_model_values(object, object_type = "estimate", value_name = "xi")
    mu <- extract_model_values(object, object_type = "estimate", value_name = "mu")
    out <- mu + sum(xi * colMeans(x))
    return(out)
}

#' Akaike's An Information Criterion
#'
#' @description Extract the AIC from an estimated model.
#' @param object an object of class \dQuote{tsarma.estimate}.
#' @param ... not currently used.
#' @param k the penalty per parameter to be used; the default k = 2 is the
#' classical AIC.
#' @return A numeric value.
#' @aliases AIC
#' @method AIC tsarma.estimate
#' @rdname AIC
#' @export
#'
#'
AIC.tsarma.estimate <- function(object, ..., k = 2)
{
    out <- -2.0 * as.numeric(logLik(object)) + k * object$npars
    return(out)
}

#' Bayesian Information Criterion
#'
#' @description Extract the BIC from an estimated model.
#' @param object an object of class \dQuote{tsarma.estimate}.
#' @param ... not currently used.
#' @return A numeric value.
#' @aliases BIC
#' @method BIC tsarma.estimate
#' @rdname BIC
#' @export
#'
#'
BIC.tsarma.estimate <- function(object, ...)
{
    out <- -2 * as.numeric(logLik(object)) + object$npars * log(nobs(object))
    return(out)
}


#' Extract the number of observations from an estimated model
#'
#' @description Extract the number of ‘observations’ from an estimated model.
#' This is principally intended to be used in computing BIC and used in other
#' tidy methods
#' @param object an object of class \dQuote{tsarma.estimate}.
#' @param ... not currently used.
#' @return A numeric value.
#' @aliases nobs
#' @method nobs tsarma.estimate
#' @rdname nobs
#' @export
#'
#'
nobs.tsarma.estimate <- function(object, ...)
{
    return(object$nobs)
}

#' Model Filtering
#'
#' @description Filters new data based on an already estimated model.
#' @param object an object of class \dQuote{tsarma.estimate} or \dQuote{tsarna.spec}.
#' @param y an xts vector of new values to filter.
#' @param newxreg model regressors with the same number of rows as y. This can be either
#' a numeric or xts matrix. Only needed if the model was estimated with regressors in the
#' conditional mean equation.
#' @param ... additional arguments for future expansion options.
#' @return A \dQuote{tsarma.estimate} object with updated information if the input
#' object was also of class \dQuote{tsarma.estimate}, else will return an object
#' of class \dQuote{tsarma.filter} for input of class \dQuote{tsarma.spec}.
#' @details The method filters new data and updates the object with this new information,
#' appending y and xreg datasets so that is can be called recursively as new
#' data arrives. For a \dQuote{tsarma.spec} object as input, the existing data
#' will be filtered with the parameter in the specification (considered fixed), but
#' the returned object will not contain information about parameter uncertainty
#' such as the hessian or scores.
#' @aliases tsfilter
#' @method tsfilter tsarma.estimate
#' @rdname tsfilter
#' @export
#'
#'
tsfilter.tsarma.estimate <- function(object, y = NULL, newxreg = NULL, ...)
{
    out <- .filter_arma_estimate(object, y = y, newxreg = newxreg, ...)
    return(out)
}

#' @aliases tsfilter
#' @method tsfilter tsarma.spec
#' @rdname tsfilter
#' @export
#'
tsfilter.tsarma.spec <- function(object, y = NULL, newxreg = NULL, ...)
{
    out <- .filter_arma_spec(object, y = y, newxreg = newxreg, ...)
    return(out)
}


#' Model Prediction
#'
#' @description Prediction function for class \dQuote{tsarma.estimate}.
#' @param object an object of class \dQuote{tsarma.estimate}.
#' @param h the forecast horizon.
#' @param newxreg regressors rows equal to h. This can be either
#' a numeric or xts matrix. Only needed if the model was estimated with regressors in the
#' conditional mean equation.
#' @param bootstrap whether to generate a bootstrap distribution for the conditional
#' volatility using re-sampling of the empirical innovations as in the paper by Pascual et al (2006).
#' @param nsim the number of simulations to use for generating the simulated
#' predictive distribution.
#' @param innov an optional matrix of dimensions nsim by h of innovations (see
#' innov_type below for types supported).
#' @param innov_type there are 3 options available. If type is \dQuote{q}, then
#' these represent quantiles and will be transformed to the appropriate distribution
#' innovations used in the model. If type is \dQuote{z}, these represent standardized
#' innovations which will be scaled by the model standard deviation. Finally, if
#' type is \dQuote{r}, then these are zero mean non-scaled innovations for which
#' the sigma has already been infused into the distribution (i.e. either static
#' or from a GARCH model output).
#' @param innov_init an optional vector of initialization values for the
#' standardized innovations. This allows the seeding of the initial innovations
#' with user supplied values (useful when simulating forward from an existing
#' model for the purpose of continuing the modeled series from some fixed point).
#' These should be of the same type as defined by innov_type.
#' @param forc_dates an optional vector of forecast dates equal to h. If NULL will use the
#' implied periodicity of the data to generate a regular sequence of dates after the
#' last available date in the data.
#' @param series_init an optional vector of values (y) to initialize the forecast.
#' If NULL, will use the last available values from the model. This must
#' be equal to the max of the ARMA order.
#' @param ... additional arguments for future expansion options.
#' @return A \dQuote{tsarma.predict} object.
#' @details The bootstrap method considered here, is based on re-sampling innovations
#' from the empirical distribution of the fitted ARMA model to generate future
#' realizations of the series.
#' @aliases predict
#' @method predict tsarma.estimate
#' @rdname predict
#' @export
#'
#'
predict.tsarma.estimate <- function(object, h = 1, newxreg = NULL, bootstrap = FALSE, nsim = 1000, innov = NULL, innov_type = "r", innov_init = NULL, forc_dates = NULL, series_init = NULL, ...)
{
    model <- object$spec$model$model
    p <- .predict_arma(object = object, h = h, newxreg = newxreg, bootstrap = bootstrap, nsim = nsim, innov = innov, innov_type = innov_type,
                       innov_init = innov_init, forc_dates = forc_dates, series_init = series_init)
    return(p)
}

#' Probability Integral Transform (PIT)
#'
#' @description Calculates and returns the conditional probability integral
#' transform given the data and estimated density
#' @param object an object.
#' @param ... additional parameters passed to the method.
#' @return An xts vector of the conditional probabilities.
#' @aliases pit
#' @method pit tsarma.estimate
#' @rdname pit
#' @export
#'
#
pit.tsarma.estimate <- function(object, ...)
{
    parameter <- NULL
    dist <- extract_model_values(object, object_type = "estimate", value_name = "distribution")
    distribution <- object$spec$distribution
    sigma <- extract_model_values(object, object_type = "estimate", value_name = "sigma")
    mu <- extract_model_values(object, object_type = "estimate", value_name = "mu")
    mu <- rep(mu, length(sigma))
    r <- as.numeric(extract_model_values(object, object_type = "estimate", value_name = "y"))
    p <- pdist(distribution, q = r, mu = mu, sigma = sigma, skew = dist[1], shape = dist[2], lambda = dist[3])
    p <- xts(p, object$spec$target$index)
    return(p)
}

#' Half Life
#'
#' @description Calculates and returns the half-life of an ARMA model.
#' @param object an object of class \dQuote{tsarma.estimate}.
#' @param ... not currently used.
#' @details The half life is defined as the period it
#' takes a series  to reach half its long-term average values. While this is know
#' for an AR(1) model, it is less clear what this should be for higher orders of
#' ARMA(p,q) models. Therefore, the function first converts the the ARMA model into
#' an infinite (truncated) MA model and then uses a spline approximation to find
#' the value for which the psi-weights (MA) equal 0.5.
#' @aliases halflife
#' @method halflife tsarma.estimate
#' @rdname halflife
#' @export
#'
#
halflife.tsarma.estimate <- function(object, ...)
{
    phi <- extract_model_values(object, object_type = "estimate", value_name = "phi")
    theta <- extract_model_values(object, object_type = "estimate", value_name = "theta")
    psi <- ARMAtoMA(phi, theta, lag.max = 100)
    if (tail(psi,1) > 0.5) {
        psi <- ARMAtoMA(phi, theta, lag.max = 500)
    }
    psi <- abs(psi)
    f <- approxfun(y = 0:length(psi), x = c(1, psi), method = "linear")
    out <- f(0.5)
    return(out)
}


#' Model Simulation
#'
#' @description Simulates paths of a ARMA model.
#' @param object an object of class \dQuote{tsarma.spec} or \dQuote{tsarma.estimate}.
#' @param h the number of time steps to simulate paths for.
#' @param nsim the number of sample paths to generate.
#' @param seed an integer that will be used in a call to set.seed before simulating.
#' @param series_init the seed value for initializing the mean equation recursion.
#' If NULL, the unconditional mean is used based on the supplied parameters.
#' This should be a vector and assumes all sample paths are seeded the same way.
#' @param innov an optional matrix of dimensions nsim by h of innovations (see
#' innov_type below for types supported).
#' @param innov_type there are 3 options available. If type is \dQuote{q}, then
#' these represent quantiles and will be transformed to the appropriate distribution
#' innovations used in the model. If type is \dQuote{z}, these represent standardized
#' innovations which will be scaled by the model standard deviation. Finally, if
#' type is \dQuote{r}, then these are zero mean non-scaled innovations for which
#' the sigma has already been infused into the distribution (i.e. either static
#' or from a GARCH model output).
#' @param innov_init an optional vector of initialization values for the
#' standardized innovations. This allows the seeding of the initial innovations
#' with user supplied values (useful when simulating forward from an existing
#' model for the purpose of continuing the modeled series from some fixed point).
#' These should be of the same type as defined by innov_type.
#' @param xreg an optional vector of length h representing any mean regressors
#' to use in the simulation.
#' @param ... not currently used.
#' @return An object of class \dQuote{tsarma.simulate} with slots for the
#' simulated series and innovations, both of which are of class
#' \dQuote{tsmodel.distribution}.
#' @aliases simulate
#' @method simulate tsarma.spec
#' @rdname simulate
#' @export
#'
#'
simulate.tsarma.spec <- function(object, nsim = 1000, seed = NULL, h = 1, innov = NULL, innov_init = NULL, innov_type = "r", series_init = NULL, xreg = NULL, ...)
{
    out <- .simulate_arma(object, h = h, seed = seed, nsim = nsim, innov = innov, innov_init = innov_init, series_init = series_init, xreg = xreg, type = "spec")
    class(out) <- "tsarma.simulate"
    return(out)
}

#' @aliases simulate
#' @method simulate tsarma.estimate
#' @rdname simulate
#' @export
#'
#'
simulate.tsarma.estimate <- function(object, nsim = 1000, seed = NULL, h = 1, innov = NULL, innov_init = NULL, series_init = NULL, xreg = NULL, ...)
{
    out <- .simulate_arma(object, h = h, nsim = nsim, seed = seed, innov = innov, innov_init = innov_init, series_init = series_init, xreg = xreg, type = "estimate")
    class(out) <- "tsarma.simulate"
    return(out)
}


#' ARMA Model Impulse Response Function
#'
#' @description Calculates the impulse response function of an ARMA model.
#' @param phi the AR coefficients in the difference equation notation.
#' @param theta the MA coefficients in the difference equation notation.
#' @param h the forecast horizon. This may need to be tweaked to ensure that
#' the output is long enough to show the decay to zero. If only an MA model is
#' passed (phi is zero length), then h should be set to the MA order.
#' @param sigma the innovations standard deviation (see details).
#' @details The impulse response function (irf) of a time series model
#' measures the changes in the future responses of all variables in the system
#' when a variable is shocked by an impulse. For an ARMA model, this simply
#' transforms the ARMA part to an infinite MA representation. Finally, the responses
#' are scaled by the lower triangular factor in the Cholesky factorization of the
#' innovation variance (which is the standard deviation in this case).
#' For standardized shocks, where the variance is 1, the scaling is unity.
#' @return A data.table with the horizon and responses, starting at period 0.
#' @aliases arma_irf
#' @rdname arma_irf
#' @export
#'
#'
arma_irf <- function(phi = numeric(), theta = numeric(), sigma = 1, h = 40)
{
    ma_rep <- c(1, ARMAtoMA(ar = phi, ma = theta, lag.max = h))
    v <- as.numeric(chol(sigma^2))
    ma_rep <- v * ma_rep
    tab <- data.table(horizon = seq(0, length(ma_rep) - 1, by = 1), response = ma_rep)
    out <- list(table = tab, phi = phi, theta = theta, sigma = sigma)
    class(out) <- "arma.irf"
    return(out)
}


arma_roots <- function(phi = numeric(), theta = numeric())
{
    if (length(phi) == 0 & length(theta) == 0) return(NULL)
    radian <- 57.29577951
    roots_ar <- re_roots_ar <- im_roots_ar <- amp_ar <- atan_ar <- degree_ar <- NULL
    roots_ma <- re_roots_ma <- im_roots_ma <- amp_ma <- atan_ma <- degree_ma <- NULL
    ar_table <- ma_table <- NULL
    if (length(phi) > 0) {
        n_ar <- length(phi)
        roots_ar <- polyroot(c(1, -phi))
        re_roots_ar <- Re(roots_ar)
        im_roots_ar <- Im(roots_ar)
        amp_ar <- apply(cbind(re_roots_ar, im_roots_ar), 1, FUN = function(x) sqrt(x[1]^2 + x[2]^2))
        atan_ar <- apply(cbind(amp_ar, re_roots_ar), 1 , FUN = function(x) atan2(x[1], x[2]))
        degree_ar <-  atan_ar * radian
        ar_table <- data.table(parameter = paste0("phi_",1:n_ar), equation = "AR", root = roots_ar, amplitude = amp_ar, atan = atan_ar, degree = degree_ar)
    }
    if (length(theta) > 0) {
        n_ma <- length(theta)
        roots_ma <- polyroot(c(1, theta))
        re_roots_ma <- Re(roots_ma)
        im_roots_ma <- Im(roots_ma)
        amp_ma <- apply(cbind(re_roots_ma, im_roots_ma), 1, FUN = function(x) sqrt(x[1]^2 + x[2]^2))
        atan_ma <- apply(cbind(amp_ma, re_roots_ma), 1 , FUN = function(x) atan2(x[1], x[2]))
        degree_ma <-  atan_ma * radian
        ma_table <- data.table(parameter = paste0("theta_",1:n_ma), equation = "MA", root = roots_ma, amplitude = amp_ma, atan = atan_ma, degree = degree_ma)
    }
    out <- list(table = rbind(ar_table, ma_table), phi = phi, theta = theta)
    class(out) <- "arma.roots"
    return(out)
}

#' Estimated Model Plots
#'
#' @description Plot method for \dQuote{tsarma.estimate} class.
#' @param x an object of class \dQuote{tsarma.estimate}.
#' @param y not used.
#' @param ... not used.
#' @method plot tsarma.estimate
#' @rdname plot.tsarma.estimate
#' @export
#'
#
plot.tsarma.estimate <- function(x, y = NULL, ...)
{
    group <- NULL
    dpars <- extract_model_values(x, object_type = "estimate", "distribution")
    distribution <- x$spec$distribution
    dist_print <- distribution_abb(distribution)
    layout(matrix(c(1,1,1,1,2,2,2,2,3,3,4,4), ncol = 4, nrow = 3, byrow = T))
    par(mar = c(3, 4, 3, 2))
    plot(as.zoo(x$spec$target$y), ylab = "", xlab = "", main = "Actual vs Fitted")
    lines(as.zoo(fitted(x)), col = 2)
    grid()
    legend("topleft",c("Actual","Fitted"), col = c("black","red"), lty = 1, bty = "n")
    hist(residuals(x, standardize = T), breaks = "fd", xlab = "", main = "Histogram (Std. Residuals)", probability = T)
    box()
    curve(ddist(distribution = distribution, x, mu = 0, sigma = 1, skew = dpars[1], shape = dpars[2], lambda = dpars[3]), col = "steelblue", add = T)
    legend("topleft", paste0(dist_print,"(0,1) Density"), col = "steelblue", lty = 1, bty = "n")
    par(mar = c(5, 4, 3, 2))
    mod_roots <- arma_roots(x$parmatrix[group == "phi"]$value, x$parmatrix[group == "theta"]$value)
    plot(mod_roots)
    mod_ird <- arma_irf(x$parmatrix[group == "phi"]$value, x$parmatrix[group == "theta"]$value, h = 25)
    mod_ird <- mod_ird$table
    plot(mod_ird$horizon, mod_ird$response, ylab = "response", xlab = "horizon", main = "Impulse Response\n (unit variance shock)", type = "l")
    grid()
}

#' ARMA Inverse Roots Plot
#'
#' @description Plot method for \dQuote{arma.roots} class.
#' @param x an object of class \dQuote{arma.roots}.
#' @param y not used.
#' @param ... not used.
#' @method plot arma.roots
#' @rdname plot.arma.roots
#' @export
#'
#
plot.arma.roots <- function(x, y = NULL, ...)
{
    equation <- NULL
    x <- x$table
    ar_roots <- x[equation == "AR"]
    ma_roots <- x[equation == "MA"]

    if (nrow(ar_roots) == 0) {
        ar_real_roots <- NULL
        ar_imag_roots <- NULL
        ar_roots <- NULL
    } else {
        ar_real_roots <- Re(ar_roots$root)
        ar_imag_roots <- Im(ar_roots$root)
        ar_roots <- ar_roots$root
    }

    if (nrow(ma_roots) == 0) {
        ma_real_roots <- NULL
        ma_imag_roots <- NULL
        ma_roots <- NULL
    } else {
        ma_real_roots <- Re(ma_roots$root)
        ma_imag_roots <- Im(ma_roots$root)
        ma_roots <- ma_roots$root

    }
    limits <- 1/max(abs(c(ar_real_roots, ma_real_roots)), 1.5, abs(c(ar_imag_roots, ma_imag_roots)))
    if (!is.null(ar_roots)) {
        plot(1/x[equation == "AR"]$root, xlim = c(-1, 1), ylim = c(-1, 1), xlab = "", ylab = "", pch = 23, col = "steelblue")
        if (!is.null(ma_roots)) {
            points(1/x[equation == "MA"]$root, pch = 21, col = "red")
            legend("topleft", c("AR", "MA"), col = c("steelblue","red"), pch = c(23, 21), bty = "n", ncol = 2)
        } else {
            legend("topleft", c("AR"), col = c("steelblue"), pch = c(23), bty = "n", ncol = 1)
        }
        cpoints <- (2 * pi/360) * (0:360)
        lines(sin(cpoints), cos(cpoints), col = "darkgrey")
        abline(h = 0, col = "darkgrey")
        abline(v = 0, col = "darkgrey")
        title("Inverse Roots and Unit Circle\n", xlab = "Real Part", ylab = "Imaginary Part")
    } else{
        if (!is.null(ma_roots)) {
            plot(1/x[equation == "MA"]$root, xlim = c(-1, 1), ylim = c(-1, 1), xlab = "", ylab = "", pch = 21, col = "red")
            cpoints <- (2 * pi/360) * (0:360)
            lines(sin(cpoints), cos(cpoints), col = "darkgrey")
            abline(h = 0, col = "darkgrey")
            abline(v = 0, col = "darkgrey")
            title("Inverse Roots and Unit Circle\n", xlab = "Real Part", ylab = "Imaginary Part")
            legend("topleft", c("MA"), col = c("red"), pch = c(21), bty = "n", xpd = TRUE, inset = c(0, -0.05), ncol = 1)
        }
    }
    return(invisible(x))
}

#' ARMA Impulse Response Function Plot
#'
#' @description Plot method for \dQuote{arma.irf} class.
#' @param x an object of class \dQuote{arma.irf}.
#' @param y not used.
#' @param ... not used.
#' @method plot arma.irf
#' @rdname plot.arma.irf
#' @export
#'
#
plot.arma.irf <- function(x, y = NULL, ...)
{
    y <- x$table
    v <- x$sigma
    plot(y$horizon, y$response, ylab = "response", xlab = "horizon", main = paste0("Impulse Response\n (variance = ",signif(v^2, 3),")"), type = "l")
    points(y$horizon, y$response)
    abline(h = 0, col = "grey")
    return(invisible(x))
}
