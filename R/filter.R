.filter_arma_estimate <- function(object, y = NULL, newxreg = NULL, ...)
{
    parameter <- group <- NULL
    if (is.null(y)) return(object)
    if (!is.xts(y)) stop("\ny must be an xts vector")
    if (any(index(y) < max(object$spec$target$index))) stop("\ntime index of new y vector must be greater than existing index of y")
    y_new <- .merge_data(object$spec$target$y, y)
    new_n <- NROW(y_new) - NROW(object$spec$target$y)
    n <- NROW(y)
    maxpq <- max(object$spec$model$order)
    model <- c(maxpq, object$spec$model$order)
    x_orig <- xts(object$spec$xreg$xreg, object$spec$target$index)
    mu <- object$parmatrix[parameter == "mu"]$value
    phi <- object$parmatrix[group == "phi"]$value
    theta <- object$parmatrix[group == "theta"]$value
    xi <- object$parmatrix[group == "xi"]$value
    x_new <- .process_filter_regressors(old_regressors = x_orig, new_regressors = newxreg, new_index = index(y), new_n = new_n,
                                        include_regressors = object$spec$xreg$include_xreg, maxpq = maxpq)
    if (maxpq > 0) initstate <- as.numeric(tail(fitted(object), maxpq)) else initstate <- 0
    x <- as.numeric(x_new %*% xi)
    epsilon <- rep(0, new_n)
    if (maxpq > 0) {
        epsilon <- c(rep(0, maxpq), epsilon)
        epsilon[1:maxpq] <- as.numeric(tail(residuals(object), maxpq))
    }
    y_new <- tail(as.numeric(y_new), maxpq + new_n)
    fitted <- .armafilter(y = y_new, epsilon = epsilon, x = x, initstate = initstate, mu = mu, phi = phi, theta = theta, model = model)
    fitted <- c(object$fitted, tail(fitted, new_n))
    object$spec$target$y_orig <- c(object$spec$target$y_orig, tail(object$spec$target$y_orig, new_n))
    y_new <- xts(tail(y_new, new_n), tail(index(y), new_n))
    object$spec$target$y <- rbind(object$spec$target$y, y_new)
    object$spec$target$y_orig <- c(object$spec$target$y_orig, as.numeric(y_new))
    object$spec$target$index <- c(object$spec$target$index, index(y_new))
    object$spec$xreg$xreg <- rbind(object$spec$xreg$xreg, tail(x_new, new_n))
    object$residuals <- object$spec$target$y - fitted
    object$fitted <- fitted
    object$nobs <- length(fitted)
    return(object)
}

.filter_arma_spec <- function(object, y = NULL, newxreg = NULL, ...)
{
    if (!is.null(y)) {
        if (!is.xts(y)) stop("\ny must be an xts vector")
        y <- .merge_data(object$target$y, y)
        if (!is.null(newxreg) & object$xreg$include_xreg) {
            if (is.null(newxreg)) stop("\nnewxreg cannot be NULL when y supplied and xreg was included in original specification")
            if (!is.xts(newxreg)) stop("\nnewxreg must be an xts object")
            if (!all.equal(index(y), index(newxreg))) stop("\nindex y must be the same as the index of xreg.")
            x_orig <- xts(object$xreg$xreg, object$target$index)
            xreg <- .merge_data(x_orig, newxreg)
        } else {
            xreg <- matrix(0, ncol = 1, nrow = NROW(y))
        }
        good <- rep(1, NROW(y))
        if (any(is.na(y))) {
            good[which(is.na(y))] <- 0
        }
        object$target$y_orig <- as.numeric(y)
        object$target$index <- index(y)
        object$target$y <- y
        object$target$good <- good
        object$xreg$xreg <- coredata(xreg)
    }
    parameter <- group <- NULL
    maxpq <- max(object$model$order)
    model <- c(maxpq, object$model$order)
    parmatrix <- object$parmatrix
    active_p <- sum(parmatrix$estimate)
    n <- NROW(object$target$y_orig)
    x <- object$xreg$xreg
    mu <- parmatrix[parameter == "mu"]$value
    phi <- parmatrix[group == "phi"]$value
    theta <- parmatrix[group == "theta"]$value
    xi <- parmatrix[group == "xi"]$value
    x <- as.numeric(x %*% xi)
    epsilon <- rep(0, n)
    y <- object$target$y_orig
    fitted <- .armafilter2(y = y, epsilon = epsilon, x = x, mu = mu, phi = phi, theta = theta, model = model)
    residuals <- y - fitted
    parmatrix[,estimate := 0]
    logl <- sum(-ddist(object$distribution, y, mu = fitted, sigma = parmatrix[parameter == "sigma"]$value,
                       skew = parmatrix[parameter == "skew"]$value, shape = parmatrix[parameter == "shape"]$value,
                       lambda = parmatrix[parameter == "lambda"]$value, log = TRUE))
    out <- list(parmatrix = object$parmatrix,
                parameter_scale = rep(1, active_p),
                loglik = logl,
                npars = 0,
                nobs = n,
                fitted = fitted,
                residuals = residuals,
                spec = object)
    class(out) <- c("tsarma.estimate","tsarma.filter")
    return(out)
}

