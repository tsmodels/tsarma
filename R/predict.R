.predict_arma <- function(object, h = 1, newxreg = NULL, bootstrap = FALSE, nsim = 1000, innov = NULL, innov_type = "r",
                          forc_dates = NULL, series_init = NULL, ...)
{
    model <- object$spec$model
    include_xreg <- object$spec$xreg$include_xreg
    y <- extract_model_values(object, "estimate", "y")
    xreg <- extract_model_values(object, "estimate", "xreg")
    mu <- extract_model_values(object, "estimate", "mu")
    xi <- extract_model_values(object, "estimate", "xi")
    phi <- extract_model_values(object, "estimate", "phi")
    theta <- extract_model_values(object, "estimate", "theta")
    sigma <- extract_model_values(object, "estimate", "sigma")
    maxpq <- max(model$order)
    epsilon <- as.numeric(residuals(object))
    order <- model$order
    if (maxpq == 0) maxpq <- 1
    forc_dates <- .forecast_dates(forc_dates = forc_dates, h = h, sampling = object$spec$target$sampling,
                                  last_index = tail(object$spec$target$index, 1))
    x <- .process_prediction_regressors(old_regressors = coredata(xreg), new_regressors = newxreg,
                                        xi = xi, h = h, include_regressors = include_xreg,
                                        maxpq = maxpq, regressor_argument = "newxreg")

    series_forc <- rep(0, h + maxpq)
    y_f <- c(tail(as.numeric(y), maxpq), rep(0, h))
    epsilon <- c(tail(as.numeric(epsilon), maxpq), rep(0, h))
    constant <- as.numeric(x) + mu
    for (i in (maxpq + 1):(maxpq + h)) {
        y_f[i] <- constant[i]
        if (order[1] > 0) {
            for (j in 1:order[1]) {
                y_f[i] <- y_f[i] + phi[j] * (y_f[i - j] - constant[i - j])
            }
        }
        if (order[2] > 0) {
            for (j in 1:order[2]) {
                if ((i - maxpq) > j) {
                    s <- 0
                } else {
                    s <- theta[j] * epsilon[i - j]
                }
                y_f[i] <- y_f[i] + s
            }
        }

    }
    series_forc <- tail(y_f, h)
    psi_weights <- c(1, ARMAtoMA(phi, theta, lag.max = h + 1))
    psi2_cumsum <- cumsum(psi_weights^2)
    sigma_forc <- sqrt((sigma^2) * psi2_cumsum[1:h])
    # simulated distribution
    innov_init <- as.numeric(tail(residuals(object), maxpq))
    series_init <- as.numeric(tail(y, maxpq))
    if (bootstrap) {
        innov <- matrix(sample(as.numeric(residuals(object)), nsim * h, replace = TRUE), nrow = nsim, ncol = h)
        sim <- simulate(object, h = h, nsim = nsim, innov = innov, innov_init = innov_init, series_init = series_init, innov_type = "r")
    } else {
        sim <- simulate(object, h = h, nsim = nsim, innov = NULL, innov_init = innov_init, series_init = series_init)
    }
    sim <- sim$simulated
    colnames(sim) <- as.character(forc_dates)
    class(sim) <- "tsmodel.distribution"
    L <- list(original_series = object$spec$target$y, distribution = sim, spec = object$spec, sigma = xts(sigma_forc, forc_dates), mean = xts(series_forc, forc_dates))
    class(L) <- c("tsarma.predict","tsmodel.predict")
    return(L)
}
