.simulate_arma <- function(object, h = 1, seed = NULL, nsim = 1000, innov = NULL, innov_type = "r", innov_init = NULL, series_init = NULL, xreg = NULL, type = "estimate", ...)
{
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }

    if (is(object, "tsarma.spec")) {
        model <- object$model
        include_xreg <- object$xreg$include_xreg
        distribution <- object$distribution
        sampling <- object$target$sampling
    } else if (is(object, "tsarma.estimate")) {
        model <- object$spec$model
        include_xreg <- object$spec$xreg$include_xreg
        distribution <- object$spec$distribution
        sampling <- object$spec$target$sampling
    }
    y <- extract_model_values(object, object_type = type, "y")
    old_xreg <- extract_model_values(object, object_type = type, "xreg")
    mu <- extract_model_values(object, object_type = type, "mu")
    xi <- extract_model_values(object, object_type = type, "xi")
    phi <- extract_model_values(object, object_type = type, "phi")
    theta <- extract_model_values(object, object_type = type, "theta")
    sigma <- extract_model_values(object, object_type = type, "sigma")
    dist <- extract_model_values(object, object_type = type, "distribution")

    maxpq <- max(model$order)
    order <- model$order
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    sim_dates <- .forecast_dates(forc_dates = NULL, h = h, sampling = sampling, last_index = tail(index(y), 1))

    epsilon <- .process_innovations(innov = innov, innov_type = innov_type,
                                    innov_init = innov_init, nsim = nsim, h = h,
                                    distribution = distribution, sigma = sigma,
                                    dpars = dist, maxpq = maxpq)

    x <- .process_prediction_regressors(old_regressors = coredata(old_xreg),
                                        new_regressors = xreg, xi = xi, h = h,
                                        include_regressors = include_xreg,
                                        maxpq = maxpq,
                                        regressor_argument = "xreg")
    constant <- as.numeric(x) + mu
    if (!is.null(series_init) & maxpq > 0) {
        if (length(series_init) != maxpq) stop(paste0("\nseries_init must be of length max(order) : ", maxpq))
        series_sim[,1:maxpq] <- matrix(series_init, nrow = nsim, ncol = maxpq, byrow = TRUE)
    } else {
        if (maxpq > 0) {
            series_sim[,1:maxpq] <- matrix(tail(as.numeric(y), maxpq), nrow = nsim, ncol = maxpq, byrow = TRUE)
        }
    }
    for (i in (maxpq + 1):(h + maxpq)) {
        series_sim[,i] <- constant[i]
        if (order[1] > 0) {
            for (j in 1:order[1]) {
                series_sim[,i] <- series_sim[,i] + phi[j] * (series_sim[,i - j] - constant[i - j])
            }
        }
        if (order[2] > 0) {
            for (j in 1:order[2]) {
                series_sim[,i] <- series_sim[,i] + theta[j] * epsilon[,i - j]
            }
        }
        series_sim[,i] <- series_sim[,i] + epsilon[,i]
    }
    if (maxpq > 0) series_sim <- series_sim[,-c(1:maxpq), drop = FALSE]
    colnames(series_sim) <- as.character(sim_dates)
    class(series_sim) <- "tsmodel.distribution"
    L <- list(simulated = series_sim, error = epsilon)
    class(L) <- "tsarma.simulate"
    return(L)
}
