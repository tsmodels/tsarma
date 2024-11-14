sampling_frequency <- function(x)
{
    if (is(x, "Date") || length(grep("POSIX", class(x))) > 0) {
        dates <- x
    }
    else {
        dates <- index(x)
    }
    u <- min(diff(dates))
    count <- attr(u, "units")
    if (count == "days") {
        u <- round(u)
        daily <- c(1, 2, 3)
        weekly <- c(4, 5, 6, 7)
        monthly <- c(27, 28, 29, 30, 31, 32)
        yearly <- 355:370
        if (u %in% daily) {
            period <- "days"
            attr(period, "date_class") <- "Date"
        }
        else if (u %in% weekly) {
            period <- "weeks"
            attr(period, "date_class") <- "Date"
        }
        else if (u %in% monthly) {
            period <- "months"
            attr(period, "date_class") <- "Date"
        }
        else if (u %in% yearly) {
            period <- "years"
            attr(period, "date_class") <- "Date"
        }
        else {
            period <- "unknown"
            attr(period, "date_class") <- "POSIXct"
        }
    }
    else if (count == "hours") {
        period <- paste0(u, " hours")
        attr(period, "date_class") <- "POSIXct"
    }
    else if (count == "mins") {
        period <- paste0(u, " mins")
        attr(period, "date_class") <- "POSIXct"
    }
    else if (count == "secs") {
        period <- paste0(u, " secs")
        attr(period, "date_class") <- "POSIXct"
    }
    else {
        period <- "unknown"
        attr(period, "date_class") <- "POSIXct"
    }
    if (period == "unknown")
        warning("\ncould not determine sampling frequency")
    return(period)
}

initialize_data <- function(y)
{
    n <- NROW(y)
    good <- rep(1, NROW(y))
    if (any(is.na(y))) {
        good[which(is.na(y))] <- 0
    }
    sampling <- sampling_frequency(index(y))
    spec <- list()
    spec$target$y_orig <- as.numeric(y)
    spec$target$index <- index(y)
    spec$target$sampling <- sampling
    spec$target$y <- y
    spec$target$good <- good
    return(spec)
}

.parameters_arma <- function(y, constant = FALSE, order = c(0,0), xreg = NULL, distribution = "norm", ...)
{
    parameter <- value <- group <- lower <- upper <- estimate <- NULL
    y <- as.numeric(y)
    maxpq <- max(order)
    n <- NROW(y)
    if (constant) {
        mu <- mean(y, na.rm = TRUE)
    } else {
        mu <- 0.0
    }
    var_y <- var(y, na.rm = T)
    parmatrix <- data.table("parameter" = "mu", value = mu,
                            lower =  -1.0 * abs(mu) * 100,
                            upper = abs(mu) * 100,
                            estimate = ifelse(constant, 1, 0),
                            scale = 1, group = "mu", equation = "[M]",
                            symbol = "\\mu")
    if (order[1] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "phi1", value = 0,
                                      lower = -1, upper = 1, estimate = 0,
                                      scale = 1, group = "phi",
                                      equation = "[AR]", symbol = "\\phi_1"))
    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("phi",1:order[1]),
                                      value = 0.01, lower = -1, upper = 1,
                                      estimate = 1, scale = 1, group = "phi",
                                      equation = "[AR]",
                                      symbol = paste0("\\phi_",1:order[1])))
    }
    if (order[2] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "theta1", value = 0, lower = -1,
                                      upper = 1, estimate = 0, scale = 1,
                                      group = "theta", equation = "[MA]",
                                      symbol = paste0("\\theta_",1)))

    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("theta",1:order[2]),
                                      value = 0.01, lower = -1, upper = 1,
                                      estimate = 1, scale = 1, group = "theta",
                                      equation = "[MA]",
                                      symbol = paste0("\\theta_",1:order[2])))

    }
    if (!is.null(xreg)) {
        m <- NCOL(xreg)
        xmatrix <- data.table("parameter" = paste0("xi",1:m), value = 1,
                              lower = -100, upper = 100, estimate = 1, scale = 1,
                              group = "xi", equation = "[M]",
                              symbol = paste0("\\xi_",1:m))
    } else {
        xmatrix <- data.table("parameter" = paste0("xi",1), value = 0, lower = -100,
                              upper = 100, estimate = 0, scale = 1, group = "xi",
                              equation = "[M]", symbol = paste0("\\xi_",1))
    }
    parmatrix <- rbind(parmatrix, xmatrix)
    sigma_matrix <- data.table("parameter" = "sigma", value = sqrt(var_y),
                               lower = 1e-12, upper = 10 * sqrt(var_y), estimate = 1, scale = 1,
                               group = "sigma", equation = "[V]",
                               symbol = paste0("\\sigma"))
    parmatrix <- rbind(parmatrix, sigma_matrix)
    dmatrix <- distribution_parameters(distribution)
    parmatrix <- rbind(parmatrix, dmatrix)
    parmatrix[,estimate := as.integer(estimate)]
    # initialize values
    start_pars <- .arma_start_pars(y, constant = constant, order = order, xreg = xreg, method = "arima")
    if (!is.null(start_pars)) {
        parmatrix[group == "mu", value := start_pars$mu]
        parmatrix[group == "mu", lower := -10 * abs(start_pars$mu)]
        parmatrix[group == "mu", upper :=  10 * abs(start_pars$mu)]
        parmatrix[group == "phi", value := start_pars$phi]
        parmatrix[group == "theta", value := start_pars$theta]
        parmatrix[group == "xi", value := start_pars$xi]
        parmatrix[group == "xi", lower := -10 * abs(start_pars$xi)]
        parmatrix[group == "xi", upper :=  10 * abs(start_pars$xi)]
        parmatrix[group == "sigma", value := start_pars$sigma]
        parmatrix[group == "sigma", upper :=  10 * start_pars$sigma]
    }
    return(parmatrix)
}

.arma_start_pars <- function(y, constant = TRUE, order = c(0,0), xreg = NULL, method = c("arima", "lm"))
{
    if (!is.null(xreg)) {
        colnames(xreg) <- paste0("xi",1:ncol(xreg))
    }
    init_pars <- arma_initialization(y, order = order, constant = constant)
    if (method == "arima") {
        f_pars <- c()
        if (order[1] > 0) f_pars <- c(f_pars, init_pars$phi0)
        if (order[2] > 0) f_pars <- c(f_pars, init_pars$theta0)
        if (constant) f_pars <- c(f_pars, as.numeric(NA))
        if (!is.null(xreg)) f_pars <- c(f_pars, rep(as.numeric(NA), ncol(xreg)))
        mod <- try(arima(y, order = c(order[1], 0, order[2]), xreg = xreg, include.mean = constant, fixed = f_pars, method = "CSS-ML", SSinit = "Rossignol2011", transform.pars = FALSE), silent =  TRUE)
        if (inherits(mod, 'try-error')) {
            return(.arma_start_pars(y, constant, order, xreg = xreg, method = "lm"))
        }
        cf <- coef(mod)
        coef_names <- names(cf)
        if (order[1] > 0) {
            phi <- unname(cf[grepl("^ar[0-9]",coef_names)])
        } else {
            phi <- 0
        }
        if (order[2] > 0) {
            theta <- unname(cf[grepl("^ma[0-9]",coef_names)])
        } else {
            theta <- 0
        }
        if (constant) {
            mu <- unname(cf[grepl("^intercept",coef_names)])
        } else {
            mu <- 0
        }
        if (!is.null(xreg)) {
            xi <- unname(cf[grepl("^xi[0-9]",coef_names)])
        } else {
            xi <- 0
        }
        sigma <- sqrt(mod$sigma2)
    } else {
        if (order[1] > 0) {
            phi <- init_pars$phi0
        } else {
            phi <- 0
        }
        if (order[2] > 0) {
            theta <- init_pars$theta0
        } else {
            theta <- 0
        }
        if (constant) {
            mu <- mean(y)
        } else {
            mu <- 0
        }
        if (!is.null(xreg)) {
            mod <- lm(y~xreg)
            cf <- coef(mod)
            cf_names <- names(cf)
            xi <- unname(tail(cf, ncol(xreg)))
            if (any(is.na(xi))) {
                xi[which(is.na(xi))] <- 0
            }
        } else {
            xi <- 0
        }
        sigma <- sd(y - mean(y))
    }

    return(list(mu = mu, phi = phi, theta = theta, xi = xi, sigma = sigma))
}



arma_initialization <- function(y, order = c(1, 1), constant = TRUE) {
    x <- as.numeric(y)
    n <- length(x)
    max_order <- max(order)
    if (max_order == 0) return(list(phi0 = 0, theta0 = 0))
    k <- round(1.1 * log(n))
    e <- as.numeric(na.omit(drop(ar.ols(x, order.max = k, aic = FALSE, demean = constant, intercept = FALSE)$resid)))
    ee <- embed(e, max_order + 1)
    x_demeaned <- y - mean(y)
    xx <- embed(x_demeaned[-(1:k)], max_order + 1)
    if (order[1] == 0) {
        coef <- lm(xx[, 1] ~ ee[, (1:order[2]) + 1] - 1)$coef
        phi <- 0
        theta <- unname(coef)
    } else if (order[2] == 0) {
        coef <- lm(xx[, 1] ~ xx[, (1:order[1]) + 1] - 1)$coef
        phi <- unname(coef)
        theta <- 0
    } else {
        coef <- lm(xx[, 1] ~ xx[, (1:order[1]) + 1] + ee[,(1:order[2]) + 1] - 1)$coef
        phi <- coef[1:order[1]]
        theta <- coef[-c(1:order[1])]
    }
    return(list(phi0  = unname(phi), theta0 = unname(theta)))
}




score_function <- function(x, env)
{
    # add one call to the fun for models which need to update data
    # (integration done in R for some values since TMB integration remains
    # challenging)
    tmp <- env$fun(x, env)
    - 1.0 * log(env$tmb$report(par = x)$ll_vector)
}

# start imports from the corpcor package ---------------------------------------
.is_positive_definite <- function(x)
{
    x <- as.matrix(x)
    eval <- eigen(x, only.values = TRUE, symmetric = TRUE)$values
    tol <- max(dim(x)) * max(abs(eval)) * .Machine$double.eps
    if (sum(eval > tol) == length(eval)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

.make_positive_definite <- function(x, tol)
{
    x <- as.matrix(x)
    d <- dim(x)[1]
    if (dim(x)[2] != d) stop("Input matrix is not square!")
    es <- eigen(x, symmetric = TRUE)
    esv <- es$values
    if (missing(tol)) tol <- d * max(abs(esv)) * .Machine$double.eps
    delta <- 2 * tol
    tau <- pmax(0, delta - esv)
    dm <- es$vectors %*% diag(tau, d) %*% t(es$vectors)
    return(x + dm)
}
# end imports from the corpcor package -----------------------------------------

.lag_vector <- function(x, n_lag = 1, remove_na = FALSE, pad = NA)
{
    # has NAs
    x <- as.matrix(x)
    n <- NROW(x)
    d <- NCOL(x)
    if (d == 1) x <- matrix(x, ncol = 1)
    z <- apply(x, 2, FUN = function(y) .embed_vector(y, n_lag + 1)[,n_lag + 1])
    if (!remove_na) z <- rbind(matrix(pad, ncol = d, nrow = n_lag),z)
    return(z)
}

.embed_vector <- function(x, k, by = 1, ascending = FALSE)
{
    x <- matrix(x, ncol = 1)
    n <- NROW(x)
    s <- seq(1, n - k + 1, by = by)
    lens <- length(s)
    cols <- if (ascending) 1:k else k:1
    return(matrix(x[s + rep(cols, rep(lens,k)) - 1], lens))
}

# for filtering: replace old values with new and append new values
.merge_data <- function(old, new)
{
    n <- NCOL(old)
    new_data <- merge(old, new)
    inc <- which(is.na(new_data[,1]))
    new_data[inc,1:n] <- as.numeric(new_data[inc,-c(1:n)])
    new_data <- new_data[,1:n]
    colnames(new_data) <- colnames(old)
    return(new_data)
}

# generation of future dates for use in prediction
.forecast_dates <- function(forc_dates = NULL, h = 1, sampling, last_index)
{
    if (is.null(forc_dates)) {
        forc_dates = .future_dates(last_index, frequency = sampling, n = h)
    } else {
        if (length(forc_dates) != h) stop("\nforc_dates must be a vector of length h")
        if (any(forc_dates <= last_index)) stop("\nforc_dates must be stricly greater than in-sample dates and times.")
    }
    return(forc_dates)
}

# check and process regressors if present for predict
.process_prediction_regressors <- function(old_regressors, new_regressors = NULL, xi = 0, h = 1, include_regressors = FALSE, maxpq = 1, regressor_argument = "newxreg")
{
    if (include_regressors) {
        if (is.null(new_regressors)) stop(paste0("\n",regressor_argument," is NULL but model has a regressor in the variance."))
        if (!is.xts(new_regressors)) new_regressors <- as.matrix(new_regressors)
        if (NROW(new_regressors) != h) stop(paste0("\n",regressor_argument," must have h rows."))
        if (NCOL(new_regressors) != NCOL(old_regressors)) stop(paste0("\n",regressor_argument," must have the same number of columns as regressors in the model."))
        new_v <- rbind(tail(old_regressors,maxpq), coredata(new_regressors))
    } else {
        new_v <- rbind(tail(old_regressors,maxpq), matrix(0, nrow = h, ncol = NCOL(old_regressors)))
    }
    v <- as.numeric(new_v %*% xi)
    return(v)
}
.process_filter_regressors <- function(old_regressors, new_regressors = NULL, new_index, new_n, include_regressors = FALSE, maxpq = 0)
{
    n <- length(new_index)
    if (include_regressors) {
        if (is.null(new_regressors)) stop("\nnewxreg is NULL but model has a regressor in the variance.")
        if (NROW(new_regressors) != n) stop("\nnewxreg must have the same number of rows as the new y vector.")
        if (NCOL(new_regressors) != NCOL(old_regressors)) stop("\nnewxreg does not have the same number of columns as xreg in model.")
        # set the index of the new_regressors to that of the new y
        new_regressors <- xts(coredata(new_regressors), new_index)
        x_new <- coredata(.merge_data(old_regressors, new_regressors))
    } else {
        # we do this in case the filtering involves updating old data without appending new data (which is a case
        # not currently supported since it will error out but may be supported in future).
        if (new_n > 0) {
            x_new <- rbind(coredata(old_regressors), matrix(0, nrow = new_n, ncol = NCOL(old_regressors)))
        } else {
            x_new <- coredata(old_regressors)
        }
    }
    x_new <- tail(x_new, maxpq + new_n)
    rownames(x_new) <- NULL
    return(x_new)
}

.process_innovations <- function(innov, innov_type, innov_init, nsim, h, distribution, sigma, dpars, maxpq)
{
    if (!is.null(innov) & maxpq > 0) {
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        if (innov_type == "q") {
            if (any(innov < 0 | innov > 1 )) {
                stop("\ninnov must be >0 and <1 (uniform samples) for innov_type = 'q'")
            }
            if (any(innov == 0)) innov[which(innov == 0)] <- 1e-12
            if (any(innov == 1)) innov[which(innov == 1)] <- (1 - 1e-12)
            epsilon <- matrix(qdist(distribution, p = innov, mu = 0, sigma = 1, skew = dpars[1], shape = dpars[2], lambda = dpars[3]), nrow = nsim, ncol = h)
            epsilon <- epsilon * sigma
        } else if (innov_type == "z") {
            epsilon <- innov * sigma
        } else {
            epsilon <- innov
        }
        if (maxpq > 0) {
            init <- matrix(rdist(distribution, maxpq * nsim, mu = 0, sigma = sigma, skew = dpars[1], shape = dpars[2], lambda = dpars[3]), nrow = nsim, ncol = maxpq)
            epsilon <- cbind(init, epsilon)
        }
    } else {
        epsilon <- matrix(rdist(distribution, (maxpq + h) * nsim, mu = 0, sigma = sigma, skew = dpars[1], shape = dpars[2], lambda = dpars[3]), nrow = nsim, ncol = h + maxpq)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(order) : ", maxpq))
        epsilon[,1:maxpq] <- matrix(innov_init, ncol = maxpq, nrow = nsim, byrow = TRUE)
    } else {
        if (maxpq > 0) {
            epsilon[,1:maxpq] <-  matrix(rdist(distribution, maxpq * nsim, mu = 0, sigma = sigma, skew = dpars[1], shape = dpars[2], lambda = dpars[3]), nrow = nsim, ncol = maxpq, byrow = TRUE)
        }
    }
    return(epsilon)
}

.calendar_eom <- function(date, ...)
{
    if (!is(date, "Date")) date <- as.Date(date)
    # Add a month, then subtract a day:
    date.lt <- as.POSIXlt(date, format = "%Y-%m-%d", tz = tz(date))
    mon <- date.lt$mon + 2
    year <- date.lt$year
    # If month was December add a year
    year <- year + as.integer(mon == 13)
    mon[mon == 13] <- 1
    iso <- ISOdate(1900 + year, mon, 1, hour = 0, tz = tz(date))
    result <- as.POSIXct(iso) - 86400 # subtract one day
    result <- result + (as.POSIXlt(iso)$isdst - as.POSIXlt(result)$isdst)*3600
    result <- as.Date(result)
    return(result)
}

.future_dates <- function(start, frequency, n = 1)
{
    if (frequency %in% c("days", "weeks", "months","years")) {
        switch(frequency,
               "days"   = as.Date(start) %m+% days(1:n),
               "weeks"  = as.Date(start) %m+% weeks(1:n),
               "months" = .calendar_eom(as.Date(start) %m+% months(1:n)),
               "years"  = as.Date(start) %m+% years(1:n))
    } else if (grepl("secs|mins|hours|",frequency)) {
        # Add one extra point and eliminate first one
        seq(as.POSIXct(start), length.out = n + 1, by = frequency)[-1]
    } else{
        as.Date(start) + (1:n)
    }
}
