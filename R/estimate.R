#--------------------------------
# arma model
arma_fun <- function(pars, env)
{
    estimate <- parameter <- value <- NULL
    parmatrix <- env$parmatrix
    parmatrix[estimate == 1, value := pars]
    llh <- env$tmb$fn(pars)
    if (!is.finite(llh) | is.na(llh)) {
        llh <- env$llh + 0.2 * abs(env$llh)
        env$llh <- llh
        return(llh)
    } else {
        env$llh <- llh
        return(llh)
    }
}

arma_grad <- function(pars, env)
{
    env$tmb$gr(pars)
}

arma_hess <- function(pars, env)
{
    env$tmb$he(pars)
}

.tmb_initializa_arma <- function(spec, ...)
{
    include <- value <- group <- .N <- parameter <- phi <- theta <- `.` <- NULL
    parmatrix <- spec$parmatrix
    parmatrix[,include := 1]
    parmatrix[estimate == 0, include := as.numeric(NA)]
    map <- lapply(split(parmatrix[,list(P = 1:.N * include), by = "group"], by = "group", keep.by = FALSE, drop = T), function(x) as.factor(x$P))
    parameters <- lapply(split(parmatrix[,list(value, group)], by = "group", keep.by = FALSE), function(x) as.numeric(x$value))
    cmodel <- spec$model_options
    # create a constraint matrix to hold the gradient of the inequalities
    cmatrix <- parmatrix[,.(parameter, estimate, group)]
    cmatrix[,phi := 0]
    cmatrix[,theta := 0]
    # augment data with max(p,q) vectors
    x <- spec$xreg$xreg
    y <- as.numeric(spec$target$y)
    pscale <- parmatrix$scale
    data <- list(y = y, x = x, pscale = pscale, cmodel = cmodel, model = "arma")
    if (sum(spec$model$order) == 0) {
        ineq_fun <- NULL
        ineq_jac <- NULL
    } else if (spec$model$order[1] > 0 & spec$model$order[2] == 0) {
        ineq_fun <- ar_ineq
        ineq_jac <- ar_ineq_jac
    } else if (spec$model$order[1] == 0 & spec$model$order[2] > 0) {
        ineq_fun <- ma_ineq
        ineq_jac <- ma_ineq_jac
    } else {
        ineq_fun <- arma_ineq
        ineq_jac <- arma_ineq_jac
    }
    return(list(fun = arma_fun, grad = arma_grad, hess = arma_hess, data = data, parameters = parameters,
                map = map, cmatrix = cmatrix, ineq_fun = ineq_fun, ineq_jac = ineq_jac))
}

.estimate_arma <- function(object, solver, control, ...)
{
    parameter <- group <- value <- NULL
    solver <- match.arg(solver, choices = "nloptr")
    L <- .tmb_initializa_arma(object)
    tmb <- MakeADFun(data = L$data, parameters = L$parameters, atomic = TRUE, map = L$map, silent = TRUE, DLL = "tsarma_TMBExports")
    env <- new.env()
    env$fun <- L$fun
    env$grad <- L$grad
    env$hess <- L$hess
    env$tmb <- tmb
    env$llh <- 1
    env$model <- "arma"
    env$parmatrix <- object$parmatrix
    env$cmatrix <- L$cmatrix
    init_pars <- object$parmatrix[estimate == 1]$value * 1/object$parmatrix[estimate == 1]$scale
    lower <- object$parmatrix[estimate == 1]$lower * 1/object$parmatrix[estimate == 1]$scale
    upper <- object$parmatrix[estimate == 1]$upper * 1/object$parmatrix[estimate == 1]$scale
    sol <- nloptr(x0 = init_pars, eval_f = env$fun, eval_grad_f = env$grad, eval_g_ineq = L$ineq_fun,
                  eval_jac_g_ineq = L$ineq_jac, lb = lower, ub = upper, env = env, opts = control)
    pmatrix <- copy(env$parmatrix)
    pmatrix[estimate == 1, value := sol$solution]
    optimal_pars <- sol$solution
    # check hessian and use for scaling
    H <- tmb$he(optimal_pars)
    object$parmatrix <- pmatrix
    scaled_solution <- .estimate_arma_scaled(optimal_pars, H, object, control)
    D <- solve(diag(scaled_solution$par_scale, length(scaled_solution$par_scale), length(scaled_solution$par_scale)))
    pars <- scaled_solution$solution$solution * scaled_solution$solution$par_scale
    solve_conditions <- solver_conditions(scaled_solution$solution$solution, scaled_solution$env$fun, scaled_solution$env$grad, scaled_solution$env$hess, scaled_solution$lower, scaled_solution$upper, scaled_solution$env)
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    spec <- object
    spec$parmatrix <- NULL
    out <- list(parmatrix = parmatrix, scaled_hessian = scaled_solution$hessian,
                scaled_scores = scaled_solution$scores,
                parameter_scale = scaled_solution$par_scale,
                conditions = solve_conditions,
                loglik = scaled_solution$solution$objective,
                npars = sum(parmatrix$estimate),
                nobs = NROW(spec$target$y),
                fitted = scaled_solution$fitted,
                residuals = scaled_solution$residuals,
                spec = spec)
    class(out) <- "tsarma.estimate"
    return(out)
}

.estimate_arma_scaled <- function(pars, H, object, control)
{
    parameter <- value <- NULL
    test <- .is_positive_definite(H)
    if (!test) H <- .make_positive_definite(H)
    par_scale <- sqrt(1/diag(H))
    # stage 2: scaled estimation
    scaled_object <- object
    # make a copy of parmatrix rather than leave the pointer in place
    scaled_object$parmatrix <- copy(object$parmatrix)
    scaled_object$parmatrix[estimate == 1, scale := par_scale]
    scaled_object$parmatrix[estimate == 1, value := pars / par_scale]
    scaled_init_pars <- pars *  1/par_scale
    scaled_lower <- scaled_object$parmatrix[estimate == 1]$lower * 1/par_scale
    scaled_upper <-  scaled_object$parmatrix[estimate == 1]$upper * 1/par_scale
    scaled_L <- .tmb_initializa_arma(scaled_object)
    scaled_tmb <- MakeADFun(data = scaled_L$data, parameters = scaled_L$parameters, map = scaled_L$map, silent = TRUE, DLL = "tsarma_TMBExports")
    scaled_env <- new.env()
    scaled_env$fun <- scaled_L$fun
    scaled_env$grad <- scaled_L$grad
    scaled_env$hess <- scaled_L$hess
    scaled_env$tmb <- scaled_tmb
    scaled_env$llh <- 1
    scaled_env$model <- "arma"
    scaled_env$distribution <- object$distribution
    scaled_env$parmatrix <- scaled_object$parmatrix
    scaled_env$cmatrix <- scaled_L$cmatrix
    scaled_sol <- nloptr(x0 = scaled_init_pars, eval_f = scaled_env$fun, eval_grad_f = scaled_env$grad,
                         eval_g_ineq = scaled_L$ineq_fun, eval_jac_g_ineq = scaled_L$ineq_jac,
                         lb = scaled_lower, ub = scaled_upper, env = scaled_env, opts = control)
    scaled_sol$par_scale <- par_scale
    hessian <- scaled_tmb$he(scaled_sol$solution)
    scores <- jacobian(score_function, scaled_sol$solution, env = scaled_env)
    # fitted, residuals
    m <- object$model_options[1]
    f <- scaled_tmb$report(scaled_sol$solution)$fitted
    r <- scaled_tmb$report(scaled_sol$solution)$epsilon
    return(list(solution = scaled_sol, env = scaled_env, lower = scaled_lower, upper = scaled_upper,
                hessian = hessian, scores = scores, par_scale = par_scale, fitted = f, residuals = r))
}
