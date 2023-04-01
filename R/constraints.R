# arma ineq and jacobian ------------------------------------------------------
ar_ineq <- function(pars, env)
{
    group <- value <- NULL
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    phi <- parmatrix[group == "phi"]$value * parmatrix[group == "phi"]$scale
    return(ar_conditions(phi))
}
ma_ineq <- function(pars, env)
{
    group <- value <- NULL
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    theta <- parmatrix[group == "theta"]$value * parmatrix[group == "theta"]$scale
    return(ma_conditions(theta))
}

arma_ineq <- function(pars, env)
{
    group <- value <- NULL
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    phi <- parmatrix[group == "phi"]$value * parmatrix[group == "phi"]$scale
    theta <- parmatrix[group == "theta"]$value * parmatrix[group == "theta"]$scale
    return(c(ar_conditions(phi), ma_conditions(theta)))
}

ar_conditions <- function(pars)
{
    if (all(pars == 0)) return(-1e10)
    1 - min(Mod(polyroot(c(1, -pars))))
}

ma_conditions <- function(pars)
{
    if (all(pars == 0)) return(-1e10)
    1 - min(Mod(polyroot(c(1, pars))))
}


arma_ineq_phi <- function(pars)
{
    return(ar_conditions(pars))
}

arma_ineq_theta <- function(theta)
{
    return(ma_conditions(theta))
}

ar_ineq_jac <- function(pars, env)
{
    group <- value <- NULL
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    phi <- parmatrix[group == "phi"]$value * parmatrix[group == "phi"]$scale
    phij <- grad(func = arma_ineq_phi, x = phi)
    cmatrix <- env$cmatrix
    cmatrix[group == "phi", phi := phij]
    out <- matrix(cmatrix[estimate == 1]$phi, nrow = 1)
    return(out)
}

ma_ineq_jac <- function(pars, env)
{
    group <- value <- NULL
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    theta <- parmatrix[group == "theta"]$value * parmatrix[group == "theta"]$scale
    thetaj <- grad(func = arma_ineq_theta, x = theta)
    cmatrix <- env$cmatrix
    cmatrix[group == "theta", theta := thetaj]
    out <- matrix(cmatrix[estimate == 1]$theta, nrow = 1)
    return(out)
}

arma_ineq_jac <- function(pars, env)
{
    group <- value <- NULL
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    phi <- parmatrix[group == "phi"]$value * parmatrix[group == "phi"]$scale
    theta <- parmatrix[group == "theta"]$value * parmatrix[group == "theta"]$scale
    phij <- grad(func = arma_ineq_phi, x = phi)
    thetaj <- grad(func = arma_ineq_theta, x = theta)
    cmatrix <- env$cmatrix
    cmatrix[group == "phi", phi := phij]
    cmatrix[group == "theta", theta := thetaj]
    out <- rbind(cmatrix[estimate == 1]$phi,
                 cmatrix[estimate == 1]$theta)
    return(out)
}
