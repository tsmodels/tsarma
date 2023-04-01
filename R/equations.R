.distribution_abb <- function(distribution)
{
    switch(distribution,
           "norm" = "N\\left(0,\\sigma\\right)",
           "snorm" = "SN\\left(0,\\sigma,\\zeta\\right)",
           "std" = "T\\left(0,\\sigma,\\nu\\right)",
           "sstd" = "ST\\left(0,\\sigma,\\zeta, \\nu\\right)",
           "ged" = "GED\\left(0,\\sigma,\\nu\\right)",
           "sged" = "SGED\\left(0,\\sigma,\\zeta, \\nu\\right)",
           "jsu" = "JSU\\left(0,\\sigma,\\zeta, \\nu\\right)",
           "nig" = "NIG\\left(0,\\sigma,\\zeta, \\nu\\right)",
           "ghst" = "GHST\\left(0,\\sigma,\\zeta, \\nu\\right)",
           "ghyp" = "GH\\left(0,\\sigma,\\zeta, \\nu,\\lambda\\right)")
}

.equation_regressors <- function(xreg = NULL)
{
    eq <- NULL
    if (!is.null(xreg)) {
        n <- NCOL(xreg)
        if (n == 1) {
            eq <- "\\xi_1 x_{1,t}"
        } else {
            eq <- paste0("\\sum_{j=1}^",n," \\xi_j x_{j,t}")
        }
    }
    return(list(regressor = eq))
}


# garch equations
.equation_arma <- function(order, constant = TRUE, xreg = NULL, distribution = "norm")
{
    regressor_equations <- .equation_regressors(xreg)
    # constant
    if (constant) {
        eq_constant <- "\\mu_t = \\mu"
        if (!is.null(xreg)) {
            eq_constant <- paste0(eq_constant," + ", regressor_equations$regressor)
        }
    } else {
        if (!is.null(xreg)) {
            eq_constant <- paste0("\\mu_t = ", regressor_equations$regressor)
        }
        eq_constant <- "\\mu_t = 0"
    }
    eq_mu <- "\\mu_t"

    if (!is.null(xreg)) {
        eq_constant <- paste0(eq_constant," + ", regressor_equations$regressor)
    }
    eq_distribution <- .distribution_abb(distribution)
    if (order[1] > 0) {
        if (order[1] > 1) {
            eq_phi <- paste0("\\sum_{j=1}^",order[1],"\\phi_j \\left(y_{t-j}-\\mu_{t-j}\\right)")
        } else {
            eq_phi <- paste0("\\phi_1\\left(y_{t-1}-\\mu_{t-1}\\right)")
        }
    } else {
        eq_phi <- NULL
    }
    if (order[2] > 0) {
        if (order[2] > 1) {
            eq_theta <- paste0("\\sum_{j=1}^",order[2],"\\theta_j \\varepsilon_{t-j}")
        } else {
            eq_theta <- paste0("\\theta_1\\varepsilon_{t-1}")
        }
    } else {
        eq_theta <- NULL
    }
    eq_arma <- paste0("\\hat y_t = ",eq_mu)
    if (!is.null(eq_phi)) eq_arma <- paste0(eq_arma," + ", eq_phi)
    if (!is.null(eq_theta)) eq_arma <- paste0(eq_arma," + ", eq_theta)
    eq_arma <- paste0(eq_arma, " + ", "\\varepsilon_t")
    eq_distribution <- paste0("\\varepsilon_t \\sim ", eq_distribution)
    out <- list(eq_constant = eq_constant, eq_arma = eq_arma, eq_distribution = eq_distribution)
    return(out)
}
