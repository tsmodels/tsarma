# validation function for currently implemented distriibutions
valid_distributions <- function()
{
    c("norm", "std", "snorm", "sstd", "ged", "sged", "nig", "ghyp", "jsu", "ghst")
}
distribution_parameters <- function(distribution)
{
    tmp <- NULL
    if (distribution == "norm") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "ged") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 2, lower = 0.1, upper = 50, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "std") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 4, lower = 2.01, upper = 100, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "snorm") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0.5, lower = 0.1, upper = 10, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "sged") {
        tmp <- rbind(
            data.table(parameter = "skew", value = 1, lower = 0.01, upper = 30, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
            data.table(parameter = "shape", value = 2, lower = 0.1, upper = 60, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
            data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "sstd") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 1, lower = 0.01, upper = 30, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 4, lower = 2.01, upper = 60, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "nig") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0.2, lower = -0.99, upper = 0.99, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 0.4, lower = 0.01, upper = 25, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "ghyp") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0.2, lower = -0.99, upper = 0.99, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 2, lower = 0.25, upper = 25, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "jsu") {
        # johnson has 2 shape parameters. The second one we model with the "skew"
        # representation in rugarch
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0, lower = -20, upper = 20, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 1, lower = 0.1, upper = 10, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "ghst") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0.1, lower = -80, upper = 80, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 8.5, lower = 4.1, upper = 25, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
}

distribution_class <- function(distribution)
{
    switch(distribution,
           "norm" = 1,
           "std" = 2,
           "snorm" = 3,
           "sstd" = 4,
           "ged" = 5,
           "sged" = 6,
           "nig" = 7,
           "ghyp" = 8,
           "jsu" = 9,
           "ghst" = 10
    )
}

distribution_abb <- function(distribution)
{
    switch(distribution,
           "norm" = "N",
           "std" = "T",
           "snorm" = "SN",
           "sstd" = "ST",
           "ged" = "GED",
           "sged" = "SGED",
           "nig" = "NIG",
           "ghyp" = "GH",
           "jsu" = "JSU",
           "ghst" = "GHST"
    )
}
