extract_model_values <- function(object, object_type, value_name, ...)
{
    group <- NULL
    if (object_type == "spec") {
        value <- switch(value_name,
                        "y" = object$target$y,
                        "xreg" = xts(object$xreg$xreg, object$target$index),
                        "xi" = object$parmatrix[group == "xi"]$value,
                        "mu" = object$parmatrix[group == "mu"]$value,
                        "phi" = object$parmatrix[group == "phi"]$value,
                        "theta" = object$parmatrix[group == "theta"]$value,
                        "sigma" = object$parmatrix[group == "sigma"]$value,
                        "distribution" = object$parmatrix[group == "distribution"]$value)
    } else if (object_type == "estimate") {
        value <- switch(value_name,
                        "y" = object$spec$target$y,
                        "xreg" = xts(object$spec$xreg$xreg, object$spec$target$index),
                        "xi" = object$parmatrix[group == "xi"]$value,
                        "mu" = object$parmatrix[group == "mu"]$value,
                        "phi" = object$parmatrix[group == "phi"]$value,
                        "theta" = object$parmatrix[group == "theta"]$value,
                        "sigma" = object$parmatrix[group == "sigma"]$value,
                        "distribution" = object$parmatrix[group == "distribution"]$value)
    } else if (object_type == "env") {

    }
    return(value)
}
