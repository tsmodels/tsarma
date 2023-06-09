% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{tsequation.tsarma.estimate}
\alias{tsequation.tsarma.estimate}
\alias{tsequation}
\title{Model Equation (LaTeX)}
\usage{
\method{tsequation}{tsarma.estimate}(object, ...)
}
\arguments{
\item{object}{an object of class \dQuote{tsarma.estimate}.}

\item{...}{not currently used.}
}
\value{
A list of equations in LaTeX which can be used in documents.
}
\description{
Generates a list of model equations in LaTeX.
}
\details{
This method is called in the summary when the format output option
chosen is \dQuote{flextable}.
}
