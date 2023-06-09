% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{summary.tsarma.estimate}
\alias{summary.tsarma.estimate}
\alias{summary}
\title{ARMA Model Estimation Summary}
\usage{
\method{summary}{tsarma.estimate}(object, digits = 4, vcov_type = "H", ...)
}
\arguments{
\item{object}{an object of class \dQuote{tsarma.estimate}.}

\item{digits}{integer, used for number formatting. Optionally, to avoid
scientific notation, set \sQuote{options(scipen=999)}.}

\item{vcov_type}{the type of standard errors based on the vcov estimate (see \code{\link{vcov}}).}

\item{...}{not currently used.}
}
\value{
A list with summary information of class \dQuote{tsarma.summary}.
}
\description{
Summary method for class \dQuote{tsarma.estimate}
}
