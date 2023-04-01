% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{print.summary.tsarma}
\alias{print.summary.tsarma}
\title{Model Estimation Summary Print method}
\usage{
\method{print}{summary.tsarma}(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  include.symbols = TRUE,
  include.equation = TRUE,
  include.statistics = TRUE,
  table.caption = paste0(toupper(x$model), " Model Summary"),
  format = c("console", "flextable"),
  ...
)
}
\arguments{
\item{x}{an object of class \dQuote{summary.tsarma}.}

\item{digits}{integer, used for number formatting. Optionally, to avoid
scientific notation, set \sQuote{options(scipen=999)}.}

\item{signif.stars}{logical. If TRUE, ‘significance stars’ are printed for each coefficient.}

\item{include.symbols}{logical. If TRUE, replaces parameter names with their symbols (if they exist).}

\item{include.equation}{logical. If TRUE, adds a section with the symbolic model equation.}

\item{include.statistics}{logical. If TRUE, adds a section with summary statistics on the model.}

\item{table.caption}{an optional string for the table caption.}

\item{format}{either prints to \dQuote{console} or prints and returns a \dQuote{flextable} object.}

\item{...}{additional arguments passed to flextable print method.}
}
\value{
Invisibly returns a flextable object if output was set to \dQuote{flextable} else
the original summary object if output was set to \dQuote{console}.
}
\description{
Print method for class \dQuote{summary.tsarma}
}