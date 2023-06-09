% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/backtest.R
\name{tsbacktest.tsarma.spec}
\alias{tsbacktest.tsarma.spec}
\alias{tsbacktest}
\title{Walk Forward Rolling Backtest}
\usage{
\method{tsbacktest}{tsarma.spec}(
  object,
  start = floor(length(object$target$y_orig))/2,
  end = length(object$target$y_orig),
  h = 1,
  estimate_every = 1,
  rolling = FALSE,
  trace = FALSE,
  ...
)
}
\arguments{
\item{object}{an object of class \\dQuote{tsarma.spec}.}

\item{start}{numeric data index from which to start the backtest.}

\item{end}{numeric data index on which to end the backtest. The backtest will
end 1 period before that date in order to have at least 1 out of sample value
to compare against.}

\item{h}{forecast horizon. As the expanding window approaches the \\dQuote{end},
the horizon will automatically shrink to the number of available out of sample
periods.}

\item{estimate_every}{number of periods at which the model is re-estimated
(defaults to 1).}

\item{rolling}{this indicates whether forecasts are made only on the estimation
date (FALSE) or whether to filter the data 1 period at a time and forecast
from the filtered data (TRUE).}

\item{trace}{whether to show the progress bar. The user is expected to have
set up appropriate handlers for this using the \\dQuote{progressr} package.}

\item{...}{not currently used.}
}
\value{
A list which includes a data.table having the following columns:
\itemize{
\item estimation_date: the date at which the model was estimated.
\item filter_date: the date on which a prediction was generated. For rolling
prediction this means that an estimated model was filtered for new data prior
to re-predicting.
\item horizon: the forecast horizon of the prediction.
\item size: the length of the data used in estimation.
\item forecast_date: the date corresponding to the forecast.
\item mu: the conditional mean prediction.
\item sigma: the conditional volatility prediction.
\item skew: the distribution skew parameter (non-time varying hence constant
across each estimation window).
\item shape: the distribution shape parameter (non-time varying hence constant
across each estimation window).
\item shape: the distribution lambda parameter (non-time varying hence constant
across each estimation window).
\item actual: the actual observation corresponding to the forecast date.
}
Additional slots in the list include the distribution used and other information
relating to the backtest setup.
}
\description{
Generates an expanding window walk forward backtest with option
for rolling the forecast by filtering (see details).
}
\details{
The rolling option allows to re-estimate the data every n periods
whilst filtering the data 1-step ahead between re-estimation dates so that overlapping
forecasts are generated.
}
\note{
The function can use parallel functionality as long as the user has
set up a \\code{\\link[future]{plan}} using the future package.
}
