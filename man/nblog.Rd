% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nblog.R
\name{dnblog}
\alias{dnblog}
\alias{dconv_nblog}
\alias{dtrunc_dnblog}
\alias{pnblog}
\alias{pconv_nblog}
\alias{rnblog}
\alias{fit_nblog}
\title{Negative-binomial log-normal distribution functions}
\usage{
dnblog(
  x,
  gcn = rep(1, length(x)),
  mu,
  sig,
  odisp,
  log = FALSE,
  constant = "G25_K51"
)

dconv_nblog(x, q, mu, sig, odisp, log = FALSE, adjust = 0)

dtrunc_dnblog(
  x,
  gcn = rep(1, length(x)),
  mu,
  sig,
  odisp,
  trunc = 0,
  constant = "G25_K51",
  log = FALSE,
  debug = FALSE
)

pnblog(
  q,
  gcn = rep(1, length(q)),
  mu,
  sig,
  odisp,
  log.p = FALSE,
  lower.tail = FALSE,
  constant = "G25_K51"
)

pconv_nblog(x, q, mu, sig, odisp, log = FALSE, adjust = 0)

rnblog(n, mu, sig, odisp)

fit_nblog(
  x,
  gcn = rep(1, length(x)),
  trunc = 0,
  start.value,
  constant = "G25_K51",
  min.odisp = 0.001,
  min.sig = 1e-05,
  ...
)
}
\arguments{
\item{x, q}{the observed species abundances; in 16S profiling, this is the observed 16S read count}

\item{gcn}{the 16S copy number of OTUs; it should be a numeric vector of the same length as \code{x}}

\item{mu}{the log mean of the log-normal distribution}

\item{sig}{the log standard deviation of the log-normal distribution}

\item{odisp}{the over-dispersion parameter \code{r} or \code{size} in the negative-binomial distribution that represents sampling bias}

\item{log, log.p}{logical; if \code{TRUE}, the log-transformed value will be returned; defaults to \code{FALSE}}

\item{constant}{the node and weight values for Gauss-Kronrod quadrature formula; see \code{Details} for more information}

\item{trunc}{the truncation cutoff}

\item{n}{the number of random numbers to be generated}
}
\value{
\code{dnblog},\code{dtrunc_dnblog}, \code{pnblog}, \code{dconv_nblog} and \code{pconv_nblog} return a numeric vector the same length as \code{x} or \code{q}

\code{rnblog} return a numeric vector of length \code{n}

\code{fit_nblog} return a \code{mle2} class object of the fitted model
}
\description{
The probability mass and random number generator functions for the negative-binomial log-normal distribution
}
\details{
\code{dnblog} and \code{dtrunc_dnblog} are the (truncated) probability mass function; \code{pnblog} is the cumulative distribution function; \code{rnblog} is the random number generator function.

\code{dconv_nblog} and \code{pconv_nblog} are the integrand functions for \code{dnblog} and \code{pnblog}, respectively; despite their role as internal functions, it is exported for user-created functions.

In negative-binomial log-normal distribution, the sampling effort \code{eta}, or the mean of the negative-binomial distribution, has perfect collinearity with the \code{mu} in the log-normal distribution. Therefore, \code{eta} is degenerated and fixed to 1 internally (before GCN correction).

\code{constant} can be either a character value naming the specific element to use in the data list \code{quadgk_constant}, or a user-supplied list of numeric vectors containing the node and weight values.

The quantile function \code{qnblog} is not implemented and thus \code{rnblog} only generates random values under the full distribution, making it not suitable for RAD prediction.

For predicting RAD, use \code{representative_nblog} for approximated quantiles.
}
