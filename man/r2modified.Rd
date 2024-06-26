% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lack_of_fit_test.R
\name{r2modified}
\alias{r2modified}
\title{Goodness of fit by R^2}
\usage{
r2modified(
  obs,
  pred,
  weight = rep(1, length(obs)),
  log = FALSE,
  normalize = FALSE
)
}
\arguments{
\item{obs}{the observed values}

\item{pred}{the predicted values}

\item{weight}{the weight of each value; defaults to equal weights; these weights are used to balance data points from datasets of different sizes.}

\item{log}{if \code{TRUE}, the R^2 is calculated on a log-log scale; defaults to \code{FALSE}.}

\item{normalize}{if \code{TRUE}, the R^2 is normalized between 0 and 1; defaults to \code{FALSE}.}
}
\value{
\code{r2modified} returns a single numeric value <=1
}
\description{
Calculate the R^2 of a 1 to 1 fit
}
\details{
This function is adapted from the modified coefficient of determination described by Shoemaker et al (2017). Instead of coefficient of determination R^2 in linear regression, \code{r2modified} calculates the R^2 of 1 to 1 fit; as a result, it can be negative.
}
\references{
Shoemaker, William R., Kenneth J. Locey, and Jay T. Lennon. "A macroecological theory of microbial biodiversity." Nature ecology & evolution 1.5 (2017): 0107.
}
