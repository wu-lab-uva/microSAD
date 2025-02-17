% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/octave_fit_test.R
\name{r2m_rad}
\alias{r2m_rad}
\title{Goodness-of-fit for octave fit by RAD-RAD comparison}
\usage{
r2m_rad(octave.x, octaves, fitted.values)
}
\arguments{
\item{octave.x}{the number of species in each octave (e.g., $Data$species from the gambin fit object)}

\item{octaves}{the octave index matching the \code{octave.x} argument (e.g., $Data$octave from the gambin fit object)}

\item{fitted.values}{the expected proportion of species in each octave, calculated from the model fit (e.g., $fitted.values from the gambin fit object)}
}
\value{
\code{r2m_rad} returns a numeric value
}
\description{
Check goodness-of-fit for octave fit in RAD-RAD comparison format to allow comparison with sad-based models
}
\details{
The octaves start with '0', which contains all the species with 1 individuals. Octave 'i' contains species with \link{2^i,2^(i+1)-1} individuals.

\code{fitted.values} generate the empirical cumulative distribution function from which an RAD can be sampled like in representative_* functions.
}
