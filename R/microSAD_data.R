#' Node and weight values for Gauss-Kronrod quadrature formula in quadgk
#'
#' A named list of 2 elements, with each element is a named list of 4 numeric vectors
#'
#' @format A named list with 2 elements, with each element is a named list of 4 numeric vectors
#'
#' \describe{
#'   \item{G7_K15}{The node and weight values for a G7_K15 Gauss-Kronrod quadrature formula}
#'   \item{G25_K51}{The node and weight values for a G25_K51 Gauss-Kronrod quadrature formula}
#'   \item{$K_node}{The node values for Kronrod rule}
#'   \item{$K_weight}{The weight values for Kronrod rule}
#'   \item{$G_node}{The node values for Gauss rule}
#'   \item{$G_weight}{The weight values for Gauss rule}
#' }
#' @source All values are cited from Pavel Holoborodko, 2011
#' \url{https://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights}
"quadgk_constant"