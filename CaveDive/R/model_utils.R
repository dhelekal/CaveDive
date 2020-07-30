#' @export
constant.rate <- function (s, N) return (1 / N)
#' @export
constant.rate.int <- function (t, s, N) {
    return(s / N)
}
#' @export
constant.rate.int_inv <- function(t, s, N) return(s * N)

#' @export
logistic.rate <- function(s, K, rate) return(1/((K-1)*(1-exp(s*rate))) + 1/((exp(-s*rate)-1))) 
#' @export
logistic.rate.int <- function(t, s, K, rate) {
    Q <- ((K-2)/((K-1)*rate))
    out <- Q*( log( (1 - exp((t+s)*rate)) / (exp((-rate)*(t+s)) - 1) ) - log( (1 - exp(t*rate)) / (exp((-rate)*t) - 1) ) )
    return(out)
}