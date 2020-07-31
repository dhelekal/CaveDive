#' @export
constant.rate <- function (s, N) return (1 / N)
#' @export
constant.rate.int <- function (t, s, N) {
    return(s / N)
}
#' @export
constant.rate.int_inv <- function(t, s, N) return(s * N)

#' @export
logistic.rate <- function(s, K, rate, t0) return(1/((K-1)*(1-exp((s+t0)*rate))) + 1/((exp((-(s+t0)*rate)-1)))) 
#' @export
logistic.rate.int <- function(t, s, K, rate, t0) {

    func1 <- function(x) (-1/((K-1)*rate))*log(exp(-rate*x)-1)
    func2 <- function(x) (1/rate)*log(1-exp(rate*x))

    #func3 <- function(x) (1/rate)*log((1-exp(rate*x)) / ((exp(-rate*x)-1)^(1/(K-1))))

    if (t+t0 < 0 && t+s+t0 < 0) {
        out <- (func1(t+t0+s)-func1(t+t0)) + (func2(t+t0+s)-func2(t+t0))
        out <- func3(t+s+t0) - func3(t+t0)
    } else{
        out <- Inf
    }
    return(out)
}