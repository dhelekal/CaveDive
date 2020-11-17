tol<-1e-6

#' @export
constant.rate <- function (s, N) return (1 / N)
#' @export
constant.rate.int <- function (t, s, N) {
    return(s / N)
}
#' @export
constant.rate.int_inv <- function(t, s, N) return(s * N)

#' @export
sat.rate <- function(s, K, rate, t0) return((1/K)*((1+rate*(t0+s)^2)/(rate*(t0+s)^2)))

#' @export
sat.rate.int <- function(t, s, K, rate, t0) {

    func1 <- function(x) ((1/K)*(x-(1/rate)*(1/x)))

    if (t+t0 < -tol && t+s+t0 < -tol) {
        out <- (func1(t+t0+s)-func1(t+t0))
    } else{
        out <- Inf
    }
    return(out)
}

#' @export
sat.rate.int_inv <- function(t, F, K, rate, t0) {

    return( (-1/rate)*log(exp(-rate*K*F)*(exp(-rate*(t+t0))-1) + 1) - (t+t0)) 

}

#' Numerically invert rates
#'
#' @param func function to invert
#' @param N function value
#' @param dfun analytic derivative for quasi newtonian scheme
#' @return value of inverse function
#' @export
inv_rates <- function(func, N, dfun=NULL){   
    f <- function(x) func(x) - N
    lo <- 0
    hi <- 1e0
    if (sign(f(lo)) == sign(f(hi))){

        sg <- sign(f(lo))

        min_hi <- 1
        while(sign(f(hi*(10^(min_hi)))) == sg && min_hi < 20){
            min_hi <- min_hi+1
        }

        if (sign(f(hi*(10^(min_hi)))) == sg) {
            warning(paste0("Cannot find initial bisection search interval. ", "min_lo: ", 0, " min_hi: ", min_hi))
            return(NA)
        }


        hi <- hi*(10^(min_hi))
    }

    b <- bisect(f, lo, hi, maxiter = 50)

    out <- b$root

    if (b$f.root > 1e-10) {
        n <- newtonRaphson(f,out,dfun=dfun)
        if (n$f.root >  1e-10) {
            warning(paste0("Function root suspiciously large, please revise. f.root: ", n$f.root))
        }
        out <-n$root
    }
    return(out)
}
