library(CaveDive)

set.seed(1)
    
    sam <- runif(100, 0, 10)

    tmax <- max(sam)

    sam <- sam - tmax
    sam <- sam[order(-sam)]
    
    colours <- trunc(runif(100, 1, 4))

    N <- 1000
    
    A <- c(1, 1.3)
    
    K <- c(1000,2000)

    div_times <- c(-15, -30, -Inf)
    div_cols <- c(3, 2, 1)

    rates <- list(function (s) logistic.rate(s, K[1], A[1], div_times[1]), function (s) logistic.rate(s, K[2], A[2], div_times[2]), function (s) constant.rate(s, N))
    rate.ints <- list(function(t,s) logistic.rate.int(t, s, K[1], A[1], div_times[1]), function(t,s) logistic.rate.int(t, s, K[2], A[2],div_times[2]), function(t,s) constant.rate.int(t,s,N))

    co <- structured_coal.simulate(sam, colours, div_times, div_cols, rates, rate.ints)
    print(co)