library(CaveDive)

set.seed(1)
    
    sam <- runif(100, 0, 10)

    tmax <- max(sam)

    sam <- sam - tmax
    sam <- sam[order(-sam)]
    
    colours <- trunc(runif(100, 1, 4))

    N <- 100
    
    A <- c(1.1, 0.3)
    
    K <- c(100,150)

    div_times <- c(-25, -40, -Inf)
    div_cols <- c(1, 2, 3)

    rates <- list(function (s) half_log.rate(s, K[1], A[1], div_times[1]), function (s) half_log.rate(s, K[2], A[2], div_times[2]), function (s) constant.rate(s, N))
    rate.ints <- list(function(t,s) half_log.rate.int(t, s, K[1], A[1], div_times[1]), function(t,s) half_log.rate.int(t, s, K[2], A[2],div_times[2]), function(t,s) constant.rate.int(t,s,N))

    co <- structured_coal.simulate(sam, colours, div_times, div_cols, rates, rate.ints)
    print(co)

    tree_str <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, co$div_from)
