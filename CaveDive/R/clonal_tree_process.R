#' @export
clonal_tree_process.simulate_params <- function(n_exp, n_tips, concentration, sampling_times, K_mean, K_sd, sampling_scale=c(0,20), r_mean=0, r_sd_mult=1, time_mean_sd_mult=4) {
    
    sam <- runif(n_tips, sampling_scale[1], sampling_scale[2])
    sam <- sam - max(sam)

    colouring <- c()
    expansion_probs <- c()

    max_it <- 40

    it <- 0
    while(length(unique(colouring)) < n_exp) {
        if (it > max_it) {
            warning("Maximum sampling iterations exceeded.")
            return (NA)
        }
        expansion_probs <- rdirichlet(1, rep(concentration, n_exp))
        colouring <- rmultinom(n_tips, 1, expansion_probs)
        colouring <- sapply(c(1:n_tips), function (i) which(colouring[,i]>0))
        it <- it + 1
    }


    N <- rlnorm(1, meanlog = K_mean, sdlog = K_sd)
    K <- c()
    A <- c()
    div_times <- c()

    if (n_exp > 1) {
         clade_sizes <- sapply(c(1:n_exp), function(i) length(which(colouring==i)))
         K <- rlnorm(n_exp-1, meanlog = K_mean, sdlog = K_sd)
         A <- rlnorm(n_exp-1, meanlog = r_mean, sdlog = r_sd)
     
         mean_time_to_MRCA <- function(n, N) N*2*(1-1/n)
         var_time_to_MRCA <- function(n,N) 4*N*sum(sapply(c(2:n), function(j) 1/((j**2)*((j-1)**2))))
     
         div_times <- rep(-1, n_exp-1)
         
         it <- 0
         for (i in c(1:(n_exp-1))) {
            
            ### centre at expected TMRCA - time_mean_sd_mult * SD
            time_var <- r_sd_mult*var_time_to_MRCA(clade_sizes[i],K[i])
            time_mean <- mean_time_to_MRCA(clade_sizes[i],K[i])-time_mean_sd_mult*sqrt(var_time_to_MRCA(clade_sizes[i],K[i]))
            time_shape <- (time_mean**2)/(time_var**2)    
            time_rate <- time_mean/(time_var**2)


            ##rejection criterion: within 1SD of expected TMRCA
            while (div_times[i] < min(sam) && 
                    div_times[i] < (mean_time_to_MRCA(clade_sizes[i],K[i]) - 1*sqrt(var_time_to_MRCA(clade_sizes[i],K[i])))) {
                if (it > max_it) {
                    warning("Maximum sampling iterations exceeded.")
                    return (NA)
                }

                div_times[i] <- min(sampling_times)-rgamma(1, shape=time_shape, scale = 1/time_rate)
                it <- it + 1
            }
         }
        div_ord <- order(-div_times)
        div_times <- div_times[div_ord]
        A <- A[div_ord]
        K <- K[div_ord]
        expansion_probs <- expansion_probs[div_ord]

        div_ord <- c(div_ord, n_exp)
        colouring <- sapply(colouring, function(i) div_ord[i])

    } 

    div_times <- c(div_times, -Inf)
    div_cols <- c(1:(n_exp))

    return(list(n_exp=n_exp, N=N, K=K, A=A, tip_colours=colouring, div_times=div_times, div_cols=div_cols, exp_probs=expansion_probs))
} 

#' @export'
clonal_tree_process.simulate_tree <- function(n_exp=n_exp, N=N, K=K, A=A, tip_colours=colouring, div_times=div_times, div_cols=div_cols, exp_probs=expansion_probs) {

    rates <- lapply(c(1:(n_exp-1)), function(i) return(function (s) sat.rate(s, K[i], A[i], div_times[i])))
    rates[[n_exp]] <- function (s) constant.rate(s, N)

    rate.ints <- lapply(c(1:(n_exp-1)), function(i) return(function (t, s) sat.rate.int(t, s, K[i], A[i], div_times[i])))
    rate.ints[[n_exp]] <- function(t, s) constant.rate.int(t, s, N)

    co <- structured_coal.simulate(sampling_times, colouring, div_times, div_cols, rates, rate.ints)

    full_lh <- co$log_lh + 
               sum(sapply((div_cols), function(i) length(which(colouring==i)) * log(expansion_probs[i])))

    return(co=co, full_lh=full_lh)
}