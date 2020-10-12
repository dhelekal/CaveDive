library(ggplot2)

func <- function(x, r, K) return(K*r*(20-x)^2/(1+r*(20-x)^2))

k = 10
alphas = c(0.01,0.1,1,5)

XX <- c()
YY <- c()
var <- c()

for (a in alphas) {
    x <- seq(0, 20, by=0.1)
    v <- rep(paste0("r: ", a, ", K: ",k), length(x))
    y <- sapply(x, function (i) func(i,a,k))

    XX <- c(XX, x)
    YY <- c(YY, y)
    var <- c(var, v)
}

df<- data.frame(t=XX, alpha =YY, Parameters=var)

pdf("alpha_plots.pdf", width=5, height=5)
plt<- ggplot(df, aes(x = t, y = alpha)) + 
  geom_line(aes(color = Parameters, linetype = Parameters)) +
  theme_bw() + theme(aspect.ratio=1)
plot(plt)
dev.off()
