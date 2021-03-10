library(CaveDive)
library(rmutil)
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)
library(gridExtra)

priors <- standard_priors(expansion_rate=1, 
                    N_mean_log=3, 
                    N_sd_log=3, 
                    t_mid_rate=5, 
                    K_sd_log=1, 
                    exp_time_nu=1/2, 
                    exp_time_kappa=1/2)

set.seed(3)

sam <- runif(99,-5, 0)
sam <- sam - max(sam)
tip_cols <- c(rep(1,33), rep(2,33), rep(3,33))
tip_cols <- tip_cols[order(-sam)]
sam <- sam[order(-sam)]

divs <- c(-30,-60, -Inf)
N <- 100
K <- c(50,150)
t_mid <-c(20, 5)
A <- sapply(t_mid, function(x) (1/x)**2)

sim <- expansions_simulate(priors, sam, 1, 
                            given=list(n_exp=2,N=N,K=K,t_mid=t_mid,div_times=divs, tip_colours=tip_cols, div_from=c(3,3)))

co <- sim$co
params <- sim$params
phy.txt <- build_coal_tree.structured(sam, co$times, params$tip_colours, co$colours, params$div_times, params$div_cols, co$div_from, include_div_nodes=TRUE)
phy <- read.tree(text=phy.txt$full)

x_lower <- min(co$times)
x_upper <- 0

rates <- lapply(c(1:2), function(i) return(function (s) sat.rate(s, K[i], A[i], divs[i])))
rates[[3]] <- function (s) constant.rate(s, N)

tree_plt <- plot_structured_tree(phy, 3) + scale_color_viridis(discrete=TRUE, option="plasma") +
            coord_flip() +
            scale_x_reverse() +
            theme(legend.position="right",
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.line.x = element_blank())

pop_fs <- data.frame(x=unlist(lapply(c(1:3), function (i) seq(from=-x_upper, to=-x_lower, length.out=500))),
                     col=unlist(lapply(c(1:3), function (i) rep(i, 500))))

pop_fs$y <- sapply(c(1:nrow(pop_fs)), function(i) 1/rates[[pop_fs$col[i]]](pop_fs$x[i]))

func_plt <- ggplot(pop_fs, aes(x=x,y=y, fill=factor(col), color=factor(col))) + 
            geom_area(alpha=.3, position = position_dodge()) +
            geom_line(size=2) +
            ylim(0, 200) + xlim(-x_upper,-x_lower)
func_plt <- func_plt + scale_fill_viridis(discrete=TRUE, option="plasma") + 
   scale_color_viridis(discrete=TRUE, option="plasma") + 
   scale_y_continuous(position = "right", name="Neg") +
   coord_flip() +
   theme_bw() +
   xlab("Time Before Present") +
   theme(axis.title.x = element_blank(),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.line.x = element_blank(),
         legend.position="none", 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(color = 'black'))

p <- arrangeGrob(
    grobs=list(func_plt,tree_plt),
    nrow=1,
    ncol=2,
    widths = c(5,5),
    heights = 5)

png("fig1.png",width=1600,height=900)
plot(p)
dev.off()
