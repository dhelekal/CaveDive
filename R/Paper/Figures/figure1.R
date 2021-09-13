library(CaveDive)
library(rmutil)
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)
library(gridExtra)
library(RColorBrewer)
library(egg)

set.seed(3)

tips <- 200
priors <- standard_priors(expansion_rate=1, 
                    N_mean_log=4, 
                    N_sd_log=1/2, 
                    t_mid_rate=5, 
                    K_sd_log=1, 
                    exp_time_nu=1/2, 
                    exp_time_kappa=1/4)
given <- list(n_exp=3)

sam <- runif(200,-10, 0)
sam <- sam - max(sam)
sam <- sam[order(-sam)]

sim <- expansions_simulate(priors, sam, 2, given=given)

co <- sim$co
params <- sim$params

set.seed(3)
phy.txt.div <- build_coal_tree.structured(sam, co$times, params$tip_colours, co$colours, params$div_times, params$div_cols, co$div_from)
phy <- read.tree(text=phy.txt.div$full)

divs <- params$div_times
N <- params$N
K <- params$K
t_mid <-params$t_mid
A <- sapply(t_mid, function(x) (1/x)**2)

x_lower <- min(co$times)
x_upper <- 0

rates <- lapply(c(1:3), function(i) return(function (s) sat.rate(s, K[i], A[i], divs[i])))
rates[[4]] <- function (s) constant.rate(s, N)

tree_plt <- plot_structured_tree(phy, 4) + scale_color_brewer(palette="Dark2") +
            scale_fill_brewer(palette="Dark2") +
            coord_flip() +
            scale_x_reverse() +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.line.x = element_blank(), 
                  legend.position = "none", 
                  text = element_text(size=20))

layers <- lapply(tree_plt$layers, function(x) if(class(x$geom)[1] == "GeomPoint") NULL else x)
layers <- layers[!sapply(layers, is.null)]
tree_plt$layers <- layers

pop_fs <- data.frame(x=unlist(lapply(c(1:4), function (i) seq(from=-x_upper, to=-x_lower, length.out=500))),
                     col=unlist(lapply(c(1:4), function (i) rep(i, 500))))

pop_fs$y <- sapply(c(1:nrow(pop_fs)), function(i) 1/rates[[pop_fs$col[i]]](pop_fs$x[i]))

func_plt <- ggplot(pop_fs, aes(x=x,y=y, fill=factor(col), color=factor(col))) + 
            geom_area(alpha=.3, position = position_dodge()) +
            geom_line(size=2) +
            ylim(0, 200) + xlim(-x_upper,-x_lower)
func_plt <- func_plt + scale_fill_brewer(palette="Dark2")  + 
   scale_color_brewer(palette="Dark2")  + 
   scale_y_continuous(position = "right", name="Effective population size") +
   coord_flip() +
   theme_bw() +
   xlab("Time (Years)") +
   theme(legend.position="none", 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(color = 'black'), 
         text = element_text(size=20))

png("fig1.png",width=1280,height=720)
ggarrange(
    func_plt,tree_plt,
    widths = c(2,2),
    heights= c(2))
dev.off()
