# ----------------------------------------------------------------------------
# Revision | Integrated UMAP and galaxy density plot
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/2.3.umap&GalaxyPlot.R
# Input    : 1.Data/revision/2.datIntegratedWithBatchEffectsRemoval.RDS
# Output   : 3.Output/revision/2.Integration/formalUMAP.pdf
#            3.Output/revision/2.Integration/galaxyPlot.png
# Upstream : 31_rev_batch_correction_simspec.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)
library(ggh4x)
library(ggrepel)
library(extrafont)
library(ggforce)
library(viridis)
library(ggdark)

dat <- readRDS("1.Data/revision/2.datIntegratedWithBatchEffectsRemoval.RDS")

DimPlot(dat, group.by = "newAnno") + 
  coord_fixed()

df <- data.frame(x = dat[["umap"]]@cell.embeddings[,1],
                 y = dat[["umap"]]@cell.embeddings[,2],
                 spl = dat$samples,
                 lab = dat$newAnno)

dfLab <- df %>% group_by(lab) %>% 
  summarise(x = mean(x), y = mean(y))

# formal UMAP ----

# theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
# circle <- cbind(cos(theta), sin(theta))
# 
# aux <- function(x, one, two, prob=0.8) {
#   if(nrow(x) <= 2) {
#     return(NULL)
#   }
#   sigma <- var(cbind(x[,one], x[,two]))
#   mu <- c(mean(x[,one]), mean(x[,two]))
#   ed <- sqrt(qchisq(prob, df = 2))
#   data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
# }
# 
# ell <- plyr::ddply(df, "lab", aux, one="x", two="y", prob = 0.2)

ggplot(df, aes(x , y, col = lab, fill = lab)) +
  geom_point(shape=21, alpha = 0.6, size = 3) +
  ggrepel::geom_text_repel(data = dfLab, aes(label = lab), 
            col = "black",
            size = 10) +
  # geom_mark_ellipse(data=ell, 
  #                   aes(x=X1, y=X2, label=lab, group=lab),
  #                   inherit.aes = F,
  #                   color = NA,
  #                   label.fill = "transparent",
  #                   label.fontsize = 25,
  #                   label.buffer = unit(15, 'mm'),
  #                   con.cap = 0) +
  theme_void()+ 
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_bw() +
  coord_fixed() +
  theme(legend.position="none",
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=2),
        axis.title = element_text(size = 25, face = "bold"),
        axis.ticks =  element_line(size = 2),
        axis.text = element_text(size = 20, face = "bold"),
        text = element_text(family = "sans", size = 15))

ggsave("3.Output/revision/2.Integration/formalUMAP.pdf", # png
       width = 10.5, height = 9, units = "in")

# Galaxy plot ----
get_density <- function(x, y, ...) {
  density_out <- MASS::kde2d(x, y, ...)
  int_x <- findInterval(x, density_out$x)
  int_y <- findInterval(y, density_out$y)
  comb_int <- cbind(int_x, int_y)
  return(density_out$z[comb_int])
}

ggplot(df, aes(x=x, y=y)) + 
  stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=1), contour=FALSE) + 
  geom_point(shape = ".") +
  stat_density2d(geom="tile", 
                 aes(fill=..density..^0.25,     
                     alpha=ifelse(..density..^0.25<0.1,1,0)), 
                 contour=FALSE) + 
  scale_fill_viridis(option="magma") +
  # geom_mark_ellipse(data=ell, aes(x=X1, y=X2, label=lab, col=lab),
  #                   inherit.aes = F, size = 2, 
  #                   label.fontsize = 19,
  #                   label.buffer = unit(2, 'mm')) +
  xlim(-9.5, 12) +
  ylim(-9, 9) +
  dark_theme_void() +
  coord_fixed() +
  facet_wrap(~spl, ncol = 2) +
  theme(legend.position = "none")

ggsave("3.Output/revision/2.Integration/galaxyPlot.png",
       width = 8.5, height = 7.3, units = "in")

# bar chart ----
data.frame(prop.table(table(dat$reps, dat$newAnno),
                      margin = 2)) %>%  
  # mutate(Var3 = ifelse(Var2 %in% Var2[grep("U", Var2)], "Na", "Cl")) %>% 
  # mutate(Var4 = ifelse(Var2 %in% Var2[grep("C", Var2)], "Fx", "Ifn")) %>% 
  ggplot(aes(x = Freq, y = Var1, 
             fill = Var2)) +
  geom_bar(stat = "identity", alpha = 0.7, 
           position = position_fill(reverse = TRUE)) +
  scale_x_continuous(labels = scales::percent) +
  ylab(NULL) +
  xlab("Percentage") +
  # facet_grid(Var3 + Var4 ~ ., scales = "free_y", space = "free_y") +
  theme_bw()+
  scale_fill_discrete(name = "Cell type") +
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) +
  theme(legend.position="bottom",
        text = element_text(family = "sans", size = 15, face = "bold"))

ggsave("3.Output/revision/2.Integration/barCellDistr.png", 
       width = 8, height = 6,
       units = "in", bg = "white")

