# ----------------------------------------------------------------------------
# Data mining | Circos plot of significant cytokine-receptor pairs
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/6.3.3.gse186898celltalkerCircls.R
# Input    : 3.Output/revision/
# Output   : 3.Output/revision/
# Upstream : 41_gse186898_celltalker.R
# Note     : ggraph
#            adjusted P < 0.05 highlighted with solid black edges
# ----------------------------------------------------------------------------

require(tidyverse)
require(igraph)
require(ggraph)

dat <- read.csv("3.Output/revision/6.Chemokins/GSE186898/3.ligandReceptor/Overall_lr.csv",
                row.names = 1) %>% 
  mutate(allInteraction = paste(cell_type1, "ligand", ligand, 
                                cell_type2, "receptpor", receptpor, 
                                sep ="_"))


edge <- dat[c("ligand", "receptpor")] %>% 
  rename(from = 1, to = 2)

node <- rbind(dat[c("ligand", "cell_type1")] %>% 
  rename(gene = 1, type = 2) %>% 
  mutate(type = paste(type, "ligand", sep = "_")),
  dat[c("receptpor", "cell_type2")] %>% 
    rename(gene = 1, type = 2) %>% 
    mutate(type = paste(type, "receptpor", sep = "_")))

cmbDf <- edge %>% 
  left_join(node %>% rename("from" = 1), by = "from") %>% 
  left_join(node %>% rename("to" = 1), by = "to") %>% 
  rename("fromG" = 1,
         "toG" = 2,
         "fromT" = 3,
         "toT" = 4) %>% # T for type, G for gene
  mutate(from = paste(fromT, fromG, sep = "_"),
         to = paste(toT, toG, sep = "_"),
         type =  paste(fromT, toT, sep = "_"),
         allInteraction = paste(from, to, sep = "_"))

cmbDf$type <- sapply(cmbDf$type, function(x) paste(
  sort(strsplit(x, split = "_")[[1]]), collapse = "_"))


# re-fomate data from graph generation ----
d1 <- data.frame(from="origin", to = unique(c(cmbDf$fromT, cmbDf$toT)))
d2 <- unique(rbind(cmbDf[c("fromT", "from")] %>% rename(from = 1, to = 2),
                   cmbDf[c("toT", "to")] %>% rename(from = 1, to = 2)))

edges <- rbind(d1, d2)

vertices  <-  data.frame(
  name = unique(c(as.character(edges$from), as.character(edges$to)))
) 

vertices$id <- NA
myleaves <- which(is.na( match(vertices$name, edges$from) ))

nleaves <- length(myleaves)
vertices$id[ myleaves ] <- seq(1:nleaves)


vertices <- vertices %>% separate(name, into = c("cell", "type", "gene"),
                                  sep = "_", 
                                  remove = F) %>% 
  arrange(cell)


mygraph <- graph_from_data_frame(edges, vertices = vertices)

from = match(cmbDf$from, vertices$name)
to = match(cmbDf$to, vertices$name)

connDf <- data.frame(from = from, 
                     to = to,
                     fromG = vertices$name[from],
                     toG = vertices$name[to])

connDf$allInteraction = paste(connDf$fromG, connDf$toG, sep = "_")

connDf <- merge(cmbDf[c("type", "allInteraction")], connDf, by = "allInteraction") %>% 
  unique() %>% 
  merge(dat, by = "allInteraction", all = T) %>% 
  drop_na()

# figure ----
# colors_use <- rainbow(n=length(unique(connDf$type)))

p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) 

for (i in 1:length(unique(connDf$allInteraction))) {
  type = unique(connDf$allInteraction)[i]
  p = p + geom_conn_bundle(data = get_con(from = connDf[connDf$allInteraction == type,]$from, 
                                          to = connDf[connDf$allInteraction == type,]$to),
                           aes(x = x * .6, y = y * .6),
                           alpha = ifelse(connDf[connDf$allInteraction == type,]$fdr < 0.05, 1, 0.2),
                           width = log10(connDf[connDf$allInteraction == type,]$interact_ratio),
                           tension = 0.8)
}

p + geom_node_point(aes(x = .64 * x, y = .64 * y,
                        filter = leaf, 
                        col = cell, 
                        shape = type), 
                    size= 5, alpha = 0.6, stroke = 2) +
  geom_node_text(aes(filter = leaf, label=gene, 
                     x = .7 * x, y = .7 * y,
                     angle = ifelse(atan(-(x/y))*(180/pi) < 0,
                                    90 + atan(-(x/y))*(180/pi),
                                    270 + atan(-x/y)*(180/pi)),
                     hjust = ifelse(x > 0, 0 ,1)), 
                 size=4, fontface = 2) +
  scale_shape_manual(label = c("Chemokin", "Receptor"),
                     values=c(19, 1)) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE)) +
  labs(color = "Cell type",
       shape = "") +
  theme_void() +
  theme(
    legend.position="bottom",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) + coord_fixed()

ggsave("3.Output/revision/6.Chemokins/GSE186898/3.ligandReceptor/Circles.png",
       width = 8, height = 8, units = "in", bg = "white")
