##R scripts for R plots in figure 3 onwards
library(ggtree)
library(viridis)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(forcats)

###FIG 3a
setwd("//wsl.localhost/ubuntu/home/limsingwei2/aeDNA_phylo/figures/TWISST2")
topofreq <- read.delim("topofreq_strict_nolabel.txt", col.names=c("freq", "topo"), header=FALSE)
topofreq$topo <- factor(topofreq$topo, levels=topofreq$topo)
topofreqplot <- ggplot(topofreq)+
  geom_col(aes(x=topo, y=freq, fill=topo))+
  scale_fill_viridis_d(option="turbo")+
  xlab("Topology no.")+
  ylab("Frequency")+
  theme_bw()+
  theme(legend.position="none")
topofreqplot
ggsave("topofreqplot.jpeg", topofreqplot, device="jpeg", width=27, height=10, units="cm", dpi=1000)

###FIG 3b
setwd("//wsl.localhost/ubuntu/home/limsingwei2/aeDNA_phylo/figures/TWISST2/toptopos")
toptopos_labels <- c("topo80", "topo34", "topo28", "topo88", "topo42", "topo20", "topo52", "topo36", "topo2", "topo54")
toptopos <- list(topo80, topo34, topo28, topo88, topo42, topo20, topo52, topo36, topo2, topo54)

topoplot <- function(toponame){
  tree <- read.tree(paste0(toponame,".nwk"))
  phylo <- ggtree(tree,layout="dendrogram")+
    geom_tiplab(geom = "label",  # labels not text
                label.padding = unit(0.05, "lines"), # amount of padding around the labels
                label.size = 0,
                size = 6)+
    geom_nodelab(geom="label",size=3.5)
  phylo
   }
topo80 <- topoplot("topo80")
topo34 <- topoplot("topo34")
topo28 <- topoplot("topo28")
topo88 <- topoplot("topo88")
topo42 <- topoplot("topo42")
topo20 <- topoplot("topo20")
topo52 <- topoplot("topo52")
topo36 <- topoplot("topo36")
topo2 <- topoplot("topo2")
topo54 <- topoplot("topo54")
topoall <- ggarrange(plotlist=toptopos, labels=toptopos_labels, ncol=5, nrow=2)
topoall
ggsave("topoall.jpeg", topoall, device="jpeg", width=45, height = 10, units = "cm", dpi = 1000)


###FIG4a, 4c
setwd("//wsl.localhost/ubuntu/home/limsingwei2/aeDNA_phylo/figures")
##for plotting phylogenies, just change newick file
tree <- read.tree("Bos.nwk")
phylo <- ggtree(tree,layout="rectangular", branch.length = 0.5)+
  geom_tiplab(geom = "label",  # labels not text
              label.padding = unit(0.05, "lines"), # amount of padding around the labels
              label.size = 0,
              size = 6)+
  geom_nodelab(geom="label",size=5)+
  xlim(0,5)
phylo
ggsave("Bos.lineage.phylo.jpeg", phylo, device="jpeg", width=27, height = 14, units = "cm", dpi = 1000)

###FIG 4b
diagsnps <- read.delim("diagnostic.txt", col.names=c("topo", "node", "no"), header=FALSE)
diagsnps$topo <- factor(diagsnps$topo, levels=diagsnps$topo)
diagplot <- diagsnps %>% mutate(topo=fct_relevel(topo, "topo80", "topo34", "topo28", "topo88", "topo42", "topo20", "topo52", "topo36", "topo2", "topo54")) %>%
  ggplot()+
  geom_col(aes(x=topo,y=no,fill=node))+
  scale_fill_viridis_d(option="turbo")+
  xlab("Topology no.")+
  ylab("Number of diagnostic SNPs")+
  theme_bw()
diagplot
ggsave("diagplot.jpeg", diagplot, device="jpeg", width=20, height=10, units="cm", dpi=1000)



library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggmap)
library(raster)
library(maps)
##FIG 5a
##MAP
setwd("//wsl.localhost/ubuntu/home/limsingwei2/aeDNA_phylo/figures")
theme_set(theme_bw())
register_google("Google maps API")
map_metadata <- read.delim("map_metadata.txt")

world <- c(left = -180, bottom = -90, right = 180, top = 90)
map <- get_map(world,source="google", maptype = "hybrid", color = "color", force=TRUE)
plot_map <- ggmap(map) + 
  geom_jitter(data = map_metadata, aes(x = Longitude, y = Latitude, fill=Location),color="white", size = 3, shape = 23)+
  scale_fill_viridis_d(option="inferno", end=0.7)
plot_map
ggsave("Samples.map1.sat.jpeg", plot_map, device="jpeg", width=27, height = 14, units = "cm", dpi = 1000)


FIG5b
snps_data <- read.delim("no_snps.txt")
snps_data$Status<- as.factor(snps_data$Status)
var.test(SNPs~Status, data=snps_data)
assigned <- snps_data[snps_data$Status=="Assigned",2]
not_assigned <- snps_data[snps_data$Status=="Not assigned",2]
wilcox.test(assigned, not_assigned, alternative = "two.sided")
snps_plot <- ggplot(snps_data, aes(Status,SNPs,fill=Status))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test", aes(label = ..p.signif..),label.x = 1.5,label.y=75000,size=7)+
  ylab("Number of SNPs called")+
  xlab("Lineage assignment")+
  scale_fill_viridis_d(option="magma", begin=0.3, end=0.7)+
  theme_bw()
snps_plot
ggsave("snps_plot.jpeg", snps_plot, device="jpeg", width=20, height = 10, units = "cm", dpi = 1000)


###FIG4c
##iceland map
iceland <- c(left = -45, bottom = 50, right = 5, top = 75)
map <- get_map(iceland,source="google", maptype = "hybrid", color = "color", force=TRUE)
plot_map <- ggmap(map) + 
  geom_point(data = map_metadata, aes(x = Longitude, y = Latitude, fill=Location),color="white", size = 4, shape = 23)+
  scale_fill_viridis_d(option="inferno", end=0.7)+
  geom_label_repel(data=map_metadata, aes(x=Longitude,y=Latitude,label=Lineage,color=Location), size=4, max.overlaps=50,force=20)+
  scale_color_viridis_d(option="inferno", end=0.7)
plot_map
ggsave("iceland.map.all.sat.jpeg", plot_map, device="jpeg", width=27, height = 14, units = "cm", dpi = 1000)



###FIG5d
plot_map <- ggmap(map) + 
  geom_jitter(data = map_metadata, aes(x = Longitude, y = Latitude, fill=Location),color="white", size = 3, shape = 23)+
  scale_fill_viridis_d(option="inferno", end=0.7)+
  geom_label_repel(data=map_metadata, aes(x=Longitude,y=Latitude,label=Lineage,color=Location), size=1.5, max.overlaps=50,force=10)+
  scale_color_viridis_d(option="inferno", end=0.7)
plot_map
ggsave("Samples.map.all.sat.jpeg", plot_map, device="jpeg", width=27, height = 14, units = "cm", dpi = 1000)


####FIG6
plot_map <- ggmap(map) + 
  geom_point(data = map_metadata, aes(x = Longitude, y = Latitude, fill=Location),color="white", size = 4, shape = 23)+
  scale_fill_viridis_d(option="inferno", end=0.7)
plot_map
ggsave("tibet.maponly.sat.jpeg", plot_map, device="jpeg", width=27, height = 14, units = "cm", dpi = 1000)
