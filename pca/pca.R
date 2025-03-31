library(tidyverse)
library(ggplot2)
library(viridis)
library(ggrepel)
m <- read_table("lineage.txt", col_names = TRUE)
pca_plot <- function(eigenvec, eigenval, metad) { 
  pca <- read_table(eigenvec, col_names = TRUE)
  eigenval <- scan(eigenval)
  # set names
  names(pca)[1] <- "ID"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  metadata <- read_table(metad, col_names = TRUE)
  names(metadata)[1] <- "ID"
  names(metadata)[2] <- "lineage"
  names(metadata)[3] <- "study"
  findmeta <- function(id, meta){
    meta[meta$ID==id, c(2,3)]
  } 
  dat <- data.frame(t(sapply(pca$ID, findmeta, meta=metadata)))
  pca$lineage <- unlist(dat$lineage)
  pca$study <- unlist(dat$study)
  #pca <- as_tibble(data.frame(pca, lineage, study))
  pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)
  b <- ggplot(pca, aes(PC1, PC2, col = factor(lineage), label=ID)) + geom_point(size = 2,alpha=0.5) #+ geom_label(size=0.5)
  b <- b + scale_colour_viridis_d(option="turbo")
  b <- b + coord_equal() + theme_light()
  b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
}
pca_plot("Bos.chr7.eigenvec", "Bos.chr7.eigenval", "../lineage.txt")
onlycows
saveplot <- function(no){
  vec<-paste0("Bos.chr",no,".eigenvec")
  val<-paste0("Bos.chr",no,".eigenval")
  p <- pca_plot(vec, val, "../lineage.txt")
  filename<-paste0("//wsl.localhost/ubuntu/home/limsingwei2/aeDNA_phylo/pca/pca_nohybrids/plots/Bos.chr",no,".pca.jpeg")
  ggsave(filename, p, device="jpeg", width=20, height = 14, units = "cm", dpi = 1000)
}
chrnos <- seq(1,29)
sapply(chrnos, saveplot)
sapply(c("X","Y","M"), saveplot)
