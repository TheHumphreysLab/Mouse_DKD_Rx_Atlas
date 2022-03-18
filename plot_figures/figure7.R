### umap plot and slc5a2 expression plot ######
x0<-min(pt@reductions$umap@cell.embeddings[,1])
library(RColorBrewer)
library(Seurat)
library(ggplot2)
load("pct_subclust.Rda")
x1<-x0+(max(pt@reductions$umap@cell.embeddings[,1])-min(pt@reductions$umap@cell.embeddings[,1]))/8
y0<-min(pt@reductions$umap@cell.embeddings[,2])
y1<-y0+(max(pt@reductions$umap@cell.embeddings[,2])-min(pt@reductions$umap@cell.embeddings[,2]))/8
tiff('pt_subtypes.tiff', units="in", width=5, height=7, res=300)
DimPlot(pt, label = T, cols = c("darkolivegreen1",'cornflowerblue'),label.box = T, label.size = 6, pt.size = 0.0001, raster = F)+NoLegend()+NoAxes()+
  geom_segment(aes(x=x0, xend = x1 , y=y0, yend = y0), size=1,
               arrow = arrow(length = unit(0.2,"cm"))) +
  geom_segment(aes(x=x0, xend = x0 , y=y0, yend = y1), size=0.8,
               arrow = arrow(length = unit(0.2,"cm"))) +
  xlab("UMAP_1")+theme(axis.title.x = element_text(hjust = 0.05, size = 12))+
  ylab('UMAP_2')+theme(axis.title.y = element_text(hjust = 0.05, angle = 90, size = 12))+
  ggtitle("PCT subclustering")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

s1_slc5a1<-FeaturePlot(subset(pt, ident='S1'),features = "Slc5a2")
s1_slc5a1<-s1_slc5a1$data
s1_slc5a1$direction<-ifelse(s1_slc5a1$Slc5a2>1, "pos","neg")
p <- ggplot(s1_slc5a1,aes(x=Slc5a2)) + 
  geom_density(alpha=0.5, fill='cornflowerblue')+
  #geom_density(data = s1_slc5a1[s1_slc5a1$direction=="pos",],aes(x=Slc5a2), fill='red')+
  xlim(-0.5,5)+geom_vline(aes(xintercept=0.5),color="red", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        strip.text = element_text(size = 14),axis.title = element_text(size = 14),
        legend.position = "bottom", legend.text = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 12, face = 'bold'))+
         ylab("Density")+xlab("")+
  geom_text(aes(label= "S1 cells \n 68.0% positive"), x=4, y=1.2)
s2_slc5a1<-FeaturePlot(subset(pt, ident='S2'),features = "Slc5a2")
s2_slc5a1<-s2_slc5a1$data
s2_slc5a1$direction<-ifelse(s2_slc5a1$Slc5a2>1, "pos","neg")
p2 <- ggplot(s2_slc5a1,aes(x=Slc5a2)) + 
  geom_density(alpha=0.5, fill='darkolivegreen1')+
  #geom_density(data = s2_slc5a1[s2_slc5a1$direction=="pos",],aes(x=Slc5a2), fill='red')+
  xlim(-0.5,5)+geom_vline(aes(xintercept=0.5),color="red", linetype="dashed", size=1)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        strip.text = element_text(size = 14),axis.title = element_text(size = 14),
        legend.position = "bottom", legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12, face = 'bold'))+
  xlab("Sglt2 expression")+ylab("Density")+
  geom_text(aes(label= "S2 cells \n 7.4% positive"), x=4, y=7.5)

tiff("sglt2 postive.tiff", width =4, height = 7, units = 'in', res = 300)
p / p2
dev.off()

#### violin plot for Srsf7
sglt2_unique<-SetIdent(pt, value = 'group')
sglt2_unique<-subset(sglt2_unique,idents=c("Grp1","Grp3a","Grp4a","Grp5a",'Grp7a'))
meta<-sglt2_unique@meta.data
meta$group2<-meta$group
meta$group2[meta$group2=='Grp1']<-"db/m"
meta$group2[meta$group2=='Grp3a']<-"AAV"
meta$group2[meta$group2=='Grp4a']<-"ACEi"
meta$group2[meta$group2=='Grp5a']<-"Rosi"
meta$group2[meta$group2=='Grp7a']<-"SGLT2i"
sglt2_unique@meta.data<-meta

plot_violin_pct<-function(seu_obj, feature){
vlnplot<-VlnPlot(seu_obj,features = feature, pt.size  = 0)
vlnplot<-vlnplot$data
names(vlnplot)[1]<-'feature'
vlnplot$group<-meta$group2
vlnplot$group<-factor(vlnplot$group, levels = c("db/m","AAV","ACEi","Rosi","SGLT2i"))
vlnplot$celltype<-meta$celltype2
p1<-ggplot(vlnplot, aes(x = group, y = feature, fill = group)) + 
  facet_wrap(~ celltype, nrow = 2, dir = "v") +
  geom_violin(scale = "width",show.legend = FALSE,lwd=0.1) +
  scale_fill_manual(values = c("magenta",'darkgoldenrod4','#58ACE4','#FF5E6F','#8BD76F')) +
  xlab("") +
  ylab(paste(feature,"expression")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic", size = 12),
        strip.background =element_rect(fill="lemonchiffon1"),
        axis.text.x = element_text(size=10, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 12))
g <- ggplot_gtable(ggplot_build(p1))
strip_t <- which(grepl('strip-t', g$layout$name))
strip_r <- which(grepl('strip-r', g$layout$name))
strip_both<-c(strip_t, strip_r)
fills <- c("darkolivegreen1",'cornflowerblue')
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
print(grid.draw(g))
}

plot_violin_pct(sglt2_unique, feature = "Srsf7")

png('Srsf7.png', units="in", width=2, height=4, res=1200)
plot_violin_pct(sglt2_unique, feature = "Srsf7")
dev.off()
