library(edgeR)
library(ggplot2)
library(Seurat)
library(pheatmap)
library(progress)
load('glom.seurat.Rda')
glom2<-SetIdent(glom, value = 'group')
glom2<-subset(glom2, idents = c('Grp1','Grp3a','Grp4a','Grp5a','Grp6a','Grp7a','Grp8a'))
glom2<-SetIdent(glom2, value = 'celltype')
levels(glom2)<-rev(levels(glom))
glom3<-SetIdent(glom2, value = 'group')
glom3<-subset(glom3, idents = c('Grp1','Grp3a'))
glom3<-SetIdent(glom3, value = 'celltype')
levels(glom3)
cpm_data<-read.csv('Raw_counts_glom.csv') # human glom bulk RNA-seq data
cpm_data<-cpm_data[!duplicated(cpm_data$Gene.symbol),]
rownames(cpm_data)<-cpm_data$Gene.symbol
cpm_data<-cpm_data[,-1]
glom_bulk<-cpm_data
glom_bulk<-glom_bulk[rowSums(glom_bulk>1)>3,]
design<-data.frame(DN=c(rep(1,19),rep(0,20)), LN=c(rep(0,19),rep(1,20)))
dge <- DGEList(counts=glom_bulk)
dge <- calcNormFactors(dge)
dge<-estimateDisp(dge, design = design)
fit<-glmQLFit(dge)
contrast.matrix <- makeContrasts(DN-LN, levels=design)
qlf_male.clk.t<-glmQLFTest(fit, contrast= contrast.matrix)
bulk_deg<-topTags(qlf_male.clk.t, n=31668)
bulk_deg<-bulk_deg$table
bulk_deg$logFDR<-(-log10(bulk_deg$FDR))
bulk_deg$gene<-rownames(bulk_deg)
bulk_deg2<-bulk_deg[bulk_deg$FDR<0.05,]
load('Pod_day2_deg.Rda')
pod<-cell1_deg[[1]]
pod<-pod[pod$p_val<0.05,]
load('PEC_day2_deg.Rda')
pec<-cell1_deg[[1]]
pec<-pec[pec$p_val<0.05,]
load('GEC_day2_deg.Rda')
gec<-cell1_deg[[1]]
gec<-gec[gec$p_val<0.05,]
load('MC_day2_deg.Rda')
mc<-cell1_deg[[1]]
mc<-mc[mc$p_val<0.05,]
all.genes<-c("Tmsb4x"  , "Rasl11a", "Iqgap1",    "Loxl2" ,  "Ggt5" ,
             "Trim7",  "Ntm", "Ms4a2", 
             "Sorbs2"  , "Gata3", "Celf2", "Vegfc",    "Arhgap15" ,"Ptn",  "Unc5b" ,   "Mcam" ,
             "Eln"   ,  "Mgll", "Qsox1"  , "Atp13a3",
             "Enpep"   , "Pak1"  ,"Tmem150c" ,"Sirpa",
             "Pkp4"  ,   "Dapk1", "Slc25a36",
             "Nr4a1" , "Egr1" ,
             "Insr" ,  "Foxo1" , "Ctnna3" )

cpm_data<-cpm(dge)
cpm_data<-data.frame(cpm_data)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
rownames(cpm_data)<-firstup(tolower(rownames(cpm_data)))
hm_palette = colorRampPalette(rev(c('#940063',"#BF0080",  'white',"#029902", "#008000")))(n = 255)
rownames(cpm_data)<-toupper(rownames(cpm_data))
cpm_data<-cpm_data[,39:1]
all.genes2<-toupper(all.genes)
colnames(cpm_data)<-c(paste0("LN",1:20), paste0("DN",1:19))
annotation_col = data.frame(
  group = factor(c(rep("Healthy human glom.",20), rep("Diabetic human glom.", 19)))
)
rownames(annotation_col)<-colnames(cpm_data[all.genes2,])
Var1        <- c("cyan", "dodgerblue4")
names(Var1) <- c("Healthy human glom.","Diabetic human glom.")
anno_colors <- list(group = Var1)

png(filename =  'humanDN_mapping1.png', width = 6, height = 6,units = 'in', res = 600)
pheatmap(as.matrix(cpm_data[all.genes2,]), scale = 'row', 
         cluster_rows = F, cluster_cols = F, 
         color = hm_palette, annotation_col = annotation_col, 
         gaps_col = 20, annotation_colors = anno_colors)
dev.off()
DotPlot(glom3,features =   rev(all.genes[1:20]),split.by = 'group')+RotatedAxis()+coord_flip()

DotPlot(glom3,features =   rev(all.genes[21:32]),split.by = 'group')+RotatedAxis()+coord_flip()

glom3@meta.data$group2<-glom3@meta.data$group
glom3@meta.data$group2[glom3@meta.data$group2=='Grp1']<-'db/m'
glom3@meta.data$group2[glom3@meta.data$group2=='Grp3a']<-'cpm_dataV2d'
glom3@meta.data$group2<-factor(glom3@meta.data$group2, levels = c('db/m','cpm_dataV2d'))

p <- PlotMultiGeneGroup(object = glom3, features = all.genes, splitby = 'group2')+ 
  scale_color_gradientn(colours  =  colorRampPalette(c('grey95','lemonchiffon1','#940063'))(255))+
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())

g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-t', g$layout$name))
fills<-c('#1f77b4','#ff7f0e','#2ca02c','#8c564b')
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
png(filename =  'humanDN_mapping2.png', width = 5, height = 9,units = 'in', res = 600)
grid.draw(g)
dev.off()
hm_palette = colorRampPalette(rev(c('#940063',"#BF0080",  'white',"#029902", "#008000")))(n = 255)
#### construct upset data frame ###
gene_list<-c("Tmsb4x"  , "Rasl11a", "Iqgap1",    "Loxl2" ,  "Ggt5")
treatments<-c('ACEi','Rosi','SGLT2i','ACEi+Rosi','ACEi+SGLT2i')
generate_upsetdata<-function(celltype, genes){
  load(paste0(celltype,'_day2_deg.Rda'))
  order2<-c(2,3,5,4,6)
  upsetdata<-matrix(data = NA, nrow = length(genes), ncol = length(treatments))
  for (i in 1:length(genes)){
    for (j in 1:length(order2)){
      acei<-cell1_deg[[order2[j]]]
      acei_up<-acei[acei$avg_log2FC>0 & acei$p_val<0.05,]
      acei_down<-acei[acei$avg_log2FC<0 & acei$p_val<0.05,]
      gene1<-rownames(acei_up)
      gene2<-rownames(acei_down)
      if(genes[i] %in% gene1){
        upsetdata[i,j]<-1
      } else if(genes[i] %in% gene2) {
        upsetdata[i,j]<-(-1)
      } else {
        upsetdata[i,j]<-0
      }
    }
  }
  upsetdata<-data.frame(upsetdata)
  rownames(upsetdata)<-genes
  colnames(upsetdata)<-treatments
  upsetdata
}

gene_list<-c("Tmsb4x"  , "Rasl11a", "Iqgap1",    "Loxl2" ,  "Ggt5")
gene_listset1<-generate_upsetdata(celltype = 'Pod', genes = gene_list)
gene_list<-c("Trim7",  "Ntm", "Ms4a2")
pec_upset1<-generate_upsetdata(celltype = 'PEC', genes = gene_list)
gene_list<-c("Sorbs2"  , "Gata3", "Celf2", "Vegfc",    "Arhgap15" ,"Ptn",  "Unc5b" ,   "Mcam" )
mc_upset1<-generate_upsetdata(celltype = 'MC', genes = gene_list)
gene_list<- c("Eln"   ,  "Mgll", "Qsox1"  , "Atp13a3")
gec_upset1<-generate_upsetdata(celltype = 'GEC', genes = gene_list)
gene_list<-c("Enpep"   , "Pak1"  ,"Tmem150c" ,"Sirpa")
gene_listset2<-generate_upsetdata(celltype = 'Pod', genes = gene_list)
gene_list<-c("Pkp4"  ,   "Dapk1", "Slc25a36")
pec_upset2<-generate_upsetdata(celltype = 'PEC', genes = gene_list)
gene_list<-c("Nr4a1" , "Egr1" )
mc_upset2<-generate_upsetdata(celltype = 'MC', genes = gene_list)
gene_list<- c("Insr" ,  "Foxo1" , "Ctnna3")
gec_upset2<-generate_upsetdata(celltype = 'GEC', genes = gene_list)
all_upset<-do.call('rbind', list(gene_listset1, pec_upset1, mc_upset1, gec_upset1, 
                                 gene_listset2, pec_upset2, mc_upset2, gec_upset2))
all_upset$gene<-rownames(all_upset)
library(reshape2)
all_upset<-melt(all_upset)
all_upset$gene<-factor(all_upset$gene, levels = rev(all.genes))
all_upset$variable<-factor(all_upset$variable, levels = treatments)
all_upset$value<-as.character(all_upset$value)
##### ploting ##############
cols <- c('0' = "grey80", '1' = '#940063', '-1'="#008000") 
dots1_f<-
  ggplot(all_upset, aes(y=gene, x=variable))+
  geom_point(shape=21, size=6, colour="black", aes(fill=factor(value)))+
  scale_fill_manual(values = cols,labels = c("Down", "n.s",'Up'))+
  theme_minimal()+
  labs(x="", y="")+
  theme(legend.position = "right",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45,vjust = 1, hjust=1))+
  labs(fill='')
png(filename =  'humanDN_mapping3.png', width = 4, height = 6.5,units = 'in', res = 600)
dots1_f
dev.off()

