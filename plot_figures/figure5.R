#### plot deconvolution ###
load("bisque_mousedb_full.Rda")
prop_est<-res$bulk.props
gene_expr<-data.frame(prop_est["Inj.PT",],colnames(prop_est))
gene_expr$group<-metadata$group
colnames(gene_expr)[1:2]<-c("gene","sample")
gene_expr1<-gene_expr[gene_expr$group %in% c("Grp1", "Grp3a","Grp4a", "Grp7a", "Grp8a"),]
gene_expr2<-gene_expr[gene_expr$group %in% c("Grp1", "Grp3b","Grp4b", "Grp7b", "Grp8b"),]
gene_expr1$Time<-"2 days"
gene_expr2$Time<-"2 weeks"
gene_expr<-rbind(gene_expr1, gene_expr2)
gene_expr$Disease<-as.character(gene_expr$group)
gene_expr$Disease[gene_expr$Disease=="Grp1"]<-"db/m"
gene_expr$Disease[gene_expr$Disease %in% c("Grp3a","Grp3b")]<-"AAV"
gene_expr$Disease[gene_expr$Disease %in% c("Grp4a","Grp4b")]<-"ACEi"
gene_expr$Disease[gene_expr$Disease %in% c("Grp7a","Grp7b")]<-"SGLT2i"
gene_expr$Disease[gene_expr$Disease %in% c("Grp8a","Grp8b")]<-"ACEi+SGLT2i"
gene_expr$Disease<-factor(gene_expr$Disease, levels = c("db/m","AAV","ACEi","SGLT2i","ACEi+SGLT2i"))
df2a <- data_summary(gene_expr[gene_expr$Time=="2 days",], varname = 'gene',groupnames='Disease')
df2b <- data_summary(gene_expr[gene_expr$Time=="2 weeks",], varname = 'gene',groupnames='Disease')
df2a$Time<-"2 days"
df2b$Time<-"2 weeks"
df2<-rbind(df2a, df2b)
df2$Disease <- factor(df2$Disease , levels = c("db/m","AAV","ACEi","SGLT2i","ACEi+SGLT2i"))
df2$gene<-round(df2$gene * 100, 2)
tiff('inj_PT_deconvolution.tiff', units="in", width=4, height=6, res=300)
ggplot(df2[df2$Time=="2 weeks",], aes(x = Disease, y = gene, fill=Disease)) + 
  geom_bar(color="gray30", position="dodge", stat = 'identity',width = 0.5, alpha=0.5) + 
  geom_errorbar(data=df2[df2$Time=="2 weeks",], aes(ymin=gene, ymax=gene+sd), width=.2,
                position=position_dodge(.9), color='gray30')+
  geom_quasirandom(data = gene_expr[gene_expr$Time=="2 weeks",], aes(x = Disease, y = gene),width = 0.2, size=0.5) + 
  ylab(paste("Inj.PT","fraction")) + xlab("") + 
  theme(panel.background = element_rect(fill = "white",colour = "gray30"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(size = 14,hjust = 0.5, face = 'bold'),
        axis.text = element_text(size = 12),axis.title=element_text(size=14),
        legend.text=element_text(size=10),
        legend.title = element_text(size = 10),legend.position="none", 
        strip.text = element_text(size = 12))+
  ggtitle(paste("Inj.PT","deconvolution on bulk RNA-seq"))
dev.off()

### plot injPT fraction ###
load("injPT_fraction_merged.Rda")
both_frac<-both_frac[!grepl("Rosi", both_frac$group2),]
both_frac$group2<-as.character(both_frac$group2)
both_frac$group2<-factor(both_frac$group2, levels = c("db/m",  "AAV","ACEi" , "SGLT2i","ACEi+SGLT2i"))
png(filename = 'injPT_fraction.png', width = 3, height = 6,units = 'in', res = 600)
ggplot(both_frac, aes(group2, Frac, fill=group2))+
  geom_bar(position = "dodge",  stat = "summary", fun='mean', width = 0.6)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.2, size=1, color='midnightblue', alpha=0.8)+
  geom_quasirandom(size=1,width = 0.2, color='midnightblue', alpha=0.8)+
  ylab('Percentage of cells')+xlab('')+
  scale_fill_manual(values = brewer.pal(7, 'Set2'))+ggtitle('Inj.PT cell fraction')+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        strip.text = element_text(size = 12),axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        legend.position = "none", axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = 'bold'), 
        plot.title = element_text(hjust =0.5))+ facet_wrap(~type, scales = 'free_y')

### injury scoring ###
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
load('../DEG/injpt_scoring.Rda')
load('../DEG/Inj.PT_week2_deg.Rda')
dn_deg<-cell1_deg[[1]]
dn_deg<-aa[order(dn_deg$avg_log2FC, decreasing = T),]
inj.PT <- AddModuleScore(object = inj.PT, features  = list(rownames(dn_deg)[1:50]), name  = "PTInjuryScore2", random.seed = 1)
metadata<-inj.PT@meta.data
metadata1<-metadata[metadata$group %in% c('Grp1','Grp3b','Grp4b','Grp5b','Grp6b'),]
metadata2<-metadata[metadata$group %in% c('Grp1','Grp3b','Grp4b','Grp7b','Grp8b'),]
metadata1$group2<-metadata1$group
metadata2$group2<-metadata2$group
metadata1$group[metadata1$group=='Grp1']<-'db/m'
metadata1$group[metadata1$group=='Grp3b']<-'AAV'
metadata1$group[metadata1$group=='Grp4b']<-'ACEi'
metadata1$group[metadata1$group=='Grp5b']<-'Rosi'
metadata1$group[metadata1$group=='Grp6b']<-'ACEi+Rosi'
metadata2$group[metadata2$group=='Grp1']<-'db/m'
metadata2$group[metadata2$group=='Grp3b']<-'AAV'
metadata2$group[metadata2$group=='Grp4b']<-'ACEi'
metadata2$group[metadata2$group=='Grp7b']<-'SGLT2i'
metadata2$group[metadata2$group=='Grp8b']<-'ACEi+SGLT2i'
metadata2$group<-factor(metadata2$group, levels = c('db/m','AAV','ACEi','SGLT2i','ACEi+SGLT2i'))
metadata1$group<-factor(metadata1$group, levels = c('db/m','AAV','ACEi','Rosi','ACEi+Rosi'))

p1<-ggplot(data=metadata1) +
  stat_density(aes(x=PTInjuryScore21,y=..density..,color=group, linetype=group), lwd=1,geom="line",position="identity") +
  theme_ipsum()+ 
  scale_color_manual(values = c('green','black', 'orange','dodgerblue','red'), name='Group')+
  scale_linetype_manual(values=c(1,1,6,6,6), name="Group")+xlab('PT injury score')+ylab('Density')+
  theme(legend.position = 'top', axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
p2<-ggplot(data=metadata2) +
  stat_density(aes(x=PTInjuryScore21,y=..density..,color=group, linetype=group), lwd=1,geom="line",position="identity") +
  theme_ipsum()+ 
  scale_color_manual(values = c('green','black', 'orange','dodgerblue','red'), name='Group')+
  scale_linetype_manual(values=c(1,1,6,6,6), name="Group")+xlab('PT injury score')+ylab('Density')+
  theme(legend.position = 'top', axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
tiff('PT injury score2.tiff', units="in", width=8, height=4, res=300)
p1 | p2
dev.off()

metadata$group[metadata$group=='Grp1']<-'db/m'
metadata$group[metadata$group=='Grp3b']<-'AAV'
metadata$group[metadata$group=='Grp4b']<-'ACEi'
metadata$group[metadata$group=='Grp5b']<-'Rosi'
metadata$group[metadata$group=='Grp6b']<-'ACEi+Rosi'
metadata$group[metadata$group=='Grp7b']<-'SGLT2i'
metadata$group[metadata$group=='Grp8b']<-'ACEi+SGLT2i'
metadata$group<-factor(metadata$group, levels = c('db/m','AAV','ACEi','Rosi','SGLT2i','ACEi+Rosi','ACEi+SGLT2i'))
tiff('PT injury score3.tiff', units="in", width=8.5, height=4.5, res=300)
ggplot(data=metadata) +
  stat_density(aes(x=PTInjuryScore21,y=..density..,color=group, linetype=group), lwd=1,geom="line",position="identity") +
  theme_ipsum(base_family = "Arial")+ 
  scale_color_manual(values = c('green','black', 'orange','dodgerblue','purple','red','darkred'), name='Group')+
  scale_linetype_manual(values=c(1,1,3,3,3,6,6), name="Group")+
  xlab('PT injury score')+ylab('Density')+
  theme(legend.position = 'right', axis.title.x = element_text(size=18), 
        axis.text = element_text(size=16),legend.title = element_blank(),
        axis.title.y = element_text(size=18),legend.text = element_text(size = 16))
dev.off()
