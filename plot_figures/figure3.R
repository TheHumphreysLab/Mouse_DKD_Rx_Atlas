PlotGeneGroup2<-function(object, feature, splitby ){
  if (is.null(levels(object@meta.data[,splitby]))){
    object@meta.data[,splitby] <-factor(object@meta.data[,splitby], levels = names(table(object@meta.data[,splitby])))
  }
  dataplot<-DotPlot(object, features = feature, split.by =  splitby, cols = material.heat(n=length((levels(object@meta.data[,splitby])))))
  dataplot<-dataplot$data
  dataplot$avg.exp<-scale(dataplot$avg.exp)
  dataplot$Cluster<-gsub( "_.*$", "", dataplot$id )
  dataplot$Disease<-gsub( ".*_", "", dataplot$id )
  dataplot$Disease<-factor(dataplot$Disease, levels = levels(object@meta.data[,splitby]))
  dataplot$Cluster<-factor(dataplot$Cluster, levels = levels(object))
  colnames(dataplot)[1:2]<-c('Avg.Exp', 'Pct.Exp')
  dotplot<-ggplot(dataplot, aes(y = Cluster, x = Disease)) +  geom_tile(fill="white", color="white") +
    geom_point(aes( colour=Avg.Exp, size =Pct.Exp))  +  scale_color_gradientn(colours  =  colorRampPalette(c('grey80','lemonchiffon1','indianred1','darkred'))(255)
    )+ scale_size(range = c(0, 10))+
    theme(axis.line = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(size = 16,hjust = 0.5, face = 'bold'),
          axis.text = element_text(size = 12),axis.title=element_text(size=8),legend.text=element_text(size=8),
          legend.title = element_text(size = 8),legend.position="right")+ylab("")+xlab("")+ggtitle(feature)
}

PlotMultiGeneGroup<-function(object, features, splitby){
  pb <- progress_bar$new(
    format = "  Ploting [:bar] :percent eta: :eta",
    clear = FALSE, total = length(features), width = 100)
  features=rev(features)
  plot_list<-list()
  for(i in 1:length(features)){
    pp<-PlotGeneGroup2(object = object, feature = features[i], splitby = splitby)
    plot_list[[i]]<-pp$data
    pb$tick()
    Sys.sleep(1 / length(features))
  }
  all_data<-do.call('rbind', plot_list)
  dotplot<-ggplot(all_data, aes(x = Disease, y = features.plot)) +  geom_tile(fill="white", color="white") +
    geom_point(aes( colour=Avg.Exp, size =Pct.Exp), alpha=0.9)  +  scale_color_gradientn(colours  =  colorRampPalette(c('grey80','lemonchiffon1','indianred1','darkred'))(255)
    )+ scale_size(range = c(0, 10))+
    theme(axis.line = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(size = 16,hjust = 0.5, face = 'bold'),
          axis.text = element_text(size = 12),axis.title=element_text(size=8),legend.text=element_text(size=8),
          legend.title = element_text(size = 8),legend.position="right",strip.text = element_text(size = 14,colour = 'yellow',face = 'bold'))+ylab("")+xlab("")+ggtitle('')+facet_wrap(~Cluster, ncol = length(levels(object)))
  print(dotplot)
}


up_genes<-c('Opcml', 'Epha6', 'Nkain3','Ankrd1','Cyp24a1',"Cyp27b1" ,'Cyp2d9','Cyp4a10',
            'Havcr1','Cxcl1','Zan','Lmntd1','Klk1','Syt2','Nwd2','Etv4','Ltc4s',
            'Gm28888','Adamts3','Sh3rf3','Pak3','Grem1','Lrrc31','Lrrc36','Tceal9',
            'Btnl9','Sgk1','Tmem45a','Acta2','Fgf2','Tagln','Apoe','Ccr5')

down_genes<-c('Mgat5b', 'Aifm3','Sptssb','Fam169a','Dnase1','Spink1','Slc22a7','Ect2',
              'Snca','Igf1','Nr2e3','Nbl1','Ano4','Dusp15','Nsg2','Unc5d','Grin3a','Sall3',
              'S100g','Scube3','Avpr2','Elf5','Slc26a7','Kcns3','Bmp8a','Avpr1a','Lepr','Jam2',
              'Lrfn2','Efcc1','Ren1','Cenpa','Coro1a','Clec9a')

dn_groups<-subset(janssen_dn, idents=c("Grp1", "Grp2","Grp3a", "Grp3b"))

p <- PlotMultiGeneGroup(object = dn_groups, features = up_genes, splitby = 'group')
g1 <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-t', g$layout$name))
fills<-c('#00288A','#DD001F','#84D000','#00CB47','#947F00','#006234','#FFC393',
         '#560060','#00F1EA','#74FEAE','#0086A5','#FF8DFF','#FF97B6','#C1249A',
         '#FF9E00','#CE7554','#7D807E')
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

png(filename =  'dotplot_mutliple1.png', width = 18, height = 8,units = 'in', res = 600)
grid.draw(g1)
dev.off()



p <- PlotMultiGeneGroup(object = aa, features = down_genes, splitby = 'group')+ 
  scale_color_gradientn(colours  =  colorRampPalette(c('grey80','lightcyan','deepskyblue','blue4'))(255))

g2 <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-t', g$layout$name))
fills<-c('#00288A','#DD001F','#84D000','#00CB47','#947F00','#006234','#FFC393',
         '#560060','#00F1EA','#74FEAE','#0086A5','#FF8DFF','#FF97B6','#C1249A',
         '#FF9E00','#CE7554','#7D807E')
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
png(filename =  'dotplot_mutliple2.png', width = 18, height = 8,units = 'in', res = 600)
grid.draw(g2)
dev.off()

