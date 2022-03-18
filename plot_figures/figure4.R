generate_deg_mtx<-function(dir_day2,dir_week2){
  load(dir_day2)
  pod_day2<-cell1_deg
  load(dir_week2)
  pod_week2<-cell1_deg
  pod_day2<-pod_day2[1:6]
  names(pod_day2)<-c('dn_d2','acei_d2','rosi_d2','acei_rosi_d2','sglt_d2', 'acei_sglt_d2')
  pod_week2<-pod_week2[1:6]
  names(pod_week2)<-c('dn_w2','acei_w2','rosi_w2','acei_rosi_w2','sglt_w2', 'acei_sglt_w2')
  pod<-c(pod_day2, pod_week2)
  for (i in 1:12){
    aa<-pod[[i]]
    aa<-aa[aa$p_val_adj<0.05,]
    pod[[names(pod)[i]]]<-aa
  }
  treatments<-names(pod)[c(-1,-7)]
  deg_save<-matrix(data = NA,nrow = 10, ncol = 10)
  for (i in 1:length(treatments)){
    if (i<6) {
      day2<-pod[[1]]
      total1<-dim(day2)[1]
      rx<-pod[[treatments[i]]]
      total_tx<-dim(rx)[1]
      dn_up_day2<-rownames(day2[day2$avg_log2FC>0,])
      dn_down_day2<-rownames(day2[day2$avg_log2FC<0,])
      rx_up_day2<-rownames(rx[rx$avg_log2FC>0,])
      rx_down_day2<-rownames(rx[rx$avg_log2FC<0,])
      c1<-length(intersect(dn_up_day2, rx_down_day2))
      c2<-length(dn_up_day2)
      down_rx<-length(rx_down_day2)
      c3<-length(intersect(dn_down_day2, rx_up_day2))
      c4<-length(dn_down_day2)
      up_rx<-length(rx_up_day2)
      c5<-c1+c3
      deg_save[i,]<-c(treatments[i],c1,c2,down_rx,c3,c4,up_rx,c5,total1,total_tx)
    } else {
      week2<-pod[[7]]
      total2<-dim(week2)[1]
      rx<-pod[[treatments[i]]]
      total_tx<-dim(rx)[1]
      dn_up_week2<-rownames(week2[week2$avg_log2FC>0,])
      dn_down_week2<-rownames(week2[week2$avg_log2FC<0,])
      rx_up_week2<-rownames(rx[rx$avg_log2FC>0,])
      rx_down_week2<-rownames(rx[rx$avg_log2FC<0,])
      c1<-length(intersect(dn_up_week2, rx_down_week2))
      c2<-length(dn_up_week2)
      down_rx<-length(rx_down_week2)
      c3<-length(intersect(dn_down_week2, rx_up_week2))
      c4<-length(dn_down_week2)
      up_rx<-length(rx_up_week2)
      c5<-c1+c3
      deg_save[i,]<-c(treatments[i],c1,c2,down_rx,c3,c4,up_rx,c5,total2,total_tx)
    }
  }
  deg_save<-data.frame(deg_save)
  colnames(deg_save)<-c('treatment','up_save','up_dn','down_rx','down_save','down_dn','up_rx','total_save','total_dn','total_rx')
  return(deg_save)
}

setwd('/mnt/sdc/Janssen_newProcessing/DEG')
celltypes<-levels(janssen_dn)
all_deg<-list()
for (i in 1:length(celltypes)){
  cell1<-generate_deg_mtx(dir_day2 = paste0(celltypes[i],'_day2_deg.Rda'), dir_week2 = paste0(celltypes[i],'_week2_deg.Rda'))
  cell1$celltype<-celltypes[i]
  cell1$time<-c(rep('Day2',5), rep("Week2",5))
  all_deg[[i]]<-cell1
  print(celltypes[i])
}
all_deg<-do.call('rbind', all_deg)
aa<-apply(all_deg[,2:10],2, as.integer)
all_deg[,2:10]<-aa
rownames(all_deg)<-all_deg$X
all_deg$percent<-100*(all_deg$total_save/all_deg$total_dn)
data_plt<-all_deg[,c(1,11:13)]

data_plt$celltype<-factor(data_plt$celltype, levels = celltypes)
data_plt$treatment2<-data_plt$treatment
data_plt$treatment2<-gsub('_d2','',data_plt$treatment2)
data_plt$treatment2<-gsub('_w2','',data_plt$treatment2)
data_plt$treatment2[data_plt$treatment2=='acei']<-'ACEi'
data_plt$treatment2[data_plt$treatment2=='rosi']<-'Rosi'
data_plt$treatment2[data_plt$treatment2=='sglt']<-'SGLT2i'
data_plt$treatment2[data_plt$treatment2=='acei_rosi']<-'ACEi+Rosi'
data_plt$treatment2[data_plt$treatment2=='acei_sglt']<-'ACEi+SGLT2i'
data_plt$treatment2<-factor(data_plt$treatment2, levels = c('ACEi','Rosi','SGLT2i','ACEi+Rosi','ACEi+SGLT2i'))

p<-ggplot(data = data_plt, aes(time, percent,fill=treatment2))+
  geom_bar(position="dodge",stat = 'identity',width = 0.6, alpha=0.8)+facet_wrap(~celltype, ncol = 3)+
  ylab('%DEG being corrected by treatment')+xlab('Treatment timepoint')+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),axis.text = element_text(size = 10),axis.title= element_text(size = 14),
        legend.text = element_text(size = 10), legend.title = element_blank(),
        strip.text = element_text(size = 12,colour = 'yellow',face = 'bold'), strip.background = element_rect(color="black", fill='mediumvioletred'))+
  scale_fill_manual(values =  c('#58ACE4','#FF5E6F','#8BD76F','#FFAF5D','purple'))
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
fills<-c('#00288A','#DD001F','#84D000','#00CB47','#947F00','#006234','#FFC393',
         '#560060','#00F1EA','#74FEAE','#0086A5','#FF8DFF','#FF97B6','#C1249A',
         '#FF9E00','#CE7554','#7D807E')
k <- 1
for (i in stripr[-6]) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  print(k)
  k <- k+1
}
tiff(filename = 'DEG corrected by treatment.tiff', width = 10, height = 8,units = 'in', res = 300)
grid.draw(g)
dev.off()
