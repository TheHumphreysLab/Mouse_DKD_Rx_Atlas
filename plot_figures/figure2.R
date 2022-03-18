library(plotly)
library(circlize)
library(dplyr)
transform_coordinates<-function(coord_data, zoom){
  center_data<-coord_data-mean(c(min(coord_data),max(coord_data)))
  max_data<-max(center_data)
  new_data<-center_data*zoom/max_data
  return(new_data)
}

add_polar_coord<-function(dat){
  celltypes<-levels(dat$Cluster)
  new_dat<-list()
  for (i in 1:length(celltypes)){
    dat$Cluster<-as.character(dat$Cluster)
    dat1<-dat[dat$Cluster==celltypes[i],]
    dat1$x_polar<-1:nrow(dat1)
    new_dat[[i]]<-dat1
  }
  new_dat<-do.call('rbind', new_dat)
  return(new_dat)
}

get_segment<-function(dat, group){
  dat<-dat[order(dat[,group],decreasing = F), ]
  rownames(dat)<-1:nrow(dat)
  dat<-dat[!duplicated(dat[,group]),]
  dat_seg<-as.integer(rownames(dat))
  return(dat_seg)
}

get_metadata<-function(obj, color){
  metadata<-obj@meta.data
  metadata$Cluster<-obj@active.ident
  metadata$UMAP1<-as.numeric(obj@reductions$umap@cell.embeddings[,1])
  metadata$UMAP2<-as.numeric(obj@reductions$umap@cell.embeddings[,2])
  metadata$x<-transform_coordinates(metadata$UMAP1, zoom = 0.8)
  metadata$y<-transform_coordinates(metadata$UMAP2, zoom = 0.8)
  color_df<-data.frame(Cluster=levels(obj), Colors=color)
  cellnames<-rownames(metadata)
  metadata$cells<-rownames(metadata)
  metadata<-merge(metadata, color_df, by='Cluster')
  rownames(metadata)<-metadata$cells
  metadata<-metadata[cellnames,]
}

title_text <- function(x0, y0, x1, y1, text, rectArgs = NULL, textArgs = NULL) {
  center <- c(mean(c(x0, x1)), mean(c(y0, y1)))
  do.call('rect', c(list(xleft = x0, ybottom = y0, xright = x1, ytop = y1), rectArgs))
  do.call('text', c(list(x = center[1], y = center[2], labels = text), textArgs))
}

mk_marker_ct<-function(dat){
  ori_names<-rownames(dat)
  zero_ct<-dat[rowSums(dat)==0,]
  non_zero<-dat[rowSums(dat)!=0,]
  max_genes<-colnames(non_zero)[max.col(non_zero,ties.method="first")]
  non_zero<-data.frame(cells=rownames(non_zero), genes=max_genes)
  zero_ct<-data.frame(cells=rownames(zero_ct), genes='No_expr')
  all_cells<-rbind(non_zero, zero_ct)
  rownames(all_cells)<-all_cells$cells
  all_cells<-all_cells[ori_names,]
  return(all_cells)
}


col_use<-c('#00288A','#DD001F','#84D000','#00CB47','#947F00','#006234','#FC705B','#FFC393',
                    '#560060','#00F1EA','#74FEAE','#0086A5','#FF8DFF','#FF97B6','#C1249A',
                    '#FF9E00','#CE7554','#7D807E')


cc<-get_metadata(jassen_dn, color = col_use)
cc %>%
  dplyr::group_by(Cluster) %>%
  summarize(x = median(x = x),y = median(x = y)) -> centers
library(MASS)
z <- kde2d(cc$x, cc$y, n=1000)

###### plot the whole figure
circos.clear()
png(filename =  'circular10.png', width = 8, height = 6,units = 'in', res = 600)
par(bg = '#F9F2E4')
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0.01,0),"track.height" = 0.01, gap.degree =c(rep(2,13),12,rep(2,4)))
circos.initialize(sectors =  cc$Cluster, x = cc$x_polar2)

### Tract 2
circos.track(cc$Cluster, cc$x_polar2, y=cc$UMAP_2, bg.border=NA,panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter,
              CELL_META$cell.ylim[2]+ mm_y(4),
              CELL_META$sector.index,
              cex=0.5, col = 'black', facing = "bending.inside", niceFacing = T)
  circos.axis(labels.cex = 0.3, col = 'black', labels.col =  'black')
})
for(i in 1:18){
  dd<-cc[cc$Cluster==celltypes[i],]
  circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), y1 = 0, col = col_use[i],  lwd=3, sector.index = celltypes[i])
}

text(x = 0.165, y=1, labels = "Type", cex = 0.4, col = 'black')

### Tract 3
circos.track(cc$Cluster, cc$x_polar2, y=cc$UMAP_2, bg.border=NA)
col_group<-c('gray90','gray80','gray50','black','pink','red','darkseagreen1','darkgreen','lightcyan','blue','yellow','goldenrod3','plum1','purple')

for(i in 1:length(celltypes)) {
  dd<-cc[cc$Cluster==celltypes[i],]
  dat_seg<-get_segment(dd, group = 'group')
  dat_seg2<-c(dat_seg[-1]-1, nrow(dd))
  scale_factor<-max(dd$x_polar2)/nrow(dd)
  dat_seg<-scale_factor*dat_seg
  dat_seg2<-scale_factor*dat_seg2
  circos.segments(x0 = dat_seg, y0 = 0, x1 = dat_seg2, y1 = 0, col = col_group, sector.index = celltypes[i], lwd=3)
}

text(x = 0.165, y=0.97, labels = "Group", cex = 0.4, col = 'black')

### Tract 4
circos.track(cc$Cluster, cc$x_polar2, y=cc$UMAP_2, bg.border=NA)
for(i in 1:length(celltypes)) {
  dd<-cc[cc$Cluster==celltypes[i],]
  dat_seg<-get_segment(dd, group = 'MouseID')
  dat_seg2<-c(dat_seg[-1]-1, nrow(dd))
  scale_factor<-max(dd$x_polar2)/nrow(dd)
  dat_seg<-scale_factor*dat_seg
  dat_seg2<-scale_factor*dat_seg2
  circos.segments(x0 = dat_seg, y0 = 0, x1 = dat_seg2, y1 = 0, col = col_mouse, sector.index = celltypes[i], lwd=3)
}
text(x = 0.165, y=0.94, labels = "Mouse", cex = 0.4, col = 'black')

### Tract 5
circos.track(cc$Cluster, cc$x_polar2, y=cc$UMAP_2, bg.border=NA)
for(i in 1:length(celltypes)) {
  dd<-cc[cc$Cluster==celltypes[i],]
  scale_factor<-max(dd$x_polar2)/nrow(dd)
  circos.points(x = scale_factor*(min(dd$x_polar):max(dd$x_polar)), y = -1, cex = 0.08, col = dd$col_markers, sector.index = celltypes[i])
}

text(x = 0.165, y=0.91, labels = "Markers", cex = 0.4, col = 'black')


### Graph in the center
points(cc$x,cc$y, pch = 19, col = alpha(cc$Colours,0.2), cex = 0.01,pos = 1);
text(centers$x,centers$y, labels=centers$Cluster, cex = 0.8, col = 'black')
contour(z, drawlabels=F, nlevels=100,levels = c(0.1,0.2),col = '#ae9c76', add=TRUE)

#points(cc$x*0.32+1.2,cc$y*0.32+0.8, pch = 19, col = alpha(cc$col_markers,0.05), cex = 0.003);
#text(centers2$x*0.32+1.2,centers2$y*0.32+0.8, labels=centers2$genes, cex = 0.4, col = 'black',pos = 1)

###EC subtypes
ec_color<-c('#bff542','#83f78f','#EBA1A2','#D70016','#eab3fc','#83b1f7')
ec_meta<-get_metadata(ec.integrated, color = ec_color)

ec_meta %>%
  dplyr::group_by(Cluster) %>%
  summarize(x = median(x = x),y = median(x = y)) -> centers_ec

points(ec_meta$x*0.32-1.2,ec_meta$y*0.32+0.73, pch = 19, col = alpha(ec_meta$Colors,0.5), cex = 0.1);
text(centers_ec$x*0.32-1.2,centers_ec$y*0.32+0.73, labels=centers_ec$Cluster, cex = 0.6, col = 'black')

#fib subtypes
fib_color<-c('#f564df','#fcba03','#1cfc03')

fib_meta<-get_metadata(fib.integrated, color = fib_color)
fib_meta %>%
  dplyr::group_by(Cluster) %>%
  summarize(x = median(x = x),y = median(x = y)) -> centers_fib

points(fib_meta$x*0.32+1.2,fib_meta$y*0.32-0.73, pch = 19, col = alpha(fib_meta$Colors,0.5), cex = 0.1);
text(centers_fib$x*0.32+1.2,centers_fib$y*0.32-0.73, labels=centers_fib$Cluster, cex = 0.6, col = 'black')

#immune subtypes
immune_color<-c('#7ad5ff','#ff7b00','#00ffaa')

immune_meta<-get_metadata(immune, color = immune_color)
immune_meta %>%
  dplyr::group_by(Cluster) %>%
  summarize(x = median(x = x),y = median(x = y)) -> centers_immune

points(immune_meta$x*0.32-1.2,immune_meta$y*0.32-0.73, pch = 19, col = alpha(immune_meta$Colors,0.5), cex = 0.2);
text(centers_immune$x*0.32-1.2,centers_immune$y*0.32-0.73, labels=centers_immune$Cluster, cex = 0.6, col = 'black')

##Loh subtype
loh_color<-c('#c08ff7','#ac6bfa','#d1b0f5','green','#b6cf00','#e1ff00','#f57758')

loh_meta<-get_metadata(loh.integrated, color = loh_color)
loh_meta %>%
  dplyr::group_by(Cluster) %>%
  summarize(x = median(x = x),y = median(x = y)) -> centers_loh

points(loh_meta$x*0.3+1.2,loh_meta$y*0.3+0.73, pch = 19, col = alpha(loh_meta$Colors,0.5), cex = 0.1);
text(centers_loh$x*0.3+1.2,centers_loh$y*0.3+0.73, labels=centers_loh$Cluster, cex = 0.6, col = 'black')


title_text(x0 = -1.35, x1 = -1.05, y0 = -1.06, y1=-1, text = 'Immune cells',
           rectArgs = list(border='#F9F2E4',lwd=0.5),
           textArgs = list(col='black',cex = 1))
title_text(x0 = 1.05, x1 = 1.35, y0 = -1.06, y1=-1, text = 'Fibroblasts',
           rectArgs = list(border='#F9F2E4',lwd=0.5),
           textArgs = list(col='black',cex = 1))

title_text(x0 = -1.35, x1 = -1.05, y0 = 1.06, y1=1, text = 'Endothelia',
           rectArgs = list(border='#F9F2E4',lwd=0.5),
           textArgs = list(col='black',cex = 1))


title_text(x0 = 1.05, x1 = 1.35, y0 = 1.06, y1=1, text = 'TAL cells',
           rectArgs = list(border='#F9F2E4',lwd=0.5),
           textArgs = list(col='black',cex = 1))


#plot group#
#points(cc$x*0.32-1.2,cc$y*0.32+0.8, pch = 19, col = alpha(cc$group_col,0.05), cex = 0.003);
#points(cc$x*0.32-1.2,cc$y*0.32-0.8, pch = 19, col = alpha(cc$ID_col,0.05), cex = 0.003);
lgd_points = Legend(labels = names(table(cc$group3)), type = "points", 
                    title_position = "topleft", 
                    title = "Group",
                    title_gp = gpar(col='black',fontsize = 7, fontface='bold'),
                    legend_gp = gpar(col = col_group),
                    labels_gp = gpar(col='black',fontsize = 5),
                    grid_height = unit(2, "mm"),
                    grid_width = unit(2, "mm"),
                    background = col_group)
draw(lgd_points, x = unit(30, "mm"), y = unit(55, "mm"),
     just = c("right", "bottom"))

lgd_points2 = Legend(labels = cc$genes, type = "points", 
                     title_position = "topleft", 
                     title = "Markers",
                     title_gp = gpar(col='black',fontsize = 7, fontface='bold'),
                     legend_gp = gpar(col = cc$col_markers),
                     labels_gp = gpar(col='black',fontsize = 5),
                     grid_height = unit(2, "mm"),
                     grid_width = unit(2, "mm"),
                     background = cc$col_markers)
draw(lgd_points2, x = unit(190, "mm"), y = unit(50, "mm"),
     just = c("right", "bottom"))

dev.off()

