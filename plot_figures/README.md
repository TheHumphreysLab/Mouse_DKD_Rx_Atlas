# Scripts to reproduce the figures in our manuscript

## 1. figure2.R
#### The R script is for reproducing the figure2 in our manuscript. We wrapped the codes into two functions in our R package plot1cell (https://github.com/HaojiaWu/plot1cell) to help users who want to generate a similar circlize graph for their own dataset. Simply install the package and here is the example to generate the graph:
```
devtools::install_github("HaojiaWu/plot1cell")
library(plot1cell)
data_plot<-prepare_circlize_data(seurat_obj)
plot_circlize(data_plot)
add_track(data_plot, group = 'disease') ### change the 'disease' to whatever column name in your metadata
```
## 2. figure3.R
#### The dotplot graph in figure3 can also be easily generated by plot1cell:
```
library(plot1cell)
genelist<-c("gene1", "gene2"...."geneN")
complex_dotplot_multiple(seu_obj = seu, features = genelist, groupby = "disease") ### change the 'disease' to whatever column name in your metadata
```