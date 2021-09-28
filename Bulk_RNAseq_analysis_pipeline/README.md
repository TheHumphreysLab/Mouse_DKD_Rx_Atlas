## a pipeline for RNA-seq analysis using Bash shell and R
fastqc and trim_galore need to be preinstalled before running this pipeline. Here are the instructions for installing these two tools.<br />
fastqc: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt <br />
trim_galore: https://github.com/FelixKrueger/TrimGalore <br />
STAR will also need to be installed and a reference genome for mapping needs to be created. We followed the STAR manual to prepare  the files (https://github.com/alexdobin/STAR).<br />
Finally, the following R packages are required to be installed in R:<br />
```R
all_packages <- c("GenomicRanges","GenomicFeatures","GenomicAlignments","AnnotationDbi","GenomeInfoDb","plyranges")
BiocManager::install(all_packages)
```
Usage: bash run_RNAseq.sh -g REFDIR -w WORKDIR -f FASTQDIR
```sh
bash run_RNAseq.sh -g REFDIR -w WORKDIR -f FASTQDIR
```
