## This directory contains scripts for processing the count matrix from CellRanger
1. We first removed the ambient RNA contamination using CellBender. The parameters we used to run Cellbender are specified in the "cellbender" subdirectory.
2. We then perform QC on the count matrix from CellBender output. This includes steps to filter cells and genes, and remove cell doublets. We made a pipeline for the whole process in the subdirectory "data_cleanup".
