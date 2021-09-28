# Parameters used to run cellbender on each sample
## 1. Install Cellbender https://github.com/broadinstitute/CellBender
## 2. Run Cellbender on the UMI count matrix from cellranger (in cuda mode)
cellbender remove-background --input cellrangeroutput.raw_feature_bc_matrix.h5 --output cellbender_output.h5 --cuda --expected-cells 10000 --total-droplets-included 30000 --epochs 200 --z-dim 200 --z-layers 1000

