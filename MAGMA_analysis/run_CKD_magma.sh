awk 'NR==1 || $6>0.01 && $6 <0.99' CKD_overall_EA_JW_20180223_nstud23.dbgap.txt > ckd_filtered.txt

python ~/celltype_traits/fast_match.py -b ~/celltype_traits/magma_auxiliary/ref_dat/g1000_eur.bim -bcols 1,0,3 -g ckd_filtered.txt -gcols 2,8,9

magma --annotate window=10,1.5 --snp-loc ckd_filtered.txt.bed --gene-loc ~/celltype_traits/magma_auxiliary/gene_loc/NCBI38.gene.loc --out ckd.annotated_10kbup_15_down

magma --bfile /home/haojiawu/celltype_traits/magma_auxiliary/g1000_eur/g1000_eur --pval ckd_filtered.txt.pval ncol=3 --gene-annot ckd.annotated_10kbup_15_down.genes.annot --out ckd.annotated_10kbup_15_down

ckd=~/celltype_traits/Magma/example/test_CKD/ckd.annotated_10kbup_15_down.genes.raw
cell_type=/mnt/sdc/Janssen_newProcessing/DEG/MAGMA/cell_deg.txt
magma --gene-results  $ckd --set-annot  $cell_type --out ckd
