awk 'NR==1 || $6>0.01 && $6 <0.99' formatted_20180517-UACR_overall-EA-nstud_18-SumMac_400.tbl.rsid > uacr_filtered.txt

python ~/celltype_traits/fast_match.py -b ~/celltype_traits/magma_auxiliary/ref_dat/g1000_eur.bim -bcols 1,0,3 -g uacr_filtered.txt -gcols 2,8,9

magma --annotate window=10,1.5 --snp-loc uacr_filtered.txt.bed --gene-loc ~/celltype_traits/magma_auxiliary/gene_loc/NCBI38.gene.loc --out uacr.annotated_10kbup_15_down

magma --bfile /home/haojiawu/celltype_traits/magma_auxiliary/g1000_eur/g1000_eur --pval uacr_filtered.txt.pval ncol=3 --gene-annot uacr.annotated_10kbup_15_down.genes.annot --out uacr.annotated_10kbup_15_down

uacr=~/celltype_traits/Magma/example/test_UACR/int.annotated_10kbup_15_down.genes.raw
cell_type=cell_type=/mnt/sdc/Janssen_newProcessing/DEG/MAGMA/cell_deg.txt
magma --gene-results  $uacr --set-annot  $cell_type --out uacr
