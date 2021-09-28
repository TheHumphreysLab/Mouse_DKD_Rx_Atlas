awk 'NR==1 || $6>0.01 && $6 <0.99' 20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt > egfr_filtered.txt

python ~/celltype_traits/fast_match.py -b ~/celltype_traits/magma_auxiliary/ref_dat/g1000_eur.bim -bcols 1,0,3 -g ~/celltype_traits/Magma/example/test_eGFR/egfr_filtered.txt -gcols 2,8,9

magma --annotate window=10,1.5 --snp-loc egfr_filtered.txt.bed --gene-loc ~/celltype_traits/magma_auxiliary/gene_loc/NCBI38.gene.loc --out egfr.annotated_10kbup_15_down

magma --bfile /home/haojiawu/celltype_traits/magma_auxiliary/g1000_eur/g1000_eur --pval egfr_filtered.txt.pval ncol=3 --gene-annot egfr.annotated_10kbup_15_down.genes.annot --out egfr.annotated_10kbup_15_down

egfr=~/celltype_traits/Magma/example/test_eGFR/egfr.annotated_10kbup_15_down.genes.raw
cell_type=/mnt/sdc/Janssen_newProcessing/DEG/MAGMA/cell_deg.txt
magma --gene-results  $egfr --set-annot  $cell_type --out egfr
