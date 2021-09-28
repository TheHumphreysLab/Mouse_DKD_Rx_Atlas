awk 'NR==1 || $6>0.01 && $6 <0.99' BUN_overall_EA_YL_20171108_METAL1_nstud24.dbgap.txt > bun_filtered.txt

python ~/celltype_traits/fast_match.py -b ~/celltype_traits/magma_auxiliary/ref_dat/g1000_eur.bim -bcols 1,0,3 -g bun_filtered.txt -gcols 2,8,9

magma --annotate window=10,1.5 --snp-loc bun_filtered.txt.bed --gene-loc ~/celltype_traits/magma_auxiliary/gene_loc/NCBI38.gene.loc --out bun.annotated_10kbup_15_down

magma --bfile /home/haojiawu/celltype_traits/magma_auxiliary/g1000_eur/g1000_eur --pval bun_filtered.txt.pval ncol=3 --gene-annot bun.annotated_10kbup_15_down.genes.annot --out bun.annotated_10kbup_15_down

bun=~/celltype_traits/Magma/example/test_BUN/bun.annotated_10kbup_15_down.genes.raw
cell_type=/mnt/sdc/Janssen_newProcessing/DEG/MAGMA/cell_deg.txt
magma --gene-results  $bun --set-annot  $cell_type --out bun
