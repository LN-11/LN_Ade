plink --vcf snp.vcf --make-bed  --out pop --allow-extra-chr
plink --bfile pop --recode --out pop --allow-extra-chr
###########################################  pca  #####################################  
gcta64 --bfile pop --maf 0.05 --make-grm --out pop_grm  --thread-num 5
gcta64 --grm pop_grm --pca 10 --out pop_pca &
###########################################  admixture  ###############################  
for K in `seq 1 27`; do admixture --cv pop.bed $K|tee log${K}.out; done
grep -h CV log*.out
###########################################  tree  ####################################  
iqtree -s pop.phy -m MFP -nt 20 -bb 1000 --prefix pop.iqtree1
