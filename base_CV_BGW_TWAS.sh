cat ${gene_list_file} | while read gene_name ; do

BGW_dir=/pine/scr/a/b/abhattac/BGW-TWAS # tool directory
GeneExpFile=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/bgw_tcgaexp_block1.txt
geno_dir=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/blocks1
Genome_Seg_Filehead=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/geno_block1_filehead.txt
GTfield=GT # specify genotype field "GT" for genotype
wkdir=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/WorkDir1
LDdir=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/LDdir1
num_cores=8 # number of cores to be used

${BGW_dir}/bin/Step1_get_sumstat.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} --GeneExpFile ${GeneExpFile} \
--geno_dir ${geno_dir} --LDdir ${LDdir} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--GTfield ${GTfield} --num_cores ${num_cores}

${BGW_dir}/bin/Step2_prune.sh --wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--p_thresh 0.001 --max_blocks 100

${BGW_dir}/bin/Step3_EM-MCMC.sh  --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --LDdir ${LDdir} \
--N 368 --hfile ${BGW_dir}/Example/hypval.txt \
--em 5 --burnin 10000 --Nmcmc 10000 \
--PCP_thresh 0.0001 --num_cores 8

rm -rf ${wkdir}/${gene_name}_scores
rm -rf ${wkdir}/${gene_name}_EM_MCMC




BGW_dir=/pine/scr/a/b/abhattac/BGW-TWAS # tool directory
GeneExpFile=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/bgw_tcgaexp_block2.txt
geno_dir=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/blocks2
Genome_Seg_Filehead=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/geno_block2_filehead.txt
GTfield=GT # specify genotype field "GT" for genotype
wkdir=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/WorkDir2
LDdir=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/LDdir2
num_cores=8 # number of cores to be used

${BGW_dir}/bin/Step1_get_sumstat.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} --GeneExpFile ${GeneExpFile} \
--geno_dir ${geno_dir} --LDdir ${LDdir} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--GTfield ${GTfield} --num_cores ${num_cores}

${BGW_dir}/bin/Step2_prune.sh --wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--p_thresh 0.001 --max_blocks 100

${BGW_dir}/bin/Step3_EM-MCMC.sh  --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --LDdir ${LDdir} \
--N 380 --hfile ${BGW_dir}/Example/hypval.txt \
--em 5 --burnin 10000 --Nmcmc 10000 \
--PCP_thresh 0.0001 --num_cores 8

rm -rf ${wkdir}/${gene_name}_scores
rm -rf ${wkdir}/${gene_name}_EM_MCMC



BGW_dir=/pine/scr/a/b/abhattac/BGW-TWAS # tool directory
GeneExpFile=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/bgw_tcgaexp_block3.txt
geno_dir=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/blocks3
Genome_Seg_Filehead=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/geno_block3_filehead.txt
GTfield=GT # specify genotype field "GT" for genotype
wkdir=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/WorkDir3
LDdir=/pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/BGW_files/LDdir3
num_cores=8 # number of cores to be used

${BGW_dir}/bin/Step1_get_sumstat.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} --GeneExpFile ${GeneExpFile} \
--geno_dir ${geno_dir} --LDdir ${LDdir} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--GTfield ${GTfield} --num_cores ${num_cores}

${BGW_dir}/bin/Step2_prune.sh --wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--p_thresh 0.001 --max_blocks 100

${BGW_dir}/bin/Step3_EM-MCMC.sh  --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --LDdir ${LDdir} \
--N 378 --hfile ${BGW_dir}/Example/hypval.txt \
--em 5 --burnin 10000 --Nmcmc 10000 \
--PCP_thresh 0.0001 --num_cores 8

echo ${gene_name} >> /pine/scr/a/b/abhattac/TCGA_omics/TCGA_omics/omics/THESEAREDONE.txt

rm -rf ${wkdir}/${gene_name}_scores
rm -rf ${wkdir}/${gene_name}_EM_MCMC

done