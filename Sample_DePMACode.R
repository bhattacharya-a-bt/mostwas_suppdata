## Define a vector geneList of genes you wish to train
require(bigsnpr)
require(data.table)
require(MOSTWAS)

tempFolder = paste0('temp',i,'/')
if (dir.exists(tempFolder)){
  system(paste0('rm -r ',tempFolder))
}
if (!dir.exists(tempFolder)){
  dir.create(tempFolder)
}
require(bigsnpr)
require(data.table)
require(MOSTWAS)
geneLocs = fread('TCGAgenelocs.txt')
geneLocs = subset(geneLocs,chr %in% c(1:22))

snpObj = snp_attach(snp_readBed('TCGA_tot.bed',
                                backingfile = paste0(tempFolder,
                                                     'temp',i,'_SnpObj.bk')))

mediator = fread('mediators_TCGA_tot.tsv')
exp = fread('inter_mRNA_TCGA_121019.tsv')
colnames(exp) = colnames(mediator)
mediator = rbind(mediator,exp)
mediator = mediator[!duplicated(mediator$Mediator),]
medLocs = fread('TCGAmediatorlocs.txt')
medLocs$right = medLocs$pos + 1
geneLocs = fread('TCGAgenelocs.txt')
colnames(medLocs) = colnames(geneLocs)
medLocs = rbind(geneLocs,medLocs)
medLocs = medLocs[!duplicated(medLocs$geneid),]
rm(geneLocs)
covariates = fread('TCGA_covs.txt')
qtlTra = fread('TCGA_eQTL_dosageTomRNA_part_tra_full.txt')
qtMed = fread('TCGA_QTmeds_dosageToMediator_cis_full.txt')
qtlTra_parts = paste0('TCGA_eQTL_dosageTomRNA_part_tra_',c(1,2,3),'.txt')
qtMed_parts = paste0('TCGA_QTmeds_dosageToMediator_cis_part',c(1,2,3),'.txt')

for (g in geneList){
  
  print(g)
  if (!paste0(g,'.wgt.med.RData') %in% list.files('DePMAModels/')){
    DePMA(geneInt = g,
      snpObj,
      mediator,
      medLocs,
      covariates,
      cisDist = 1e6,
      qtlTra,
      qtMed,
      h2Pcutoff = 1.5,
      dimNumeric = 5,
      verbose = F,
      seed = 1218,
      sobel = F,
      nperms = 1000,
      k = 5,
      parallel = F,
      parType = 'no',
      prune = T,
      ldThresh = .5,
      cores = 5,
      qtlTra_parts,
      qtMed_parts,
      modelDir = 'DePMAModels/',
      tempFolder = tempFolder,
      R2Cutoff = -1)
    fff = paste0(tempFolder,list.files(tempFolder))
    fff = fff[!(fff %in%  paste0(tempFolder,'temp',
                                 i,'_SnpObj.bk',
                                 c('.bk','.rds')))]
    file.remove(fff)
    }
  
}
system(paste0('rm -r ',tempFolder))
