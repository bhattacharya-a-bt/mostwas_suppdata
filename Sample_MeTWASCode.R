### Define a vector geneList with the genes you wish to train
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

snpObj = snp_attach(snp_readBed('TCGA_tot.bed',
		backingfile = paste0(tempFolder,'temp',i,'_SnpObj.bk')))
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
qtlFull = fread('TCGA_QTL_medsTomRNA_top5.tsv')

for (g in geneList){

        print(g)
	if (!paste0(g,'.wgt.med.RData') %in% list.files('MeTWASModels/')){
        MeTWAS(geneInt = g,
               snpObj = snpObj,
               mediator = mediator,
               medLocs = medLocs,
               covariates = covariates,
               dimNumeric = 5,
               qtlFull = qtlFull,
               h2Pcutoff = 1,
               numMed = 5,
               seed = 1218,
               k = 5,
               cisDist = 1e6,
               parallel = F,
               prune = T,
               ldThresh = .5,
               cores = 5,
               verbose = F,
               R2Cutoff = -1,
               modelDir = 'MeTWASModels/',
               tempFolder = tempFolder)
        fff = paste0(tempFolder,list.files(tempFolder))
        fff = fff[!(fff %in%  paste0(tempFolder,'temp',
                                       i,'_SnpObj.bk',
                                c('.bk','.rds')))]
        file.remove(fff)}


        }
system(paste0('rm -r ',tempFolder))
