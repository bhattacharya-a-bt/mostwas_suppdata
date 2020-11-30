setwd('/pine/scr/a/b/abhattac/TCGA_omics/sims/')

sim_params = expand.grid(c(0.05,0.2),
                         c(0.05,.2,0.5),
                         c(.05,.25),
                         c(1),
                         c(.2,.5,.8))
ppp = sim_params[combo,]
nsims = 1e3
p.c.cis = as.numeric(ppp[1])
p.c.trans = as.numeric(ppp[1])
h2.e.c = as.numeric(ppp[2])
h2.e.t = as.numeric(ppp[3])
nqtl = as.numeric(400)
nM = as.numeric(ppp[4])
ngwas = 1500
h2p = as.numeric(ppp[5])



require(data.table)
require(MOSTWAS)

ref = fread('snps_sims.tsv')
locs = fread('snplocs_sims.tsv')

cis.ref = subset(ref,
                 ref$SNP %in% subset(locs,chr == 'chr6')$snpid)
trans.ref = subset(ref,
                 ref$SNP %in% subset(locs,chr == 'chr14')$snpid)

#### SIMULATE GENOTYPES

sim_geno <- function(geno,n){

  ## estimate LD
  G = as.matrix(geno[,-1])
  require(pbapply)
  mafs = pbapply(G,1,function(x) mean(x)/2)
  G = pbapply(G,1,scale)
  #LD = (t(G) %*% G) / nrow(G) + diag(ncol(G)) * 0.1
  LD = (t(G) %*% G) / nrow(G) + diag(ncol(G))
  L = chol(LD)
  Z = t(L %*% t(matrix(rnorm(ncol(G) * n),ncol = ncol(G))))
  Z = as.matrix(scale(Z))

  return(Z)

}


sim_trait <- function(g,h2g){

  n = length(g)

  if (h2g > 0){

    s2g = var(g)
    s2e = s2g * ((1/h2g)-1)
    e = rnorm(n,mean = 0, sd = sqrt(s2e))
    y = g + e

  }  else {

    e = rnorm(n)
    y = e

  }

  y = as.numeric(scale(y))
  return(y)

}



sim_beta <- function(p.causal, eqtl_h2, n_snps){

  # number of QTLs
  n_qtls = floor(p.causal * n_snps)

  # select which SNPs are causal
  c_qtls = sample(1:n_snps,n_qtls)
  b_qtls = rep(0,n_snps)

  # sample effects from normal prior
  b_qtls[c_qtls] = rnorm(n_qtls,
                         mean = 0,
                         sd = sqrt(eqtl_h2/n_qtls))

  return(b_qtls)

}


sim_eqtl = function(geno.cis,
                    geno.trans,
                    nqtl,
                    b_qtls.cis,
                    b_qtls.trans,
                    eqtl_h2.cis,
                    eqtl_h2.trans,
                    p.causal.trans,
                    numMed){

  ## CIS-HERTIABLE PORTION
  Z_qtl.cis = sim_geno(geno.cis, nqtl)
  n.cis = nrow(Z_qtl.cis)
  p.cis = ncol(Z_qtl.cis)
  # GRM AND LD
  A.cis = (Z_qtl.cis %*% t(Z_qtl.cis)) / p.cis
  LD_qtl.cis = (t(Z_qtl.cis) %*% Z_qtl.cis) / n.cis

  # sim gene expression
  gexpr.cis = sim_trait(Z_qtl.cis %*% b_qtls.cis, eqtl_h2.cis)

  mediator.exp = matrix(nrow = length(gexpr.cis),
                        ncol = numMed)
  Z_qtl.trans = sim_geno(geno.trans, nqtl)
  n.trans = nrow(Z_qtl.trans)
  p.trans = ncol(Z_qtl.trans)
  # GRM AND LD
  A.trans = (Z_qtl.trans %*% t(Z_qtl.trans)) / p.trans
  LD_qtl.trans = (t(Z_qtl.trans) %*% Z_qtl.trans) / n.trans
  # sim gene expression
  for (i in 1:numMed){
    mediator.exp[,i] = sim_trait(Z_qtl.trans %*% b_qtls.trans[,i],
                                 eqtl_h2.trans)
  }
  w_med = rnorm(numMed,mean = 0,sd = sqrt(eqtl_h2.trans/numMed))
  gexpr.trans = c(mediator.exp %*% w_med) +
    rnorm(nrow(mediator.exp),
          mean = 0,
          sd = 1 - eqtl_h2.trans)
  fin.qtls.trans = b_qtls.trans %*% w_med
  gexpr = gexpr.cis + gexpr.trans

  return(list(Z.cis = Z_qtl.cis,
              Z.trans = Z_qtl.trans,
              exp = gexpr,
              med = mediator.exp,
              fin.qtls.trans = fin.qtls.trans))


}

trainCis <- function(pheno,
                     k = 10,
                     snpCur,
                     lab = 'Cis'){

  train = caret::createFolds(y = pheno,
                             k=k,
                             returnTrain = T)
  test = caret::createFolds(y = pheno,
                            k = k,
                            returnTrain = F)

  pred.blup = pred.enet = vector(mode = 'numeric',length = length(pheno))
  for (i in 1:k){

    blup = rrBLUP::mixed.solve(y = pheno[train[[i]]],
                               Z = t(snpCur[,train[[i]]]))
    pred.blup[test[[i]]] = as.numeric(t(snpCur[,test[[i]]]) %*% blup$u)
    enet = glmnet::cv.glmnet(y = pheno[train[[i]]],
                             x = t(snpCur[,train[[i]]]),
                             nfolds = 5)
    pred.enet[test[[i]]] = as.numeric(predict(enet,
                                              newx = t(snpCur[,test[[i]]]),
                                              s = 'lambda.min'))

  }

  model = ifelse(adjR2(pheno,pred.blup) >= adjR2(pheno,pred.enet),
                 'LMM','Elastic net')

  fin.model.enet = glmnet::cv.glmnet(y = pheno,
                                     x = snpCur,
                                     nfolds = 5)
  mod.df.enet = data.frame(SNP = paste0(lab,1:ncol(snpCur)),
                           Effect = as.numeric(coef(fin.model.enet,s = 'lambda.min'))[-1])
  mod.df.enet = subset(mod.df.enet,Effect != 0)

  fin.model.blup = rrBLUP::mixed.solve(y = pheno,
                                       Z = snpCur)
  mod.df.blup = data.frame(SNP = paste0(lab,1:ncol(snpCur)),
                           Effect = as.numeric(fin.model.blup$u))

  if (model == 'Elastic net' & nrow(mod.df.enet) != 0){

    return(list(Model = mod.df.enet,
                Predicted = pred.enet,
                CVR2 = adjR2(pheno,pred.enet)))
  }

  if (model == 'LMM' | nrow(mod.df.enet) == 0){

    return(list(Model = mod.df.blup,
                Predicted = pred.blup,
                CVR2 = adjR2(pheno,pred.blup)))
  }
}

trainMeTWAS_simple <- function(pheno,
                               Z.cis,
                               Z.trans,
                               mediator,
                               k = 10){

  originalCis = trainCis(pheno,
                         k = k,
                         Z.cis,
                         lab = 'Cis')

  medTrainList = apply(mediator,2,trainCis,snpCur = Z.trans,lab = 'Trans')
  names(medTrainList) = paste0('Med',1:length(medTrainList))

  fe.R2 = 0
  if (length(medTrainList) > 0){
    medTrainList = medTrainList[as.numeric(which(sapply(medTrainList,
                                                        function(x) x[3]) >= .01))]
    if (length(medTrainList) > 0){
      fixedEffects = as.data.frame(matrix(ncol = length(medTrainList),
                                          nrow = nrow(mediator)))
      colnames(fixedEffects) = names(medTrainList)
      for (i in 1:ncol(fixedEffects)){
        fixedEffects[,i] = medTrainList[[i]][2]
      }
      fixedEffects = as.data.frame(apply(fixedEffects,2,scale))
      fixedEffects$pheno = pheno
      lmCVFit <- lm(pheno~.,
                    data = fixedEffects)

      fe.R2 = fe.R2 + adjR2(as.numeric(predict(lmCVFit)),pheno)

      trans.mod.df = as.data.frame(abind::abind(lapply(1:length(medTrainList),
                                                       amplifyTrans,
                                                       medTrainList = medTrainList,
                                                       lmCaretObj = lmCVFit),
                                                along = 1))
      trans.mod.df$Effect = as.numeric(as.character(trans.mod.df$Effect))
      trans.mod.df = subset(trans.mod.df,SNP != 'Intercept')
      rownames(trans.mod.df) = NULL
      pheno = pheno - as.numeric(predict(lmCVFit))
    }
  }


  cisGenoMod = trainCis(pheno,
                    k = k,
                    Z.cis,
                    lab = 'Cis')

  cisGenoMod$Model$Mediator = 'Cis'
  if (exists('trans.mod.df')){
    cisGenoMod$Model = rbind(cisGenoMod$Model,trans.mod.df)
  }
  cisGenoMod$CVR2 = cisGenoMod$CVR2 + fe.R2
  cisGenoMod$CVR2.cis.first = originalCis$CVR2
  cisGenoMod$Model.cisOnly = originalCis$Model

  return(cisGenoMod)


}

regress <- function(Z,pheno){

  extractLM = function(x,pheno){

    a = summary(lm(pheno~x))
    return(coef(a)[2,])

  }

  require(pbapply)
  gwas = t(as.data.frame(apply(Z,2,extractLM,pheno = pheno)))
  colnames(gwas) = c('Beta','SE','t','P')
  gwas = as.data.frame(gwas)
  return(gwas[,-3])

}

trainDePMA_simple <- function(pheno,
                              Z.cis,
                              Z.trans,
                              mediator,
                              k = 3){

  originalCis = trainCis(pheno,
                         k = k,
                         Z.cis,
                         lab = 'Cis')

  colnames(Z.cis) = paste0('Cis',1:ncol(Z.cis))
  colnames(Z.trans) = paste0('Trans',1:ncol(Z.trans))
  colnames(mediator) = paste0('Med',1:ncol(mediator))
  transeQTL = regress(Z.trans,pheno)
  transeQTL$SNP = rownames(transeQTL)
  transeQTL$Gene = 'Gene'

  eqtm = apply(mediator,2,regress,Z=Z.trans)
  eqtm = as.data.frame(abind::abind(eqtm,along=1))
  eqtm$SNP = rep(transeQTL$SNP,ncol(mediator))
  eqtm$Med = rep(colnames(mediator),each = ncol(Z.trans))
  transeQTL = subset(transeQTL, P < 0.05)
  eqtm = subset(eqtm, P < .05)

  mediatedSNPs = intersect(transeQTL$SNP,eqtm$SNP)
  Z.trans.this = as.matrix(Z.trans[,mediatedSNPs])
  medTest = pbapply(Z.trans.this,2,MOSTWAS::permuteTME,
              expression = pheno, mediators = mediator,
              nperms = 1000,nc = 1,covs = NULL)
  TME = sapply(medTest,function(x) x[[1]])
  TME.P = sapply(medTest,function(x) x[[2]])
  p_weights = p.adjust(TME.P,'BH')
  include.trans = mediatedSNPs[p_weights < 0.10]
  enet = glmnet::cv.glmnet(y = pheno,
                           x = cbind(Z.cis,Z.trans[,include.trans]),
                           nfolds = 5)
  blup = rrBLUP::mixed.solve(y = pheno,
                             Z = cbind(Z.cis,Z.trans[,include.trans]))
  names = colnames(cbind(Z.cis,Z.trans[,include.trans]))
  Zc = Z.cis
  Zt = Z.trans
  p = pheno
  m = mediator

  train = caret::createFolds(y = pheno,
                             k=k,
                             returnTrain = T)
  test = caret::createFolds(y = pheno,
                            k = k,
                            returnTrain = F)
  pred.enet = pred.blup = vector('numeric',length(pheno))
  for (i in 1:k){
    Z.cis = Zc[train[[i]],]
    Z.trans = Zt[train[[i]],]
    pheno = p[train[[i]]]
    mediator = as.matrix(m[train[[i]],])
    colnames(Z.cis) = paste0('Cis',1:ncol(Z.cis))
    colnames(Z.trans) = paste0('Trans',1:ncol(Z.trans))
    colnames(mediator) = paste0('Med',1:ncol(mediator))
    transeQTL = regress(Z.trans,pheno)
    transeQTL$SNP = rownames(transeQTL)
    transeQTL$Gene = 'Gene'

    eqtm = apply(mediator,2,regress,Z=Z.trans)
    eqtm = as.data.frame(abind::abind(eqtm,along=1))
    eqtm$SNP = rep(transeQTL$SNP,ncol(mediator))
    eqtm$Med = rep(colnames(mediator),each = ncol(Z.trans))
    transeQTL = subset(transeQTL, P < 0.05)
    eqtm = subset(eqtm, P < .05)

    mediatedSNPs = intersect(transeQTL$SNP,eqtm$SNP)
    Z.trans.this = as.matrix(Z.trans[,mediatedSNPs])
    medTest = pbapply(Z.trans.this,2,MOSTWAS::permuteTME,
                      expression = pheno, mediators = mediator,
                      nperms = 1000,nc = 1,covs = NULL)
    TME = sapply(medTest,function(x) x[[1]])
    TME.P = sapply(medTest,function(x) x[[2]])
    p_weights = p.adjust(TME.P,'BH')
    include.trans = mediatedSNPs[p_weights < 0.10]

    this.enet = glmnet::cv.glmnet(y = pheno,
                             x = cbind(Z.cis,Z.trans[,include.trans]),
                             nfolds = 5)
    this.blup = rrBLUP::mixed.solve(y = pheno,
                               Z = cbind(Z.cis,Z.trans[,include.trans]))

    pred.enet[test[[i]]] = predict(this.enet,newx = cbind(Zc[test[[i]],],
                                              Zt[test[[i]],
                                                      include.trans]),
                        s = 'lambda.min')
    pred.blup[test[[i]]] = as.numeric(cbind(Zc[test[[i]],],
                                            Zt[test[[i]],
                                               include.trans]) %*% this.blup$u)

  }

  r2.blup = adjR2(p,pred.blup)
  r2.enet = ifelse(sd(pred.enet) == 0,0,adjR2(p,pred.enet))

  if (r2.blup < r2.enet){
    Model = data.frame(SNP = names,
                       Effect = as.numeric(coef(enet,
                                                s='lambda.min'))[-1])
    Predicted = pred.enet
    Model = subset(Model,Effect!=0)
    if (nrow(Model) <= 1){
      r2.enet = -1
    }

  }
  if (r2.blup >= r2.enet){
    Model = data.frame(SNP = names,
                       Effect = blup$u)
    Predicted = pred.blup
  }

  return(list(Model = Model,
              Predicted = Predicted,
              CVR2 = max(r2.blup,r2.enet),
              CVR2.cis.first = originalCis$CVR2,
              Model.cisOnly = originalCis$Model))

}


sim_gwas <- function(geno.cis,
                     geno.trans,
                     ngwas,
                     b_qtls.cis,
                     b_qtls.trans,
                     var_explained,
                     null = F){

  Z_gwas.cis = sim_geno(geno.cis, ngwas)
  colnames(Z_gwas.cis) = paste0('Cis',1:ncol(Z_gwas.cis))
  Z_gwas.trans = sim_geno(geno.trans, ngwas)
  colnames(Z_gwas.trans) = paste0('Trans',1:ncol(Z_gwas.trans))

  #var_explained only reflects that due to genetics
  gwas_expr = cbind(Z_gwas.cis,Z_gwas.trans) %*% c(b_qtls.cis,b_qtls.trans)

  if (var_explained > 0){
    alpha = rnorm(1)
  }  else {alpha = 0}
  
  if (null){ alpha = 0  }

  y = sim_trait(gwas_expr * alpha, var_explained)
  Z_gwas = cbind(Z_gwas.cis,Z_gwas.trans)
  gwas = regress(Z_gwas,y)
  return(list(gwas = gwas, alpha = alpha, Z_gwas = Z_gwas,y=y))
}


compute_depma = function(gwas, mod, LD,Z_gwas,y,impute = T){

  mod = subset(mod,SNP != '')
  if (!impute){
  gwas = subset(gwas,rownames(gwas) %in% as.character(mod$SNP))
  gwas$Z = gwas$Beta/gwas$SE
  coef = as.matrix(mod$Effect)
  score = t(coef) %*% as.matrix(gwas$Z)
  within_var = as.numeric(t(coef) %*% LD[as.character(mod$SNP),as.character(mod$SNP)] %*% coef)
  Z = score/within_var
  p = 2 * pnorm(-1*abs(Z))

  return(list(score = score,
              within_var = within_var,
              tot.Z = Z,
              p = p))}

  if (impute){
    Z_gwas = Z_gwas[,as.character(mod$SNP)]
    if (nrow(mod) == 1){imp = Z_gwas * mod$Effect
    } else { imp = c(Z_gwas %*% mod$Effect) }
    p = coef(summary(lm(y~imp)))[2,4]
    return(list(p = p))
  }


}

compute_metwas = function(gwas, mod, LD,Z_gwas,y,impute = T){

  mod = subset(mod,SNP != '')
  gwas = subset(gwas,rownames(gwas) %in% as.character(mod$SNP))
  gwas$Z = gwas$Beta/gwas$SE
  require(dplyr)
  a = as.data.frame(mod %>% dplyr::group_by(SNP) %>% dplyr::summarize(sum(Effect)))
  colnames(a) = c('SNP','Effect')
  mod = a[match(rownames(gwas),a$SNP),]

  if (!impute){
  coef = as.matrix(mod$Effect)
  score = t(coef) %*% as.matrix(gwas$Z)
  within_var = as.numeric(t(coef) %*% LD[as.character(mod$SNP),as.character(mod$SNP)] %*% coef)
  Z = score/within_var
  p = 2 * pnorm(-1*abs(Z))

  return(list(score = score,
              within_var = within_var,
              tot.Z = Z,
              p = p))}

  if (impute){
      Z_gwas = Z_gwas[,as.character(mod$SNP)]
      if (nrow(mod) == 1){imp = Z_gwas * mod$Effect
      } else { imp = c(Z_gwas %*% mod$Effect) }
      p = coef(summary(lm(y~imp)))[2,4]
      return(list(p = p))
  }
  }

cis.r2 = metwas.r2 = depma.r2 = vector('numeric',nsims)
twas_sig.metwas = twas_sig.depma = twas_sig.cis = vector('numeric',nsims)
twas_sig.metwas0 = twas_sig.depma0 = twas_sig.cis0 = vector('numeric',nsims)
for (j in 1:nsims){
  print(paste0(j,'/',nsims))
  G = rbind(cis.ref,trans.ref)[,-1]
  require(pbapply)
  mafs = pbapply(G,1,function(x) mean(x)/2)
  G = pbapply(G,1,scale)
  #LD = (t(G) %*% G) / nrow(G) + diag(ncol(G)) * 0.1
  LD = (t(G) %*% G) / nrow(G) + diag(ncol(G))
  colnames(LD) = rownames(LD) =
    c(paste0('Cis',1:nrow(cis.ref)),
      paste0('Trans',1:nrow(trans.ref)))


  b_eqtls.cis = sim_beta(p.causal = p.c.cis,
                       eqtl_h2 = h2.e.c,
                       n_snps = nrow(cis.ref))
  b_eqtls.trans = replicate(nM,sim_beta(p.causal = p.c.trans,
                                            eqtl_h2 = h2.e.t,
                                            n_snps = nrow(trans.ref)))
  exp.sim = sim_eqtl(cis.ref,
                     trans.ref,
                     nqtl,
                     b_eqtls.cis,
                     b_eqtls.trans,
                     eqtl_h2.cis = h2.e.c,
                     eqtl_h2.trans = h2.e.t,
                     p.causal.trans = p.c.trans,
                     numMed = nM)
  metwas = trainMeTWAS_simple(pheno = exp.sim$exp,
                            Z.cis = exp.sim$Z.cis,
                            Z.trans = exp.sim$Z.trans,
                            mediator = exp.sim$med,
                            k = 3)
  depma = trainDePMA_simple(pheno = exp.sim$exp,
                          Z.cis = exp.sim$Z.cis,
                          Z.trans = exp.sim$Z.trans,
                          mediator = as.matrix(exp.sim$med),
                          k = 3)
  cis.r2[j] = min(max(metwas$CVR2.cis.first,
                  depma$CVR2.cis.first),h2.e.c)
  metwas.r2[j] = min(metwas$CVR2,h2.e.c+h2.e.t)
  depma.r2[j] = min(depma$CVR2,h2.e.c+h2.e.t)
  sim.gwas.null = sim_gwas(cis.ref,
                           trans.ref,
                           ngwas,
                           b_qtls.cis = rep(0,length(b_eqtls.cis)),
                           b_qtls.trans = rep(0,length(exp.sim$fin.qtls.trans)),
                           h2p)
  twas_sig.cis[j] = compute_depma(sim.gwas$gwas,metwas$Model.cisOnly,LD,
                                  y = sim.gwas$y,Z_gwas = sim.gwas$Z_gwas,impute = T)$p
  twas_sig.metwas[j] = compute_metwas(sim.gwas$gwas,metwas$Model,LD,
                                      y = sim.gwas$y,Z_gwas = sim.gwas$Z_gwas,impute = T)$p
  twas_sig.depma[j] = compute_depma(sim.gwas$gwas,depma$Model,LD,
                                    y = sim.gwas$y,Z_gwas = sim.gwas$Z_gwas,impute = T)$p

  twas_sig.cis0[j] = compute_depma(sim.gwas.null$gwas,metwas$Model.cisOnly,LD,
                                   y = sim.gwas.null$y,Z_gwas = sim.gwas.null$Z_gwas,impute = T)$p
  twas_sig.metwas0[j] = compute_metwas(sim.gwas.null$gwas,metwas$Model,LD,
                                       y = sim.gwas.null$y,Z_gwas = sim.gwas.null$Z_gwas,impute = T)$p
  twas_sig.depma0[j] = compute_depma(sim.gwas.null$gwas,depma$Model,LD,
                                     y = sim.gwas.null$y,Z_gwas = sim.gwas.null$Z_gwas,impute = T)$p

}
summary = data.frame(p.causal.Cis = p.c.cis,
                     p.causal.Trans = p.c.trans,
                     h2e.cis = h2.e.c,
                     h2e.trans = h2.e.t,
                     nqtl = nqtl,
                     ngwas = ngwas,
                     numMed = nM,
                     h2p = h2p,
                     Mean.cisR2 = mean(cis.r2),
                     SE.cisR2 = sd(cis.r2),
                     Mean.metwasR2 = mean(metwas.r2),
                     SE.metwasR2 = sd(metwas.r2),
                     Mean.depmaR2 = mean(depma.r2),
                     SE.depmaR2 = sd(depma.r2),
                     power.metwas = mean(twas_sig.metwas <  2.5e-6),
                     power.depma = mean(twas_sig.depma <  2.5e-6),
                     power.cis = mean(twas_sig.cis <  2.5e-6),
                     power.metwas0 = mean(twas_sig.metwas0 <  2.5e-6),
                     power.depma0 = mean(twas_sig.depma0 <  2.5e-6),
                     power.cis0 = mean(twas_sig.cis0 <  2.5e-6))
fwrite(summary,'sims_power_r2.tsv',
       append=T,row.names = F,quote=F,sep='\t')

