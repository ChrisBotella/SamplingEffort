library(glmnet)

RepoDir = '~/LOF article/'
saveDir = '~/LOF article/'

setwd(RepoDir)
source('_functions.R')

originalVes = c('etp','chbio_12','chbio_1','chbio_5','alti','slope','awc_top','bs_top','clc')
predictor_formula = " ~ etp + I(etp^2) + I(chbio_12-etp) + I((chbio_12-etp)^2) + chbio_1 + I(chbio_1^2) + chbio_5 + I(chbio_5^2) + alti + I(alti^2) + slope + I(slope^2) + awc_top + I(awc_top^2) + bs_top + I(bs_top^2) + spht + slope:I(chbio_12-etp)"

nSpecies = 170
scoreThresh = .85

setwd(saveDir)
BG = readRDS(paste('BackgroundPts_lof',nSpecies,'_score',scoreThresh,sep=""))
occ = readRDS(paste('occ_lof',nSpecies,'_score',scoreThresh,sep=""))


colnames(occ)[colnames(occ)=="squareId"]="q_hash"
colnames(occ)[colnames(occ)=="glc19SpId"]="taxa"
occo = occ[,c('q_hash',originalVes,'spht','taxa')]
BGo = BG[,c('q_hash',originalVes,'spht')]

Data = bind.background.pts(occo,BGo,Bias_killer=300,equalSizeCells=T)


#model = lof_glmnet(Data,Data$pseudo_y,sub('~(.*)','\\1' , predictor_formula),weights=Data$w,lambdaMinRatio=1e-8,nlambda=200)
#setwd(saveDir)
#saveRDS(model,paste('model_lof',nSpecies,'_score',scoreThresh,sep=''))

matrixFormu = as.formula(paste(' ~ q_hash + taxa * (',sub('~(.*)','\\1' , predictor_formula),')'))
begin = dim(occo)[1]
Informations = data.frame(nli=seq(begin,dim(Data)[1],by=200000),size=NA,elapsed=NA)
for(nli in ){
  print(nli)
  data = Data[1:nli,]
  data$spht = factor(data$spht)
  data$q_hash = factor(data$q_hash)
  deb=Sys.time()
  SparseDes = sparse.model.matrix(matrixFormu,data)
  dur = Sys.time()-deb
  print('OK')
  a = object.size(SparseDes)
  print(format(a,units='Mb'))

  Informations$size[k] = format(a,units='Mb')
  Informations$elapsed[k] = paste(dur,attr(dur,'units'))
  setwd(saveDir)
  write.table(Informations,'infos.csv',sep=";",row.names=F,col.names=T)
  gc(reset=T)
}


