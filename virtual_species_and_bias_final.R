require(raster)
require(MASS)
require(ggplot2)
require(data.table)
require(plyr)
require(glmnet)

user = "admbotella/Documents/"

RepoDir = paste('C:/Users/',user,"/pCloud local/boulot/Github/SamplingEffort/",sep="")
dataDir = paste('C:/Users/',user,"/pCloud local/boulot/data/",sep="")
saveDir = paste("C:/Users/",user,"/pCloud local/boulot/data/personalWorker/files/",sep="")


setwd(RepoDir)
source('_functions.R')

# Script Sections must be executed in order.
# However, each sections save its output as a checkpoint
# and thus one can restart the script from any section 
# if the previous steps were executed in a previous run

#####
# Define species and occurrences
#####

setwd(saveDir)
occ = read.csv('PL_complete.csv',sep=";",header=T,stringAsFactor=F)               

lonMin = 1.5
lonMax = 8
latMin = 41
latMax = 45
ext = extent(lonMin,lonMax,latMin,latMax)

Min = 50

nSpecies = 50

cd = occ$Longitude>=lonMin & occ$Longitude<=lonMax & occ$Latitude>=latMin & occ$Latitude<=latMax
occo = occ[cd,,drop=F]

counts = table(occo$glc19SpId)
counts = counts[order(counts,decreasing=T)]
counts = counts[1:nSpecies]

setwd(saveDir)
r= readRDS('queries_counts_raster.R')
r_Buf = readRDS('basisGrid')
squareSize=res(r_Buf)[1]

#####
# Generate virtual species niches
#####

enviroVariables = c('alti','chbio_12')

setwd(saveDir)
r2 = readRDS(paste('effortMap_real_size101_H',H,sep=""))
print('effortMap loaded')
vals = as.data.frame(rasterToPoints(r2))
colnames(vals) = c('Longitude','Latitude',"obs")
vals$ve = NA


set.seed(32)

for(ev in enviroVariables){
  print(ev)
  ### Load environmental raster
  setwd(paste(dataDir,'0_mydata/',ev,'/',sep=""))
  r_ev = raster(paste(ev,'.tif',sep=""))
  r_ev = crop(r_ev,ext)
  print('extract enviro values over grid')
  vals$ve = extract(r_ev,vals[,c('Longitude','Latitude')])
  print('generate species params')
  quantiles = quantile(vals$ve,c(.1,.9),na.rm = T)
  
  # optimums positions of species
  mus = runif(nSpecies,quantiles[1],quantiles[2])
  
  # standard deviations of species
  sigs = rgamma(nSpecies,shape=3,scale=50)
  
  # Compute loglinear reparametrization
  df_esp = data.frame(sp_id=1:nSpecies,
                      I.ve = 0,
                      I.ve.2 = 0)
  
  for(i in 1:nSpecies){
    Beta = nauralToLoglinearParam( mus[i] , sigs[i])
    df_esp$I.ve[i] = Beta[2]
    df_esp$I.ve.2[i] = Beta[3]
  }
  
  setwd(saveDir)
  saveRDS(df_esp,paste('df_esp_simu_realistic_',ev,sep=""))
}




#####
# Generate sampling effort with varying roughness
#####


setwd(saveDir)
r= readRDS('queries_counts_raster.R')
r_Buf = readRDS('basisGrid')

refMax = r@data@max
size = 101
Hs = c(20,50,80,100,-1)
for(H in Hs){
  if(H==(-1)){
    print('H: 20 and effort is averaged per sampling cell')
    # Average original sampling effort over LOF quadrats
    setwd(saveDir)
    r_true_effort = readRDS('effortMap_real_size101_H20')
    vals = as.data.frame(rasterToPoints(r_true_effort))
    colnames(vals)[1:2]=c('Longitude','Latitude')
    vals$q_hash = extract(r_Buf,vals[,c('Longitude','Latitude')])
    effort = aggregate(list(obs=vals$layer),by=list(q_hash=vals$q_hash),mean)
    vals= merge(vals,effort,by="q_hash",all.x=T)
    vals = vals[,-which(colnames(vals)=="layer")]
    colnames(vals)[colnames(vals)=="obs"]="layer"
    r2 = rasterFromXYZ(vals[,c('Longitude','Latitude','layer')],crs=CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
  }else{
    print(paste('H:',H))
    W = kde2d(0,0,h=H,n=c(size,size),
              lims=c(-(size-1)/2,(size-1)/2,-(size-1)/2,(size-1)/2))$z
    
    r2 = focal( r , w = W ,pad=T,padValue=0)
    r2 = refMax * r2 / r2@data@max
  }
  setwd(saveDir)
  saveRDS(r2,paste('effortMap_real_size',size,'_H',H,sep=""))
  print('effortMap saved')
  
  vals = as.data.frame(rasterToPoints(r2))
  
  p = ggplot()+ geom_tile(data=vals,aes(x=x,y=y,fill=layer),alpha=.95)
  p = p +scale_fill_gradient(low = "darkorchid4",high="goldenrod")+xlab('Longitude')+ylab('Latitude')
  p = p +theme_bw()+ theme(legend.title = element_text(size=18),
                           legend.text = element_text(size=15),
                           axis.title.x = element_text(size=18),
                           axis.title.y = element_text(size=18),
                           axis.text.x = element_text(size=13),
                           axis.text.y = element_text(size=13))
  
  setwd(saveDir)
  png(paste('effortMap_simu_real_size',size,'_H',H,'.png',sep=""),width=1500,height=1150)
  print(p)
  dev.off()
  
}


######
# Generate species occurrences
######

lonMin = 1.5;lonMax = 8;latMin = 41;latMax = 45
ext = extent(lonMin,lonMax,latMin,latMax)
set.seed(32)

occToSimu = expand.grid(ev=enviroVariables,H=c(-1,20,50,80,100))
size=101
for(i in 1:dim(occToSimu)[1]){
  
  H = occToSimu$H[i]
  ev = occToSimu$ev[i]
  
  print(paste('Load sampling effort with bandwidth:',H))
  setwd(saveDir)
  r = readRDS(paste('effortMap_real_size',size,'_H',H,sep=""))
  vals = as.data.frame(rasterToPoints(r))
  colnames(vals)=c('Longitude','Latitude','obs')
  
  print('load species params')
  df_esp = readRDS(paste('df_esp_simu_realistic_',ev,sep=""))
  nSpecies = dim(df_esp)[1]
  
  print(paste('Load enviro. variable: ',ev))
  setwd(paste(dataDir,'0_mydata/',ev,'/',sep=""))
  r_ev = raster(paste(ev,'.tif',sep=""))
  r_ev = crop(r_ev,ext)
  print('extract enviro values over grid')
  vals$ve = extract(r_ev,vals[,c('Longitude','Latitude')])
  
  ### Simulate species 
  print('simulate occurrences')
  XMat = matrix(NA,dim(vals)[1],2)
  XMat[!is.na(vals$ve),] = model.matrix(~ ve + I(ve^2) - 1 , data=vals)
  
  for(i in 1:nSpecies){
    print(paste('esp',i))
    # Load species abundance fonction parameters
    par = as.numeric(df_esp[i,c('I.ve','I.ve.2')])
    # Vector of values of species abundance function over grid points 
    vals$esp = exp( XMat %*% par )
    # Vector of values of observed points intensity over grid points
    vals$int = vals$obs * vals$esp 
    # Continuously draw points from intensity values (vals$int) gridded (Longitude/Latitude)
    df = DrawPts( vals , counts[i] ) 
    
    df$taxa=as.character(i)
    if(i==1){DF = df}else{DF=rbind(DF,df)}
  }
  
  print('save occurrences')
  DF = DF[!is.na(DF$Longitude) | !is.na(DF$Latitude),,drop=F]
  setwd(saveDir)
  write.table(DF,paste('occurrences_ev_',ev,'_H_',H,'.csv',sep=""),sep=";",row.names=F,col.names=T)
}
  
#####
# Remove squares with less than Min occurrences 
# Count occurrences per cell
#####

lonMin = 1.5;lonMax = 8;latMin = 41;latMax = 45
ext = extent(lonMin,lonMax,latMin,latMax)
set.seed(32)

occToSimu = expand.grid(ev=enviroVariables,H=c(-1,20,50,80,100))
size=101
Min = 50
squareSize = .1
setwd(saveDir)
r_Buf = readRDS('basisGrid')

for(i in 1:dim(occToSimu)[1]){
  
  H = occToSimu$H[i]
  ev = occToSimu$ev[i]
  
  
  print(paste('Load enviro. variable: ',ev))
  setwd(paste(dataDir,'0_mydata/',ev,'/',sep=""))
  r_ev = raster(paste(ev,'.tif',sep=""))
  r_ev = crop(r_ev,ext)
  
  print('Load occurrences')
  setwd(saveDir)
  DF = read.csv(paste('occurrences_ev_',ev,'_H_',H,'.csv',sep=""),sep=";",header=T,stringsAsFactors = F)
  
  DF$q_hash = extract( r_Buf, DF[,c('Longitude','Latitude')]  )
  tab = table(DF$q_hash)
  tab = tab[order(tab,decreasing = T)] 
  #barplot(tab)
  #print(sum(tab>Min)) # Number of kept cells for fitting LOF
  indicesToKeep = as.numeric(names(tab)[tab>Min])
  r_lofCells = r_Buf
  r_lofCells[!r_lofCells[]%in%indicesToKeep] = NA
  
  occ_kept = DF[DF$q_hash%in%indicesToKeep,]
  print(paste('Number of occurrences for final fit:',sum(occ_kept$q_hash%in%indicesToKeep)))
  
  # Get env. variable for each occurrence
  occ_kept$ve = extract(r_ev,occ_kept[,c('Longitude','Latitude')])
  occ_kept = occ_kept[!is.na(occ_kept$ve),,drop=F]
  
  setwd(saveDir)
  saveRDS(r_lofCells,paste('r_lofCells50_size',squareSize,"_Min",Min,'_H',H,'_ev_',ev,sep=""))
  saveRDS(occ_kept,paste('occ_lof50_size',squareSize,"_Min",Min,'_H',H,'_ev_',ev,sep=""))
}

######
# Draw background points uniformly
# over sampling cells
######

MinimumBgPtsPerCell = 10

enviroVariables= c('alti','chbio_12')

setwd(saveDir)
r_Buf = readRDS('basisGrid')
ext = extent(r_Buf)
cellsIds = unique(getValues(r_Buf))
cellsIds = cellsIds[!is.na(cellsIds)]
squareSize = res(r_Buf)[1]
nTmp = 10000
LOF_background = data.frame(x = NA , y = NA , q_hash = NA)
LOF_background = LOF_background[-1,,drop=F]
minNPerSquareLOF = 0
while(minNPerSquareLOF<MinimumBgPtsPerCell){
  tmp = data.frame(x = runif(nTmp,ext[1],ext[2]) , y = runif(nTmp,ext[3],ext[4]) , q_hash = NA)
  # Filter for LOF
  tmp$q_hash = extract(r_Buf,tmp[,1:2])
  cdLOF = !is.na(tmp$q_hash)
  LOF_background = rbind(LOF_background,tmp[cdLOF,,drop=F])
  
  emptyCells = setdiff(cellsIds,unique(LOF_background$q_hash))
  if(length(emptyCells)>0){
    minNPerSquareLOF = 0
  }else{minNPerSquareLOF = min(table(LOF_background$q_hash))}
  
  flush.console()
  cat('     \r  Mini. n° points/LOF square:',minNPerSquareLOF,', tot. n° LOF pts:',dim(LOF_background)[1],'        \r')
}

colnames(LOF_background)[1:2] = c('Longitude','Latitude')

occToSimu = expand.grid(ev=c('alti','chbio_12'),H=c(-1,20,50,80,100))
for(i in 1:dim(occToSimu)[1]){
  tmp = LOF_background
  
  ev = occToSimu$ev[i]
  H = occToSimu$H[i]
  
  print(ev)
  print(H)
  
  setwd(paste(dataDir,'0_mydata/',ev,'/',sep=""))
  r_ev = raster(paste(ev,'.tif',sep=""))
  r_ev = crop(r_ev,ext)
  tmp$ve = extract(r_ev,tmp[,c('Longitude','Latitude')])
  
  setwd(saveDir)
  occ = readRDS(paste('occ_lof50_size',squareSize,"_Min",Min,'_H',H,'_ev_',ev,sep=""))
  tmp = tmp[tmp$q_hash%in%unique(occ$q_hash),,drop=F]
  nona = complete.cases(tmp)
  tmp = tmp[nona,]
  
  setwd(saveDir)
  saveRDS(tmp,paste('BG_lof50_size',squareSize,"_Min",Min,'_ev_',ev,'_H',H,sep=""))
}

#####
# Fit LOF
#####

predictor_formula = " ~ ve + I(ve^2) "
squareSize = .1
Min = 50

occToSimu = expand.grid(ev=c('alti','chbio_12'),H=c(-1,20,50,80,100))

for(i in 1:dim(occToSimu)[1]){
  ev = occToSimu$ev[i]
  H = occToSimu$H[i]
  print(paste('H:',H,' ev:',ev))
  setwd(saveDir)
  occ = readRDS(paste('occ_lof50_size',squareSize,"_Min",Min,'_H',H,'_ev_',ev,sep=""))
  bg = readRDS(paste('BG_lof50_size',squareSize,"_Min",Min,'_ev_',ev,'_H',H,sep=""))
  
  OCC = occ[,c('q_hash','ve','taxa')]
  BG = bg[,c('q_hash','ve')]
  
  # Assemble design matrix
  Data = bind.background.pts(OCC,BG,Bias_killer=1000)
  # Fit model
  lof.mod = lof_glmnet(Data,Data$pseudo_y,sub('~(.*)','\\1' , predictor_formula),weights=Data$w,lambdaMinRatio=1e-12,nlambda=300)
  
  setwd(saveDir)
  saveRDS(lof.mod,paste('model_lof50_size',squareSize,"_Min",Min,'_ev_',ev,'_H',H,sep=""))
}

#####
# Compute sampling effort R2
#####

df = expand.grid(H=c(-1,20,50,80,100),ev=c('alti','chbio_12'),R2=NA,R2_squaresAvg=NA,R2_log=NA,R2_log_squaresAvg=NA)
setwd(saveDir)
for(i in 1:dim(df)[1]){
  print(i)
  H = df$H[i]
  ev = df$ev[i]
  print(paste('H:',H,' ev:',ev))
  
  setwd(saveDir)
  lof.mod = try(readRDS(paste('model_lof50_size0.1_Min50_ev_',ev,'_H',H,sep="")))
  
  if(!is.character(lof.mod)){
    r_lofCells = readRDS( paste('r_lofCells50_size0.1_Min50_H',H,'_ev_',ev,sep=""))
    
    OCC = readRDS(paste('occ_lof50_size0.1_Min50_H',H,'_ev_',ev,sep=""))
    
    r_true_effort = readRDS(paste('effortMap_real_size',size,'_H',H,sep=""))
    vals = as.data.frame(rasterToPoints(r_true_effort))
    colnames(vals)[1:3]=c('Longitude','Latitude','obs')
    
    
    nl = dim(lof.mod$beta)[2]
    coefficients = get.treat.coef(lof.mod,nl)$coefficients
    qCoefs = coefficients[regexpr("q_hash",names(coefficients))>0]
    cells = data.frame(q_hash=paste('q_hash',unique(OCC$q_hash),sep=""))
    coefs = data.frame(coef= qCoefs, q_hash = names(qCoefs))
    cells = merge(cells,coefs,by="q_hash",all.x=T)
    cells$coef[is.na(cells$coef)] = 0
    
    coos = rasterToPoints(r_lofCells)
    colnames(coos)[1:3]=c('Longitude','Latitude','q_hash')
    coos=as.data.frame(coos)
    coos$q_hash = paste('q_hash',coos$q_hash,sep="")
    
    
    cells = merge(cells,coos,by="q_hash",all.x=T)
    
    
    cells$logPred = cells$coef
    cells$Pred = exp(cells$coef)
    
    # correlation of prediction with true effort averaged over quadrats
    vals$q_hash = extract(r_lofCells,vals[,c('Longitude','Latitude')])
    tmp = vals[!is.na(vals$q_hash),,drop=F]
    tmp$q_hash = paste('q_hash',tmp$q_hash,sep="")
    realEffort = aggregate(list(sampEffort=tmp$obs),by=list(q_hash=tmp$q_hash),mean)
    
    cells = merge(cells,realEffort,by="q_hash",all.x=T)
    cells$logTrue = log(cells$sampEffort)
    
    df$R2_log_squaresAvg[i] = cor(cells$logTrue,cells$logPred)^2
    
    df$R2_squaresAvg[i] = cor(cells$sampEffort,cells$Pred)^2
    
    # Plot 
    X= cells$logTrue
    Y = cells$logPred
    #png(paste('realistic_simu_trueQavg_vs_fitted_effort_H',H,'_ev_',ev,'.png',sep=""))
    #plot(X,Y,xlab="Log10-True sampling effort (averaged over quadrats)",ylab="Log10-Predicted sampling effort")
    #text(x=min(X)+0.3*(max(X)-min(X)),y=min(Y)+0.8*(max(Y)-min(Y)),labels=paste('R2=',df$R2_squaresAvg[i]),cex=1)
    #dev.off()
    
    # Correlation of prediction with true effort at fine resolution
    fineRes = merge(tmp,cells[,c('q_hash','logPred','Pred')],by="q_hash",all.x=T)
    X = log(fineRes$obs)
    Y = fineRes$logPred
    df$R2_log[i] = cor(X,Y)^2
    
    # Plot
    #png(paste('realistic_simu_true_vs_fitted_effort_H_',df$H[i],df$veSuffix[i],'.png',sep=""))
    #plot(X,Y,xlab="Log10-True sampling effort",ylab="Log10-Predicted sampling effort")
    #text(x=min(X)+0.3*(max(X)-min(X)),y=min(Y)+0.8*(max(Y)-min(Y)),labels=paste('R2=',df$R2[i]),cex=1)
    #dev.off()  
    
    df$R2[i] = cor(fineRes$obs,fineRes$Pred)^2
    
  }
}

write.table(df,'R2_sampEffort_lof50.csv',sep=";",row.names=F,col.names=T)

setwd(saveDir)
df = read.csv('R2_sampEffort_lof50.csv',sep=";",header=T,stringsAsFactors = F)
df$label = paste('x:',df$ev,' H:',df$H)
df$label = factor(df$label,levels=df$label)

# Log R2
p = ggplot()+geom_point(data=df,aes(x=label,y=R2_log_squaresAvg),colour="black",size=3)+geom_line(data=df,aes(x=label,y=R2_log_squaresAvg,group=ev),colour="black")
p = p + geom_point(data=df,aes(x=label,y=R2_log),colour="red",size=3) + geom_line(data=df,aes(x=label,y=R2_log,group=ev),colour="red")
p = p + theme_bw()+theme(axis.text.x = element_text(angle=-90,vjust=.5))
p = p + scale_y_continuous(limits=c(0,1))
p = p + xlab('Simulation scenario')
p = p + ylab('R2(log-true effort,log-fit effort)')
print(p)

png('effort_R2_log_z.png',height=500,width=700)
print(p)
dev.off()

# R2
p = ggplot()+geom_point(data=df,aes(x=label,y=R2_squaresAvg),colour="black",size=3)+geom_line(data=df,aes(x=label,y=R2_squaresAvg,group=ev),colour="black")
p = p + geom_point(data=df,aes(x=label,y=R2),colour="red",size=3) + geom_line(data=df,aes(x=label,y=R2,group=ev),colour="red")
p = p + theme_bw()+theme(axis.text.x = element_text(size=19,angle=-90),
                         axis.text.y=element_text(size=18),
                         axis.title.y=element_text(size=17),
                         axis.title.x = element_text(size=22))
p = p + scale_y_continuous(limits=c(0,1))
p = p + xlab('Simulation scenario')
p = p + ylab('R2(true effort,est. effort)')
print(p)

png('effort_R2_z.png',height=500,width=700)
print(p)
dev.off()


#####
# Compute species R2
##### 

df=expand.grid(H=c(-1,20,50,80,100),
               ev=c('alti','chbio_12'))
R2s = expand.grid(H=c(-1,20,50,80,100),
                 ev=c('alti','chbio_12'),
                 species=1:50,
                 R2_log_x=NA,
                 R2_log_z=NA,
                 R2_x=NA,
                 R2_z=NA)
setwd(saveDir)
r_effort = readRDS('effortMap_real_size101_H20')
vals = as.data.frame(rasterToPoints(r_effort))
colnames(vals)[1:3]=c('Longitude','Latitude','obs')
regVes=data.frame(init=rep(NA,1000))
for(ev in unique(df$ev)){
  setwd(paste(dataDir,'0_mydata/',ev,'/',sep=""))
  r_ev = raster(paste(ev,'.tif',sep=""))
  r_ev = crop(r_ev,extent(r_effort))
  eval(parse(text=paste('vals$',ev,'= extract(r_ev,vals[,c("Longitude","Latitude")])',sep="")))
  minVe = eval(parse(text=paste('min(vals$',ev,',na.rm=T)',sep="")))
  maxVe = eval(parse(text=paste('max(vals$',ev,',na.rm=T)',sep="")))
  eval(parse(text=paste('regVes$',ev,' = seq(minVe,maxVe,length.out=1000)',sep="")))
}
setwd(saveDir)
for(i in 1:dim(df)[1]){
  print(i)
  H = df$H[i]
  ev = df$ev[i]
  print(paste('H:',H,' ev:',ev))
  
  df_esp = readRDS(paste('df_esp_simu_realistic_',ev,sep=""))
  
  lof.mod = readRDS(paste('model_lof50_size0.1_Min50_ev_',ev,'_H',H,sep=""))
  
  nl = dim(lof.mod$beta)[2]
  coefficients = get.treat.coef(lof.mod,nl)$coefficients
  spCoefs = coefficients[regexpr("ve",names(coefficients))>0]
  
  for(j in 1:dim(df_esp)[1]){
    beta_true = c(df_esp$I.ve[j],df_esp$I.ve.2[j])
    
    spRefLin = spCoefs['ve']
    spRefQua = spCoefs['I(ve^2)']
    cdSp = regexpr(paste('taxa',j,':',sep=""),names(spCoefs))>0
    if(sum(cdSp)){
      spContrLin = spCoefs[cdSp & regexpr(':ve',names(spCoefs))>0]
      spContrQua = spCoefs[cdSp & regexpr('ve\\^2',names(spCoefs))>0]
      beta_fit = c(spRefLin+spContrLin,
               spRefQua+spContrQua)
    }else{
      beta_fit = c(spRefLin,
               spRefQua)
    }
    
    cd = R2s$H==H & R2s$ev==ev & R2s$species==j
    
    # Compute spatial performance metrics
    vecEv = eval(parse(text=paste('vals$',ev,sep="")))
    X = matrix(c(vecEv,vecEv^2),length(vecEv),2)
    trueLogInt = as.vector(X %*% beta_true)
    trueLogInt=trueLogInt[!is.na(trueLogInt)]
    fitLogInt = as.vector(X %*% beta_fit)
    fitLogInt = fitLogInt[!is.na(fitLogInt)]
    
    R2s$R2_log_z[cd] = cor( trueLogInt ,fitLogInt )^2
    R2s$R2_z[cd]= cor(exp(trueLogInt) , exp(fitLogInt))^2 
       
    # Compute environmental performance metrics
    vecEv = eval(parse(text=paste('regVes$',ev,sep="")))
    X = matrix(c(vecEv,vecEv^2),length(vecEv),2)
    trueLogInt = as.vector(X %*% beta_true)
    trueLogInt=trueLogInt[!is.na(trueLogInt)]
    fitLogInt = as.vector(X %*% beta_fit)
    fitLogInt = fitLogInt[!is.na(fitLogInt)]
    
    R2s$R2_log_x[cd] = cor( trueLogInt ,fitLogInt )^2
    R2s$R2_x[cd]= cor(exp(trueLogInt) , exp(fitLogInt))^2
  }
}


write.table(R2s,'R2_species_lof50.csv',sep=";",row.names=F,col.names=T)

R2s$label = paste('x:',R2s$ev,' H:',R2s$H,sep="")
R2s$label = factor(R2s$label,levels=unique(R2s$label))

### Plots
setwd(saveDir)

# R2_x
p = ggplot(R2s[!is.na(R2s$R2_x),],aes(x=label,y=R2_x))+geom_violin(colour="black",width=1.)
p = p + geom_boxplot(width=.1)
p = p + stat_summary(fun.y=mean, geom="point",colour="black",size=2)
p = p + scale_y_continuous(limits=c(0,1))
p = p + xlab('Simulation scenario')
p = p + ylab('R2(enviro. density,est. enviro. density)')
p = p + theme_bw() + theme(axis.text.x = element_text(size=19,angle=-90),
                           axis.text.y=element_text(size=18),
                           axis.title.y=element_text(size=17),
                           axis.title.x = element_text(size=22))
png('species_R2_x.png',height=500,width=700)
print(p)
dev.off()

# R2_z
p = ggplot(R2s[!is.na(R2s$R2_z),],aes(x=label,y=R2_z))+geom_violin(colour="black",width=1.)
p = p + geom_boxplot(width=.1)
p = p + stat_summary(fun.y=mean, geom="point",colour="black",size=2)
p = p + scale_y_continuous(limits=c(0,1))
p = p + xlab('Simulation scenario')
p = p + theme_bw() + theme(axis.text.x = element_text(size=22,angle=-90),
                           axis.text.y=element_text(size=18),
                           axis.title.y=element_text(size=22),
                           axis.title.x = element_text(size=22))

png('species_R2_z.png',height=500,width=700)
print(p)
dev.off()


# R2_log_x
p = ggplot(R2s[!is.na(R2s$R2_log_x),],aes(x=label,y=R2_log_x))+geom_violin(colour="black",width=1.)
p = p + geom_boxplot(width=.1)
p = p + stat_summary(fun.y=mean, geom="point",colour="black",size=2)
p = p + theme_bw()+ theme(axis.text.x = element_text(angle=-90))
print(p)

# R2_log_z
p = ggplot(R2s[!is.na(R2s$R2_log_z),],aes(x=label,y=R2_log_z))+geom_violin(colour="black",width=1.)
p = p + geom_boxplot(width=.1)
p = p + stat_summary(fun.y=mean, geom="point",colour="black",size=2)
p = p + theme_bw()+ theme(axis.text.x = element_text(angle=-90))
print(p)

#####
# Compute corr( Err(log s_x) , Err(Log \hat{s_x})) along enviro
#####

df = expand.grid(H=c(-1,20,50,80,100),ev=c('alti','chbio_12'),R2=NA,R2_squaresAvg=NA,R2_log=NA,uu=NA)
df$ev = as.character(df$ev)

# get points and define environmental breaks
setwd(saveDir)
r_true_effort = readRDS('effortMap_real_size101_H20')
vals = as.data.frame(rasterToPoints(r_true_effort))
colnames(vals)[1:3]=c('Longitude','Latitude','true')
nW = 30
enviroSeq = data.frame(init=rep(NA,nW))
steps = rep(NA,length(unique(df$ev)))
names(steps) = unique(df$ev)
for(ev in unique(df$ev)){
  setwd(paste(dataDir,'0_mydata/',ev,'/',sep=""))
  r_ev = raster(paste(ev,'.tif',sep=""))
  r_ev = crop(r_ev,extent(r_effort))
  eval(parse(text=paste('vals$',ev,'= extract(r_ev,vals[,c("Longitude","Latitude")])',sep="")))
  enviroSeq[,ev] = seq(min(vals[,ev],na.rm=T),max(vals[,ev],na.rm=T),length.out=nW)
  steps[ev] = enviroSeq[2,ev] - enviroSeq[1,ev]
}
enviroSeq = enviroSeq[,-which(colnames(enviroSeq)=="init")]

setwd(saveDir)
for(i in 1:dim(df)[1]){
  print(i)
  H = df$H[i]
  ev = df$ev[i]
  print(paste('H:',H,' ev:',ev))
  
  df_esp = readRDS(paste('df_esp_simu_realistic_',ev,sep=""))
  
  # Get fitted coef per cell
  lof.mod = readRDS(paste('model_lof50_size0.1_Min50_ev_',ev,'_H',H,sep=""))
  r_lofCells = readRDS( paste('r_lofCells50_size0.1_Min50_H',H,'_ev_',ev,sep=""))
  OCC = readRDS(paste('occ_lof50_size0.1_Min50_H',H,'_ev_',ev,sep=""))
  nl = dim(lof.mod$beta)[2]
  coefficients = get.treat.coef(lof.mod,nl)$coefficients
  qCoefs = coefficients[regexpr("q_hash",names(coefficients))>0]
  cells = data.frame(q_hash=paste('q_hash',unique(OCC$q_hash),sep=""))
  coefs = data.frame(coef= qCoefs, q_hash = names(qCoefs))
  cells = merge(cells,coefs,by="q_hash",all.x=T)
  cells$coef[is.na(cells$coef)] = 0
  
  # Get true effort for all points
  tmp = vals
  tmp$true = rasterToPoints(r_true_effort)[,3]
  if(sum(tmp$true==0)>0){tmp$true[tmp$true==0] = min(tmp$true[tmp$true>0])}
  tmp = tmp[!is.na(tmp[,ev]),]
  
  # Get cells ids for all points
  tmp$q_hash = extract(r_lofCells,tmp[,c('Longitude','Latitude')])
  tmp = tmp[!is.na(tmp$q_hash),]
  tmp$q_hash = paste('q_hash',tmp$q_hash,sep="")
  tmp = merge(tmp,cells[,c('q_hash','coef')],by="q_hash",all.x=T)
  
  enviro = enviroSeq
  enviro[,paste('Int_s_logsDiff',ev,sep="")] = NA
  enviro[,paste('Int_s_',ev,sep="")] = NA
  enviro[,paste('mu_x_',ev,sep="")] = NA 
  enviro[,paste('Sum_n_lamb_',ev,sep="")] = NA 
  enviro[,paste('Sum_n_lamb_betaDiff_w_',ev,sep="")] = NA 
  for(j in 1:nW){
    print(j)
    w = enviro[j,ev]
    cdW = tmp[,ev]>=(w-steps[ev]/2) & tmp[,ev]<(w+steps[ev]/2)
    
    # Number of points for w proportional to the area where environment is w (mu_x) 
    enviro[j,paste('mu_x_',ev,sep="")] = sum(cdW)
    # 
    enviro[j,paste('Int_s_',ev,sep="")] = sum(tmp$true[cdW],na.rm=T)
    # 
    enviro[j,paste('Int_s_logsDiff',ev,sep="")] = sum(tmp$true[cdW]*(log(tmp$true[cdW])-tmp$coef[cdW]),na.rm=T)
    #
    enviro[j,paste('Sum_n_lamb_',ev,sep="")] = 0
    enviro[j,paste('Sum_n_lamb_betaDiff_w_',ev,sep="")] = 0
    for(k in 1:dim(df_esp)[1]){
      beta_true = c(df_esp$I.ve[k],df_esp$I.ve.2[k])
      spRefLin = spCoefs['ve']
      spRefQua = spCoefs['I(ve^2)']
      cdSp = regexpr(paste('taxa',k,':',sep=""),names(spCoefs))>0
      if(sum(cdSp)){
        spContrLin = spCoefs[cdSp & regexpr(':ve',names(spCoefs))>0]
        spContrQua = spCoefs[cdSp & regexpr('ve\\^2',names(spCoefs))>0]
        beta_fit = c(spRefLin+spContrLin,
                     spRefQua+spContrQua)
      }else{
        beta_fit = c(spRefLin,
                     spRefQua)
      }
      
      mu = - beta_true[1]/ (2*beta_true[2])
      sig = sqrt(-1/(2*beta_true[2]))
      Lambda_k_w = exp(-(w-mu)^2/(2*sig^2))
      
      n_k = sum(OCC$taxa==k)
      logBetaDiff_w = sum(c(w,w^2)*(beta_true-beta_fit))
      enviro[j,paste('Sum_n_lamb_',ev,sep="")] = enviro[j,paste('Sum_n_lamb_',ev,sep="")] + Lambda_k_w * n_k 
      enviro[j,paste('Sum_n_lamb_betaDiff_w_',ev,sep="")] = enviro[j,paste('Sum_n_lamb_betaDiff_w_',ev,sep="")] + Lambda_k_w * n_k *logBetaDiff_w
    }
  }
  
  write.table(enviro,paste('sx_table_H',H,'_ev_',ev,'.csv',sep=""),sep=";",row.names=F,col.names=T)
  
  Err_logs = enviro[,paste('Sum_n_lamb_',ev,sep="")] * enviro[,paste('Int_s_logsDiff',ev,sep="")] /enviro[,paste('mu_x_',ev,sep="")] 
  
  Err_logsp = enviro[,paste('Int_s_',ev,sep="")] * enviro[,paste('Sum_n_lamb_betaDiff_w_',ev,sep="")] /enviro[,paste('mu_x_',ev,sep="")] 
  
}


enviro[,paste('Int_s_',ev,sep="")] / enviro[,paste('mu_x_',ev,sep="")]

xxx=seq(0,3000,length.out=50)
for(k in 1:10){
  beta_true = c(df_esp$I.ve[k],df_esp$I.ve.2[k])
  mu = - beta_true[1]/ (2*beta_true[2])
  sig = sqrt(-1/(2*beta_true[2]))
  print(paste('mu =', mu))
  print(paste('sig = ',sig))
  print(max(exp(-(xxx-mu)^2/(2*sig^2))))
  
  print(max(exp(sum(c(xxx,xxx^2)*beta_true))))
}


plot(xxx,exp(-(xxx-mu)^2/(2*sig^2)),type="l")


plot(xxx,-(xxx-mu)^2/(2*sig^2),type="l")

y = sapply(xxx,function(w){sum(c(w,w^2)*beta_true)})
plot(xxx,y)
