# Created on 20/03/2019
# Simulation experiment of article 
# "A novel method for estimating relative sampling effort across space from large amounts of geolocated species occurrences coming from citizen-science data"
library(data.table)
library(raster)
library(ggplot2)
library(glmnet)

dir = "C:/Users/Christophe/pCloud local/0_These/Github/SamplingEffort/"
setwd(dir)
source('_functions.R')

SaveModels = T
expeName = "test"
runName = "_0"
SaveDir = dir

setwd(SaveDir)


#####
# Define simulation domain & var. enviro.  & observation
##### 

resolution = 0.005
deltaBorder = 1e-7 # Points are strictly generated inside the domain limits
# by avoiding machine numerical approximations

xmin = 0;xmax=10;ymin=0;ymax=1
# Create points grid over spatial domain
vals = expand.grid(Longitude=seq(xmin+(resolution/2)+deltaBorder,xmax-(resolution/2)-deltaBorder,resolution),
                   Latitude=seq(ymin+(resolution/2)+deltaBorder,ymax-(resolution/2)-deltaBorder,resolution),
                   esp=NA)
# linear environmental gradient, increase from west to east
vals$axe3 = vals$Longitude - 5

# 1 on left, 0 on right
vals$mid = mid(vals$axe3)

cd = vals$Latitude==vals$Latitude[1]

# SigmoFast
vals$SigmoFast = SigmoFast(vals$axe3)

vals$SigmoMedium = SigmoMedium(vals$axe3)
# scalpHill
vals$scalpHill = scalpHill(vals$axe3)
# valley
vals$valley = valley(vals$axe3)
# cutHurts
vals$cutHurts = cutHurts(vals$axe3)
# cutNice 
vals$cutNice = cutNice(vals$axe3)

# Background points
N_0 = 50000# Number of background points

# Draw uniform background points on spatial domain
dom = data.frame(Longitude=runif(N_0,xmin+deltaBorder,xmax-deltaBorder),
                 Latitude=runif(N_0,ymin+deltaBorder,ymax-deltaBorder))
dom$axe3 = dom$Longitude - 5

#####
# Design of virtual species
#####

# Define groups of species
# Each group may contain several "TG species" and a main species (the last one)  
groups = list("Monosp"=1,
              "Bisp"=c(2,3),
              "Trisp"=c(2,3,1),
              "Hexasp"=c(2,3,4,5,6,1))

# optimums positions of species
mus = c(0,-2.5,2.5,-4.5,4.5,2)

# standard deviations
sigs = c(1.6,1.6,1.6,.8,.8,1)
# Compute loglinear reparametrization
df_esp = data.frame(sp_id=1:length(mus),
                    I.axe3 = 0,
                    I.axe3.2. = 0)

for(i in 1:dim(df_esp)[1]){
  Beta = nauralToLoglinearParam( mus[i] , sigs[i])
  df_esp$I.axe3[i] = Beta[2]
  df_esp$I.axe3.2.[i] = Beta[3]
}
setwd(SaveDir)
write.table(df_esp,"df_esp.csv",sep=";",row.names=F,col.names=T)

gr = "Bisp"
x = seq(-5,5,0.1)
yy = rep(0,length(x))
for(i in groups[gr][[1]]){
  tmp = gaussDens(x,mus[i],sigs[i])
  yy = yy + tmp/sum(tmp)
}
plot(x,yy,ylim=c(0,max(yy)),lwd=2,type="l")

#####
# Species abundance model formula
#####
f= "I( axe3 ) + I( axe3 ^2)"

#####
# 1) Simulation and Fit
#####

## Define the experimental planning into "Sched" data.frame

nrep = 2 # Number of repetitions per (sampling,species) model configurations
# True effort
trueEffort = as.character(
  expand.grid(1:nrep,
              c('mid','cutNice','cutHurts'))[,2])  

# Pool of virtual species to simulate
spGroup = c( rep('Bisp', length(trueEffort)) )  
# Class of estimators for sampling effort
effortClass = c( rep('histo_glmnet',length(trueEffort))) 
# mesh 
mesh = c(rep(10,length(trueEffort)))
# nb of points to draw 
nnn = c(rep(400000,length(trueEffort)))

Sched = data.frame( id=1:length(spGroup),
                    spGroup=spGroup,
                    effortClass=effortClass,
                    trueEffort = trueEffort,
                    mesh=mesh,
                    n = nnn )
setwd(SaveDir)
write.table(Sched,paste('Schedule_',expeName,runName,".csv",sep="") , sep=";",row.names=F,col.names=T)

# variables matrix
XMat = matrix(NA,dim(vals)[1],2)
XMat[!is.na(vals$axe3),] = model.matrix( as.formula(paste('~',f)) , vals )[,2:3]

for(j in Sched$id){
  n_i = Sched$n[j] 
  gr = as.character(Sched$spGroup[j]) # Name of species group 
  pool = groups[gr][[1]]  # Vector of species identifiers
  espece = pool[length(pool)] # Focal species
  type = Sched$trueEffort[j] # Name of the sampling effort type
  # Must be present in columns of var:vals
  
  NI = round( n_i / length(pool) ) # Number of points to generate per species
  # Simulate species occurrences points
  for(i in 1:length(pool)){
    print(paste('esp',pool[i]))
    # Load species abundance fonction parameters
    par = as.numeric(df_esp[df_esp$sp_id== pool[i],c('I.axe3','I.axe3.2.')])
    # Vector of values of species abundance function over grid points 
    vals$esp = exp( XMat %*% par )
    # Vector of values of observed points intensity over grid points
    vals$int = vals[,as.character(type)] * vals$esp 
    # Continuously draw points from intensity values (vals$int) gridded (Longitude/Latitude)
    df = DrawPts( vals , NI ) 
    
    df$taxa=as.character(pool[i])
    if(i==1){DF = df}else{DF=rbind(DF,df)}
  }
  DF$axe3 = DF$Longitude - 5

  # msh : number of quadrats along Longitude axis
  msh = Sched$mesh[j]
  # q : size of quadrats for sampling effort model
  q = 10 / msh
  # Get the quadrat identifier of all points
  # It is then used as a factor in glmnet
  DF = get_q_hash_aniso(DF,q,1,0,0)
  dom = get_q_hash_aniso(dom,q,1,0,0)
  # Duplicate and attach background points (dom) to occurrences (DF)
  df = bind.background.pts(DF,bg=dom,Bias_killer=200)
  
  # Fit model with GLMNET
  deb=Sys.time()
  Momo = lof_glmnet(df,y = df$pseudo_y,occupationString = f,
                    weights = df$w,
                    lambdaMinRatio = 1e-10,
                    nlambda=200)
  print(Sys.time()-deb)

  if(SaveModels){
    setwd(SaveDir)
    saveRDS(Momo,paste('From_expe_',expeName,runName,'_model_',Sched$id[j],sep=""))
  }
}

#####
# 2) Functions summary statistics
#####

step2 = 0.01
evalPts = seq(-5+1e-6,5-1e-6,step2)

toread = paste("From_expe_",expeName,runName,sep="")
schedName = paste("Schedule_",expeName,runName,".csv",sep="") # CHANGE
setwd(SaveDir)
Sched = read.csv(schedName,sep=";",header=T)


res = data.frame(id = Sched$id , coefficients = NA)

for(id in Sched$id){
  model = readRDS(paste(toread,'_model_',id,sep=""))
  res$coefficients[res$id==id] = list(model$coefficients)
}


x = data.frame(axe3=evalPts,Longitude=evalPts+5,Latitude=0.5)
evalEstimates = data.frame(id=1,axe3=0,EstName="no",value=NA)
evalEstimates = evalEstimates[-1,]
for(id in Sched$id){
  msh = Sched$mesh[Sched$id==id]
  Coefs = res$coefficients[res$id==id][[1]]
  
  # Effort ---- Histo
  q = 10 / msh
  
  x_q = get_q_hash_aniso(x,q,10,0,0)
  x_q$coefName = paste('q_hash',x_q$q_hash,sep="")
  CoefsDf = data.frame(coefName=names(Coefs),coef=Coefs)
  x_q = merge(x_q,CoefsDf,by='coefName',all.x=T)
  if(sum(is.na(x_q$coef))>0){x_q$coef[is.na(x_q$coef)]=0}
  x_q$value =   exp(x_q$coef)
  
  evalEstimates = rbind(evalEstimates,
                        data.frame(id=id,axe3=x_q$axe3,EstName="effort",value=x_q$value) )
  # Species
  sps = groups[as.character(Sched$spGroup[Sched$id==id])][[1]]
  X = model.matrix(~ I(axe3) + I(axe3^2),x)[,2:3]
  Coefs_ = Coefs[c('I(axe3)','I(axe3^2)')]
  EvalSp = exp( X %*% matrix(Coefs_,2,1) )
  sp = min(sps)
  evalEstimates= rbind(evalEstimates,
                       data.frame(id=id,axe3=evalPts,EstName=paste('taxa',sp,sep=""),value=EvalSp))
  sps = sps[!sps==sp]
  if(length(sps)>0){
    for(sp in sps){
      Coefs__ = Coefs[paste('taxa',sp,':',c('I(axe3)','I(axe3^2)'),sep="")]
      EvalSp = exp( X %*% matrix(Coefs_+Coefs__,2,1) )
      evalEstimates= rbind(evalEstimates,
                           data.frame(id=id,axe3=evalPts,EstName=paste("taxa",sp,sep=""),value=EvalSp))
    }
  }
}
setwd(SaveDir)
write.table(evalEstimates,paste("evalEstimates_",runName,".csv",sep=""),sep=";",row.names=F,col.names=T)


TAB = merge(Sched[,c('id','spGroup','effortClass','trueEffort','mesh','n')],evalEstimates,by="id")


TABO = aggregate(list(value=evalEstimates$value),
                 by=list(spGroup=TAB$spGroup,effortClass=TAB$effortClass,
                         trueEffort=TAB$trueEffort,mesh=TAB$mesh,n=TAB$n,
                         axe3=evalEstimates$axe3,EstName=evalEstimates$EstName),
                 FUN=mean)
TABO$stat='Average Estimator'

sdTop= function(x){mean(x)+2*sd(x)}

TMP = aggregate(list(value=evalEstimates$value),
                by=list(spGroup=TAB$spGroup,effortClass=TAB$effortClass,
                        trueEffort=TAB$trueEffort,mesh=TAB$mesh,n=TAB$n,
                        axe3=evalEstimates$axe3,EstName=evalEstimates$EstName),
                FUN=sdTop )
TMP$stat = 'Avg. Est.+2*st. dev.'

TABO = rbind(TABO,TMP)

sdAbo=function(x){mean(x)-2*sd(x)}
TMP = aggregate(list(value=evalEstimates$value),
                by=list(spGroup=TAB$spGroup,effortClass=TAB$effortClass,
                        trueEffort=TAB$trueEffort,mesh=TAB$mesh,n=TAB$n,
                        axe3=evalEstimates$axe3,EstName=evalEstimates$EstName),
                FUN=sdAbo )
TMP$stat = 'Avg. Est.-2*st. dev.'

TABO = rbind(TABO,TMP)

setwd(SaveDir)
write.table(TABO,paste("CurveStats",runName,".csv",sep=""),sep=";",row.names=F,col.names=T)

##### 
# 3) PLOT all fitted estimates:
# Sampling effort (vs Longitude) 
# and species intensities (vs enviro. variable) 
######

setwd(SaveDir)
#TABO = read.csv("CurveStats.csv",sep=";",header=T)
TABO = read.csv(paste("CurveStats",runName,".csv",sep=""),sep=";",header=T)

TABO$Longitude = TABO$axe3+5

### SAVE PLOTS
plots = aggregate( list(na=TABO$axe3),by=list(effortClass=TABO$effortClass,
                                              trueEffort=TABO$trueEffort,
                                              spGroup=TABO$spGroup,
                                              mesh=TABO$mesh,
                                              n=TABO$n,
                                              EstName=TABO$EstName),FUN=mean)


for(i in 1:dim(plots)[1]){
  
  class = as.character(plots$effortClass[i])
  eff = as.character(plots$trueEffort[i])
  group = as.character(plots$spGroup[i])
  msh = plots$mesh[i]
  n = plots$n[i]
  est = as.character(plots$EstName[i])
  try(dev.off())
  save.summary.functionals.estimate.graphs.V2(runName,
                                              resTable=TABO,
                                              resTableResolution=step2,
                                              saveDirectory=SaveDir,
                                              class,
                                              eff,
                                              group,
                                              msh,
                                              est,
                                              n,
                                              grid = vals[vals$Latitude==vals$Latitude[1],])
  
}



#####
# Single model Plot: Estimated vs True effort 
#####

## Load fitted model and simulation model
setwd(SaveDir)
toread = paste("From_expe_",expeName,runName,sep="")
schedName = paste("Schedule_",expeName,runName,".csv",sep="") # CHANGE
Sched = read.csv(schedName,sep=";",header=T)
id= 1
Momo = readRDS(paste(toread,'_model_',id,sep=""))


stepo = 0.01
trueDomain = data.frame(Longitude=seq(1e-6,10-1e-6,stepo),Latitude=0.5)
q = 1
trueDomain = get_q_hash_aniso(trueDomain,step_x = q,step_y = 1,x_0 = 0,y_0=0 )
unique(trueDomain$q_hash)
trueDomain$coefName = paste('q_hash',trueDomain$q_hash,sep="")

Coefs = data.frame(coefName=names(Momo$coefficients),coef = Momo$coefficients)
trueDomain = merge(trueDomain,Coefs,by='coefName',all.x=T)
trueDomain$coef[is.na(trueDomain$coef)]= 0
trueDomain$val = exp(trueDomain$coef)
trueDomain$val = trueDomain$val/sum(trueDomain$val*stepo)

trueDomain$type="est"
trueDomain2 = data.frame(Longitude=trueDomain$Longitude,val=NA,type="true")
cd = Sched$id==id
if(Sched$trueEffort[cd]=="cutNice"){
  trueDomain2$val= cutNice(trueDomain$Longitude-5)
}else if(Sched$trueEffort[cd]=="mid"){
  trueDomain2$val= mid(trueDomain$Longitude-5)
}else if(Sched$trueEffort[cd]=="cutHurts"){
  trueDomain2$val= cutHurts(trueDomain$Longitude-5)
}else if(Sched$trueEffort[cd]=="Sigmo"){
  trueDomain2$val= SigmoFast(trueDomain$Longitude-5)
}
trueDomain2$val = trueDomain2$val/sum(trueDomain2$val * stepo)

trueDomainF = rbind(trueDomain[,c('Longitude','val','type')],trueDomain2)
trueDomainF$type = factor(trueDomainF$type,levels=c('true','est'))
library(ggplot2)
p = ggplot(trueDomainF,aes(x=Longitude,y=val))
p = p+geom_line(aes(colour=trueDomainF$type ,size=trueDomainF$type))
p = p+scale_y_continuous(limits=c(0,0.5))
p = p+scale_size_manual('curve',values=c(3,1))  
p = p+scale_colour_manual('curve',values=c('goldenrod','blue'))
p = p+theme_bw()
print(p)

