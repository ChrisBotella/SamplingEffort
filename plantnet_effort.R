

library(glmnet)
library(raster)
library(rgdal)
library(rgeos)

user = "Christophe"

RepoDir = paste('C:/Users/',user,"/pCloud local/0_These/Github/SamplingEffort/",sep="")
saveDir = paste('C:/Users/',user,"/pCloud local/0_These/data/LOF data/19_03_20 for article/",sep="")
dataDir = paste('C:/Users/',user,"/pCloud local/0_These/data/",sep="")

setwd(RepoDir)
source('_functions.R')

originalVes = c('etp','chbio_12','chbio_1','chbio_5','alti','slope','awc_top','bs_top','clc')
predictor_formula = " ~ etp + I(etp^2) + I(chbio_12-etp) + I((chbio_12-etp)^2) + chbio_1 + I(chbio_1^2) + chbio_5 + I(chbio_5^2) + alti + I(alti^2) + slope + I(slope^2) + awc_top + I(awc_top^2) + bs_top + I(bs_top^2) + spht + slope:I(chbio_12-etp)"

nSpecies = 170
scoreThresh = .85

######
# Extract species occurrences
# from Pl@ntNet dataset
# Requires C.Botella environmental database
# Download zip : 
######

setwd('P:/Partage GLC19/')
OCC= read.csv("PL_complete.csv",sep=";",header=T)
OCC = OCC[OCC$FirstResPLv2Score>scoreThresh,,drop=F]

# remove species with less than MinOcc occurrences
length(unique(OCC$glc19SpId))
counts = table(OCC$glc19SpId)
counts = counts[order(counts,decreasing=T)]
spToKeep = names(counts)[1:nSpecies]
occ = OCC[OCC$glc19SpId%in% spToKeep,]

# Get environmental variables 
occ = get_variables(originalVes,occ,dpath = dataDir,fix=F)
occ$spht = factor(get_spht(occ$clc))
# remove NAs
tokeep = complete.cases(occ[,colnames(occ)%in%originalVes])
occ = occ[tokeep,,drop=F]

print(paste('Number of occurrences :',dim(occ)[1]))

######
# Create LOF grid
# remove empty or scarce squares
######

squareSize = 4000
Min= 5

## Make raster of squares including all the metropolitan French territory  
# get france polygon 
setwd(RepoDir)
france = getData('GADM',country="FRA",level="0")
# Project to lambert 93 
frProj = spTransform(france,CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

plot(frProj)

frBuff= gBuffer(frProj,width=squareSize,joinStyle="ROUND") # create France polygon with 4km buffer area around borders
ext = extent(frBuff)
r = raster(ext= ext,resolution=c(squareSize,squareSize),crs = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
#r[] = runif(ncell(r),0,1)
r_frBuf = rasterize(frBuff, r)  
# Values of the raster = Indices of squares  
r_frBuf[!is.na(r_frBuf[])] = 1:(sum(!is.na(getValues(r_frBuf))))
plot(r_frBuf)

# Number of cells
sum(!is.na(getValues(r_frBuf)))

## Remove squares with less than Min=5 occurrences 
# Count occurrences per cell
pts = SpatialPoints(occ[,c('Longitude','Latitude')],proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
occProj = spTransform(pts,CRSobj = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

occ = cbind(occ,as.data.frame(occProj@coords))
colnames(occ)[(dim(occ)[2]-1):dim(occ)[2]]= c('x_lamb93','y_lamb93')

occ$squareId = extract( r_frBuf, occ[,c('x_lamb93','y_lamb93')]  )

tab = table(occ$squareId)
tab = tab[order(tab,decreasing = T)] 
barplot(tab)
sum(tab>Min)

indicesToKeep = as.numeric(names(tab)[tab>Min])

r_lofCells = r_frBuf
r_lofCells[!r_lofCells[]%in%indicesToKeep] = NA
plot(r_lofCells)

occTot = occ
occ = occTot[occTot$squareId%in%indicesToKeep,]
sum(occ$squareId%in%indicesToKeep) # Number of occurrences for fitting LOF

cellsIds = unique(occ$squareId)
length(cellsIds)

setwd(saveDir)
saveRDS(occ,paste('occ_lof',nSpecies,'_score',scoreThresh,sep=""))

######
# Prepare background points LOF 
######

# We uniformly draw radom points inside the square extent of the raster and
# (i) MAXENT: Keep those falling inside the france raster (4 km buffer) for MAXENT, until at least 3 per cell
# (i) LOF: Only keep those falling inside LOF squares (until at least 5 per square)

nTmp = 30000
ext = extent(r_frBuf)
LOF_background = data.frame(x_lamb93 = NA , y_lamb93 = NA , q_hash = NA)
LOF_background = LOF_background[-1,,drop=F]
minNPerSquareLOF = 0
while(minNPerSquareLOF<4){
  tmp = data.frame(x_lamb93 = runif(nTmp,ext[1],ext[2]) , y_lamb93 = runif(nTmp,ext[3],ext[4]) , q_hash = NA)
  # Filter for LOF
  tmp$q_hash = extract(r_lofCells,tmp[,1:2])
  cdLOF = !is.na(tmp$q_hash)
  LOF_background = rbind(LOF_background,tmp[cdLOF,,drop=F])
  
  emptyCells = setdiff(cellsIds,unique(LOF_background$q_hash))
  if(length(emptyCells)>0){
    minNPerSquareLOF = 0
  }else{minNPerSquareLOF = min(table(LOF_background$q_hash))}
    
  flush.console()
  cat('     \r  Mini. n° points/LOF square:',minNPerSquareLOF,', tot. n° LOF pts:',dim(LOF_background)[1],'        \r')
}

## Get environmental variables
# LOF
pts = SpatialPoints(LOF_background[,1:2],proj4string = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") )
pts = spTransform(pts, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
LOF_background = cbind(LOF_background,data.frame(pts@coords))
colnames(LOF_background)[4:5] = c('Longitude','Latitude')
LOF_background = get_variables(originalVes,LOF_background,dpath = dataDir,fix=F)
LOF_background$spht = factor(get_spht(LOF_background$clc))
nona = complete.cases(LOF_background)
LOF_background = LOF_background[nona,]


#####
# Optimal drawing of background points
##### 

# We use an algorithm to draw background points so that
# (i) the whole number of points is enought for their 90% quantile enveloppe
# (in the features space) match the one of the whole geographic domain
# (ii) the concentration of points per unit of volume in the features
# space is the same in all cells

# Parameters
n = 20
delta = 2
eigRatio = 0.95
Alpha = 0.92
bootstrapSamples = 5
envelopeVolumeError = 0.05
beginWith = 20000
step = 10000

### Algorithm

## 1) Draw a large set of points representing the whole geographic domain 
# => n points per cell
resX= (extent(r_lofCells)[2]-extent(r_lofCells)[1])/ncol(r_lofCells)
resY = (extent(r_lofCells)[4]-extent(r_lofCells)[3])/nrow(r_lofCells)
cellsCoo = coordinates(r_lofCells)
cellsCoo = cellsCoo[!is.na(getValues(r_lofCells)),]
q_hashes = expand.grid(1:n,cellsIds)[,2]
x=rep(-9999,length(cellsIds)*n)
y=rep(-9999,length(cellsIds)*n)
toReDraw = rep(T,length(cellsIds)*n)
Df = data.frame(q_hash = q_hashes,x_lamb93=NA,y_lamb93=NA,Longitude=NA,Latitude=NA)
for(i in 1:length(originalVes)){
  eval(parse(text=paste('Df$',originalVes[i],'=NA',sep="")))
}

while(sum(toReDraw)>0){
  qLack = unique(q_hashes[toReDraw])
  print(paste('Missing ',sum(toReDraw),' points over into',length(qLack),' cells'))
  
  for(i in 1:length(qLack)){
    toFill = q_hashes==qLack[i] & toReDraw
    nLack = sum(x[toFill]==-9999)
    
    xmin = cellsCoo[i,1]-resX/2+delta
    xmax = cellsCoo[i,1]+resX/2-delta
    ymin = cellsCoo[i,2]-resY/2+delta
    ymax = cellsCoo[i,2]+resY/2-delta
    
    x[toFill] = runif(nLack,xmin,xmax)
    y[toFill] = runif(nLack,ymin,ymax)
    
    if(i/1000==round(i/1000)){
      flush.console()
      cat('    \r     Process...',100*i/length(qLack),'%       \r    ')
    }
  }
  Df[toReDraw,c('x_lamb93')] = x[toReDraw] 
  Df[toReDraw,c('y_lamb93')] = y[toReDraw] 
  
  
  pts = SpatialPoints(Df[toReDraw,c('x_lamb93','y_lamb93')],proj4string = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "))
  pts = spTransform(pts,CRSobj = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs '))
  Df[toReDraw,c('Longitude')] = pts@coords[,1]
  Df[toReDraw,c('Latitude')]=pts@coords[,2]
  
  
  Df[toReDraw,originalVes] = get_variables(originalVes,Df[toReDraw,c('Longitude','Latitude')],dpath = dataDir,fix=F)[,originalVes]
  # Points that having NA environmental variables must be redrawn
  # until we have all required points
  toReDraw = !complete.cases(Df)
  if(sum(toReDraw)>0){
    x[toReDraw]= -9999
    y[toReDraw]= -9999
  }
  
}
Df$spht = get_spht(Df$clc)

## 2) Approximation of the environmental volume of the geographic space 
# Get points coordinates over the PCA axis
# the first axis such that the sum of their 
# Eigen values makes a proportion of "eigRatio" of the total sum
Des = model.matrix(as.formula(predictor_formula),Df)
Des = Des[,2:dim(Des)[2]]
Des = scale(Des)
SVD = svd(Des)
d = 0
sumEig = 0
while(sumEig/sum(SVD$d)<eigRatio){
  d = d+1
  sumEig = sumEig + SVD$d[d]
}
print(d)
toAdd = SVD$u[,1:d]
colnames(toAdd)=paste('axe_',1:d,sep="")
toAdd=data.frame(toAdd)
hist(toAdd$axe_1,breaks="fd") 
# The distribution of our points over the SVD axis 
# is not that gaussian... But still...
# We approximate the total hyper-volume of the Alpha-confidence interval
# envelope of the points distribution in the environmental space  
# with an ellipsoïde
axesLength = NULL
for(di in 1:d){
  quantiles = SVD$d[di]* quantile(toAdd[,di],c(0.5-Alpha/2,.5+Alpha/2))
  axesLength = c(axesLength,quantiles[2]-quantiles[1])
}
# We now approximate the Alpha-enveloppe volume
# of the whole domain with an ellipsoïd volume
Gamma = 1/ (exp(1)*dgamma(1,shape=(d/2)+1,scale=1) ) 
Vol = prod(axesLength)* pi^(d/2) / Gamma 

#3) Compute the volume of every cell
ptPerCell = data.frame(q_hash =cellsIds, volume=rep(NA,length(cellsIds)),nPt = NA)
for(i in 1:length(cellsIds)){
  cellPts = toAdd[Df$q_hash==cellsIds[i],]
  axesLength = NULL
  for(di in 1:d){
    quantiles = SVD$d[di]* quantile(cellPts[,di],c(0.5-Alpha/2,.5+Alpha/2))
    axesLength = c(axesLength,quantiles[2]-quantiles[1])
  }
  # We now approximate the Alpha-enveloppe volume
  # of the whole domain with an ellipsoïd volume
  ptPerCell$volume[i] =  prod(axesLength)* pi^(d/2) / Gamma
  if(i/1000==round(i/1000)){
    flush.console()
    cat('    \r   Process...',100*i/length(cellsIds),'%       \r     ')
  }
}
totVolumes = sum(ptPerCell$volume)
ptPerCell$share = ptPerCell$volume/totVolumes


if(F){
  ## Facultative
  # 4) Chose total number of background 
  # so that their Alpha-enveloppe has approximately the same volume
  # as the whole domain Alpha-enveloppe
  # with tolerance in ratio : envelopeVolumeError  
  v = 0
  N = 0
  brick = data.frame(N=NA,vMean = NA,vSd = NA)
  brick = brick[-1,]
  while(((v-2*sdv)/Vol)<(1-envelopeVolumeError)){
    if(N==0){N=beginWith}else{N=N+step}
    print(N)
    
    # Compute number of points to draw per cell
    for(i in 1:length(cellsIds)){
      expected = N * ptPerCell$share[i]
      ptPerCell$nPt[i] = 1 + floor(expected-1) + rbinom(1,1,expected-floor(expected))
    }
    
    vs = NULL
    for(s in 1:bootstrapSamples){
      print(paste('sample',s))
      ids = NULL
      # Draw points in each cells
      for(i in 1:length(cellsIds)){
        ids = c(ids,sample( which(q_hashes==cellsIds[i]) , min(ptPerCell$nPt[i],n) ) )
        if(i/1000==round(i/1000)){
          flush.console()
          cat('    \r      Process...',100*i/length(cellsIds),'%       \r ')
        }
      }
      
      tmp = toAdd[ids,]
      tmpAxesLen = NULL
      for(di in 1:d){
        quantiles = SVD$d[di]* quantile(tmp[,di],c(0.5-Alpha/2,.5+Alpha/2))
        tmpAxesLen = c(tmpAxesLen,quantiles[2]-quantiles[1])
      }
      # Add total volume
      vs = c(vs,  prod(tmpAxesLen) * pi^(d/2) / Gamma )
    }
    v = mean(vs)
    sdv = sd(vs)
    brick = rbind(brick,data.frame(N=N,vMean=v,vSd=sdv))
  }
  brick$type = "mean"
  topl = rbind(brick,data.frame(N=brick$N,vMean=brick$vMean-2*brick$vSd,vSd=brick$vSd,type="Mean-2*sd"))
  topl = rbind(topl,data.frame(N=brick$N,vMean=brick$vMean+2*brick$vSd,vSd=brick$vSd,type="Mean+2*sd"))
  topl = rbind(topl,data.frame(N=brick$N,vMean=Vol,vSd=brick$vSd,type="whole domain volume"))
  ggplot(topl,aes(x=N,y=vMean,colour=type))+geom_line()+scale_y_continuous(limits=c(0,max(topl$vMean)))+theme_bw()+ylab('Approximate Alpha-enveloppe volume')+xlab('Number of points drawn')
}

# CHOSE N
N = 40000

# We draw in each cell a number of points equal to 
# n_q = 1 + floor(expected-1) + Bernoulli(expected-floor(expected))
# if expected>1 => E[n_q] = expected  , and 1 otherwise 
# Total number of background points is thus superior in expectation
# to N, but this insures there is at least one point per cell 
for(i in 1:length(cellsIds)){
  expected = N * ptPerCell$share[i]
  ptPerCell$nPt[i] = min( 1 + max(0,floor(expected-1)) + rbinom(1,1,expected-floor(expected)) , n)
}
# If it remains points to be distributed
Reminder = N - sum(ptPerCell$nPt)
ptPerCell= ptPerCell[order(ptPerCell$nPt,decreasing = T),]
for(i in which(ptPerCell$nPt<20)){
  toPut = 20-ptPerCell$nPt[i]
  if(Reminder>=toPut){
    Reminder = Reminder - toPut
    ptPerCell$nPt[i] = ptPerCell$nPt[i] + toPut
  }else{
    ptPerCell$nPt[i] = ptPerCell$nPt[i] + Reminder
    Reminder = 0
  }
}


# Draw points from Df
pool = NULL
for(i in 1:length(cellsIds)){
  pool = c(pool, sample( which(q_hashes==cellsIds[i]) , min(ptPerCell$nPt[i],n ) ))
  if(i/1000==round(i/1000)){
    flush.console()
    cat('    \r      Process...',100*i/length(cellsIds),'%       \r ')
  }
}
toKeep = (1:dim(Df)[1])%in%pool

BG = Df[toKeep,,drop=F]
# Save 
setwd(saveDir)
saveRDS(BG,paste('BackgroundPts_lof',nSpecies,'_score',scoreThresh,sep=""))

#####
# Fit LOF-GLMNET
#####

library(glmnet)
library(raster)
library(rgdal)
library(rgeos)
RepoDir = paste('C:/Users/',user,"/pCloud local/0_These/Github/SamplingEffort/",sep="")
saveDir = paste('C:/Users/',user,"/pCloud local/0_These/data/LOF data/19_03_20 for article/",sep="")

setwd(RepoDir)
source('_functions.R')

originalVes = c('etp','chbio_12','chbio_1','chbio_5','alti','slope','awc_top','bs_top','clc')
predictor_formula = " ~ etp + I(etp^2) + I(chbio_12-etp) + I((chbio_12-etp)^2) + chbio_1 + I(chbio_1^2) + chbio_5 + I(chbio_5^2) + alti + I(alti^2) + slope + I(slope^2) + awc_top + I(awc_top^2) + bs_top + I(bs_top^2) + spht + slope:I(chbio_12-etp)"

setwd(saveDir)
BG = readRDS('BackgroundPts_lof170_score0.85')
occ = readRDS('occBackgroundPts_lof170_score0.85')

colnames(occ)[colnames(occ)=="squareId"]="q_hash"
colnames(occ)[colnames(occ)=="glc19SpId"]="taxa"
occ = occ[,c('q_hash',originalVes,'spht','taxa')]
BG = BG[,c('q_hash',originalVes,'spht')]

Data = bind.background.pts(occ,BG,Bias_killer=300,equalSizeCells=T)

lof.mod = lof_glmnet(Data,Data$pseudo_y,sub('~(.*)','\\1' , predictor_formula),weights=Data$w,lambdaMinRatio=1e-8,nlambda=200)

setwd(saveDir)
saveRDS(lof.mod,'model_lof170_score.85')

