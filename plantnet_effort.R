library(glmnet)
library(raster)
library(rgdal)
library(rgeos)
library(data.table)
library(plyr)
library(ggplot2)

dir = "C:/Users/Christophe/Downloads/SamplingEffort-master/"
RepoDir=dir
dataDir =dir
saveDir =dir

#### Parameters
# Model formula 
originalVes = c('etp','chbio_12','chbio_1','chbio_5','alti','slope','awc_top','bs_top','clc')
predictor_formula = " ~ etp + I(etp^2) + I(chbio_12-etp) + I((chbio_12-etp)^2) + chbio_1 + I(chbio_1^2) + chbio_5 + I(chbio_5^2) + alti + I(alti^2) + slope + I(slope^2) + awc_top + I(awc_top^2) + bs_top + I(bs_top^2) + spht + slope:I(chbio_12-etp)"

# species occurrences selection parameters
nSpecies = 300
scoreThresh = .85
# Sampling cells grid parameters
squareSize = 4000
Min= 5
# Background points parameters
n = 6
delta = 2

# Chosen species for plot of the intensity over France 
SpeciesForIntensityPlot= "Phytolacca americana L."

setwd(dir)
source('_functions.R')


######
# Extract species occurrences
# from Pl@ntNet dataset
# Requires C.Botella environmental database
# Download zip : 
######

setwd(dir)

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

## Make raster of squares including all the metropolitan French territory  
# get france polygon 
setwd(RepoDir)
#france = getData('GADM',country="FRA",level="0")
france = readRDS('gadm36_FRA_0_sp.rds')

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

# Number of cells
print('Number of cells:',sum(!is.na(getValues(r_frBuf))))

## Remove squares with less than Min=5 occurrences 
# Count occurrences per cell
pts = SpatialPoints(occ[,c('Longitude','Latitude')],proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
occProj = spTransform(pts,CRSobj = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

occ = cbind(occ,as.data.frame(occProj@coords))
colnames(occ)[(dim(occ)[2]-1):dim(occ)[2]]= c('x_lamb93','y_lamb93')

occ$q_hash = extract( r_frBuf, occ[,c('x_lamb93','y_lamb93')]  )

tab = table(occ$q_hash)
tab = tab[order(tab,decreasing = T)] 

indicesToKeep = as.numeric(names(tab)[tab>Min])

r_lofCells = r_frBuf
r_lofCells[!r_lofCells[]%in%indicesToKeep] = NA
#plot(r_lofCells)


occTot = occ
occ = occTot[occTot$q_hash%in%indicesToKeep,]
#sum(occ$q_hash%in%indicesToKeep) # Number of occurrences for fitting LOF

cellsIds = unique(occ$q_hash)
print(paste('Number of sampling cells retained:',length(cellsIds)))

setwd(saveDir)
saveRDS(r_frBuf,'r_frBuf')
r_lofCellsName = paste('r_lofCells',nSpecies,'_score',scoreThresh,'_size',squareSize,sep="")
saveRDS(r_lofCells,r_lofCellsName)
occName = paste('occ_lof',nSpecies,'_score',scoreThresh,'_size',squareSize,sep="")
saveRDS(occ,occName)

######
# Draw n background point per cell 
######
## Draw n points per cell
Df = draw.pts.in.cells(n,r_lofCells,originalVes,delta,dataDir,miniPerCell=3)
cd = complete.cases(Df)
Df = Df[cd,]

# Save 
setwd(saveDir)
BGptsName = paste('BackgroundPts_nPerCell_lof',nSpecies,'_score',scoreThresh,'_size',squareSize,sep="")
saveRDS(Df,BGptsName)


#####
# Fit and save model
#####
rm(Df)
gc(reset=T)

setwd(saveDir)
BG = readRDS(BGptsName)
occ = readRDS(occName)

if("glc19SpId"%in%colnames(occ)){
  colnames(occ)[colnames(occ)=="glc19SpId"]="taxa"
}
occ = occ[,c('q_hash',originalVes,'spht','taxa')]
BG = BG[,c('q_hash',originalVes,'spht')]

Data = bind.background.pts(occ,BG,Bias_killer=300)

lof.mod = lof_glmnet(Data,Data$pseudo_y,sub('~(.*)','\\1' , predictor_formula),weights=Data$w,lambdaMinRatio=1e-12,nlambda=300)

setwd(saveDir)
modelName = paste('model_lof',nSpecies,'_score',scoreThresh,'_squareSize',squareSize,'_BG',n,'ptsPerCell')
saveRDS(lof.mod,modelName)

rm(Data)
gc(reset=T)


######
# Load fitted model and data 
# TO RUN before ALL plots
######


setwd(saveDir)
df = readRDS(paste(saveDir,"df_departements_FR",sep=""))
occ = readRDS(occName)
BG = readRDS(BGptsName)
lof.mod = readRDS(modelName)
r_lofCells = readRDS(r_lofCellsName)
r_frBuf = readRDS('r_frBuf')

if("squareId"%in%colnames(occ)){
  colnames(occ)[which(colnames(occ)=="squareId")] = "q_hash"
}

tabS = table(occ$glc19SpId)
tabS = tabS[order(tabS,decreasing = T)]
tabS = data.frame(glc19SpId=names(tabS),nOcc=as.numeric(tabS))

occo = occ[,c('glc19SpId','scName')]
occo$scName = as.character(occo$scName)
occo$glc19SpId = as.character(occo$glc19SpId)
occo = unique(occo)
tabS = merge(tabS,occo,all.x=T)
tabS = tabS[order(tabS$nOcc,decreasing = T),]

SpTable = tabS[,c('glc19SpId','nOcc','scName')]
setwd(saveDir)
write.table(SpTable,paste("speciesTable",nSpecies,".csv",sep=""),sep=";",row.names=F,col.names=T)

# older model version
nl = dim(lof.mod$beta)[2]
coefficients = get.treat.coef(lof.mod,nl)$coefficients
qCoefs = coefficients[regexpr("q_hash",names(coefficients))>0]
cells = data.frame(q_hash=paste('q_hash',unique(occ$q_hash),sep=""))
coefs = data.frame(coef= qCoefs, q_hash = names(qCoefs))
cells = merge(cells,coefs,by="q_hash",all.x=T)
cells$coef[is.na(cells$coef)] = 0
hist(log10(exp(cells$coef)),breaks='fd')
coos = rasterToPoints(r_lofCells)
colnames(coos)[1:3]=c('x_l93','y_l93','q_hash')
coos=as.data.frame(coos)
coos$q_hash = paste('q_hash',coos$q_hash,sep="")
cells = merge(cells,coos,by="q_hash",all.x=T)


######
# Plot effort estimate on a Map
######

#### Plot effort sur departements
tmp = log10(exp(cells$coef))
classes = cut(tmp,c(min(tmp)-.1,quantile(tmp,c(.1,.4,.6,.9)),max(tmp)+.1),dig.lab=2)

labels = data.frame(x=c(500000,650000,470000,400000,600000),
                    y=c(6300000,6250000,6500000,6340000,6650000),
                    label=c("32","11","16","40","36"))


#sum(is.na(classes))
#unique(classes)
#unique(tmp[which(is.na(classes))])
colo = colorRampPalette(c("darkorchid4","goldenrod"))(length(unique(classes)))
lev= levels(classes)
labels = paste(lev,paste('/ quantile',c('0 to .1','.1 to .4','.4 to .6','.6 to .9','.9 to 1')))

p = ggplot()+geom_polygon(data=df,aes(x=lon,y=lat,group=id),fill = NA,color = "grey20",size=.35)+theme_bw()
p = p + geom_tile(data=cells,aes(x=x_l93,y=y_l93,fill=classes),alpha=.95)+scale_fill_manual(values=colo,name='log10(samp. eff.)',labels=labels)+xlab('Longitude L93')+ylab('Latitude L93')
p = p + theme(legend.title = element_text(size=18),
              legend.text = element_text(size=15),
              axis.title.x = element_text(size=18),
              axis.title.y = element_text(size=18),
              axis.text.x = element_text(size=13),
              axis.text.y = element_text(size=13))
setwd(saveDir)
png(paste('effortMap_',modelName,'.png'),width=1000,height=750)
print(p)
dev.off()

#####
# Plot species intensity over France
#####

resCell = 4000
nInsideCell= 2

# Discrete intensity map
### Make prediction grid on all France
resInsideCell = resCell/ nInsideCell
sequence = seq( -resInsideCell*(nInsideCell-1)/2 ,
                resInsideCell*(nInsideCell-1)/2 ,
                resInsideCell)
AllCells = rasterToPoints(r_frBuf)
colnames(AllCells) = c('x_l93','y_l93','cellId')
AllCells = data.frame(AllCells)
grid = expand.grid(cellId=unique(AllCells$cellId),xdif=sequence,ydiff=sequence)
grid = merge(AllCells,grid,by='cellId')
grid$x_l93= grid$x_l93 + grid$xdif
grid$y_l93 = grid$y_l93 + grid$ydiff
grid = grid[,c('cellId','x_l93','y_l93')]
pts = SpatialPoints(grid[,c('x_l93','y_l93')],proj4string = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
pts = spTransform(pts,CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
pts = pts@coords
colnames(pts)=c('Longitude','Latitude')
grid = cbind(grid,pts)

Sp = as.character(tabS$glc19SpId[tabS$scName==SpeciesForIntensityPlot])

glc19SpIdVeco = as.character(SpTable$glc19SpId)
originalVes = c('etp','chbio_12','chbio_1','chbio_5','alti','slope','awc_top','bs_top','clc',"spht")
predMatrix = predict.lof.spLogRelativeIntensity(grid,
                                                coefficients,
                                                Intensityformula = as.formula(" ~ etp + I(etp^2) + I(chbio_12-etp) + I((chbio_12-etp)^2) + chbio_1 + I(chbio_1^2) + chbio_5 + I(chbio_5^2) + alti + I(alti^2) + slope + I(slope^2) + awc_top + I(awc_top^2) + bs_top + I(bs_top^2) + spht + slope:I(chbio_12-etp)"),
                                                originalVes = originalVes, 
                                                taxasToPredict = Sp,
                                                glc19SpIdVec= glc19SpIdVeco,
                                                extractorScript = paste(dir,'_functions.R',sep=""),
                                                dataDir=dir)
pred = as.vector(predMatrix)
tmp = log10(exp(pred))
grido = grid[!is.na(tmp),]
tmp = tmp[!is.na(tmp)]

# PLOT
df = readRDS(paste(saveDir,"df_departements_FR",sep=""))
classes = cut(tmp,c(min(tmp)-.1,quantile(tmp,c(.1,.4,.6,.9),na.rm=T),max(tmp)+.1),dig.lab=2,na.rm=T)
colo = colorRampPalette(c("darkorchid4","goldenrod"))(length(unique(classes)))
lev= levels(classes)
#colVals = sapply(classes,function(cl) colo[lev==cl])
dptsNum = data.frame(x=c(500000,650000,470000,400000,600000),
                    y=c(6300000,6230000,6520000,6340000,6635000),
                    label=c("32","11","16","40","36"))
labels = paste(lev,paste('/ quantile',c('0 to .1','.1 to .4','.4 to .6','.6 to .9','.9 to 1')))
p = ggplot()+geom_polygon(data=df,aes(x=lon,y=lat,group=id),fill = NA,color = "grey20",size=.35)
p = p + geom_tile(data=grido,aes(x=x_l93,y=y_l93,fill=classes),size=1,alpha=.8)+scale_fill_manual(values=colo,name="log10 of species intensity",labels=labels)+theme_bw()
p = p + geom_text(data=dptsNum,aes(x=x,y=y,label=label),size=6)
p = p + xlab("Longitude in Lambert 93") + ylab("Latitude in Lambert 93")
p = p + theme(legend.title = element_text(size=18),
              legend.text = element_text(size=15),
              axis.title.x = element_text(size=18),
              axis.title.y = element_text(size=18),
              axis.text.x = element_text(size=13),
              axis.text.y = element_text(size=13))
setwd(saveDir)
png(paste("SpIntens_",tabS$scName[tabS$glc19SpId==Sp],'_',modelName,".png",sep=""),width=1200,height=750)
print(p)
dev.off()

#####
# Occurrences map
#####

occi = occ[occ$glc19SpId==Sp,]
p=ggplot()+geom_polygon(data=df,aes(x=lon,y=lat,group=id),fill = NA,color = "grey20",size=.35)
p = p + geom_point(data=occi,aes(x=occi$x_lamb93,y=occi$y_lamb93),alpha=.35,pch=18,cex=3,color="blue") + theme_bw()
p = p + xlab("Longitude in Lambert 93") + ylab("Latitude in Lambert 93")
p = p + theme(legend.title = element_text(size=18),
              legend.text = element_text(size=15),
              axis.title.x = element_text(size=18),
              axis.title.y = element_text(size=18),
              axis.text.x = element_text(size=13),
              axis.text.y = element_text(size=13))
setwd(saveDir)
png(paste("SpOcc_",tabS$scName[tabS$glc19SpId==Sp],'_',modelName,".png",sep=""),width=1000,height=750)
print(p)
dev.off()

