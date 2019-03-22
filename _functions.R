# created on 29/11/18
# Functions for LOF-GLMNET and diagnostic



nauralToLoglinearParam = function(mu,sig){
  Beta_linear = mu/sig^2
  Beta_quadratic = -1/(2*sig^2)
  Beta_0 = - (log(2*pi*sig^2)/2) - (mu^2 /(2*sig^2)) 
  return(c(Beta_0,Beta_linear,Beta_quadratic))
}

# return 1D gaussian density value
gaussDens = function(x,mu,sig){
  return( exp(-(x-mu)^2/(2*sig^2))/(sig*sqrt(2*pi)) )
}

# Trace les différentes courbes catégorisées par la colonne "stat" de resTable
# ainsi que la courbe de référence
save.summary.functionals.estimate.graphs.V2 = function(runName,resTable,
                                                    resTableResolution,
                                                    saveDirectory,
                                                    class,
                                                    eff,
                                                    group,
                                                    msh,
                                                    est,
                                                    n,
                                                    grid=vals[vals$Latitude==vals$Latitude[1],],
                                                    gridResolution=resolution,
                                                    dfEspeces=df_esp){
  Truth = "True"
  others=c('Avg. Est.+2*st. dev.','Avg. Est.-2*st. dev.','Average Estimator')
  
  
  if(sum(unique(resTable$stat)%in%others)==3){
    
    cd1 = resTable$spGroup==group & resTable$effortClass==class & resTable$trueEffort==eff & resTable$mesh==msh & resTable$EstName==est & resTable$n==n
     
    if(est=="effort"){tplot = resTable[cd1,c('Longitude','value','stat'),drop=F]
    }else{tplot = resTable[cd1,c('axe3','value','stat'),drop=F]}
    
    # We scale the density function 
    for(s in unique(tplot$stat)){tplot$value[tplot$stat==s] = tplot$value[tplot$stat==s] / sum(resTableResolution * tplot$value[tplot$stat==s])}
    
    if(est=="effort"){
      tplot = rbind(tplot, data.frame(Longitude=grid$Longitude,value=grid[,eff],stat=Truth))
    }else{
      sp_id = as.numeric(substr(est,5,nchar(est)))
      M = model.matrix( ~ I(axe3)+I(axe3^2) , grid)[,2:3]
      v = matrix(as.numeric(dfEspeces[dfEspeces$sp_id==sp_id,2:3]),2,1)
      sp_Intens = exp( M  %*% v )
      sp_Intens = sp_Intens/ sum(gridResolution* sp_Intens)
      tplot = rbind(tplot, data.frame(axe3=grid$axe3,value=sp_Intens,stat=Truth))
    }
    tplot$value[tplot$stat==Truth] = tplot$value[tplot$stat==Truth] / sum(tplot$value[tplot$stat==Truth]* gridResolution)
    
    # We re-factor the curves for plotting them in the desired superposition order
    tplot$stat = factor(tplot$stat,levels=c(Truth,others))
    
    if(est=="effort"){
      p = ggplot(tplot,aes(x=Longitude,y=value)) + geom_line(aes(colour=tplot$stat,size=tplot$stat))
      p = p + scale_color_manual('',values=c("goldenrod","blue",'blue',"brown"))
      p = p+scale_size_manual('',values = 2*(1+c(3,1,1,1)))
      p = p + theme_bw()+theme(legend.text = element_text(size=35),
                               legend.title = element_text(size=30),
                               panel.grid.major=element_line(color="grey",size=1),
                               axis.title.x=element_text(size=35),
                               axis.title.y=element_text(size=35),
                               plot.title=element_text(size=40),
                               axis.text.x = element_text(size=25),
                               axis.text.y = element_text(size=25),
                               panel.border = element_blank(),
                               axis.line.x = element_line(color="black", size = 3),
                               axis.line.y = element_line(color="black", size = 3))
      p = p + ylab('Density on [0,10]')
      p = p + xlab('Longitude')
      p = p + scale_y_continuous(limits=c(0,0.45))
      p = p + scale_x_continuous(minor_breaks = seq(0 , 10, 1),breaks=seq(0,10,1))
    }else{
      p = ggplot(tplot,aes(x=axe3,y=value)) + geom_line(aes(colour=tplot$stat,size=tplot$stat))
      p = p + scale_color_manual('',values=c("goldenrod","blue",'blue',"brown"))
      p = p+scale_size_manual('',values = 2*(1+c(3,1,1,1)))
      p = p + theme_bw()+theme(legend.text = element_text(size=35),
                               legend.title = element_text(size=30),
                               panel.grid.major=element_line(color="grey",size=1),
                               axis.title.x=element_text(size=35),
                               axis.title.y=element_text(size=35),
                               plot.title=element_text(size=40),
                               axis.text.x = element_text(size=25),
                               axis.text.y = element_text(size=25),
                               panel.border = element_blank(),
                               axis.line.x = element_line(color="black", size = 3),
                               axis.line.y = element_line(color="black", size = 3))
      p = p + ylab('Density on [-5,5]')
      p = p + xlab('Environmental variable x')
      p = p + scale_y_continuous(limits=c(0,.35))
      p = p + scale_x_continuous(minor_breaks = seq(-5 , 5, 1),breaks=seq(-5,5,1))
    }
    
    
    setwd(saveDirectory)
    png(paste(runName,"EST_",est,'_TrueEff_',eff,'_estClass_',class,'_meshSize_',msh,'_nPts_',n,'.png',sep=""),width=1500,height=1500)
    print(p)
    dev.off()
    
  }else{
    print('Error: Stat names not ok')
  }
  
}



# approximate Fun with a step function over an irregular 
# 1 dimensional mesh defined by Breaks  
approx_Fun = function(Fun,vec,Breaks,resolution=1e-3,...){
  M = max(vec,na.rm=T)
  m = min(vec,na.rm=T)
  if(M>max(Breaks) | m<min(Breaks)){
    print('input out of bounds.')
  }else{
    f = rep(NA,length(Breaks)-1)
    for(i in 1:(length(Breaks)-1)){
      f[i] = mean(Fun(seq(Breaks[i],Breaks[i+1],resolution),...))
    }
    index = rep(NA,length(Breaks)-1)
    values = sapply(1:length(vec),function(j){
      index_inf = vec[j]>=Breaks[1:(length(Breaks)-1)]
      index = max(which(index_inf))
      return(f[index])
    })
    return(values)
  }
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}





# approximate Fun with a step function: Values are the mean of 
# Fun over each interval of the form [x_min + q*k , x_min + q*(k+1)[  
approx_Fun_regular = function(Fun,vec,q,x_min,resolution=1e-3,...){
  M = max(vec,na.rm=T)
  quo = (M - x_min)%/% q
  if (is.logical(all.equal(quo, (M - x_min)/q ,tolerance = 1e-7))){
    x_max =  quo + x_min
  }else{
    x_max = quo+1 + x_min
  }
  breaks = seq(x_min,x_max,q)
  f = rep(NA,length(breaks)-1)
  for(i in 1:(length(breaks)-1)){
    f[i] = mean(Fun(seq(breaks[i],breaks[i+1],resolution),...))
  }
  index = rep(NA,length(breaks)-1)
  values = sapply(1:length(vec),function(j){
    index_inf = vec[j]>=breaks[1:(length(breaks)-1)]
    index = max(which(index_inf))
    return(f[index])
  })
  return(values)
}

# approximate Fun with a step function over an irregular 
# 1 dimensional mesh defined by Breaks  
approx_Fun = function(Fun,vec,Breaks,resolution=1e-3,...){
  M = max(vec,na.rm=T)
  m = min(vec,na.rm=T)
  if(M>max(Breaks) | m<min(Breaks)){
    print('input out of bounds.')
  }else{
    f = rep(NA,length(Breaks)-1)
    for(i in 1:(length(Breaks)-1)){
      f[i] = mean(Fun(seq(Breaks[i],Breaks[i+1],resolution),...))
    }
    index = rep(NA,length(Breaks)-1)
    values = sapply(1:length(vec),function(j){
      index_inf = vec[j]>=Breaks[1:(length(Breaks)-1)]
      index = max(which(index_inf))
      return(f[index])
    })
    return(values)
  }
}

step_fun = function(x,breaks,values){
  y = sapply(1:length(x),function(i){
    values[which(x[i]<breaks[2:length(breaks)]  & x[i]>=breaks[1:(length(breaks)-1)])[1]] })
  if( is.null(y) & length(y)<length(x)){
    print('x out of breaks range')
  }else{return(y)}
}              

# Simulate an IPP from intensity raster and number of points
simu_IPP = function(r_int,N){
  e = extent(r_int)
  Max = max(getValues(r_int),na.rm=T)
  df = matrix(0,0,1)
  while(dim(df)[1]<N){
    tmp = data.frame(x=runif(2*N,e[1],e[2]),y=runif(2*N,e[3],e[4]),z=runif(N,0,Max))
    v = extract(r_int,tmp[,c('x','y')])
    tmp = tmp[!is.na(v) & tmp$z<v,]
    if(dim(df)[1]==0){df = tmp[,c('x','y'),drop=F]}else{df=rbind(df,tmp[,c('x','y')])}
  }
  if(dim(df)[1]>N){df=df[1:N,]}
  colnames(df)=c('Longitude','Latitude')
  return(df)
}

# Draw points according to a IPP over a continuous spatial domain, 
# based on rasterizable intensity values 
DrawPts = function(grid,nPts){
  # Create raster from data.frame
  r_int = rasterFromXYZ(grid[,c('Longitude','Latitude','int')]) 
  # Draw points from Poisson process (function in the bottom of script)
  df=simu_IPP(r_int,nPts) 
  return(df)
}

# Get quadrat identifier (q_hash) from points coordinates (Longitude, Latitude)
get_q_hash = function(pts,q,x_0,y_0){
  q_x = (pts$Longitude - x_0) %/% q
  q_y = (pts$Latitude - y_0) %/% q
  pts$q_hash = paste( as.character(q_x),as.character(q_y) )
  return( pts )
}

get_q_hash_aniso = function(pts,step_x,step_y,x_0,y_0){
  q_x = (pts$Longitude - x_0) %/% step_x
  q_y = (pts$Latitude - y_0) %/% step_y
  pts$q_hash = paste( as.character(q_x),as.character(q_y) )
  return( pts )
}

get_hash = function(tab,breaks){
  hashes = as.character(1:(length(breaks)-1))
  tab$q_hash = sapply(1:dim(tab)[1],function(i){
    hashes[which(tab$Longitude[i]<breaks[2:length(breaks)]  & tab$Longitude[i]>=breaks[1:(length(breaks)-1)])[1]] })
  return(tab)
}

### duplicate and attach background points to data for LOF 
# (i) Complete occurrences data.frame with given background points
# duplicated into all species
# (ii) compute weights and y
bind.background.pts = function(occ,bg,factorCols=NULL,Bias_killer=100,equalSizeCells=T){
  occ = occ[complete.cases(occ),]
  occ$taxa = factor(occ$taxa)
  taxas = levels(occ$taxa)
  bg$taxa = factor(NA,levels=taxas)
  n_0=dim(bg)[1]
  # Attribute likelihood weights
  q_hashes = unique(c(occ$q_hash,bg$q_hash))
  if(equalSizeCells){
    for(qh in q_hashes){
      inCell = bg$q_hash==qh
      nInCell = sum(inCell)
      bg$w[inCell] = (Bias_killer-1) / (n_0*Bias_killer*nInCell)
    }
  }else{
    print('Error: Non-equal cell size not implemented yet')
    return(NULL)
  }
  
  # background output value is zero
  bg$pseudo_y=0
  n=NULL
  occ$pseudo_y = NA
  occ$w = NA
  for(e in 1:length(taxas)){
    sp = occ$taxa==taxas[e]
    n = sum(sp)
    occ$w[sp] = 1/(Bias_killer*n)
    occ$pseudo_y[sp] = (Bias_killer*n)/1
  }
  # Bind occurrences and background points
  # with a repetition of background points per species 
  for(e in 1:length(taxas)){
    bg$taxa=taxas[e]
    occ=rbindlist(list(occ,bg))
    if(e/20==round(e/20)){
      flush.console()
      cat('     \r    Process...',100*e/length(taxas),'%       \r')
    }
  }
  if(length(factorCols)>0){
    for(c in factorCols){
      eval(parse(text=paste('occ$',c,'=factor(occ$',c,')',sep="")))
    }
  }
  return(occ)
}


### Low memory + scalable wrapper of sparse.model.matrix 
memFree.sparse.model.matrix= function(matrixFormu,Data,band=50000){
  quo = dim(Data)[1] %/% band
  res = dim(Data)[1] %% band
  npass = quo+as.numeric(res>0)
  store = data.frame(is=NA,js=NA,xs= NA)
  store = store[-1,]
  deb= Sys.time()
  for(i in 1:npass){
    if(i==npass & res>0){
      tmp = Data[(1+(i-1)*band):dim(Data)[1],,drop=F]
    }else{
      tmp = Data[(1+(i-1)*band):(i*band),,drop=F]
    }
    
    mTmp = sparse.model.matrix(matrixFormu,tmp)
    mTmp = as(mTmp, "dgTMatrix")
    store = rbindlist(list(store,data.frame(is=mTmp@i,js=mTmp@j,xs=mTmp@x)))
    dur = Sys.time()-deb
    gc(reset=T)
    if(i/5==round(i/5)){
      flush.console()
      cat('\r     Process...',100*i/npass,' time elapsed:',dur,attr(dur,'units'),'       \r')
    }
  }
  BGf = sparseMatrix(i=store$is+1,j=store$js+1,x=store$xs,dims=c(dim(Data)[1],dim(mTmp)[2]))
  colnames(BGf) = colnames(mTmp)
  gc(reset=T)
  return(BGf)
}


### Fit LOF (GLMNET version)
# fit LOF model parameters to data with the glmnet package
lof_glmnet=function(data,y,occupationString,weights=NULL,lambdaMinRatio=NULL,nlambda=100){
  
  glmnet.control(fdev=0,eps=1e-10)
  # Reset with
  #glmnet.control(factory = TRUE)
  
  if(is.null(weights)){weights=rep(1,dim(data)[1])}
  
  ntaxa = length(unique(data$taxa))
  # Make design matrix
  data$q_hash = factor(data$q_hash)
  contrasts(data$q_hash) = contr.sum(length(levels(data$q_hash)))
  
  if(ntaxa>1){
    matrixFormu = paste(' ~ q_hash + taxa * (',occupationString,')')
  }else{
    matrixFormu = paste(' ~ q_hash + ',occupationString)
  }
  matrixFormu = as.formula(matrixFormu)
  SparseDes = memFree.sparse.model.matrix(matrixFormu,data,band=100000)
  taxaItc= intersect( grep('taxa',colnames(SparseDes)) , grep('I(.*)',colnames(SparseDes),invert=T))
  
  # GLMNET formula
  if(ntaxa>1){
    formu = paste('pseudo_y ~ q_hash + taxa * (',occupationString,')')
  }else{
    formu = paste('pseudo_y ~ q_hash + ',occupationString)
  }
  formu = as.formula(formu )
  
  # Attribute penalties to groups of terms
  coNa = colnames(SparseDes)
  penaltyFactorsIdx = 0 * as.numeric(regexpr('(Intercept)',coNa)>0) + 1*as.numeric(regexpr('q_hash',coNa)>0) + 2*as.numeric(grepl('I\\((.*)\\)',coNa))
  penaltyFactors = rep(1,length(penaltyFactorsIdx))                                                                         
  penaltyFactors[penaltyFactorsIdx==0] = 1e-3
  penaltyFactors[penaltyFactorsIdx==1] = 1
  
  # Fit model
  if(is.null(lambdaMinRatio)){
    lof.mod = glmnet(x=SparseDes,y=y,family="poisson",weights = weights,penalty.factor=penaltyFactors,nlambda = nlambda,standardize = T)
  }else{
    lof.mod = glmnet(x=SparseDes,y=y,family="poisson",weights = weights,penalty.factor=penaltyFactors,lambda.min.ratio = lambdaMinRatio,nlambda = nlambda,standardize = T)
  }
  
  # Add q_hash contrasts detail
  lof.mod$q_hashContrasts = contrasts(data$q_hash)
  colnames(lof.mod$q_hashContrasts)=as.character(1:dim(contrasts(data$q_hash))[2])
  
  # Compute coefficients translation into "contr.treatment" 
  Contr = lof.mod$q_hashContrasts
  nQi = dim(Contr)[2]
  nQ = dim(Contr)[1]
  coefo = lof.mod$beta[,dim(lof.mod$beta)[2]]
  coefo_qc = coefo[regexpr('q_hash',names(coefo))>0]
  orderedIdx = sapply(1:nQi,function(k) which(paste('q_hash',k,sep="")==names(coefo_qc)) )
  coefo_qc=coefo_qc[orderedIdx]
  Qnames = rownames(Contr)
  coefo_q = rep(NA,nQ)
  names(coefo_q) = Qnames
  for(qi in 1:nQ){
    coefo_q[qi]= sum(Contr[qi,] * coefo_qc)
  }
  diffRef = 0-coefo_q[1]
  coefo_q = coefo_q + diffRef
  names(coefo_q)= paste('q_hash',names(coefo_q),sep="")
  newCoefo = coefo[regexpr('q_hash',names(coefo))<=0]
  newCoefo['(Intercept)'] = lof.mod$a0[length(lof.mod$a0)]-diffRef
  newCoefo = c(newCoefo,coefo_q[2:length(coefo_q)])
  lof.mod$coefficients = newCoefo
  
  return(lof.mod)
}

# Mid effort
mid = function(x){
  0.9 * as.numeric(x<0 ) + 0.1
}
# Sigmo rapide
SigmoFast = function(x){
  9* exp(-x*20)/(1+exp(-x*20)) + 1
}
# Coline ecretee 
scalpHill = function(x){
  1 + as.numeric(x<0)*3* exp((x+2)*20)/(1+exp((x+2)*20)) +  as.numeric(x>=0)*3* exp(-(x-2)*20)/(1+exp(-(x-2)*20))
}
# Vallee
valley = function(x){
  1 + as.numeric(x<0)*3* exp(-(x+2)*20)/(1+exp(-(x+2)*20)) +  as.numeric(x>=0)*3* exp((x-2)*20)/(1+exp((x-2)*20))
}
# Cuts where it hurts
cutHurts = function(x){
  1+ 5 * as.numeric( (x>=-4.5 & x<(-2.5)) | (x>=-.5 & x<1.5) | (x>=2.5 & x<4.5))
}
# Cuts Nice
cutNice = function(x){
  1+ 3 * as.numeric(x>=(-4) & x<(-3)) +  7*as.numeric(x>=0 & x<1) +  2*as.numeric(x>=3 & x<4)
}

#####
# Functions for real data experiment
#####

# catalogue des variables environnementales disponibles
load_variables=function(){
  variables=list()
  # variables extraites par la méthode 1 :
  # directement depuis le raster
  # sources données :  WorldClim 1.4,Chelsea 1.1, ETP , SRTM2010
  variables[[1]] = c( 
    "etp",   
    "bio_1",               "bio_2",               "bio_3",               "bio_4",               "bio_5",              
    "bio_6",               "bio_7",               "bio_8",               "bio_9",               "bio_10",             
    "bio_11",              "bio_12",              "bio_13",              "bio_14",              "bio_15",             
    "bio_16",              "bio_17",              "bio_18",              "bio_19",              
    "chbio_1",             "chbio_2",             "chbio_3",             "chbio_4",             "chbio_5",
    "chbio_6",             "chbio_7",             "chbio_8",             "chbio_9",             "chbio_10", 
    "chbio_11",            "chbio_12",            "chbio_13",            "chbio_14",            "chbio_15",
    "chbio_16",            "chbio_17",            "chbio_18",            "chbio_19",
    "alti" ,               "slope",               "shade",               "droute_fast",         "dmer_fast",
    "proxi_eau_fast",
    paste('axe',1:10,sep=""))
  # variables extraites par la méthode 2 : 
  # projection des coordonnées sur le système européen, 
  # substitution des valeurs du raster selon la règle dictée par attrib.csv 
  # gestion particulière de "text"
  # sources données :  ESDB v2.0 => PTRDB
  variables[[2]] =  c("awc_top",            "bs_top",              "cec_top",             "crusting",            "dgh",                
                      "dimp", "erodi",     "oc_top", "pd_top",   "vs", "text" ,'text_orga','text_autre','text_roche','text_noinfo')
  # variables extraites par la méthode 3 : 
  # projection des coordonnées sur le système français IGN
  # sources données :  BD Carthage
  variables[[3]] =  c("deau","dmer","proxi_eau")
  # variables extraites par la méthode 4 : 
  # sources données :  ROUTE500, IGN
  variables[[4]] = c("droute")
  # variables extraites par la méthode 5 :
  # extraction directe depuis le raster
  # création de variables binaires, une par catégorie.
  # sources données :  CLC2012
  variables[[5]] = c(
    "discont_ur_fab",      "indus",               "road_or_rail",        "port",                "airport",            
    "mine",                "dump",                "construction_site",   "green_urb",           "sport_facilit",      
    "arable_nonirrig",     "arable_irrig",        "rice",                "vine",                "fruit",              
    "olive",               "pasture",             "annual_crop",         "complex_culti",       "agri_and_nat",       
    "agroforestry",        "brl_forest",          "coni_forest",         "mixed_forest",        "nat_grass",          
    "moors_heat",          "sclero_veg",          "transi_wood",         "sand",                "rock",               
    "sparse",              "burnt",               "snow",                "marsh",               "peat",             
    "salt_marsh",          "saline",              "intertid_flat",       "wat_course",          "wat_body",           
    "coast_lag",           "estua",               "ocean",               "cont_urb_fab" )
  # variable extraite par la méthode 6 : 
  # extraction directe depuis le raster CLC dont les coordonnées sont en système : 
  # "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m"
  variables[[6]]=c('clc')
  return(variables)
}

# attribution des variables environnementales aux occurrences à partir d'un data.frame
# où chaque ligne est une occurrence et sont présentes les colonnes "Longitude" et "Latitude" au format WGS84
get_variables = function(variables,table,dpath='C:/Users/Christophe/hubiC/Documents/0_These/data/',d_NA=2000,fix=T){
  # variables : vecteur de chaînes de caractères, contient les variables à attribuer à chaque occurence selon  
  # la nomenclature fournie par load_variables
  
  # table : data.frame, au moins deux colonnes "Longitude" et "Latitude" (en système de coordonnées WGS84)
  # chaque ligne représente une occurence 
  
  # dpath : chaîne de caractères, adresse du repertoire contenant le dossier "0_mydata" 
  
  # d_NA : numeric , nombre de mètres maximal pour lequel on va chercher la valeur la plus proche non NA.
  method = load_variables()
  loaded_data = list(clc=0,eau=0,cours_eau=0,rd=0,text=0)
  for (var in variables){
    if(var%in%method[[1]]){
      print(var)
      # meth1
      r = raster(paste(dpath,'0_mydata/',var,'/',var,'.tif',sep=""))
      
      # extraction des valeurs des occcurrences
      print('extracting values...')
      
      table[,var] = extract( r, data.frame(table$Longitude,table$Latitude) )
      #command = "extract( r, data.frame(table$Longitude,table$Latitude) )"
      #eval(parse(text=paste('table$',var,'=',command,sep="")))
      
      if(fix & sum(is.na(table[,var]))>0){
        # correction des occurences cotières à une tolérance de 2km
        print('correcting values...')
        #vals[is.na(vals)] = get_value_fromXY(r, table[is.na(vals),] , tol_NA = T , d_to_closest_if_NA = 2000 )
        if (var =="dmer_fast"){
          #vals[is.na(vals)] = 32767
          table[is.na(table[,var]),var] = 32767
        }else{
          #vals[is.na(vals)] = get_value_fromXY_par(r,table[is.na(vals),],d_to_closest_if_NA = d_NA)
          table[is.na(table[,var]),var] = get_value_fromXY_par(r,table[is.na(table[,var]),],d_to_closest_if_NA = d_NA)
        }
      }
      #eval(parse(text=paste('table$',var,'=vals',sep="")))
      rm(r)
      gc(reset=T)
      print(paste(var,' DONE.                    '))
    }else if(var%in%method[[2]]){
      print(var)
      # meth2
      text_vars = c('text','text_orga','text_roche','text_autre','text_noinfo')
      condition_generale = loaded_data$text==0 | !(var%in%text_vars)
      if(condition_generale){
        
        if(var %in% text_vars){
          path = paste(dpath,'/0_mydata/text',sep="")
        }else{
          path = paste(dpath,'/0_mydata/',var,sep="")
        }
        x = new("GDALReadOnlyDataset", path)
        getDriver(x)
        getDriverLongName(getDriver(x))
        xx=asSGDF_GROD(x)
        r = raster(xx)
        crs(r) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m" 
        print('projecting occurences..')
        coos = data.frame(table$Longitude,table$Latitude)
        pts_occ = SpatialPointsDataFrame( coords=coos,data=coos,proj4string = CRS("+proj=longlat +datum=WGS84"))
        coords = spTransform( pts_occ , CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m") )@coords
        proj_table = data.frame( Longitude =coords[,1] , Latitude =coords[,2]  )
        rm(coords,pts_occ,coos)
        print('cropping France...')
        lims <- extent(3200000,4350000,1500000,3200000)
        r2 <- crop(r, lims)
        print('extracting values...')
        vals = extract( r2 ,  proj_table )
        if(fix & sum(is.na(vals))>0){
          print('correcting values...')
          vals[is.na(vals)] = get_value_fromXY_par(r2, proj_table[is.na(vals),] , d_to_closest_if_NA = d_NA )
        }
        rm(r,r2)
        gc(reset=T)
        setwd(path)
        attrib = read.csv('attrib.csv',sep=";",header=T)
        attrib = data.frame(attrib)
        print('got values')
        new_vals = vals
        # Remplacer les valeurs par la valeur QUANTI associée
        for (cate in attrib$VALUE){
          new_vals[vals==cate] = attrib$QUANTI[attrib$VALUE==cate]
        }
        
        if( var%in%text_vars ){
          if('text'%in%variables){
            # chaque catégorie non pertinente pour la texture se voit associée 
            # - 0 si la variable binaire est créée par ailleurs : Ainsi, la valeur 0 stocke
            # l'ensemble des valeurs non pertinentes sans affecter l'estimation
            # du paramètre de "text", car une autre variable binaire est ajustée
            # pour ces occurrences
            # - NA sinon, ce qui aboutira à éliminer cette occurrence
            table$text= new_vals
            if('text_noinfo'%in%variables){table$text[vals==0]=0}else{table$text[vals==0]=NA}
            if('text_orga'%in%variables){table$text[vals==8]=0}else{table$text[vals==8]=NA}
            if('text_roche'%in%variables){table$text[vals==7]=0}else{table$text[vals==7]=NA}
            if('text_autre'%in%variables){table$text[vals==6]=0}else{table$text[vals==6]=NA}
          }
          if('text_noinfo'%in%variables){
            table$text_noinfo=F
            table$text_noinfo[new_vals==0] = T
          }
          if('text_orga'%in%variables){
            table$text_orga=F
            table$text_orga[new_vals==8]= T
          }
          if('text_roche'%in%variables){
            table$text_roche=F
            table$text_roche[new_vals==7]=T
          }
          if('text_autre'%in%variables){
            table$text_autre=F
            table$text_autre[new_vals==6]=T
          }
          loaded_data$text=1
        }else{
          eval( parse(text= paste("table$",var,'=new_vals',sep=""))  )
          print(paste(var,' DONE.                    '))
        }
        
      }
      
    }else if(var%in%method[[3]]){
      print(var)
      # meth3
      path = paste(dpath,'0_mydata/hydro',sep="")
      
      # Projection des coordonnées occurences
      print('projecting occurences...')
      coos = data.frame(table$Longitude,table$Latitude)
      pts_occ = SpatialPointsDataFrame( coords=coos,data=coos,proj4string = CRS("+proj=longlat +datum=WGS84"))
      pts_occ=spTransform(pts_occ,CRS("+proj=lcc +lat_1=44 +lat_2=49 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"))
      
      
      # Chargement données requises si pas déjà chargées
      print('loading data...')
      if ((var%in% c('deau','proxi_eau')) & loaded_data$cours_eau==0){
        cours_eau = readOGR(dsn = path,layer="COURS_D_EAU", stringsAsFactors=FALSE)
        loaded_data$cours_eau=1
      }
      if(var%in%c('dmer','deau','proxi_eau') & loaded_data$eau==0){
        eau = readOGR(dsn = path,layer="HYDROGRAPHIE_SURFACIQUE", stringsAsFactors=FALSE)
        loaded_data$eau=1
      }
      
      # préparation données
      print('preparing data...')
      if (var %in% c('deau','proxi_eau')){
        eau_douce = eau[eau$NATURE=="Eau douce permanente",]
        liste = list( pt= pts_occ,cours_eau=cours_eau,eau_douce=eau_douce,gDistance=gDistance)
      }else if(var=='dmer'){
        mer = eau[eau$TYPE=="Pleine mer",]
        liste = list( pt= pts_occ,mer=mer,gDistance=gDistance)
      }
      # Fonction de calcul de la VE
      if(var=='deau'){
        FUNC = function(x,pt,cours_eau,eau_douce,gDistance){ min( gDistance(pt[x,],cours_eau) , gDistance(pt[x,],eau_douce) ) }
      }else if(var=='proxi_eau'){
        FUNC = function(x,pt,cours_eau,eau_douce,gDistance){ 
          d = min( gDistance(pt[x,],cours_eau) , gDistance(pt[x,],eau_douce) )
          if(d<25){return(TRUE)}else{return(FALSE)}  }
      }else if(var=='dmer'){ FUNC = function(x,pt,mer,gDistance){gDistance(pt[x,],mer)} }
      
      print('calculating variable..')
      vals = fragmented_parallel( dim(table)[1] , FUNC , list_arg = liste, n_jobs = 2 , type = 'vector')
      eval( parse(text= paste("table$",var,'=vals',sep=""))  )
      # nettoyage des variables
      rm(pts_occ)
      print(paste(var,' DONE.                    '))
    }else if(var%in%method[[4]]){
      print(var)
      # meth4
      print('projecting occurences...')
      coos = data.frame(table$Longitude,table$Latitude)
      pts_occ = SpatialPointsDataFrame( coords=coos,data=coos,proj4string = CRS("+proj=longlat +datum=WGS84"))
      pts_occ=spTransform(pts_occ,CRS("+proj=lcc +lat_1=44 +lat_2=49 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"))
      
      print('loading data...')
      path = paste(dpath,'/0_mydata/routes',sep="")
      routes = readOGR(dsn = path,layer="ROUTES", stringsAsFactors=FALSE)
      
      print('preparing data...')
      liste = list(pt=pts_occ,rd=routes,dist=gDistance)
      FUNC = function(i,pt,rd,dist){dist(pt[i,],rd)}
      
      print('calculating variable...')
      vals= fragmented_parallel(dim(table)[1],FUNC,list_arg = liste,n_jobs = 2 ,type = 'vector')
      eval( parse(text= paste("table$",var,'=vals',sep=""))  )
      rm(routes,pts_occ)
      gc(reset=T)
      print(paste(var,' DONE.                    '))
    }else if(var%in%method[[5]] & loaded_data$clc==0){
      print(variables[variables %in% method[[5]]])
      # meth5: binary category of 'clc' variable 
      r= raster( paste(dpath,'0_mydata/clc/g100_clc12_v18_5_FR_WGS84.tif',sep=''))
      proj = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m"
      crs(r) <-  proj
      coos = data.frame(table$Longitude,table$Latitude)
      pts_occ = SpatialPointsDataFrame( coords=coos,data=coos,proj4string = CRS("+proj=longlat +datum=WGS84"))
      coords = spTransform( pts_occ , proj  )@coords
      proj_table = data.frame( Longitude =coords[,1] , Latitude =coords[,2]  )
      rm(coords,pts_occ,coos)
      print('extracting values...')
      vals = extract(r, proj_table )
      if(fix & sum(is.na(vals))){
        print('correcting values...')
        vals[is.na(vals)] = get_value_fromXY_par(r,proj_table[is.na(vals),],d_to_closest_if_NA = 2000 )
      }
      print('transforming to binary variables...')
      setwd(paste(dpath,'0_mydata/clc',sep=""))
      legend  = read.csv('clc_legend.csv',sep=";",header=T)
      legend = legend[,c(5,1)]
      legend[,1] = as.character(legend[,1])
      legend = legend[ legend[,1] %in% variables ,]
      df = quanti_to_category(vals,legend)
      table = cbind(table,df)
      rm(proj_table,r)
      gc(reset=T)
      print('variable(s) DONE.')
      loaded_data$clc=1
    }else if(var%in%method[[6]]){
      print(variables[variables %in% method[[6]]])
      # meth6 : variable "clc" toutes catégories
      r= raster( paste(dpath,'0_mydata/clc/g100_clc12_v18_5_FR_WGS84.tif',sep=''))
      proj = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m"
      crs(r) <-  proj
      coos = data.frame(table$Longitude,table$Latitude)
      pts_occ = SpatialPointsDataFrame( coords=coos,data=coos,proj4string = CRS("+proj=longlat +datum=WGS84"))
      coords = spTransform( pts_occ , proj  )@coords
      proj_table = data.frame( Longitude =coords[,1] , Latitude =coords[,2]  )
      rm(coords,pts_occ,coos)
      print('extracting values...')
      vals = extract(r, proj_table )
      if(fix & sum(is.na(vals))>0){
        print('correcting values...')
        vals[is.na(vals)] = get_value_fromXY_par(r,proj_table[is.na(vals),],d_to_closest_if_NA = 2000 )
      }
      eval( parse(text= paste("table$",var,'=vals',sep=""))  )
    }else{
      print(paste(var,' est inconnue au bataillon -> Se réferrer à la liste des variables reconnues.'))
    }
  }
  return(table)
}


# Change land cover to spht := simplified plant habitat type = urban / arable / grasses / forest / other
get_spht = function(clcVec){
  setwd(paste('C:/Users/',user,'/pCloud local/0_These/data/MTAP article/data/',sep=""))
  clc.variables = readRDS('clc.variables.RData')
  names = clc.variables[[1]]
  
  spht = rep(NA,length(clcVec))
  cd = clcVec %in% clc.variables[[2]][[which(names=='arable')]] | clcVec%in%c(12,13,15,16,20) 
  spht[cd] = "cultivated"
  cd = clcVec %in% clc.variables[[2]][[which(names=='pasture')]] | clcVec %in% clc.variables[[2]][[which(names=='nat_grass')]] | clcVec %in% clc.variables[[2]][[which(names=='moors')]] | clcVec %in% clc.variables[[2]][[which(names=='sclero')]]
  spht[cd] = "grasses"
  cd = clcVec %in% clc.variables[[2]][[which(names=='brl_for')]] | clcVec %in% clc.variables[[2]][[which(names=='coni_for')]] | clcVec %in% clc.variables[[2]][[which(names=='mixed_for')]] | clcVec %in% clc.variables[[2]][[which(names=='transi_wood')]]
  spht[cd] = "forest"
  cd = clcVec %in% clc.variables[[2]][[which(names=='arti')]] | clcVec %in% clc.variables[[2]][[which(names=='semi_arti')]] | clcVec%in%c(11)
  spht[cd] = "urban"
  spht[is.na(spht)] = "other"
  return(spht)
}

