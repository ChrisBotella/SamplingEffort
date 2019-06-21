# created on 29/11/18
# Functions for LOF-GLMNET and diagnostic

#####
# Simulation
#####

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

## Draw a large set of points representing the whole geographic domain 
draw.pts.in.cells = function(n,r_lofCells,originalVes,delta,dataDir,miniPerCell=1){
  resX= (extent(r_lofCells)[2]-extent(r_lofCells)[1])/ncol(r_lofCells)
  resY = (extent(r_lofCells)[4]-extent(r_lofCells)[3])/nrow(r_lofCells)
  cellsCoo = coordinates(r_lofCells)
  cellsCoo = cellsCoo[!is.na(getValues(r_lofCells)),]
  
  cellsIds = unique(getValues(r_lofCells))
  cellsIds = cellsIds[!is.na(cellsIds)]
  
  q_hashes = expand.grid(1:n,cellsIds)[,2]
  x=rep(0,length(cellsIds)*n)
  y=rep(0,length(cellsIds)*n)
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
      nLack = sum(toFill)
      
      cd = cellsIds==qLack[i]
      xmin = cellsCoo[cd,1]-resX/2+delta
      xmax = cellsCoo[cd,1]+resX/2-delta
      ymin = cellsCoo[cd,2]-resY/2+delta
      ymax = cellsCoo[cd,2]+resY/2-delta
      
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
    qLack = unique(q_hashes[toReDraw])
    for(i in 1:length(qLack)){
      cd = q_hashes==qLack[i] & toReDraw
      if(sum(cd)<=(n-miniPerCell)){
        # We do NOT redraw points whose cells for
        # which at least one background point have
        # been drawn
        toReDraw[cd] = F
      }
    }
  }
  Df$spht = get_spht(Df$clc)
  return(Df)
}

# Get points coordinates of "Df" over their PCA axeses
# It will include all PCA axeses until the sum of 
# their Eigen values makes a proportion of "eigRatio" of the trace
principal.axes.coordinates=function(matrixFormula,Df,d=NULL,eigRatio=NULL){
  Des = model.matrix(matrixFormula,Df)
  Des = Des[,2:dim(Des)[2]]
  Des = scale(Des)
  SVD = svd(Des)
  if(is.null(d) & !is.null(eigRatio)){
    d = 0
    sumEig = 0
    while(sumEig/sum(SVD$d)<eigRatio){
      d = d+1
      sumEig = sumEig + SVD$d[d]
    }
  }else if(is.null(d) & is.null(eigRatio)){
    print('Error: either supply d or eigRatio')
    return(NULL)
  }
  PAcoos = SVD$u[,1:d]
  PAcoos = t(SVD$d[1:d] * t(PAcoos))
  colnames(PAcoos)=paste('axe_',1:d,sep="")
  PAcoos=data.frame(PAcoos)
  return(PAcoos)
}

# Approximate the total hyper-volume of the Alpha-confidence interval
# envelope of given points as an ellipsoïde, with their coordinates
# over the main inerty axeses of the ellipsoïde
hyperVolume.from.PAcoos = function(PAcoos,Alpha){
  axesLength = NULL
  d = dim(PAcoos)[2]
  for(di in 1:d){
    quantiles = quantile(PAcoos[,di],c(0.5-Alpha/2,.5+Alpha/2),na.rm=T)
    axesLength = c(axesLength,quantiles[2]-quantiles[1])
  }
  # We now approximate the Alpha-enveloppe volume
  # of the whole domain with an ellipsoïd volume
  Gamma = 1/ (exp(1)*dgamma(1,shape=(d/2)+1,scale=1) ) 
  volume = prod(axesLength)* pi^(d/2) / Gamma 
  return(volume)
}

# Mid effort
mid = function(x){
 as.numeric(x<0 ) 
}
# Sigmo rapide
SigmoFast = function(x){
  9* exp(-x*20)/(1+exp(-x*20)) + 1
}

# Sigmo rapide
SigmoMedium = function(x){
  9* exp(-x*5)/(1+exp(-x*5)) + 1
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

clc.groups = function(){
  V1 = c("arti","semi_arti","arable","pasture",
         "brl_for","coni_for","mixed_for","nat_grass",
         "moors","sclero","transi_wood","no_veg",      
         "coastal_area","ocean" )
  liste = list( V1, list(c(1,10),c(2,3,4,6),
  c(21,22),18,23,24,25,26, 27,28,29,
  c(31,32),c(37,38,39,42,30),c( 44) ))
  
  return(liste)
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
get_spht = function(clcVec,dir=NA){
  if(!is.na(dir)){
    setwd(dir)
  }
  clc.variables = clc.groups()

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



#####
# LOF core
#####


# Duplicate and attach background points to data for LOF 
# (i) Complete occurrences data.frame with given background points
# duplicated into all species
# (ii) compute weights and y
bind.background.pts = function(occ,bg,factorCols=NULL,Bias_killer=100){
  ### PRINT
  print('check if na in occ lines')
  print(paste('n lines without NAs:',sum(complete.cases(occ))))
  
  occ = occ[complete.cases(occ),]
  occ$taxa = factor(occ$taxa)
  taxas = levels(occ$taxa)
  
  q_hashes = union( as.character(unique(bg$q_hash)), as.character(unique(occ$q_hash)) )
  occ$q_hash = factor(occ$q_hash,levels=q_hashes)
  bg$q_hash = factor(bg$q_hash, levels=q_hashes)
  ### PRINT
  print('check taxas')
  print(str(taxas))
  print('check occ$taxa')
  print(str(occ$taxa))
  
  # Compute WEIGHTS of background points 
  # We assume here all cells have same area
  nCell = length(q_hashes) 
  bg$w= NA
  q_hashCol= as.character(bg$q_hash)
  for(qh in q_hashes){
    inCell = q_hashCol==qh
    nInCell = sum(inCell)
    bg$w[inCell] = (Bias_killer-1) / (Bias_killer*nCell*nInCell)
  }
  
  # background output value is zero
  bg$pseudo_y=0
  occ$pseudo_y = NA
  occ$w = NA
  TaxaCol = as.character(occ$taxa)
  for(e in 1:length(taxas)){
    sp = TaxaCol==taxas[e]
    n_e = sum(sp)
    if(is.na(n_e)){
      print(taxas[e])
      print(sum(sp))
      print(sum(TaxaCol==as.character(taxas[e])))}
    occ$w[sp] = 1/(Bias_killer*n_e)
    occ$pseudo_y[sp] = (Bias_killer*n_e)/1
  }
  # Bind occurrences and background points
  # with a repetition of background points per species 
  for(e in 1:length(taxas)){
    bg$taxa=factor(taxas[e],levels=taxas)
    occ=rbindlist(list(occ,bg),use.names=T)
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
    store = rbindlist(list(store,data.frame(is=mTmp@i+1+(i-1)*band,js=mTmp@j+1,xs=mTmp@x)))
    dur = Sys.time()-deb
    gc(reset=T)
    if(i/5==round(i/5)){
      flush.console()
      cat('\r     Process...',100*i/npass,' time elapsed:',dur,attr(dur,'units'),'       \r')
    }
  }
  BGf = sparseMatrix(i=store$is,j=store$js,x=store$xs,dims=c(dim(Data)[1],dim(mTmp)[2]))
  colnames(BGf) = colnames(mTmp)
  gc(reset=T)
  return(BGf)
}


## Get treatment coefficients of q_hash
get.treat.coef = function(lof.mod,wLambda=NULL){
  if(is.null(wLambda)){wLambda=dim(lof.mod$beta)[2]}
  Contr = lof.mod$q_hashContrasts
  nQi = dim(Contr)[2]
  nQ = dim(Contr)[1]
  coefo = lof.mod$beta[,wLambda]
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
  newCoefo['(Intercept)'] = lof.mod$a0[wLambda]-diffRef
  
  newCoefo = c(newCoefo,coefo_q[2:length(coefo_q)])
  lof.mod$coefficients = newCoefo
  return(lof.mod)
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
  print('make sparse model matrix')
  SparseDes = memFree.sparse.model.matrix(matrixFormu,data,band=100000)
  taxaItc= intersect( grep('taxa',colnames(SparseDes)) , grep('I(.*)',colnames(SparseDes),invert=T))
  
  # Attribute penalties to groups of terms
  coNa = colnames(SparseDes)
  penaltyFactorsIdx = 0 * as.numeric(regexpr('(Intercept)',coNa)>0) + 1*as.numeric(regexpr('q_hash',coNa)>0) + 2*as.numeric(grepl('I\\((.*)\\)',coNa))
  penaltyFactors = rep(1,length(penaltyFactorsIdx))
  penaltyFactors[penaltyFactorsIdx==0] = 1e-3 # Set relative penalty of intercept
  penaltyFactors[penaltyFactorsIdx==1] = .1 # Set relative penalty of sampling effort
  
  # Fit model
  print('call glmnet')
  if(is.null(lambdaMinRatio)){
    lof.mod = glmnet(x=SparseDes,y=y,family="poisson",weights = weights,penalty.factor=penaltyFactors,nlambda = nlambda)
  }else{
    lof.mod = glmnet(x=SparseDes,y=y,family="poisson",weights = weights,penalty.factor=penaltyFactors,lambda.min.ratio = lambdaMinRatio,nlambda = nlambda)
  }
  
  # Add q_hash contrasts detail
  lof.mod$q_hashContrasts = contrasts(data$q_hash)
  colnames(lof.mod$q_hashContrasts)=as.character(1:dim(contrasts(data$q_hash))[2])
  
  # Compute coefficients translation into "contr.treatment" 
  lof.mod =get.treat.coef(lof.mod) 
  return(lof.mod)
}


### predict relative species intensity 
predict.lof.spRelativeIntensity = function(coos,
                                           model,
                                           Intensityformula,
                                           originalVes,
                                           taxasToPredict,
                                           glc19SpIdVec,
                                           nOccPerSp = rep(dim(coos)[1],length(glc19SpId)),
                                           extractorScript = "C:/Users/Christophe/pCloud local/Github/Run/all_functions.R",
                                           dataDir="C:/Users/Christophe/pCloud local/data/"){
  
  # coos : numeric matrix with 2 columns named Longitude and Latitude (coordinates in the WGS84 system)
  listVE = load_variables()
  avail = NULL
  for(i in 1:length(listVE)){avail =  c(avail,listVE[[i]])}
  toGet = intersect(avail,originalVes)
  
  extra = setdiff(originalVes,avail)
  if(length(extra)!=1 | extra!="spht"){
    print('Error: Unknown environmental variables required in originalVes')
    return(NULL)
  }
  coos = get_variables(toGet,coos,dpath = dataDir,fix=F)
  if("clc"%in%toGet & length(extra)==1 & extra=="spht"){coos$spht = get_spht(coos$clc)
  }else{print('Error: spht required whereas clc (original source of it) not in originalVes')}
  
  Mat = sparse.model.matrix(Intensityformula,data = coos)
  predMatrix = matrix(NA,dim(Mat)[1],length(glc19SpIdVec))
  colnames(predMatrix) = as.character(glc19SpIdVec)
  
  for(i in 1:length(glc19SpIdVec)){
    Sp = glc19SpIdVec[i]
    coefSp = coefficients[regexpr(Sp,names(coefficients))>0 & regexpr("q_hash",names(coefficients))<=0]
    if(!is.null(coefSp)){
      names(coefSp) = sub(paste('taxa',Sp,':(.*)',sep=""),'\\1',names(coefSp))
      names(coefSp)[1] = '(Intercept)'
      coefRef = coefficients[regexpr("taxa",names(coefficients))<=0 & regexpr("q_hash",names(coefficients))<=0]
      coefSp = coefSp + coefRef
    }else{
      coefSp = coefficients[regexpr("taxa",names(coefficients))<=0 & regexpr("q_hash",names(coefficients))<=0]
    }
    
    #matrixFormu = as.formula(" ~ etp + I(etp^2) + I(chbio_12-etp) + I((chbio_12-etp)^2) + chbio_1 + I(chbio_1^2) + chbio_5 + I(chbio_5^2) + alti + I(alti^2) + slope + I(slope^2) + awc_top + I(awc_top^2) + bs_top + I(bs_top^2) + spht + slope:I(chbio_12-etp)")
    
    pred = as.numeric(Mat %*% coefSp)
    predMatrix[,i] = pred * nOccPerSp[i] / sum(pred)
  }
  
  if(i/10 == round(i/10)){
    flush.console()
    cat('    \r     Process...',100*i/length(glc19SpIdVec),'%        \r')
  }
  
  # /!\ predMatrix might have NA cells /!\ 
  return(predMatrix)
}




### Get species Intensity coefficients
# OK
get.spIntensityCoefficients = function(SpName,coefficients,spList){
  # SpName : character string, species identifier as in the model variable "taxa"
  # coefficients: named numeric vector, the vector of all estimated coefficients of the model extracted through 
  # the function 
  # spList : character vector, all identifiers of species fitted in the model
  coefSp = coefficients[regexpr(SpName,names(coefficients))>0 & regexpr("q_hash",names(coefficients))<=0]
  if(length(coefSp)>0){
    names(coefSp) = sub(paste('taxa',SpName,':(.*)',sep=""),'\\1',names(coefSp))
    names(coefSp)[1] = '(Intercept)'
    coefRef = coefficients[regexpr("taxa",names(coefficients))<=0 & regexpr("q_hash",names(coefficients))<=0]
    coefSp = coefSp + coefRef
  }else if(as.character(SpName) %in%as.character(spList)){
    coefSp = coefficients[regexpr("taxa",names(coefficients))<=0 & regexpr("q_hash",names(coefficients))<=0]
  }else{
    print('Species not in the list')
    return(NULL)
  }
  return(coefSp)
}

### predict relative species intensity 
# OK
predict.lof.spLogRelativeIntensity = function(coos,
                                              coefficients,
                                              Intensityformula,
                                              originalVes,
                                              taxasToPredict,
                                              glc19SpIdVec,
                                              nOccPerSp = rep(dim(coos)[1],length(glc19SpIdVec)),
                                              extractorScript = "C:/Users/Christophe/pCloud local/0_These/Github/Run/all_functions.R",
                                              dataDir="C:/Users/Christophe/pCloud local/0_These/data/"){
  
  # coos : numeric matrix with 2 columns named Longitude and Latitude (coordinates in the WGS84 system)
  listVE = load_variables()
  avail = NULL
  for(i in 1:length(listVE)){avail =  c(avail,listVE[[i]])}
  toGet = intersect(avail,originalVes)
  
  extra = setdiff(originalVes,avail)
  if(length(extra)!=1 || extra!="spht"){
    print('Error: Unknown environmental variables required in originalVes')
    return(NULL)
  }
  coos = get_variables(toGet,coos,dpath = dataDir,fix=F)
  if("clc"%in%toGet & length(extra)==1 & extra=="spht"){coos$spht = get_spht(coos$clc)
  }else{print('Error: spht required whereas clc (original source of it) not in originalVes')}
  
  Mat = sparse.model.matrix(Intensityformula,data = coos)
  
  Matos = matrix(NA,dim(coos)[1],dim(Mat)[2])
  Matos = data.frame(Matos)
  Matos[complete.cases(coos),] = Mat
  Mat = Matos 
  Mat = as.matrix(Mat)
  rm(Matos)
  gc(reset=T)
  predMatrix = matrix(NA,dim(Mat)[1],length(taxasToPredict))
  colnames(predMatrix) = as.character(taxasToPredict)
  
  for(i in 1:length(taxasToPredict)){
    Sp = taxasToPredict[i]
    coefSp = get.spIntensityCoefficients(Sp,coefficients,spList=glc19SpIdVec)
    pred = as.numeric(Mat %*% coefSp)
    predMatrix[,i] = pred * nOccPerSp[i] / sum(pred,na.rm=T)
    if(i/10 == round(i/10)){
      flush.console()
      cat('    \r     Process...',100*i/length(taxasToPredict),'%        \r')
    }
  }      
  
  # /!\ predMatrix might have NA cells /!\ 
  return(predMatrix)
}

#####
# MAXNET
#####

predict.maxnet.edit = function(mod,data,type="maxModel"){
  if(type=="maxModel"){
    terms <- sub("hinge\\((.*)\\):(.*):(.*)$", "hingeval(\\1,\\2,\\3)", 
                 names(mod$betas))
    terms <- sub("categorical\\((.*)\\):(.*)$", "categoricalval(\\1,'\\2')", 
                 terms)
    terms <- sub("thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)", 
                 terms)
  }else if(type=="lofModel"){
    terms = sub("spht(.*)","categoricalval(spht,'\\1')",names(mod$betas))
  }
  f <- formula(paste("~", paste(terms, collapse = " + "), "-1"))
  D = sparse.model.matrix(f,data)
  p = as.vector(D %*% mod$betas)
  return(exp(p - mean(p)))
}

#####
# Intensity map over France
#####

# OK
create.french.departements.dataFrame= function(saveRdsFileDir,dptDir='C:/Users/Christophe/pCloud local/0_These/data/zones_administratives_fr/regions administratives/'){
  library(rgdal)
  
  # Emplacement de l'archive décompressée, à remplacer par le votre
  #adresse_communes_geofla = paste("C:/Users/",user,"/pCloud local/0_These/data/zones_administratives_fr/communes_geofla/GEOFLA_2-2_COMMUNE_SHP_LAMB93_FXX_2016-06-28/GEOFLA/1_DONNEES_LIVRAISON_2016-06-00236/GEOFLA_2-2_SHP_LAMB93_FR-ED161/COMMUNE",sep="")
  #adresse_cantons_geofla = paste("C:/Users/",user,"/pCloud local/0_These/data/zones_administratives_fr/communes_geofla/GEOFLA_2-0_CANTON_SHP_LAMB93_FXX_2015-07-01/GEOFLA/1_DONNEES_LIVRAISON_2015/GEOFLA_2-0_SHP_LAMB93_FR-ED151/CANTON",sep="")
  
  # Description des données via orgInfo
  # Attention à ne pas mettre l'extension à la fin du nom
  setwd(dptDir)
  ogrInfo(dsn = paste(dptDir,"departements-20180101.shp",sep=""),layer="departements-20180101")
  reg <- readOGR(dsn = paste(dptDir,"departements-20180101.shp",sep=""),layer="departements-20180101", stringsAsFactors=FALSE)
  toRm = c('FR940','FR920','FR910','FR930','<NA>')
  reg = reg[!reg@data$nuts3%in%toRm,]
  reg = reg[reg@data$nom!="Mayotte",]
  regProj = spTransform(reg,CRSobj = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  
  data = regProj@data
  data$id = 1:dim(data)[1]
  df = data[,c('id','code_insee')]
  df$lon = NA
  df$lat = NA
  df = df[-1,]
  for(i in 1:dim(regProj)[1]){
    tmp = data[i,c('id','code_insee')]
    coos = as.data.frame(regProj@polygons[[i]]@Polygons[[1]]@coords)
    coos$id = i
    tmp = merge(tmp,coos,by="id")
    colnames(tmp)[3:4]=c('lon','lat')
    df = rbind(df,tmp)
  }
  saveRDS(df,saveRdsFileDir)
}


map.spIntensity.over.France = function(pred,coos_l93,outFileLoc,dptFileLoc){
  library(ggplot2)
  
  df = readRDS(dptFileLoc)
  
  classes = cut(pred,c(min(pred)-.1,quantile(pred,c(.1,.4,.6,.9)),max(pred)+.1),dig.lab=2)
  colo = colorRampPalette(c("darkorchid4","goldenrod"))(length(unique(classes)))
  lev= levels(classes)
  
  labels = paste(lev,paste('/ quantile',c('0 to .1','.1 to .4','.4 to .6','.6 to .9','.9 to 1')))
  p = ggplot()+geom_polygon(data=df,aes(x=lon,y=lat,group=id),fill = NA,color = "grey20",size=.35)
  p = p + geom_tile(data=coos_l93,aes(x=x_l93,y=y_l93,fill=classes),size=1,alpha=.8)+scale_fill_manual(values=colo,name="log10 of species intensity",labels=labels)+theme_bw()
  p = p + xlab("Longitude in Lambert 93") + ylab("Latitude in Lambert 93")
  p = p + theme(legend.title = element_text(size=18),
                legend.text = element_text(size=15),
                axis.title.x = element_text(size=18),
                axis.title.y = element_text(size=18),
                axis.text.x = element_text(size=13),
                axis.text.y = element_text(size=13))
  png(outFileLoc,width=1200,height=750)
  print(p)
  dev.off()
}
