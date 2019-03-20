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
bind.background.pts = function(occ,bg,Bias_killer=100){
  occ = occ[complete.cases(occ),]
  occ$taxa = factor(occ$taxa)
  taxas = levels(occ$taxa)
  bg$taxa = factor(NA,levels=taxas)
  n_0=dim(bg)[1]
  # Attribute likelihood weights
  bg$w = (Bias_killer-1) / (n_0*Bias_killer)
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
  for(e in 1:length(taxas)){bg$taxa=taxas[e];occ=rbind(occ,bg)}
  return(occ)
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
  SparseDes = sparse.model.matrix(matrixFormu,data)
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
  orderedIdx = sapply(1:nQi,function(k) which(grepl(paste('q_hash',k,sep=""),names(coefo_qc))) )
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

