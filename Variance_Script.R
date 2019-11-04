# computation of LOF Variance matrix 
require(ggplot2)
require(rgeos)
require(sp)
require(raster)

#####
# Functions
#####

integrateOver = function(f,D,tolerance=.01,verbose=F){
  # Monte Carlo integral approximation
  AD = gArea(D)  
  
  ElemBand = 100
  rep = 5
  coefOfVariation = 1
  k=0
  zz = NULL
  while(coefOfVariation>tolerance){
    k=k+1
    if(length(zz)==0){
      zz= spsample(D, rep*ElemBand,type = "random")@coords
    }else{
      zz = rbind(zz,spsample(D,rep*ElemBand,type = "random")@coords)
    }
    #print(dim(zz))
    Integrals = sapply(1:rep,function(i){mean(f(zz[(1+(i-1)*k*ElemBand):(i*k*ElemBand),,drop=F]))*AD})
    #print(mean(Integrals))
    M = mean(Integrals)
    if(M!=0){coefOfVariation=sd(Integrals)/abs(M)}else{coefOfVariation=tolerance-1}
  }
  if(verbose){
    print(paste('Number of points required to achieve required precision:',dim(zz)[1]))
    print(paste("Integral approximation:",mean(Integrals)))
    print(paste("Sample sizes of last serie:",k*ElemBand))
    print(paste("Observed st. dev. across last serie:",sd(Integrals)))
  }
  return(mean(Integrals))
}

integrateOverMultivariate = function(f,D,tolerance=.01,verbose=F){
  # when applied to serie of points
  # f must return a matrix where line i 
  # are the values of f on point i 
  # Monte Carlo integral approximation
  AD = gArea(D)  
  
  ElemBand = 100
  rep = 5
  coefOfVariation = 1
  k=0
  zz = NULL
  while(coefOfVariation>tolerance){
    k=k+1
    if(length(zz)==0){
      zz= spsample(D, rep*ElemBand,type = "random")@coords
    }else{
      zz = rbind(zz,spsample(D,rep*ElemBand,type = "random")@coords)
    }
    #print(dim(zz))
    Integrals = sapply(1:rep,function(i){
      rowSums( t(f(zz[(1+(i-1)*k*ElemBand):(i*k*ElemBand),,drop=F])) )*AD /(k*ElemBand) })
    #print(mean(Integrals))
    M = rowSums(Integrals)/dim(Integrals)[2]
    SD = sapply(1:dim(Integrals)[1],function(i){sd(Integrals[i,])})
    coefOfVariations = rep(NA,length(M))
    coefOfVariations[M==0] = 0
    coefOfVariations[M!=0] = SD[M!=0]/M[M!=0]
    coefOfVariation= max(coefOfVariations)
  }
  if(verbose){
    print(paste('Number of points required to achieve required precision:',dim(zz)[1]))
    print("Integral approximation:")
    print(M)
    print(paste("Sample sizes of last serie:",k*ElemBand))
    print("Observed st. dev. across last serie:")
    print(SD)
  }
  return(M)
}

Compute.AlphasFromN = function(n_is,mu_is,sig_is,tol=0.05,verb=F){
  N = length(mu_is)
  alphas = NULL
  for(i in 1:N){
    fTmp = function(zz){sampEff(zz) * gaussianNiche(x(zz),mu_is[i],sig_is[i],1)}
    alphas[i] = n_is[i] / integrateOver(fTmp, D,tol,verbose = verb)
  }
  return(alphas)
}

Compute.I_gamma_beta = function(CrossBetaFuns,cellPolygons,p,tol=0.005,verbose=F ){
  N = length(CrossBetaFuns)
  Q = length(cell_polygons)
  # We compute the cross information between cell j parameter and species i density parameter vector \beta_i
  # LateX:
  # I(\gamma_j,\beta_i) = \int_{c_j} x^i(z) s(z) \lambda^i(x^i(z)) dz
  I_gamma_beta = array(NA,dim=c(Q,N,p))
  for(j in 1:Q){
    for(i in 1:N){
      I_gamma_beta[j,i,] = integrateOverMultivariate( CrossBetaFuns[[i]] , cellPolygons[[j]] ,tolerance = globTol,verbose=verbose )
    }
  }
  return(I_gamma_beta)
}

Compute.I_gamma_alpha = function(SampEffLambdaFun,cellPolygons,tol=0.05,verb=F){
  Q = length(cellPolygons)
  N= length(SampEffLambdaFun)
  # We build the matrix of dimensions Q * N
  # where the component of row j and col i is
  # LateX:
  # I(\gamma_j,\alpha_i) = \int_{cell_j} s(z) \lambda^i(x(z))dz
  # It is the cross information of cell j parameter and species i intercept 
  I_gamma_alpha = matrix(NA,Q,N)
  for(j in 1:Q){
    for(i in 1:N){
      I_gamma_alpha[j,i] = integrateOver( SampEffLambdaFun[[i]], cellPolygons[[j]], tolerance = tol, verbose=verb )
    }
  }
  return(I_gamma_alpha)
}

Compute.I_betas = function(D,sampEffFun,LambdaFun,x_featuresFun,p,tol=0.05,verb=F){
  # We compute the information matrix of density parameters for each species
  # LateX:
  # I(\beta_i) = \int_D x(z)x^T(z) s(z) \lambda^i(x(z)) dz
  I_betas = list()
  
  AD = gArea(D)  
  ElemBand = 100
  rep = 5
  for(i in 1:N){
    if(verb){print(paste('Compute I_betas[[',i,']]',sep=""))}
    coefOfVariation = 1
    k=0
    zz = NULL
    while(coefOfVariation>tol){
      k=k+1
      if(length(zz)==0){
        zz= spsample(D, rep*ElemBand,type = "random")@coords
      }else{
        zz = rbind(zz,spsample(D,rep*ElemBand,type = "random")@coords)
      }
      
      Matrices = list()
      for(l in 1:rep){
        cd = (1+(l-1)*k*ElemBand):(l*k*ElemBand)
        Matrices[[l]] = AD * t(x_featuresFun(zz[cd,]) * sampEffFun(zz[cd,]) * as.numeric(LambdaFun[[i]](x(zz[cd,])))) %*% x_featuresFun(zz[cd,]) / (k*ElemBand)
      }
      
      Mv = NULL
      SDv = NULL
      for(m in 1:(p^2)){
        Mv[m] = mean( sapply( 1:rep, function(l) as.vector(Matrices[[l]])[m] ) )
        SDv[m] = sd( sapply( 1:rep, function(l) as.vector(Matrices[[l]])[m] ) )
      }
      
      coefOfVariations = rep(NA,length(Mv))
      coefOfVariations[Mv==0] = 0
      coefOfVariations[Mv!=0] = SDv[Mv!=0]/Mv[Mv!=0]
      coefOfVariation= max(coefOfVariations)
    }
    
    I_betas[[i]] = matrix( Mv , p , p)
  }
  return(I_betas)
}

build.Information.Matrix = function(I_gammas,I_alphas,I_gamma_alpha,I_betas,I_beta_alpha,I_gamma_beta){
  Q= length(I_gammas)
  N = length(I_alphas)
  p = dim(I_gamma_beta)[3]
  
  names = c(paste("gamma_",1:Q,sep=""),  sapply(1:N,function(i) paste(c('alpha_',rep('beta_',p)),i,sep="") ) )  
  
  I = matrix(0,length(names),length(names))
  colnames(I) = names
  rownames(I) = names
  
  diag(I[paste("gamma_",1:Q,sep=""),paste("gamma_",1:Q,sep="")]) <- I_gammas
  
  for(i in 1:N){
    # I(\beta_i)
    cd1 = names==paste("beta_",i,sep="")
    I[cd1,cd1] = I_betas[[i]]
    
    # I(\alpha_i)
    cd2 = names==paste("alpha_",i,sep="")
    I[cd2,cd2] = I_alphas[i]
    
    # I(\beta_i,\alpha_i)
    I[cd1,cd2] = I_beta_alpha[,i]
    I[cd2,cd1] = t(I_beta_alpha[,i])
    for(j in 1:Q){
      cd3 = names==paste('gamma_',j,sep="")
      I[cd1,cd3] = I_gamma_beta[j,i,]
      I[cd3,cd1] = t(I_gamma_beta[j,i,])
      
      I[cd2,cd3] = I_gamma_alpha[j,i]
      I[cd3,cd2] = I_gamma_alpha[j,i]
    }
  }
  return(I)
}


#####
# Define intensity components
#####

# Define the spatial domain
D = Polygon(matrix(c(0,10,10,0,0,0,0,10,10,0),5,2))
D = SpatialPolygons( list(Polygons( list(D), "D")) )
#plot(D)

# Environmental variable 
x = function(z){if(length(dim(z))==2){z[,1,drop=F]-5}else{z[1]-5}}
#z = D@polygons[[1]]@Polygons[[1]]@coords
#x(z)

# environmental features (same for all species here)
x_features = function(z){if(length(dim(z))==2){cbind( z[,1,drop=F]-5 , (z[,1,drop=F]-5)^2 )}else{c(z[1]-5, (z[1]-5)^2 )}}
p = length(x_features(spsample(D,1,type="random")@coords))

## Define the sampling effort function 
# Constant
cutNice = function(x){
  1+ 3 * as.numeric(x>=(-4) & x<(-3)) +  7*as.numeric(x>=0 & x<1) +  2*as.numeric(x>=3 & x<4)
}
sampEff = function(z){cutNice(x(z))}


## We define sampling cells with a raster 
Q=10
cells = raster(nrows=1,ncols=Q, ext= extent(D@bbox) )
cells[] = sampEff( rasterToPoints(cells)[,c(1,2)] )
plot(cells) # Plot sampling effort

# Define species 
mus = c(-2.5,2.5)
sigs = c(1.6,1.6)
gaussianNiche = function(x,mu,sig,alpha){
  return(alpha * exp(- (x-mu)^2/(2*sig^2)))
}

globTol = .003
saveDir="C:/Users/Christophe/pCloud local/0_These/data/LOF data/19_11_04 variance test/"

######
# Effect of the number of occurrences
######

nFocs= c(100,250,500,750,1000,1250,1500)
nStds = nFocs
compteur=1
for(nFoc in nFocs){
  for(nStd in nStds){
    n_is = c(nFoc,nStd)
    print(paste('Step ',compteur ,' over ',length(nFocs)*length(nStds)))
    # Species densities 
    N = length(mus)
    alphas = NULL
    for(i in 1:N){
      fTmp = function(zz){sampEff(zz) * gaussianNiche(x(zz),mus[i],sigs[i],1)}
      alphas[i] = n_is[i] / integrateOver(fTmp, D,globTol,verbose = F)
    }
    
    ## We build species intensity functions as elements of a list of size N 
    Lambda = list()
    for(i in 1:N){
      Lambda[[i]]= eval(parse(text=paste( "function(x){gaussianNiche(x,mu=mus[",i,"],sig=sigs[",i,"],alpha=alphas[",i,"])}" ,sep="")))
    }
    
    ## We build species occurrence intensity functions as elements of a list of size N 
    # written for species i
    # z -> s(z) \lambda_i(x(z))
    SampEffLambda = list()
    for(i in 1:N){
      SampEffLambda[[i]]= eval(parse(text=paste( "function(z){sampEff(z)*Lambda[[",i,"]](x(z))}" ,sep="")))
    }
    
    ## We build the multivariate function
    # z -> x(z) s(z) \lambda_i(x(z))
    toCrossBeta = list()
    for(i in 1:N){
      toCrossBeta[[i]]= eval(parse(text=paste( "function(z){x_features(z)*sampEff(z)*as.numeric(Lambda[[",i,"]](x(z)))}" ,sep="")))
    }
    
    #####
    # Compute Fisher Information Matrix elements 
    #####
    
    # make cells polygons 
    cell_polygons = list()
    for(j in 1:Q){ 
      cells[]=F
      cells[j]=T
      # Convert to sp object
      spdf <- as(cells, "SpatialPolygonsDataFrame")
      cell_polygons[[j]] <- spdf[rowSums(spdf@data) == 1, ]
    }
    
    # We build the matrix of dimensions Q * N
    # where the component of row j and col i is
    # LateX:
    # I(\gamma_j,\alpha_i) = \int_{cell_j} s(z) \lambda^i(x(z))dz
    # It is the cross information of cell j parameter and species i intercept 
    I_gamma_alpha = matrix(NA,Q,N)
    for(j in 1:Q){
      for(i in 1:N){
        I_gamma_alpha[j,i] = integrateOver( SampEffLambda[[i]] , cell_polygons[[j]] ,tolerance = globTol,verbose=F )
      }
    }
    
    ## We compute the information associated with each sampling cell parameter
    # LateX:
    # I(\gamma_j) = \sum_{i=1}^N \int_{cell_j} s(z) \lambda^i(x(z)) dz 
    I_gammas = rowSums(I_gamma_alpha)
    ## We compute the information associated with species intercept 
    # = expected number of species occurrences
    # LateX:
    # I(\alpha_i) = \int_D s(z) \lambda^i(x(z)) dz  
    #I_alphas = rowSums(t(I_gamma_alpha))
    I_alphas = n_is
    
    # We compute the cross information between cell j parameter and species i density parameter vector \beta_i
    # LateX:
    # I(\gamma_j,\beta_i) = \int_{c_j} x^i(z) s(z) \lambda^i(x^i(z)) dz
    I_gamma_beta = array(NA,dim=c(Q,N,p))
    
    for(j in 1:Q){
      for(i in 1:N){
        I_gamma_beta[j,i,] = integrateOverMultivariate( toCrossBeta[[i]] , cell_polygons[[j]] ,tolerance = globTol,verbose=F )
      }
    }
    
    # We now just have to sum the previous ones over cells to get the 
    # cross information of species i intercept \alpha_i and species i density parameter vector \beta_i 
    # LateX:
    # I(\beta_i,\alpha_i) = \int_D x^i(z) s(z) \lambda^i(x^i(z)) dz
    I_beta_alpha = matrix(NA,p,N)
    for(i in 1:N){
      I_beta_alpha[,i] = rowSums(t(I_gamma_beta[,i,]))
    }
    
    # We compute the information matrix of density parameters for each species
    # LateX:
    # I(\beta_i) = \int_D x(z)x^T(z) s(z) \lambda^i(x(z)) dz
    I_betas = list()
    
    AD = gArea(D)  
    ElemBand = 100
    rep = 5
    for(i in 1:N){
      coefOfVariation = 1
      k=0
      zz = NULL
      while(coefOfVariation>globTol){
        k=k+1
        if(length(zz)==0){
          zz= spsample(D, rep*ElemBand,type = "random")@coords
        }else{
          zz = rbind(zz,spsample(D,rep*ElemBand,type = "random")@coords)
        }
        
        Matrices = list()
        for(l in 1:rep){
          cd = (1+(l-1)*k*ElemBand):(l*k*ElemBand)
          Matrices[[l]] = AD * t(x_features(zz[cd,]) * sampEff(zz[cd,]) * as.numeric(Lambda[[i]](x(zz[cd,])))) %*% x_features(zz[cd,]) / (k*ElemBand)
        }
        
        Mv = NULL
        SDv = NULL
        for(m in 1:(p^2)){
          Mv[m] = mean( sapply( 1:rep, function(l) as.vector(Matrices[[l]])[m] ) )
          SDv[m] = sd( sapply( 1:rep, function(l) as.vector(Matrices[[l]])[m] ) )
        }
        
        coefOfVariations = rep(NA,length(Mv))
        coefOfVariations[Mv==0] = 0
        coefOfVariations[Mv!=0] = SDv[Mv!=0]/Mv[Mv!=0]
        coefOfVariation= max(coefOfVariations)
      }
      
      I_betas[[i]] = matrix( Mv , p , p)
    }
    
    setwd(saveDir)
    saveRDS(I_gammas,paste("Q_",Q,"_nFoc_",nFoc,"_nStd_",nStd,"_I_gammas",sep=""))
    saveRDS(I_alphas,paste("Q_",Q,"_nFoc_",nFoc,"_nStd_",nStd,"_I_alphas",sep=""))  
    saveRDS(I_gamma_alpha,paste("Q_",Q,"_nFoc_",nFoc,"_nStd_",nStd,"_I_gamma_alpha",sep=""))  
    saveRDS(I_gamma_beta,paste("Q_",Q,"_nFoc_",nFoc,"_nStd_",nStd,"_I_gamma_beta",sep=""))  
    saveRDS(I_betas,paste("Q_",Q,"_nFoc_",nFoc,"_nStd_",nStd,"_I_betas",sep=""))  
    compteur = compteur+1
  }
}

#####
# Plot variance vs number of occurrences
#####

nOccu = nFocs
Q = 10
paramsToPlot = c('beta_1_1','beta_1_2','beta_2_1','beta_2_2')
toplot = data.frame(expand.grid(nOcc = nOccu,
                                curve=c("nFoc=100 & nStd grows","nStd=100 & nFoc grows"),
                                param=paramsToPlot,
                                variance=NA))
setwd(saveDir)
for(curve in c("nFoc=100 & nStd grows","nStd=100 & nFoc grows")){
  
  for(n in nOccu){
    if(regexpr('nFoc=',curve)>0){nFoc=100;nStd=n}else{nFoc=n;nStd=100}
    I_gammas=readRDS(paste("Q_",Q,"_nFoc_",nFoc,"_nStd_",nStd,"_I_gammas",sep=""))
    I_alphas=readRDS(paste("Q_",Q,"_nFoc_",nFoc,"_nStd_",nStd,"_I_alphas",sep=""))  
    I_gamma_alpha=readRDS(paste("Q_",Q,"_nFoc_",nFoc,"_nStd_",nStd,"_I_gamma_alpha",sep=""))  
    I_gamma_beta=readRDS(paste("Q_",Q,"_nFoc_",nFoc,"_nStd_",nStd,"_I_gamma_beta",sep=""))  
    I_betas=readRDS(paste("Q_",Q,"_nFoc_",nFoc,"_nStd_",nStd,"_I_betas",sep=""))  
    
    I = build.Information.Matrix(I_gammas,I_alphas,I_gamma_alpha,I_betas,I_beta_alpha,I_gamma_beta)
    
    Eigen = svd(I[-1,-1])
    conditionNumber = Eigen$d[1]/Eigen$d[length(Eigen$d)]
    print(paste("Condition number:",conditionNumber))
    if(conditionNumber<1e6){
      print('Ok')
      Sigma = solve(I[-1,-1])
    }
    
    Names = colnames(Sigma)
    cd = toplot$nOcc==n & toplot$curve==curve
    toplot$variance[cd & toplot$param=='beta_1_1'] = Sigma[Names=='beta_1',Names=='beta_1'][1,1]
    toplot$variance[cd & toplot$param=='beta_1_2'] = Sigma[Names=='beta_1',Names=='beta_1'][2,2]
    toplot$variance[cd & toplot$param=='beta_2_1'] = Sigma[Names=='beta_2',Names=='beta_2'][1,1]
    toplot$variance[cd & toplot$param=='beta_2_2'] = Sigma[Names=='beta_2',Names=='beta_2'][2,2]
  }
}

write.table(toplot,"toplot_Q_10_Bisp_cutNice.csv",sep=";",row.names=F,col.names=T)
tmp = toplot
toplot = toplot[toplot$nOcc!=1500,]
for(param in paramsToPlot){
  cd = toplot$param==param
  pl=ggplot()+geom_line(data=toplot[cd,],aes(x=nOcc,y=variance,group=curve,colour=curve),size=2)+theme_bw()+scale_y_continuous(limits=c(0,max(toplot$variance[cd])))+ggtitle(paste(param," variance vs number of occurrences"))+xlab("Number of occurrences")+ylab('Estimation variance')
  png(paste("Q_10_",param,'_Bisp_cutNice.png',sep=""),width=600,height=400)
  print(pl)
  dev.off()
}

#####
# Effect of the number of cells
##### 

globTol=0.003
n_is = c(1000,1000)
Qs=c(4,10,20,30,50,100)
compteur=1
for(Q in Qs){
  ## We define sampling cells with a raster
  cells = raster(nrows=1,ncols=Q, ext= extent(D@bbox) )

  # make cells polygons 
  cell_polygons = list()
  for(j in 1:Q){ 
    cells[]=F
    cells[j]=T
    # Convert to sp object
    spdf <- as(cells, "SpatialPolygonsDataFrame")
    cell_polygons[[j]] <- spdf[rowSums(spdf@data) == 1, ]
  }
  
  print(paste('Step ',compteur ,' over ',length(Qs)))
  ### Species densities 
  # Compute species intensity multiplying constants 
  alphas = Compute.AlphasFromN(n_is,mus,sigs,tol = globTol)
  
  ## We build species intensity functions as elements of a list of size N 
  Lambda = list()
  for(i in 1:N){
    Lambda[[i]]= eval(parse(text=paste( "function(x){gaussianNiche(x,mu=mus[",i,"],sig=sigs[",i,"],alpha=alphas[",i,"])}" ,sep="")))
  }
  
  ## We build species occurrence intensity functions as elements of a list of size N 
  # written for species i
  # z -> s(z) \lambda_i(x(z))
  SampEffLambda = list()
  for(i in 1:N){
    SampEffLambda[[i]]= eval(parse(text=paste( "function(z){sampEff(z)*Lambda[[",i,"]](x(z))}" ,sep="")))
  }
  
  ## We build the multivariate function
  # z -> x(z) s(z) \lambda_i(x(z))
  toCrossBeta = list()
  for(i in 1:N){
    toCrossBeta[[i]]= eval(parse(text=paste( "function(z){x_features(z)*sampEff(z)*as.numeric(Lambda[[",i,"]](x(z)))}" ,sep="")))
  }
  
  # We build the matrix of dimensions Q * N
  # where the component of row j and col i is
  # LateX:
  # I(\gamma_j,\alpha_i) = \int_{cell_j} s(z) \lambda^i(x(z))dz
  # It is the cross information of cell j parameter and species i intercept 
  I_gamma_alpha = Compute.I_gamma_alpha(SampEffLambda,cell_polygons,tol = globTol)
  
  ## We compute the information associated with each sampling cell parameter
  # LateX:
  # I(\gamma_j) = \sum_{i=1}^N \int_{cell_j} s(z) \lambda^i(x(z)) dz 
  I_gammas = rowSums(I_gamma_alpha)
  ## We compute the information associated with species intercept 
  # = expected number of species occurrences
  # LateX:
  # I(\alpha_i) = \int_D s(z) \lambda^i(x(z)) dz  
  #I_alphas = rowSums(t(I_gamma_alpha))
  I_alphas = n_is
  
  # We compute the cross information between cell j parameter and species i density parameter vector \beta_i
  # LateX:
  # I(\gamma_j,\beta_i) = \int_{c_j} x^i(z) s(z) \lambda^i(x^i(z)) dz
  I_gamma_beta = Compute.I_gamma_beta(CrossBetaFuns = toCrossBeta,cellPolygons = cell_polygons,p=p,tol = globTol,verbose = F)
  
  # We now just have to sum the previous ones over cells to get the 
  # cross information of species i intercept \alpha_i and species i density parameter vector \beta_i 
  # LateX:
  # I(\beta_i,\alpha_i) = \int_D x^i(z) s(z) \lambda^i(x^i(z)) dz
  I_beta_alpha = matrix(NA,p,N)
  for(i in 1:N){I_beta_alpha[,i] = rowSums(t(I_gamma_beta[,i,]))}
  
  # We compute the information matrix of density parameters for each species
  # LateX:
  # I(\beta_i) = \int_D x(z)x^T(z) s(z) \lambda^i(x(z)) dz
  I_betas = Compute.I_betas(D,sampEff,Lambda,x_features,p=p,tol=globTol,verb=T)
  
  setwd(saveDir)
  saveRDS(I_gammas,paste("_cellXP_Q_",Q,"_nFoc_",n_is[1],"_nStd_",n_is[2],"_I_gammas",sep=""))
  saveRDS(I_alphas,paste("_cellXP_Q_",Q,"_nFoc_",n_is[1],"_nStd_",n_is[2],"_I_alphas",sep=""))  
  saveRDS(I_gamma_alpha,paste("_cellXP_Q_",Q,"_nFoc_",n_is[1],"_nStd_",n_is[2],"_I_gamma_alpha",sep=""))  
  saveRDS(I_gamma_beta,paste("_cellXP_Q_",Q,"_nFoc_",n_is[1],"_nStd_",n_is[2],"_I_gamma_beta",sep=""))  
  saveRDS(I_betas,paste("_cellXP_Q_",Q,"_nFoc_",n_is[1],"_nStd_",n_is[2],"_I_betas",sep=""))  
  compteur = compteur+1
}

#####
# Compute density error metric
#####



#####
# Plot
#####
 

#####
# Bin
#####
# totalOccurrenceIntensity 
# is the function written as
# z -> s(z) \sum_{i=1}^N \lambda^i(x(z))
#totalOccurrenceIntensity = function(z){sampEff(z)*rowSums(sapply(1:N,function(i) Lambda[[i]](x(z)))) }

