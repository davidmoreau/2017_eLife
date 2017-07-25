#################################################################################################################
### High-intensity Training Enhances Executive Function in Children in a Randomized, Placebo-Controlled Trial ###
### eLife, 2017
### David Moreau
### Centre for Brain Research, University of Auckland, New Zealand
### d.moreau@auckland.ac.nz
#################################################################################################################


### Setup ###

# Clear cache
rm(list = ls())

# Load packages (+ install if required)
packages <- c("car", "BayesFactor", "LaplacesDemon", "ggplot2", "gridExtra", "dplyr", 
              "rjags", "lsr", "psych", "pwr") 
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(packages, suppressPackageStartupMessages(require), warn.conflicts=F, quietly=T, character.only=T)


### EFA ###

cor(na.omit(efa.data))
fit <- factanal(na.omit(efa.data), 2, 
                rotation = "promax", scores="regression") #the model
print(fit, digits=2, cutoff=.1, sort=TRUE) 
nrow(na.omit(efa.data)) # N for efa

# Note: The chi-square statistic and p-value in factanal are testing the hypothesis that the model fits the data perfectly. 
# When the p value is low we can reject this hypothesis - so in this case, the 2-factor model does not fit the data perfectly
# Uniquenesses are the variance in each item that is not explained by the two factors.
fit2 <- fa.promax(efa.data, factors=2, digits=2, sort=T)
load <- fit$loadings[,1:2] 
plot(load,type="l",main="Screeplot") # set up plot 
text(load,labels=names(mydata),cex=.7) # add variable names


### PCA ###
# PCA analysis confirms the two-factor model
# (overall executive function construct with CC and WM factors
pca.data <- hit %>%
  select(A, B, C, D, E, F)
pca <- prcomp(na.omit(pca.data))
print(pca)
plot(pca, type="l", main="Screeplot") #plot variances
sumr <- summary(pca); sumr
plot(sumr[3,], type="l") #plot cummulative proportions


### Additional ###

# Split HR_rest on hit.wide
hit.wide$HR_split_quant <- ifelse(hit.wide$HR_rest_1 >= quantile(hit.wide$HR_rest_1, na.rm=T)[3], "high", "low") 
high.HR <- hit.wide %>%
  filter(HR_split_quant == "high")

high.HR$Condition <- as.factor(high.HR$Condition)
t.test(high.HR$HR_gain ~ high.HR$Condition) #p-value = 0.01128
table(high.HR$Condition)
exp(ttest.tstat(t=2.5663, n1=75, n2=66, rscale = 0.707)[['bf']]) #BF=3.531438 in favor or alt. 
exp(ttest.tstat(t=2.5663, n1=72, n2=81, rscale = 0.707)[['bf']]) #BF=3.471319 in favor or alt. 
describeBy(high.HR$HR_gain, group = high.HR$Condition)


# Check if HR_gain predicts cognitive gains in high.HR sample
model_ATT <- lm(-high.HR$ATT_gain ~ high.HR$HR_gain)
summary(model_ATT)
plot(-high.HR$ATT_gain ~ high.HR$HR_gain)
abline(model_ATT, col = "blue")

model_WM <- lm(-high.HR$WM_gain ~ high.HR$HR_gain)
summary(model_WM)
plot(-high.HR$WM_gain ~ high.HR$HR_gain)
abline(model_WM, col = "blue")

### Parameters ###
# Parts of script from:
# Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Tutorial with R, JAGS, and Stan. Academic Press, Elsevier
# Check that required packages are installed:
want = c("parallel","rjags","runjags","compute.es")
have = want %in% rownames(installed.packages())
if ( any(!have) ) { install.packages( want[!have] ) }

# Load rjags
try( library(rjags) )
# Load runjags
try( library(runjags) )
try( runjags.options( inits.warning=FALSE , rng.warning=FALSE ) )

# set default number of chains and parallelness for MCMC:
library(parallel) # for detectCores().
nCores = detectCores() 
if ( !is.finite(nCores) ) { nCores = 1 } 
if ( nCores > 4 ) { 
  nChainsDefault = 4  # because JAGS has only 4 rng's.
  runjagsMethodDefault = "parallel"
}
if ( nCores == 4 ) { 
  nChainsDefault = 3  # save 1 core for other processes.
  runjagsMethodDefault = "parallel"
}
if ( nCores < 4 ) { 
  nChainsDefault = 3 
  runjagsMethodDefault = "rjags" # NOT parallel
}


openGraph = function( width=7 , height=7 , mag=1.0 , ... ) {
  if ( .Platform$OS.type != "windows" ) {
    tryInfo = try( X11( width=width*mag , height=height*mag , type="cairo" , 
                        ... ) )
    if ( class(tryInfo)=="try-error" ) {
      lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      X11( width=width*mag , height=height*mag , type="cairo" , ... )
    }
  } else { # Windows OS
    tryInfo = try( windows( width=width*mag , height=height*mag , ... ) )
    if ( class(tryInfo)=="try-error" ) {
      lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      windows( width=width*mag , height=height*mag , ... )    
    }
  }
}

saveGraph = function( file="saveGraphOutput" , type="pdf" , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    if ( any( type == c("png","jpeg","jpg","tiff","bmp")) ) {
      sptype = type
      if ( type == "jpg" ) { sptype = "jpeg" }
      savePlot( file=paste0(file,".",type) , type=sptype , ... )     
    }
    if ( type == "pdf" ) {
      dev.copy2pdf(file=paste0(file,".",type) , ... )
    }
    if ( type == "eps" ) {
      dev.copy2eps(file=paste0(file,".",type) , ... )
    }
  } else { # Windows OS
    file=paste0(file,".",type) 
    savePlot( file=file , type=type , ... )
  }
}


HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

HDIofICDF = function( ICDFname , credMass=0.95 , tol=1e-8 , ... ) {
  incredMass =  1.0 - credMass
  intervalWidth = function( lowTailPr , ICDFname , credMass , ... ) {
    ICDFname( credMass + lowTailPr , ... ) - ICDFname( lowTailPr , ... )
  }
  optInfo = optimize( intervalWidth , c( 0 , incredMass ) , ICDFname=ICDFname ,
                      credMass=credMass , tol=tol , ... )
  HDIlowTailPr = optInfo$minimum
  return( c( ICDFname( HDIlowTailPr , ... ) ,
             ICDFname( credMass + HDIlowTailPr , ... ) ) )
}

HDIofGrid = function( probMassVec , credMass=0.95 ) {
  sortedProbMass = sort( probMassVec , decreasing=TRUE )
  HDIheightIdx = min( which( cumsum( sortedProbMass ) >= credMass ) )
  HDIheight = sortedProbMass[ HDIheightIdx ]
  HDImass = sum( probMassVec[ probMassVec >= HDIheight ] )
  return( list( indices = which( probMassVec >= HDIheight ) ,
                mass = HDImass , height = HDIheight ) )
}

# Function(s) for plotting properties of mcmc coda objects.

DbdaAcfPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat = NULL
  yMat = NULL
  for ( cIdx in 1:nChain ) {
    acfInfo = acf(codaObject[,c(parName)][[cIdx]],plot=FALSE) 
    xMat = cbind(xMat,acfInfo$lag)
    yMat = cbind(yMat,acfInfo$acf)
  }
  matplot( xMat , yMat , type="o" , pch=20 , col=plColors , ylim=c(0,1) ,
           main="" , xlab="Lag" , ylab="Autocorrelation" )
  abline(h=0,lty="dashed")
  EffChnLngth = effectiveSize(codaObject[,c(parName)])
  text( x=max(xMat) , y=max(yMat) , adj=c(1.0,1.0) , cex=1.25 ,
        labels=paste("ESS =",round(EffChnLngth,1)) )
}

DbdaDensPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject) # or nchain(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat = NULL
  yMat = NULL
  hdiLims = NULL
  for ( cIdx in 1:nChain ) {
    densInfo = density(codaObject[,c(parName)][[cIdx]]) 
    xMat = cbind(xMat,densInfo$x)
    yMat = cbind(yMat,densInfo$y)
    hdiLims = cbind(hdiLims,HDIofMCMC(codaObject[,c(parName)][[cIdx]]))
  }
  matplot( xMat , yMat , type="l" , col=plColors , 
           main="" , xlab="Param. Value" , ylab="Density" )
  abline(h=0)
  points( hdiLims[1,] , rep(0,nChain) , col=plColors , pch="|" )
  points( hdiLims[2,] , rep(0,nChain) , col=plColors , pch="|" )
  text( mean(hdiLims) , 0 , "" , adj=c(0.5,-0.2) )
  EffChnLngth = effectiveSize(codaObject[,c(parName)])
  MCSE = sd(as.matrix(codaObject[,c(parName)]))/sqrt(EffChnLngth) 
  text( max(xMat) , max(yMat) , adj=c(1.0,1.0) , cex=1.25 ,
        paste("MCSE =\n",signif(MCSE,3)) )
}

diagMCMC = function( codaObject , parName=varnames(codaObject)[1] ,
                     saveName=NULL , saveType="jpg" ) {
  DBDAplColors = c("darkolivegreen3", "black", "#FF8F00", "#4F94CD")
  openGraph(height=5,width=7)
  par( mar=0.5+c(3,4,1,0) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
       cex.lab=1.5 )
  layout(matrix(1:4,nrow=2))
  # traceplot and gelman.plot are from CODA package:
  require(coda)
  coda::traceplot( codaObject[,c(parName)] , main="" , ylab="Param. Value" ,
                   col=DBDAplColors ) 
  tryVal = try(
    coda::gelman.plot( codaObject[,c(parName)] , main="" , auto.layout=FALSE , 
                       col=DBDAplColors )
  )  
  # if it runs, gelman.plot returns a list with finite shrink values:
  if ( class(tryVal)=="try-error" ) {
    plot.new() 
    print(paste0("Warning",parName))
  } else { 
    if ( class(tryVal)=="list" & !is.finite(tryVal$shrink[1]) ) {
      plot.new() 
      print(paste0("Warning",parName))
    }
  }
  DbdaAcfPlot(codaObject,parName,plColors=DBDAplColors)
  DbdaDensPlot(codaObject,parName,plColors=DBDAplColors)
  mtext( text=parName , outer=TRUE , adj=c(0.5,0.5) , cex=2.0 )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,"Diag",parName), type=saveType)
  }
}

diagStanFit = function( stanFit , parName ,
                        saveName=NULL , saveType="jpg" ) {
  codaFit = mcmc.list( lapply( 1:ncol(stanFit) , 
                               function(x) { mcmc(as.array(stanFit)[,x,]) } ) )
  DBDAplColors = c("darkolivegreen3", "black", "#FF8F00", "#4F94CD")
  openGraph(height=5,width=7)
  par( mar=0.5+c(3,4,1,0) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , cex.lab=1.5 )
  layout(matrix(1:4,nrow=2))
  # traceplot is from rstan package
  require(rstan)
  traceplot(stanFit,pars=parName,nrow=1,ncol=1)#,main="",ylab="Param. Value",col=DBDAplColors) 
  # gelman.plot are from CODA package:
  require(coda)
  tryVal = try(
    coda::gelman.plot( codaObject[,c(parName)] , main="" , auto.layout=FALSE , 
                       col=DBDAplColors )
  )
  # if it runs, gelman.plot returns a list with finite shrink values:
  if ( class(tryVal)=="try-error" ) {
    plot.new() 
    print(paste0("Warning",parName))
  } else { 
    if ( class(tryVal)=="list" & !is.finite(tryVal$shrink[1]) ) {
      plot.new() 
      print(paste0("Warning",parName))
    }
  }
  DbdaAcfPlot(codaFit,parName,plColors=DBDAplColors)
  DbdaDensPlot(codaFit,parName,plColors=DBDAplColors)
  mtext( text=parName , outer=TRUE , adj=c(0.5,0.5) , cex=2.0 )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,"Diag",parName), type=saveType)
  }
}

# Functions for summarizing and plotting distribution of a large sample; 
# typically applied to MCMC posterior.

normalize = function( v ){ return( v / sum(v) ) }

require(coda) # loaded by rjags, but redundancy doesn't hurt

summarizePost = function( paramSampleVec , 
                          compVal=NULL , ROPE=NULL , credMass=0.95 ) {
  meanParam = mean( paramSampleVec )
  medianParam = median( paramSampleVec )
  dres = density( paramSampleVec )
  modeParam = dres$x[which.max(dres$y)]
  mcmcEffSz = round( effectiveSize( paramSampleVec ) , 1 )
  names(mcmcEffSz) = NULL
  hdiLim = HDIofMCMC( paramSampleVec , credMass=credMass )
  if ( !is.null(compVal) ) {
    pcgtCompVal = ( 100 * sum( paramSampleVec > compVal ) 
                    / length( paramSampleVec ) )
  } else {
    compVal=NA
    pcgtCompVal=NA
  }
  if ( !is.null(ROPE) ) {
    pcltRope = ( 100 * sum( paramSampleVec < ROPE[1] ) 
                 / length( paramSampleVec ) )
    pcgtRope = ( 100 * sum( paramSampleVec > ROPE[2] ) 
                 / length( paramSampleVec ) )
    pcinRope = 100-(pcltRope+pcgtRope)
  } else { 
    ROPE = c(NA,NA)
    pcltRope=NA 
    pcgtRope=NA 
    pcinRope=NA 
  }  
  return( c( Mean=meanParam , Median=medianParam , Mode=modeParam , 
             ESS=mcmcEffSz ,
             HDImass=credMass , HDIlow=hdiLim[1] , HDIhigh=hdiLim[2] , 
             CompVal=compVal , PcntGtCompVal=pcgtCompVal , 
             ROPElow=ROPE[1] , ROPEhigh=ROPE[2] ,
             PcntLtROPE=pcltRope , PcntInROPE=pcinRope , PcntGtROPE=pcgtRope ) )
}

plotPost = function( paramSampleVec , cenTend=c("mode","median","mean")[1] , 
                     compVal=NULL, ROPE=NULL, credMass=0.95, HDItextPlace=0.7, 
                     xlab=NULL , xlim=NULL , yaxt=NULL , ylab=NULL , 
                     main=NULL , cex=NULL , cex.lab=NULL ,
                     col=NULL , border=NULL , showCurve=FALSE , breaks=NULL , 
                     ... ) {
  # Override defaults of hist function, if not specified by user:
  # (additional arguments "..." are passed to the hist function)
  if ( is.null(xlab) ) xlab="Param. Val."
  if ( is.null(cex.lab) ) cex.lab=1.5
  if ( is.null(cex) ) cex=1.4
  if ( is.null(xlim) ) xlim=range( c( compVal , ROPE , paramSampleVec ) )
  if ( is.null(main) ) main=""
  if ( is.null(yaxt) ) yaxt="n"
  if ( is.null(ylab) ) ylab=""
  if ( is.null(col) ) col="darkolivegreen3"
  if ( is.null(border) ) border="white"
  
  # convert coda object to matrix:
  if ( class(paramSampleVec) == "mcmc.list" ) {
    paramSampleVec = as.matrix(paramSampleVec)
  }
  
  summaryColNames = c("ESS","mean","median","mode",
                      "hdiMass","hdiLow","hdiHigh",
                      "compVal","pGtCompVal",
                      "ROPElow","ROPEhigh","pLtROPE","pInROPE","pGtROPE")
  postSummary = matrix( NA , nrow=1 , ncol=length(summaryColNames) , 
                        dimnames=list( c( xlab ) , summaryColNames ) )
  
  # require(coda) # for effectiveSize function
  postSummary[,"ESS"] = effectiveSize(paramSampleVec)
  
  postSummary[,"mean"] = mean(paramSampleVec)
  postSummary[,"median"] = median(paramSampleVec)
  mcmcDensity = density(paramSampleVec)
  postSummary[,"mode"] = mcmcDensity$x[which.max(mcmcDensity$y)]
  
  HDI = HDIofMCMC( paramSampleVec , credMass )
  postSummary[,"hdiMass"]=credMass
  postSummary[,"hdiLow"]=HDI[1]
  postSummary[,"hdiHigh"]=HDI[2]
  
  # Plot histogram.
  cvCol = "#FF8F00"
  ropeCol = "#4F94CD"
  if ( is.null(breaks) ) {
    if ( max(paramSampleVec) > min(paramSampleVec) ) {
      breaks = c( seq( from=min(paramSampleVec) , to=max(paramSampleVec) ,
                       by=(HDI[2]-HDI[1])/18 ) , max(paramSampleVec) )
    } else {
      breaks=c(min(paramSampleVec)-1.0E-6,max(paramSampleVec)+1.0E-6)
      border="darkolivegreen3"
    }
  }
  if ( !showCurve ) {
    par(xpd=NA)
    histinfo = hist( paramSampleVec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                     freq=F , border=border , col=col ,
                     xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                     breaks=breaks , ... )
  }
  if ( showCurve ) {
    par(xpd=NA)
    histinfo = hist( paramSampleVec , plot=F )
    densCurve = density( paramSampleVec , adjust=2 )
    plot( densCurve$x , densCurve$y , type="l" , lwd=5 , col=col , bty="n" ,
          xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
          main=main , cex=cex , cex.lab=cex.lab , ... )
  }
  cenTendHt = 0.9*max(histinfo$density)
  cvHt = 0.7*max(histinfo$density)
  ROPEtextHt = 0.55*max(histinfo$density)
  # Display central tendency:
  mn = mean(paramSampleVec)
  med = median(paramSampleVec)
  mcmcDensity = density(paramSampleVec)
  mo = mcmcDensity$x[which.max(mcmcDensity$y)]
  if ( cenTend=="mode" ){ 
    text( mo , cenTendHt ,
          bquote(mode==.(signif(mo,3))) , adj=c(.5,0) , cex=cex )
  }
  if ( cenTend=="median" ){ 
    text( med , cenTendHt ,
          bquote(median==.(signif(med,3))) , adj=c(.5,0) , cex=cex , col=cvCol )
  }
  if ( cenTend=="mean" ){ 
    text( mn , cenTendHt ,
          bquote(mean==.(signif(mn,3))) , adj=c(.5,0) , cex=cex )
  }
  # Display the comparison value.
  if ( !is.null( compVal ) ) {
    pGtCompVal = sum( paramSampleVec > compVal ) / length( paramSampleVec ) 
    pLtCompVal = 1 - pGtCompVal
    lines( c(compVal,compVal) , c(0.96*cvHt,0) , 
           lty="dotted" , lwd=2 , col=cvCol )
    text( compVal , cvHt ,
          #bquote( .(round(100*pLtCompVal,1)) * "% < " *
          #.(signif(compVal,3)) * " < " * 
          #.(round(100*pGtCompVal,1)) * "%" ) ,
          adj=c(pLtCompVal,0) , cex=0.8*cex , col=cvCol )
    postSummary[,"compVal"] = compVal
    postSummary[,"pGtCompVal"] = pGtCompVal
  }
  # Display the ROPE.
  if ( !is.null( ROPE ) ) {
    pInROPE = ( sum( paramSampleVec > ROPE[1] & paramSampleVec < ROPE[2] )
                / length( paramSampleVec ) )
    pGtROPE = ( sum( paramSampleVec >= ROPE[2] ) / length( paramSampleVec ) )
    pLtROPE = ( sum( paramSampleVec <= ROPE[1] ) / length( paramSampleVec ) )
    lines( c(ROPE[1],ROPE[1]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol )
    lines( c(ROPE[2],ROPE[2]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol)
    text( mean(ROPE) , ROPEtextHt ,
          #bquote( .(round(100*pLtROPE,1)) * "% < " * .(ROPE[1]) * " < " * 
          # .(round(100*pInROPE,1)) * "% < " * .(ROPE[2]) * " < " * 
          # .(round(100*pGtROPE,1)) * "%" ) ,
          adj=c(pLtROPE+.5*pInROPE,0) , cex=1 , col=ropeCol )
    
    postSummary[,"ROPElow"]=ROPE[1] 
    postSummary[,"ROPEhigh"]=ROPE[2] 
    postSummary[,"pLtROPE"]=pLtROPE
    postSummary[,"pInROPE"]=pInROPE
    postSummary[,"pGtROPE"]=pGtROPE
  }
  # Display the HDI.
  lines( HDI , c(0,0) , lwd=4 , lend=1 )
  text( mean(HDI) , 0 , #bquote(.(100*credMass) * "% HDI" ) ,
        adj=c(.5,-1.7) , cex=cex )
  text( HDI[1] , 0 , #bquote(.(signif(HDI[1],3))) ,
        adj=c(HDItextPlace,-0.5) , cex=cex )
  text( HDI[2] , 0 , #bquote(.(signif(HDI[2],3))) ,
        adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
  par(xpd=F)
  #
  return( postSummary )
}

# Shape parameters from central tendency and scale:

betaABfromMeanKappa = function( mean , kappa ) {
  if ( mean <=0 | mean >= 1) stop("must have 0 < mean < 1")
  if ( kappa <=0 ) stop("kappa must be > 0")
  a = mean * kappa
  b = ( 1.0 - mean ) * kappa
  return( list( a=a , b=b ) )
}

betaABfromModeKappa = function( mode , kappa ) {
  if ( mode <=0 | mode >= 1) stop("must have 0 < mode < 1")
  if ( kappa <=2 ) stop("kappa must be > 2 for mode parameterization")
  a = mode * ( kappa - 2 ) + 1
  b = ( 1.0 - mode ) * ( kappa - 2 ) + 1
  return( list( a=a , b=b ) )
}

betaABfromMeanSD = function( mean , sd ) {
  if ( mean <=0 | mean >= 1) stop("must have 0 < mean < 1")
  if ( sd <= 0 ) stop("sd must be > 0")
  kappa = mean*(1-mean)/sd^2 - 1
  if ( kappa <= 0 ) stop("invalid combination of mean and sd")
  a = mean * kappa
  b = ( 1.0 - mean ) * kappa
  return( list( a=a , b=b ) )
}

gammaShRaFromMeanSD = function( mean , sd ) {
  if ( mean <=0 ) stop("mean must be > 0")
  if ( sd <=0 ) stop("sd must be > 0")
  shape = mean^2/sd^2
  rate = mean/sd^2
  return( list( shape=shape , rate=rate ) )
}

gammaShRaFromModeSD = function( mode , sd ) {
  if ( mode <=0 ) stop("mode must be > 0")
  if ( sd <=0 ) stop("sd must be > 0")
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}


### Low-level MCMC ###

genMCMC = function( datFrm, yName="y" , xName="x" ,
                    saveName=NULL , numSavedSteps=50000 , thinSteps=1 ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  y = as.numeric(datFrm[,yName])
  x = as.numeric(as.factor(datFrm[,xName]))
  xLevels = levels(as.factor(datFrm[,xName]))

  if ( any( x!=1 & x!=2 ) ) { stop("All x values must be 1 or 2.") }
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( length(x) != length(y) ) { stop("x and y must be same length.") }
  Ntotal = length(y)

  dataList = list(
    y = y ,
    x = x ,
    Ntotal = Ntotal ,
    meanY = mean(y) ,
    sdY = sd(y)
  )

  # Model
  modelString = "
  model {
  for ( i in 1:Ntotal ) {
  y[i] ~ dt( mu[x[i]] , 1/sigma[x[i]]^2 , nu )
  }
  for ( j in 1:2 ) { # 2 groups
  mu[j] ~ dnorm( meanY , 1/(100*sdY)^2 )
  sigma[j] ~ dunif( sdY/1000 , sdY*1000 )
  }
  nu ~ dexp(1/30.0)
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  
  # Init. chains
  # Initial values of MCMC chains based on data:
  mu = c( mean(y[x==1]) , mean(y[x==2]) )
  sigma = c( sd(y[x==1]) , sd(y[x==2]) )
  initsList = list( mu = mu , sigma = sigma , nu = 5 )

  # Run chains
  parameters = c( "mu" , "sigma" , "nu" )     # The parameters to be monitored
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

smryMCMC = function(  codaSamples , RopeMuDiff=NULL , RopeSdDiff=NULL , 
                      RopeEff=NULL , saveName=NULL ) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  summaryInfo = rbind( summaryInfo , 
                       "mu[1]" = summarizePost( mcmcMat[,"mu[1]"] ) )
  summaryInfo = rbind( summaryInfo , 
                       "mu[2]" = summarizePost( mcmcMat[,"mu[2]"] ) )
  summaryInfo = rbind( summaryInfo , 
                       "muDiff" = summarizePost( 
                         mcmcMat[,"mu[2]"] - mcmcMat[,"mu[1]"] , 
                         compVal=0.0 , ROPE=RopeMuDiff ) )
  summaryInfo = rbind( summaryInfo , 
                       "sigma[1]" = summarizePost( mcmcMat[,"sigma[1]"] ) )
  summaryInfo = rbind( summaryInfo , 
                       "sigma[2]" = summarizePost( mcmcMat[,"sigma[2]"] ) )
  summaryInfo = rbind( summaryInfo , 
                       "sigmaDiff" = summarizePost( 
                         mcmcMat[,"sigma[2]"] - mcmcMat[,"sigma[1]"] , 
                         compVal=0.0 , ROPE=RopeSdDiff ) )
  summaryInfo = rbind( summaryInfo , 
                       "nu" = summarizePost( mcmcMat[,"nu"] ) )
  summaryInfo = rbind( summaryInfo , 
                       "log10(nu)" = summarizePost( log10(mcmcMat[,"nu"]) ) )
  summaryInfo = rbind( summaryInfo , 
                       "effSz" = summarizePost( 
                         ( mcmcMat[,"mu[2]"]-mcmcMat[,"mu[1]"] ) 
                         / sqrt((mcmcMat[,"sigma[1]"]^2+mcmcMat[,"sigma[2]"]^2)/2) ,
                         compVal=0.0 , ROPE=RopeEff ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

plotMCMC = function( codaSamples , datFrm , yName="y" , xName="x" ,
                     RopeMuDiff=NULL , RopeSdDiff=NULL , RopeEff=NULL , 
                     showCurve=TRUE, #Usually false 
                     pairsPlot=FALSE ,
                     saveName=NULL , saveType="jpg" ) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  mu1 = mcmcMat[,"mu[1]"]
  mu2 = mcmcMat[,"mu[2]"]
  sigma1 = mcmcMat[,"sigma[1]"]
  sigma2 = mcmcMat[,"sigma[2]"]
  nu = mcmcMat[,"nu"]

  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph(width=7,height=7)
    nPtToPlot = 1000
    plotIdx = floor(seq(1,length(mu1),by=length(mu1)/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( mu1 , mu2 , sigma1 , sigma2 , log10(nu) )[plotIdx,] ,
           lower.panel=panel.cor , col="darkolivegreen3" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }

  # Set up window and layout:
  openGraph(width=6.0,height=8.0)
  layout( matrix( c(4,5,7,8,3,1,2,6,9,10) , nrow=5, byrow=FALSE ) )
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  # Select thinned steps in chain for plotting of posterior predictive curves:
  nCurvesToPlot = 20
  stepIdxVec = seq( 1 , chainLength , floor(chainLength/nCurvesToPlot) )
  # Compute limits for plots of data with posterior pred. distributions
  y = as.numeric(datFrm[,yName])
  x = as.numeric(as.factor(datFrm[,xName]))
  xLevels = levels(as.factor(datFrm[,xName]))
  y1 = y[x==1]
  y2 = y[x==2]
  xLim = c( min(y)-0.1*(max(y)-min(y)) , max(y)+0.1*(max(y)-min(y)) )
  xBreaks = seq( xLim[1] , xLim[2] , 
                 length=ceiling((xLim[2]-xLim[1])/(mean(c(sd(y1),sd(y2)))/4)) )
  histInfo1 = hist(y1,breaks=xBreaks,plot=FALSE)
  histInfo2 = hist(y2,breaks=xBreaks,plot=FALSE)
  yMax = 1.2 * max( c( histInfo1$density , histInfo2$density ) )
  xVec = seq( xLim[1] , xLim[2] , length=501 )

  # Plot data y1 and smattering of posterior predictive curves:
  histInfo = hist( y1 , prob=TRUE , xlim=xLim , ylim=c(0,yMax) , breaks=xBreaks,
                   col="red2" , border="white" , xlab="y" , ylab="" , 
                   yaxt="n" , cex.lab=1.5 , 
                   main=paste("Data for",xLevels[1],"w. Post. Pred.") )
  for ( stepIdx in 1:length(stepIdxVec) ) {
    lines(xVec, dt( (xVec-mu1[stepIdxVec[stepIdx]])/sigma1[stepIdxVec[stepIdx]], 
                    df=nu[stepIdxVec[stepIdx]] )/sigma1[stepIdxVec[stepIdx]] , 
          type="l" , col="darkolivegreen3" , lwd=1 )
  }
  text( max(xVec) , yMax , bquote(N[1]==.(length(y1))) , adj=c(1.1,1.1) )

  # Plot data y2 and smattering of posterior predictive curves:
  histInfo = hist( y2 , prob=TRUE , xlim=xLim , ylim=c(0,yMax) , breaks=xBreaks,
                   col="red2" , border="white" , xlab="y" , ylab="" , 
                   yaxt="n" , cex.lab=1.5 , 
                   main=paste("Data for",xLevels[2],"w. Post. Pred.") )
  for ( stepIdx in 1:length(stepIdxVec) ) {
    lines(xVec, dt( (xVec-mu2[stepIdxVec[stepIdx]])/sigma2[stepIdxVec[stepIdx]], 
                    df=nu[stepIdxVec[stepIdx]] )/sigma2[stepIdxVec[stepIdx]] , 
          type="l" , col="darkolivegreen3" , lwd=1 )
  }
  text( max(xVec) , yMax , bquote(N[2]==.(length(y2))) , adj=c(1.1,1.1) )

  # Plot posterior distribution of parameter nu:
  histInfo = plotPost( log10(nu) , col="darkolivegreen3" , # breaks=30 ,
                       showCurve=showCurve ,
                       xlab=bquote("log10("*nu*")") , cex.lab = 1.75 , 
                       cenTend="mode" ,
                       main="Normality" ) #  (<0.7 suggests kurtosis)
  
  # Plot posterior distribution of parameters mu1, mu2, and their difference:
  xlim = range( c( mu1 , mu2 ) )
  histInfo = plotPost( mu1 ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                       xlab=bquote(mu[1]) , main=paste(xLevels[1],"Mean") , 
                       col="darkolivegreen3" )
  histInfo = plotPost( mu2 ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                       xlab=bquote(mu[2]) , main=paste(xLevels[2],"Mean") , 
                       col="darkolivegreen3" )
  histInfo = plotPost( mu2-mu1 , compVal=0 ,  showCurve=showCurve ,
                       xlab=bquote(mu[2] - mu[1]) , cex.lab = 1.75 , 
                       ROPE=RopeMuDiff,
                       main="Difference of Means" , col="darkolivegreen3" )

  # Plot posterior distribution of param's sigma1, sigma2, and their difference:
  xlim=range( c( sigma1 , sigma2 ) )
  histInfo = plotPost( sigma1 ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                       xlab=bquote(sigma[1]) , 
                       main=paste(xLevels[1],"Scale") , 
                       col="darkolivegreen3" , cenTend="mode" )
  histInfo = plotPost( sigma2 ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                       xlab=bquote(sigma[2]) , 
                       main=paste(xLevels[2],"Scale") , 
                       col="darkolivegreen3" , cenTend="mode" )
  histInfo = plotPost( sigma2 - sigma1 , 
                       compVal=0 ,  showCurve=showCurve ,
                       xlab=bquote(sigma[2] - sigma[1]) , cex.lab = 1.75 , 
                       ROPE=RopeSdDiff ,
                       main="Difference of Scales" , col="darkolivegreen3" , 
                       cenTend="mode" )

  # Plot effect size. 
  effectSize = ( mu2 - mu1 ) / sqrt( ( sigma1^2 + sigma2^2 ) / 2 )
  histInfo = plotPost( effectSize , compVal=0 ,  ROPE=RopeEff ,
                       showCurve=showCurve,
                       xlab="Effect Size",
                       cenTend="" , cex.lab=1.0 , main="" ,
                       col="darkolivegreen3" )

  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Post",sep=""), type=saveType)
  }
}


### High-level MCMC ###

## CC
graphics.off() 
rm(list=ls()) 

myDataFrame = read.csv(file="BDNF_CC.csv")
yName="Score"
xName="Group"
fileNameRoot = "BDNF_CC_" 
RopeMuDiff=c(-0.5,0.5) ; RopeSdDiff=c(-0.5,0.5) ; RopeEff=c(-0.1,0.1)

graphFileType = "eps" 

source("BDNF_Jags_Low.R")

mcmcCoda = genMCMC( datFrm=myDataFrame , yName=yName , xName=xName ,
                    numSavedSteps=50000 , saveName=fileNameRoot )

parameterNames = varnames(mcmcCoda) 
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}

summaryInfo = smryMCMC( mcmcCoda , RopeMuDiff=RopeMuDiff , 
                        RopeSdDiff=RopeSdDiff , RopeEff=RopeEff ,
                        saveName=fileNameRoot )
show(summaryInfo)

plotMCMC( mcmcCoda , datFrm=myDataFrame , yName=yName , xName=xName , 
          RopeMuDiff=RopeMuDiff , RopeSdDiff=RopeSdDiff , RopeEff=RopeEff ,
          pairsPlot=TRUE , saveName=fileNameRoot , saveType=graphFileType )

## WM
graphics.off() 
rm(list=ls()) 

myDataFrame = read.csv(file="BDNF_WM.csv")
yName="Score"
xName="Group"
fileNameRoot = "BDNF_WM_" 
RopeMuDiff=c(-0.5,0.5) ; RopeSdDiff=c(-0.5,0.5) ; RopeEff=c(-0.1,0.1)

graphFileType = "eps" 

source("BDNF_Jags_Low.R")

mcmcCoda = genMCMC( datFrm=myDataFrame , yName=yName , xName=xName ,
                    numSavedSteps=50000 , saveName=fileNameRoot )

parameterNames = varnames(mcmcCoda) 
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}

summaryInfo = smryMCMC( mcmcCoda , RopeMuDiff=RopeMuDiff , 
                        RopeSdDiff=RopeSdDiff , RopeEff=RopeEff ,
                        saveName=fileNameRoot )
show(summaryInfo)

plotMCMC( mcmcCoda , datFrm=myDataFrame , yName=yName , xName=xName , 
          RopeMuDiff=RopeMuDiff , RopeSdDiff=RopeSdDiff , RopeEff=RopeEff ,
          pairsPlot=TRUE , saveName=fileNameRoot , saveType=graphFileType )


