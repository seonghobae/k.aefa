# Kwangwoon Automated Exploratory Factor Analysis (K.AEFA)
# Seongho Bae (seongho@kw.ac.kr)
# last: 2015. 3. 13. 10am


##############
# aefa frontend #
##############

## usage
## k.aefa(brand.imputated, estimator="ML", variables = c(...))
## k.aefa(brand.imputated, estimator="ML", variables = attributes(brand.imputated)$names)
## attributes(brand.imputated)$names

Sys.setlocale('LC_CTYPE', 'ko_KR.UTF-8')
Sys.setlocale('LC_ALL', 'ko_KR.UTF-8')

message("\n Check Packages: Processing......")

# check packages
try(update.packages(ask = F, repos = 'http://cran.nexr.com'))

if(!require(depmixS4)) {
  try(install.packages("depmixS4", dependencies = TRUE), silent=TRUE); require(depmixS4)
}

if(!require(Rsolnp)) {
  try(install.packages("Rsolnp", dependencies = TRUE), silent=TRUE); require(Rsolnp)
}

if(!require(Cairo)) {
  try(install.packages("Cairo", dependencies = TRUE), silent=TRUE); require(Cairo)
}
if(!require(cairoDevice)) {
  try(install.packages("cairoDevice", dependencies = TRUE), silent=TRUE); require(cairoDevice)
}
if(!require(stringr)) {
  try(install.packages("stringr", dependencies = TRUE), silent=TRUE); require(stringr)
}

if(!require(SQUAREM)) {
  try(install.packages("SQUAREM", dependencies = TRUE), silent=TRUE); require(SQUAREM)
}

if(!require(rrcovNA)) {
  try(install.packages("rrcovNA", dependencies = TRUE), silent=TRUE); require(rrcovNA)
}

# check packages
if(!require(FAiR)) {
  message("\n Installing Packages: Processing......")
  try(install.packages("FAiR", dependencies = TRUE), silent = T)
}

# check packages
if(!require(bfa)) {
  try(install.packages("bfa", dependencies = TRUE), silent = T)
}

# check packages
if(!require(mirt)) {
  if(Sys.info()["sysname"] == "Linux" | Sys.info()["sysname"] == "Darwin"){
    try(install.packages("devtools", dependencies = TRUE), silent = T)
    try(library('devtools'), silent = T)
    try(install_github('philchalmers/mirt'), silent = T)
    try(install_github('philchalmers/mirtCAT'), silent = T)
  }
  else {
    try(install.packages("mirt", dependencies = TRUE, repos = 'http://cran.nexr.com'), silent = T)
  }
  
}

if(!require(latticeExtra)) {
  try(install.packages('latticeExtra', dependencies = TRUE), silent = T)
}

if(!require(plyr)) {
  try(install.packages('plyr', dependencies = TRUE), silent = T)
}

if(!require(multilevel)) {
  try(install.packages('multilevel', dependencies = TRUE), silent = T)
}

if(!require(nlme)) {
  try(install.packages('nlme', dependencies = TRUE), silent = T)
}

if(!require(lsr)) {
  try(install.packages('lsr', dependencies = TRUE), silent = T)
}

if(!require(meta)) {
  try(install.packages('meta', dependencies = TRUE), silent = T)
}

if(!require(metafor)) {
  try(install.packages('metafor', dependencies = TRUE), silent = T)
}

if(!require(rmeta)) {
  try(install.packages('rmeta', dependencies = TRUE), silent = T)
}
if(!require(psychometric)) {
  try(install.packages('psychometric', dependencies = TRUE), silent = T)
}
if(!require(pracma)) {
  try(install.packages('pracma', dependencies = TRUE), silent = T)
}
if(!require(rsm)) {
  try(install.packages('rsm', dependencies = TRUE), silent = T)
}
if(!require(car)) {
  try(install.packages('car', dependencies = TRUE), silent = T)
}
if(!require(TAM)) {
  try(install.packages('TAM', dependencies = TRUE), silent = T)
}
if(!require(GPArotation)) {
  try(install.packages('GPArotation', dependencies = TRUE), silent = T)
}  
if(!require(lavaan)) {
  try(install.packages("lavaan", dependencies = TRUE), silent = T)
}
if(!require(semTools)) {
  try(install.packages("semTools", dependencies = TRUE), silent = T)
}
if(!require(psych)) {
  try(install.packages("psych", dependencies = TRUE), silent = T)
}


#update.packages(ask = F, dependencies = TRUE, checkBuilt=TRUE)
message(" Checking Packages: OK")

message("\n Loading Packages: Processing......")
try(require(depmixS4), silent = T)
try(require(Rsolnp), silent = T)
try(require(Cairo), silent = T)
try(require(cairoDevice), silent = T)
try(require(stringr), silent = T)
try(require(SQUAREM), silent = T)
try(require(psychometric), silent = T)
try(require(psych), silent = T)
try(require(FAiR), silent = T)
try(require(bfa), silent = T)
try(require(mirt), silent = T)
try(require(latticeExtra), silent = T)
try(require(pracma), silent = T)
try(require(multilevel), silent = T)
try(require(nlme), silent = T)
try(require(lsr), silent = T)
try(require(rsm), silent = T)
try(require(car), silent = T)
try(require(TAM), silent = T)
try(require(GPArotation), silent = T)
try(require(lavaan), silent = T)
try(require(semTools), silent = T)



# try(mirtCluster(remove=TRUE), silent=TRUE)
# if(Sys.info()["sysname"] == "Linux"){
# try(mirtCluster(), silent=F)
# } else {
#   try(mirtCluster(as.numeric(parallel::detectCores() * 2)), silent=F)
# }


message(" Loading Packages: OK")

## pre-defined function for Full-information item factor analysis
findM2 <- function(mirtModel, increase = 15000, iterations = 1000, ...){
  
  # for treating unstable M2() function values
  quadpts <- 0 # initial quadpts
  for(i in 1:iterations){
    
    quadpts <- quadpts + increase
    
    message('quadpts: ', paste0(quadpts), ' / iteration: ', paste0(i))
    fitStats <- M2(mirtModel, QMC = TRUE, quadpts = quadpts, ...)
    
    if(i > 1){
      # decision-making
      if(abs((fitStats$RMSEA_5 - fitStats.old$RMSEA_5)) < .001 && 
           abs((fitStats$TLI - fitStats.old$TLI)) < .001 && 
           abs((fitStats$CFI - fitStats.old$CFI)) < .001)
        return(fitStats)
    }
    
    fitStats.old <- fitStats
  }
  
}

## Kwangwoon AEFA
k.aefa <- function(dataset, request_factors = ..., estimator = ..., variables = ..., fscore = ..., group = ..., demographic_range = ..., irtmodel = ..., se = ..., test = ..., control = ..., ...) {
  
  if(!exists(as.character(substitute(request_factors)))) {
    minimum_factors <- request_factors
    print(request_factors)
    message('\n')
  } else {
    request_factors = 1
    minimum_factors = 1
  }
  
  
  #if(!request_factors) {request_factors = 100}
  
  message("\nKwangwoon Automated Exploratory Factor Analysis [k.aefa]\nSeongho Bae(seongho@kw.ac.kr)\n\n\nChecking some Packages: Processing......")
  
  # make worksheet to doing data imputation
  workingdata <- subset(dataset, select = variables)
  
  if(estimator == "irt" | estimator == "IRT" | estimator == "oa" | estimator == "mcmc" | estimator == "MCMC" | estimator == "bfa" | estimator == "BFA" | estimator == "BCFA" | estimator == "mirt" | estimator == "WLSMV" | estimator == "wlsmv" | estimator == "sem" | estimator == "lavaan"  | estimator == "multiplegroup" | estimator == "multigroup" | estimator == "multilevel" | estimator == "mixedmirt" | estimator == "WLSM" | estimator == "wlsm" | estimator == "wlsmvs" | estimator == "WLSMVS" | estimator == "TAM"  | estimator == "tam") {
    fa_covdata <- workingdata
    fa_covdata <- subset(fa_covdata, rowSums(is.na(fa_covdata)) != ncol(fa_covdata))
  }
  
  else if(!sum(is.na(workingdata)==1)==0) {
    message("\n\nWarning: \nI was found missing data.\nI'll fill out missing data now! \nIt may be wait for longtime.\n\nProcessing......")
    
    corrected <- k.imputation(workingdata, ...)
    fa_covdata <- corrected
    
  } else {fa_covdata <- workingdata}
  
  message("\n\n Raw data Processing......")
  print(apply(fa_covdata, 2, table))
  
  
  #############################################
  #                                           #
  # Let's Start Here!                         #
  # This is Bayesian Couplar Factor Analysis  #
  # (Maybe be good at the Econometric Data)   #
  #                                           #
  #############################################
  
  # Bayesian FA
  if(estimator == "mcmc" | estimator == "MCMC" | estimator == "bcfa" | estimator == "bfa" | estimator == "BFA" | estimator == "BCFA"){ 
    
    
    message("\n\nFactor number(s):", paste0(request_factors))
    
    message("\n\nProcessing... please wait ");
    bfa.result <- bfa_copula(~., data=fa_covdata, num.factor=request_factors, px=TRUE, nsim=1000, nburn=100, thin=4, loading.prior="gdp", keep.scores=T, print.status=250)
    if(request_factors==1){
      
      bfa.rotated <- mean(bfa.result)$loadings
      bcfa.loadings <- data.frame(bfa.rotated)
      row.names(bcfa.loadings) <- colnames(fa_covdata)
      
      h2 <- bfa.rotated^2
      bcfa.loadings$h2 <- data.frame(h2)
      
      message("\nLoadings\n")
      print(bcfa.loadings)
      
      message("\nEigenvalues (geominQ rotation)\n")
      
      Eigenvalues <- colSums(bfa.rotated^2)
      print(Eigenvalues)
      
      message("\nVariance Explained (geominQ rotation)\n")
      
      Explained <- Eigenvalues / sum(h2)
      print(Explained)
      
      message("\nTotal Variance Explained (geominQ rotation)\n")
      print(sum(Explained))
      
      return(bcfa.loadings)
      
    } else {
      message('Rotating Factor Solution now')
      bfa.rotated <- geominQ(mean(bfa.result)$loadings, maxit=999999)
      
      bcfa.loadings <- data.frame(bfa.rotated$loadings)
      row.names(bcfa.loadings) <- colnames(fa_covdata)
      
      bcfa.phi <- data.frame(bfa.rotated$Phi)
      row.names(bcfa.phi) <- colnames(bcfa.loadings)
      
      
      #h2 <- abs(rowSums(apply(bcfa.loadings, c(1,2), cumprod)))
      h2 <- rowSums((bfa.rotated$loadings^2))
      bcfa.loadings$h2 <- data.frame(h2)
      
      message("\nLoadings (geominQ rotation)\n")
      print(bcfa.loadings)
      
      message("\nEigenvalues (geominQ rotation)\n")
      
      Eigenvalues <- colSums(bfa.rotated$loadings^2)
      print(Eigenvalues)
      
      message("\nVariance Explained (geominQ rotation)\n")
      
      Explained <- Eigenvalues / sum(h2)
      print(Explained)
      
      message("\nTotal Variance Explained (geominQ rotation)\n")
      print(sum(Explained))
      
      message("\nfactor correlation (geominQ rotation)\n")
      print(bcfa.phi)
      
      #plot(get_coda(bfa.result))
      
      return(bcfa.loadings)
    }
    
  }
  
  ################################################
  #                                              #
  # This is Unidimensional Item Response Theory  #
  #                                              #
  ################################################
  
  # UIRT
  else if(estimator == "irt" | estimator == "IRT") {
    
    
    for(i in request_factors:100) {
      message("\n\nFactor number(s): ", paste0(i))
      result <- irt.fa(fa_covdata, rotate='geominQ', nfactors=i, correct=TRUE, plot=FALSE, maxit=100000)
      print(result, sort=TRUE ,short=FALSE)
      #return(result)
      message("\n\n\n\n\n")
      
      # heywood case check
      
      
      # Model fit check
      
      if(is.na(result$fa$RMSEA[2])) {
        result$fa$RMSEA[2] <- result$fa$RMSEA[1]
        message('RMSEA 90% Confidence Interval estimation were failed')
      }
      
      if(result$fa$TLI[1] >= .9 && result$fa$RMSEA[2] < .05) {
        message("\n\n Maximum Factor number archieved: ", paste0(i))
        message("\n\n Perfect!")
        return(result)
        break()
      } else if(result$fa$TLI[1] >= .9 | result$fa$RMSEA[2] < .08) {
        message("\n\n Maximum Factor number archieved: ", paste0(i))
        
      } else if(result$fa$TLI[1] >= .99 | result$fa$TLI[1] == "-Inf") {
        message("\n\n Warning: This model is something was wrong!")
        i=i-1
        result <- irt.fa(fa_covdata, nfactors=i, correct=TRUE, plot=TRUE)
        print(result, sort=TRUE)
        message("\n\n Use this result. Thanks a lot!")
        break()
      }else{
        message("\n\n I will extract more factors. Current factor number(s): ", paste0(i))}
      
    }
  }
  
  ##################################################
  #                                                #
  # This is Multidiminsional Item Response Theory  #
  #                                                #
  ##################################################
  
  # MIRT (MML -- tam) -- alpha version (need to implement bcfa code)
  else if(estimator == "tam" | estimator == "TAM"){ 
    
    request_factors <- 2
    message("\n\nProcessing... please wait ");
    for(i in request_factors:ncol(fa_covdata)){
      
      message("\n\nFactor number(s): ", paste0(i))
      
      if(i > request_factors){
        mod1a_old <- mod1a
      }
      try(mod1a <- tam.fa(resp=fa_covdata, irtmodel="efa", nfactors=i, control=list(maxiter=999999, maxit = 1e+5, QMC = T, Msteps = 2000, snodes = 1000, increment.factor=2, fac.oldxsi=.5), ...), silent = T)
      if(i > request_factors){
        if(mod1a_old$ic$aBIC <= mod1a$ic$aBIC | mod1a_old$ic$AICc <= mod1a$ic$AICc | colSums(data.frame(mod1a$EAP.rel) < .7) > 0){
          rotation <- geominQ(mod1a_old$B.stand, maxit=999999)
          paste0(rotation)
          
          message("rotation: geominQ(mod$B.stand, maxit=999999)")
          
          return(mod1a_old)
        }
      }
    }
  }
  
  
  # MIRT (WLSMV) -- alpha version (model fit evaluation is first!)
  else if(estimator == "WLSMV" | estimator == "wlsmv" | estimator == "sem" | estimator == "lavaan" | estimator == "MLR" | estimator == "MLM" | estimator == "MLMV" | estimator == "MLMVS" | estimator == "WLSM" | estimator == "wlsm" | estimator == "wlsmvs"  | estimator == "WLSMVS" | estimator == "PML" | estimator == "pml"){
    
    if(estimator == "WLSMV" | estimator == "wlsmv" | estimator == "sem" | estimator == "lavaan"){
      factoring_method <- 'WLSMV'
    } else if(estimator == "WLSM" | estimator == "wlsm" | estimator == "wlsmvs"  | estimator == "WLSMVS") {
      factoring_method <- estimator
    }
    else {
      factoring_method <- estimator
    }
    
    # Ignore RMSEA evaluation if small items where evaluation of model fits
    if(ncol(fa_covdata) < 15){
      RMSEA_Criterion <- .05
      #fa_covdata <- k.rescale(fa_covdata, 2)
    } else {
      RMSEA_Criterion <- .05
    }
    
    for(i in 1:100) {
      
      if(estimator == "MLR" | estimator == "MLM" | estimator == "MLMV" | estimator == "MLMVS"){
        message("\n\n\n Exploratory Factor Analysis with Lavaan \n")
      } else {
        message("\n\n\n Exploratory Factor Analysis using Limited-information Multidimensional Item Response Theory with Lavaan (PML or WLSMV Estimator; Mean- and Variance-adjusted pairwise maximum likelihood or Weighted Least Squares for categorical data method) \n")
        message("\n Please be Wait! \n")
      }
      
      
      message("\n\n\nFactor number(s): ", paste0(i))
      
      if(i > 1) {
        unrotated_old <- unrotated
        Model_Fit_old <- Model_Fit
        remove(unrotated); remove(Model_Fit)
      }
      
      # calculating factor solution
      message('\nCalculating factor solution...\n')
      if(estimator == "MLR" | estimator == "MLM" | estimator == "MLMV" | estimator == "MLMVS"){
        unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method)        
      } else {
        
        #unrotated <- efaUnrotate(fa_covdata, nf=i, estimator="WLSMV", ordered=names(fa_covdata), verbose=TRUE)
        unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method, missing='default', ordered=names(fa_covdata), start=FALSE, verbose=T, zero.add = c(.01, .01), ...)
        
      }
      message('Done...\n')
      
      message('\nCalculating Model fit Indices...\n')
      if(estimator == 'pml' | estimator == 'PML'){
        Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
        Model_Fit <- round(Model_Fit, 4)
        Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
        message('Done...\n')
      } else {
        Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
        Model_Fit <- round(Model_Fit, 4)
        Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
        message('Done...\n')
      }
      
      
      message("\n\nResults of Model Fit Indices")
      colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
      print(Model_Fit)      
      
      message("\n\n\n\n\n")
      
      
      # Model fit check
      #if(Model_Fit[1,1] >= .99999 | Model_Fit[2,1] >= .99999 | !sum(is.na(R_SQUARE)==1)==0 | as.character(Model_Fit[1,1]) == "ERROR" | as.character(Model_Fit[2,1]) == "ERROR" | as.character(Model_Fit[3,1]) == "ERROR") {
      if(Model_Fit[1,1] >= 1 | Model_Fit[2,1] >= 1 | as.character(Model_Fit[1,1]) == "ERROR" | as.character(Model_Fit[2,1]) == "ERROR" | as.character(Model_Fit[3,1]) == "ERROR") {
        
        message("\n\n Warning: This model is something was wrong!")
        
        i = i-1
        
        message("\n\n\nFactor number(s): ", paste0(i))
        
        # calculating factor solution
        message('\nCalculating factor solution...\n')
        if(estimator == "MLR" | estimator == "MLM" | estimator == "MLMV" | estimator == "MLMVS"){
          unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method)        
        } else {
          
          #unrotated <- efaUnrotate(fa_covdata, nf=i, estimator="WLSMV", ordered=names(fa_covdata), verbose=TRUE)
          unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method, missing='default', ordered=names(fa_covdata), start=FALSE, verbose=T, zero.add = c(.01, .01), ...)
          
        }
        message('Done...\n')
        
        
        
        if(i==1){
          summary(unrotated, std=TRUE)
          
          if(estimator == 'pml' | estimator == 'PML'){
            Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
            Model_Fit <- round(Model_Fit, 4)
            Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
            message('Done...\n')
          } else {
            Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
            Model_Fit <- round(Model_Fit, 4)
            Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
            message('Done...\n')
          }
          
          message("\n\nResults of Model Fit Indices")
          colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
          print(Model_Fit)      
          
          message("\n\n Use this result. Thanks a lot!")
          print(reliability(unrotated))
          return(unrotated)
        } else {
          #unrotated <- unrotated_old
          
          # estimate R-square (h2)
          message('\nCalculating Rotated factor solution...\n')
          rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
          message('Done...\n')
          
          message('\nCalculating Communality...\n')
          if(sum(is.na(data.frame(inspect(unrotated, "rsquare")))==1)==0){
            R_SQUARE <- data.frame(inspect(unrotated, "rsquare"))
            COMMUNALITY_METHOD <- estimator
          } else {
            # SMC
            message('\nCalculating SMC for h2...\n')
            loading <- data.frame(rotated@loading)
            R_SQUARE <- rowSums((loading^2))
            COMMUNALITY_METHOD <- 'SMC'
          }
          message('Done...\n')
          
          message("\nLoadings\n")
          #rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
          h2 <- R_SQUARE
          loading <- data.frame(rotated@loading)
          loading$h2 <- data.frame(h2)
          print(loading)
          
          message("\nCommunality\n")
          print(COMMUNALITY_METHOD)
          
          message("\nEigenvalues\n")
          Eigenvalues <- abs(colSums(attributes(rotated)$loading^2))
          print(Eigenvalues)
          
          
          message("\nVariance Explained\n")
          Explained <- Eigenvalues / sum(h2)
          print(Explained)
          
          message("\nTotal Variance Explained\n")
          print(sum(Explained))
          
          message("\nSummary results\n")
          print(summary(rotated))
          
          
          if(estimator == 'pml' | estimator == 'PML'){
            Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
            Model_Fit <- round(Model_Fit, 4)
#             Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
            message('Done...\n')
          } else {
            Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
            Model_Fit <- round(Model_Fit, 4)
            Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
            message('Done...\n')
          }
          
          message("\n\n Model Fit Indices")
          colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
          print(Model_Fit)      
          
          # if one of Eigenvalues under 1, decreasing factor numbers 'm' - Step 1
          if(sum(Eigenvalues < 1) > 0) {
            message("\n\n Warning: This model is can not be readable!")
            
            
            i = i-1
            
            message("\n\n\nFactor number(s): ", paste0(i))
            
            # calculating factor solution
            message('\nCalculating factor solution...\n')
            if(estimator == "MLR" | estimator == "MLM" | estimator == "MLMV" | estimator == "MLMVS"){
              unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method)        
            } else {
              
              #unrotated <- efaUnrotate(fa_covdata, nf=i, estimator="WLSMV", ordered=names(fa_covdata), verbose=TRUE)
              unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method, missing='default', ordered=names(fa_covdata), start=FALSE, verbose=T, zero.add = c(.01, .01), ...)
              
            }
            message('Done...\n')
            
            
            
            if(i==1){
              summary(unrotated, std=TRUE)
              
              if(estimator == 'pml' | estimator == 'PML'){
                Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                Model_Fit <- round(Model_Fit, 4)
#                 Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                message('Done...\n')
              } else {
                Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                Model_Fit <- round(Model_Fit, 4)
                Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                message('Done...\n')
              }
              
              message("\n\nResults of Model Fit Indices")
              colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
              print(Model_Fit)      
              
              message("\n\n Use this result. Thanks a lot!")
              print(reliability(unrotated))
              return(unrotated)
            } else {
              #unrotated <- unrotated_old
              
              # estimate R-square (h2)
              message('\nCalculating Rotated factor solution...\n')
              rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
              message('Done...\n')
              
              message('\nCalculating Communality...\n')
              if(sum(is.na(data.frame(inspect(unrotated, "rsquare")))==1)==0){
                R_SQUARE <- data.frame(inspect(unrotated, "rsquare"))
                COMMUNALITY_METHOD <- estimator
              } else {
                # SMC
                message('\nCalculating SMC for h2...\n')
                loading <- data.frame(rotated@loading)
                R_SQUARE <- rowSums((loading^2))
                COMMUNALITY_METHOD <- 'SMC'
              }
              message('Done...\n')
              
              message("\nLoadings\n")
              #rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
              h2 <- R_SQUARE
              loading <- data.frame(rotated@loading)
              loading$h2 <- data.frame(h2)
              print(loading)
              
              message("\nCommunality\n")
              print(COMMUNALITY_METHOD)
              
              message("\nEigenvalues\n")
              Eigenvalues <- abs(colSums(attributes(rotated)$loading^2))
              print(Eigenvalues)
              
              
              message("\nVariance Explained\n")
              Explained <- Eigenvalues / sum(h2)
              print(Explained)
              
              message("\nTotal Variance Explained\n")
              print(sum(Explained))
              
              message("\nSummary results\n")
              print(summary(rotated))
              
              
              if(estimator == 'pml' | estimator == 'PML'){
                Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                Model_Fit <- round(Model_Fit, 4)
#                 Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                message('Done...\n')
              } else {
                Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                Model_Fit <- round(Model_Fit, 4)
                Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                message('Done...\n')
              }
              
              message("\n\n Model Fit Indices")
              colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
              print(Model_Fit)      
              
              # if one of Eigenvalues under 1, decreasing factor numbers 'm' - Step 2
              if(sum(Eigenvalues < 1) > 0) {
                message("\n\n Warning: This model is can not be readable!")
                
                
                i = i-1
                
                message("\n\n\nFactor number(s): ", paste0(i))
                
                # calculating factor solution
                message('\nCalculating factor solution...\n')
                if(estimator == "MLR" | estimator == "MLM" | estimator == "MLMV" | estimator == "MLMVS"){
                  unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method)        
                } else {
                  
                  #unrotated <- efaUnrotate(fa_covdata, nf=i, estimator="WLSMV", ordered=names(fa_covdata), verbose=TRUE)
                  unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method, missing='default', ordered=names(fa_covdata), start=FALSE, verbose=T, zero.add = c(.01, .01), ...)
                  
                }
                message('Done...\n')
                
                
                
                if(i==1){
                  summary(unrotated, std=TRUE)
                  
                  if(estimator == 'pml' | estimator == 'PML'){
                    Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                    Model_Fit <- round(Model_Fit, 4)
#                     Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                    message('Done...\n')
                  } else {
                    Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                    Model_Fit <- round(Model_Fit, 4)
                    Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                    message('Done...\n')
                  }
                  
                  message("\n\nResults of Model Fit Indices")
                  colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
                  print(Model_Fit)      
                  
                  message("\n\n Use this result. Thanks a lot!")
                  print(reliability(unrotated))
                  return(unrotated)
                } else {
                  #unrotated <- unrotated_old
                  
                  # estimate R-square (h2)
                  message('\nCalculating Rotated factor solution...\n')
                  rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
                  message('Done...\n')
                  
                  message('\nCalculating Communality...\n')
                  if(sum(is.na(data.frame(inspect(unrotated, "rsquare")))==1)==0){
                    R_SQUARE <- data.frame(inspect(unrotated, "rsquare"))
                    COMMUNALITY_METHOD <- estimator
                  } else {
                    # SMC
                    message('\nCalculating SMC for h2...\n')
                    loading <- data.frame(rotated@loading)
                    R_SQUARE <- rowSums((loading^2))
                    COMMUNALITY_METHOD <- 'SMC'
                  }
                  message('Done...\n')
                  
                  message("\nLoadings\n")
                  #rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
                  h2 <- R_SQUARE
                  loading <- data.frame(rotated@loading)
                  loading$h2 <- data.frame(h2)
                  print(loading)
                  
                  message("\nCommunality\n")
                  print(COMMUNALITY_METHOD)
                  
                  message("\nEigenvalues\n")
                  Eigenvalues <- abs(colSums(attributes(rotated)$loading^2))
                  print(Eigenvalues)
                  
                  
                  message("\nVariance Explained\n")
                  Explained <- Eigenvalues / sum(h2)
                  print(Explained)
                  
                  message("\nTotal Variance Explained\n")
                  print(sum(Explained))
                  
                  message("\nSummary results\n")
                  print(summary(rotated))
                  
                  
                  if(estimator == 'pml' | estimator == 'PML'){
                    Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                    Model_Fit <- round(Model_Fit, 4)
#                     Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                    message('Done...\n')
                  } else {
                    Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                    Model_Fit <- round(Model_Fit, 4)
                    Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                    message('Done...\n')
                  }
                  
                  message("\n\n Model Fit Indices")
                  colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
                  print(Model_Fit)      
                  
                  
                  message("\n\n Use this result. Thanks a lot!")
                  return(rotated)
                }
                
                message("\n\n Use this result. Thanks a lot!")
                return(rotated)
              }
              
            }
            
            message("\n\n Use this result. Thanks a lot!")
            return(rotated)
          }
          
          message("\n\n Use this result. Thanks a lot!")
          return(rotated)
        }
        
        
      }
      
      else if(Model_Fit[1,1] >= .9 && Model_Fit[2,1] >= .9 && Model_Fit[4,1] < RMSEA_Criterion) {
        message("\n\n Maximum Factor number archieved: ", paste0(i))
        message("\n\n Statistically Perfect! This solution is final solution.")
        
        if(i==1){
          summary(unrotated, std=TRUE)
          print(reliability(unrotated))
          return(unrotated)
        } else {
          
          # estimate R-square (h2)
          message('\nCalculating Rotated factor solution...\n')
          rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
          message('Done...\n')
          
          message('\nCalculating communality...\n')
          if(sum(is.na(data.frame(inspect(unrotated, "rsquare")))==1)==0){
            R_SQUARE <- data.frame(inspect(unrotated, "rsquare"))
            COMMUNALITY_METHOD <- estimator
          } else {
            # SMC
            message('\nCalculating SMC for h2...\n')
            loading <- data.frame(rotated@loading)
            R_SQUARE <- rowSums((loading^2))
            COMMUNALITY_METHOD <- 'SMC'
          }
          message('Done...\n')
          
          message("\nLoadings\n")
          #rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
          h2 <- R_SQUARE
          loading <- data.frame(rotated@loading)
          loading$h2 <- data.frame(h2)
          print(loading)
          
          message("\nCommunality\n")
          print(COMMUNALITY_METHOD)
          
          message("\nEigenvalues\n")
          Eigenvalues <- abs(colSums(attributes(rotated)$loading^2))
          print(Eigenvalues)
          
          
          message("\nVariance Explained\n")
          Explained <- Eigenvalues / sum(h2)
          print(Explained)
          
          message("\nTotal Variance Explained\n")
          print(sum(Explained))
          
          message("\nSummary results\n")
          print(summary(rotated))
          
          
          if(estimator == 'pml' | estimator == 'PML'){
            Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
            Model_Fit <- round(Model_Fit, 4)
#             Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
            message('Done...\n')
          } else {
            Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
            Model_Fit <- round(Model_Fit, 4)
            Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
            message('Done...\n')
          }
          
          message("\n\n Model Fit Indices")
          colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
          print(Model_Fit)      
          
          # if one of Eigenvalues under 1, decreasing factor numbers 'm' - Step 1
          if(sum(Eigenvalues < 1) > 0) {
            message("\n\n Warning: This model is can not be readable!")
            
            
            i = i-1
            
            message("\n\n\nFactor number(s): ", paste0(i))
            
            # calculating factor solution
            message('\nCalculating factor solution...\n')
            if(estimator == "MLR" | estimator == "MLM" | estimator == "MLMV" | estimator == "MLMVS"){
              unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method)        
            } else {
              
              #unrotated <- efaUnrotate(fa_covdata, nf=i, estimator="WLSMV", ordered=names(fa_covdata), verbose=TRUE)
              unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method, missing='default', ordered=names(fa_covdata), start=FALSE, verbose=T, zero.add = c(.01, .01), ...)
              
            }
            message('Done...\n')
            
            
            
            if(i==1){
              summary(unrotated, std=TRUE)
              
              if(estimator == 'pml' | estimator == 'PML'){
                Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                Model_Fit <- round(Model_Fit, 4)
#                 Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                message('Done...\n')
              } else {
                Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                Model_Fit <- round(Model_Fit, 4)
                Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                message('Done...\n')
              }
              
              message("\n\nResults of Model Fit Indices")
              colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
              print(Model_Fit)      
              
              message("\n\n Use this result. Thanks a lot!")
              print(reliability(unrotated))
              return(unrotated)
            } else {
              #unrotated <- unrotated_old
              
              # estimate R-square (h2)
              message('\nCalculating Rotated factor solution...\n')
              rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
              message('Done...\n')
              
              message('\nCalculating Communality...\n')
              if(sum(is.na(data.frame(inspect(unrotated, "rsquare")))==1)==0){
                R_SQUARE <- data.frame(inspect(unrotated, "rsquare"))
                COMMUNALITY_METHOD <- estimator
              } else {
                # SMC
                message('\nCalculating SMC for h2...\n')
                loading <- data.frame(rotated@loading)
                R_SQUARE <- rowSums((loading^2))
                COMMUNALITY_METHOD <- 'SMC'
              }
              message('Done...\n')
              
              message("\nLoadings\n")
              #rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
              h2 <- R_SQUARE
              loading <- data.frame(rotated@loading)
              loading$h2 <- data.frame(h2)
              print(loading)
              
              message("\nCommunality\n")
              print(COMMUNALITY_METHOD)
              
              message("\nEigenvalues\n")
              Eigenvalues <- abs(colSums(attributes(rotated)$loading^2))
              print(Eigenvalues)
              
              
              message("\nVariance Explained\n")
              Explained <- Eigenvalues / sum(h2)
              print(Explained)
              
              message("\nTotal Variance Explained\n")
              print(sum(Explained))
              
              message("\nSummary results\n")
              print(summary(rotated))
              
              
              if(estimator == 'pml' | estimator == 'PML'){
                Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                Model_Fit <- round(Model_Fit, 4)
#                 Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                message('Done...\n')
              } else {
                Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                Model_Fit <- round(Model_Fit, 4)
                Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                message('Done...\n')
              }
              
              message("\n\n Model Fit Indices")
              colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
              print(Model_Fit)      
              
              # if one of Eigenvalues under 1, decreasing factor numbers 'm' - Step 2
              if(sum(Eigenvalues < 1) > 0) {
                message("\n\n Warning: This model is can not be readable!")
                
                
                i = i-1
                
                message("\n\n\nFactor number(s): ", paste0(i))
                
                # calculating factor solution
                message('\nCalculating factor solution...\n')
                if(estimator == "MLR" | estimator == "MLM" | estimator == "MLMV" | estimator == "MLMVS"){
                  unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method)        
                } else {
                  
                  #unrotated <- efaUnrotate(fa_covdata, nf=i, estimator="WLSMV", ordered=names(fa_covdata), verbose=TRUE)
                  unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method, missing='default', ordered=names(fa_covdata), start=FALSE, verbose=T, zero.add = c(.01, .01), ...)
                  
                }
                message('Done...\n')
                
                
                
                if(i==1){
                  summary(unrotated, std=TRUE)
                  
                  if(estimator == 'pml' | estimator == 'PML'){
                    Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                    Model_Fit <- round(Model_Fit, 4)
#                     Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                    message('Done...\n')
                  } else {
                    Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                    Model_Fit <- round(Model_Fit, 4)
                    Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                    message('Done...\n')
                  }
                  
                  message("\n\nResults of Model Fit Indices")
                  colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
                  print(Model_Fit)      
                  
                  message("\n\n Use this result. Thanks a lot!")
                  print(reliability(unrotated))
                  return(unrotated)
                } else {
                  #unrotated <- unrotated_old
                  
                  # estimate R-square (h2)
                  message('\nCalculating Rotated factor solution...\n')
                  rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
                  message('Done...\n')
                  
                  message('\nCalculating Communality...\n')
                  if(sum(is.na(data.frame(inspect(unrotated, "rsquare")))==1)==0){
                    R_SQUARE <- data.frame(inspect(unrotated, "rsquare"))
                    COMMUNALITY_METHOD <- estimator
                  } else {
                    # SMC
                    message('\nCalculating SMC for h2...\n')
                    loading <- data.frame(rotated@loading)
                    R_SQUARE <- rowSums((loading^2))
                    COMMUNALITY_METHOD <- 'SMC'
                  }
                  message('Done...\n')
                  
                  message("\nLoadings\n")
                  #rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
                  h2 <- R_SQUARE
                  loading <- data.frame(rotated@loading)
                  loading$h2 <- data.frame(h2)
                  print(loading)
                  
                  message("\nCommunality\n")
                  print(COMMUNALITY_METHOD)
                  
                  message("\nEigenvalues\n")
                  Eigenvalues <- abs(colSums(attributes(rotated)$loading^2))
                  print(Eigenvalues)
                  
                  
                  message("\nVariance Explained\n")
                  Explained <- Eigenvalues / sum(h2)
                  print(Explained)
                  
                  message("\nTotal Variance Explained\n")
                  print(sum(Explained))
                  
                  message("\nSummary results\n")
                  print(summary(rotated))
                  
                  
                  if(estimator == 'pml' | estimator == 'PML'){
                    Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                    Model_Fit <- round(Model_Fit, 4)
#                     Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                    message('Done...\n')
                  } else {
                    Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                    Model_Fit <- round(Model_Fit, 4)
                    Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                    message('Done...\n')
                  }
                  
                  message("\n\n Model Fit Indices")
                  colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
                  print(Model_Fit)      
                  
                  
                  message("\n\n Use this result. Thanks a lot!")
                  return(rotated)
                }
                
                message("\n\n Use this result. Thanks a lot!")
                return(rotated)
              }
              
            }
            
            message("\n\n Use this result. Thanks a lot!")
            return(rotated)
          }
          
          message("\n\n Use this result. Thanks a lot!")
          
          
          return(rotated)
        }
        
      }         
      
      else {
        # protect for overfactoring
        if(i > 1) {
          if(Model_Fit_old[1,1] >= Model_Fit[1,1] | Model_Fit_old[2,1] >= Model_Fit[2,1] | Model_Fit_old[4,1] <= Model_Fit[4,1]) {
            
            message("\n\n Warning: Overfactoring Protecton function was activated!")
            
            i = i-1
            
            message("\n\n\nFactor number(s): ", paste0(i))
            
            # calculating factor solution
            message('\nCalculating factor solution...\n')
            if(estimator == "MLR" | estimator == "MLM" | estimator == "MLMV" | estimator == "MLMVS"){
              unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method)        
            } else {
              
              #unrotated <- efaUnrotate(fa_covdata, nf=i, estimator="WLSMV", ordered=names(fa_covdata), verbose=TRUE)
              unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method, missing='default', ordered=names(fa_covdata), start=FALSE, verbose=T, ...)
              
            }
            message('Done...\n')
            
            
            
            if(i==1){
              summary(unrotated, std=TRUE)
              
              if(estimator == 'pml' | estimator == 'PML'){
                Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                Model_Fit <- round(Model_Fit, 4)
#                 Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                message('Done...\n')
              } else {
                Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                Model_Fit <- round(Model_Fit, 4)
                Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                message('Done...\n')
              }
              
              message("\n\nResults of Model Fit Indices")
              colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
              print(Model_Fit)      
              
              message("\n\n Use this result. Thanks a lot!")
              print(reliability(unrotated))
              return(unrotated)
            } else {
              #unrotated <- unrotated_old
              
              # estimate R-square (h2)
              message('\nCalculating Rotated factor solution...\n')
              rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
              message('Done...\n')
              
              message('\nCalculating communality...\n')
              if(sum(is.na(data.frame(inspect(unrotated, "rsquare")))==1)==0){
                R_SQUARE <- data.frame(inspect(unrotated, "rsquare"))
                COMMUNALITY_METHOD <- estimator
              } else {
                # SMC
                message('\nCalculating SMC for h2...\n')
                loading <- data.frame(rotated@loading)
                R_SQUARE <- rowSums((loading^2))
                COMMUNALITY_METHOD <- 'SMC'
              }
              message('Done...\n')
              
              message("\nLoadings\n")
              #rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
              h2 <- R_SQUARE
              loading <- data.frame(rotated@loading)
              loading$h2 <- data.frame(h2)
              print(loading)
              
              message("\nCommunality\n")
              print(COMMUNALITY_METHOD)
              
              message("\nEigenvalues\n")
              Eigenvalues <- abs(colSums(attributes(rotated)$loading^2))
              print(Eigenvalues)
              
              
              message("\nVariance Explained\n")
              Explained <- Eigenvalues / sum(h2)
              print(Explained)
              
              message("\nTotal Variance Explained\n")
              print(sum(Explained))
              
              message("\nSummary results\n")
              print(summary(rotated))
              
              
              if(estimator == 'pml' | estimator == 'PML'){
                Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                Model_Fit <- round(Model_Fit, 4)
#                 Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                message('Done...\n')
              } else {
                Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                Model_Fit <- round(Model_Fit, 4)
                Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                message('Done...\n')
              }
              
              message("\n\n Model Fit Indices")
              colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
              print(Model_Fit)      
              
              # if one of Eigenvalues under 1, decreasing factor numbers 'm' - Step 1
              if(sum(Eigenvalues < 1) > 0) {
                message("\n\n Warning: This model is can not be readable!")
                
                
                i = i-1
                
                message("\n\n\nFactor number(s): ", paste0(i))
                
                # calculating factor solution
                message('\nCalculating factor solution...\n')
                if(estimator == "MLR" | estimator == "MLM" | estimator == "MLMV" | estimator == "MLMVS"){
                  unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method)        
                } else {
                  
                  #unrotated <- efaUnrotate(fa_covdata, nf=i, estimator="WLSMV", ordered=names(fa_covdata), verbose=TRUE)
                  unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method, missing='default', ordered=names(fa_covdata), start=FALSE, verbose=T, zero.add = c(.01, .01), ...)
                  
                }
                message('Done...\n')
                
                
                
                if(i==1){
                  summary(unrotated, std=TRUE)
                  
                  if(estimator == 'pml' | estimator == 'PML'){
                    Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                    Model_Fit <- round(Model_Fit, 4)
#                     Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                    message('Done...\n')
                  } else {
                    Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                    Model_Fit <- round(Model_Fit, 4)
                    Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                    message('Done...\n')
                  }
                  
                  message("\n\nResults of Model Fit Indices")
                  colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
                  print(Model_Fit)      
                  
                  message("\n\n Use this result. Thanks a lot!")
                  print(reliability(unrotated))
                  return(unrotated)
                } else {
                  #unrotated <- unrotated_old
                  
                  # estimate R-square (h2)
                  message('\nCalculating Rotated factor solution...\n')
                  rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
                  message('Done...\n')
                  
                  message('\nCalculating Communality...\n')
                  if(sum(is.na(data.frame(inspect(unrotated, "rsquare")))==1)==0){
                    R_SQUARE <- data.frame(inspect(unrotated, "rsquare"))
                    COMMUNALITY_METHOD <- estimator
                  } else {
                    # SMC
                    message('\nCalculating SMC for h2...\n')
                    loading <- data.frame(rotated@loading)
                    R_SQUARE <- rowSums((loading^2))
                    COMMUNALITY_METHOD <- 'SMC'
                  }
                  message('Done...\n')
                  
                  message("\nLoadings\n")
                  #rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
                  h2 <- R_SQUARE
                  loading <- data.frame(rotated@loading)
                  loading$h2 <- data.frame(h2)
                  print(loading)
                  
                  message("\nCommunality\n")
                  print(COMMUNALITY_METHOD)
                  
                  message("\nEigenvalues\n")
                  Eigenvalues <- abs(colSums(attributes(rotated)$loading^2))
                  print(Eigenvalues)
                  
                  
                  message("\nVariance Explained\n")
                  Explained <- Eigenvalues / sum(h2)
                  print(Explained)
                  
                  message("\nTotal Variance Explained\n")
                  print(sum(Explained))
                  
                  message("\nSummary results\n")
                  print(summary(rotated))
                  
                  
                  if(estimator == 'pml' | estimator == 'PML'){
                    Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                    Model_Fit <- round(Model_Fit, 4)
#                     Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                    message('Done...\n')
                  } else {
                    Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                    Model_Fit <- round(Model_Fit, 4)
                    Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                    message('Done...\n')
                  }
                  
                  message("\n\n Model Fit Indices")
                  colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
                  print(Model_Fit)      
                  
                  # if one of Eigenvalues under 1, decreasing factor numbers 'm' - Step 2
                  if(sum(Eigenvalues < 1) > 0) {
                    message("\n\n Warning: This model is can not be readable!")
                    
                    
                    i = i-1
                    
                    message("\n\n\nFactor number(s): ", paste0(i))
                    
                    # calculating factor solution
                    message('\nCalculating factor solution...\n')
                    if(estimator == "MLR" | estimator == "MLM" | estimator == "MLMV" | estimator == "MLMVS"){
                      unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method)        
                    } else {
                      
                      #unrotated <- efaUnrotate(fa_covdata, nf=i, estimator="WLSMV", ordered=names(fa_covdata), verbose=TRUE)
                      unrotated <- efaUnrotate(fa_covdata, nf=i, estimator=factoring_method, missing='default', ordered=names(fa_covdata), start=FALSE, verbose=T, ...)
                      
                    }
                    message('Done...\n')
                    
                    
                    
                    if(i==1){
                      summary(unrotated, std=TRUE)
                      
                      if(estimator == 'pml' | estimator == 'PML'){
                        Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                        Model_Fit <- round(Model_Fit, 4)
#                         Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                        message('Done...\n')
                      } else {
                        Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                        Model_Fit <- round(Model_Fit, 4)
                        Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                        message('Done...\n')
                      }
                      
                      message("\n\nResults of Model Fit Indices")
                      colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
                      print(Model_Fit)      
                      
                      message("\n\n Use this result. Thanks a lot!")
                      print(reliability(unrotated))
                      return(unrotated)
                    } else {
                      #unrotated <- unrotated_old
                      
                      # estimate R-square (h2)
                      message('\nCalculating Rotated factor solution...\n')
                      rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
                      message('Done...\n')
                      
                      message('\nCalculating Communality...\n')
                      if(sum(is.na(data.frame(inspect(unrotated, "rsquare")))==1)==0){
                        R_SQUARE <- data.frame(inspect(unrotated, "rsquare"))
                        COMMUNALITY_METHOD <- estimator
                      } else {
                        # SMC
                        message('\nCalculating SMC for h2...\n')
                        loading <- data.frame(rotated@loading)
                        R_SQUARE <- rowSums((loading^2))
                        COMMUNALITY_METHOD <- 'SMC'
                      }
                      message('Done...\n')
                      
                      message("\nLoadings\n")
                      #rotated <- oblqRotate(unrotated, method="geomin", maxit=999999)
                      h2 <- R_SQUARE
                      loading <- data.frame(rotated@loading)
                      loading$h2 <- data.frame(h2)
                      print(loading)
                      
                      message("\nCommunality\n")
                      print(COMMUNALITY_METHOD)
                      
                      message("\nEigenvalues\n")
                      Eigenvalues <- abs(colSums(attributes(rotated)$loading^2))
                      print(Eigenvalues)
                      
                      
                      message("\nVariance Explained\n")
                      Explained <- Eigenvalues / sum(h2)
                      print(Explained)
                      
                      message("\nTotal Variance Explained\n")
                      print(sum(Explained))
                      
                      message("\nSummary results\n")
                      print(summary(rotated))
                      
                      
                      if(estimator == 'pml' | estimator == 'PML'){
                        Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue")))
                        Model_Fit <- round(Model_Fit, 4)
#                         Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "wrmr", "chisq", "df", "pvalue"))))] <- .98
                        message('Done...\n')
                      } else {
                        Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled")))
                        Model_Fit <- round(Model_Fit, 4)
                        Model_Fit[is.na(Model_Fit <- data.frame(fitmeasures(unrotated, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "wrmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))))] <- .98
                        message('Done...\n')
                      }
                      
                      message("\n\n Model Fit Indices")
                      colnames(Model_Fit) <- c('Coefficients'); rownames(Model_Fit) <- c('CFI', 'TLI', 'RMSEA', 'RMSEA Lower 5%', 'RMSEA Upper 95%', 'RMSEA < .05 p-value', 'WRMR', "chisq.scaled", "df.scaled", "pvalue.scaled")
                      print(Model_Fit)      
                      
                      
                      message("\n\n Use this result. Thanks a lot!")
                      return(rotated)
                    }
                    
                    message("\n\n Use this result. Thanks a lot!")
                    return(rotated)
                  }
                  
                }
                
                message("\n\n Use this result. Thanks a lot!")
                return(rotated)
              }
            }
          }
        } else {
          # for i = 1 message
          message("\n\n One factor solution wasn't Good model. I'll extract more factors. Current factor number(s): ", paste0(i))
        }
        
      }
      
    }
  }
  
  
  # MIRT (Full-Information Factor Analysis)
  else if(estimator == "mirt" | estimator == "MIRT" | estimator == "fifa") {
    
    # Ignore RMSEA evaluation if small items where evaluation of model fits
    if(ncol(fa_covdata) < 13){
      RMSEA_Criterion <- .05
      #fa_covdata <- k.rescale(fa_covdata, 2)
    } else {
      RMSEA_Criterion <- .05
    }
    
    rq <- 0
    
    for(i in request_factors:100) {
      message("\n\n\n Full-Information Item Factor analysis \n")
      #message("Always need to select model : model = c(2PL, 3PL, 4PL, graded)")
      message("\n\nFactor number(s): ", paste0(i))
      message("\nEstimating...\n");
      
      rq <- rq + 1
      
      if(!exists(as.character(substitute(irtmodel)))) {
        itemtype <- irtmodel
        message('Item type: ', paste0(itemtype))
        
        message('\n')
      } else {
        itemtype <- NULL
      }
      
      if(rq > 1) {
        result_old <- result
        remove(result)
      }
      
      if(rq > 1) { # for treating Rasch model
        if(length(itemtype)==0){ # for itemtype == NULL
          
        }
        else if(itemtype == 'Rasch' | irtmodel == 'Rasch'){ # Rasch model / Multidimensional Rasch model can be estimate on CFA settings
          message('The Rasch model is Unidimensional Model. I can not extract more factors. Returning one factor solution.')
          return(result_old)
        } else { # any other models
          
        }
      }
      
      
      result <- mirt::mirt(fa_covdata, i, TOL = 1e-3, accelerate = 'squarem', SE = T, itemtype=itemtype, calcNull=T, method = "MHRM", rotate = 'geominQ', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999), control = list(maxit=100000))
      #
      #message(" \n ")
      #message("\n\n IRT coefficients")
      #print(coef(result))
      #print(summary(result))
      
      # Print Item Trace Line
      if(i == 1){
        #print(plot(result, type = 'infotrace'))
        try(print(plot(result, type = 'infoSE')))
        try(print(plot(result, type = 'infotrace', facet_items = TRUE)))
        try(print(plot(result, type = 'trace')))
      }
      
      message("\n\n Model Fit Indices")
      #message(" CFI"); print(getClass(result)@CFI)
      #message(" TLI"); print(getClass(result)@TLI)
      #message(" RMSEA"); print(getClass(result)@RMSEA)
      
      # calculate and evaluation model via AIC, BIC, AICc, SABIC
      if(rq > 1) {
        anova_fit <- anova(result_old, result)
        print(anova_fit)
        message("\n\n")
#         if((anova_fit[1,1] < anova_fit[2,1]) + (anova_fit[1,2] < anova_fit[2,2]) + (anova_fit[1,3] < anova_fit[2,3]) + (anova_fit[1,4] < anova_fit[2,4]) + (anova_fit[1,5] > anova_fit[2,5]) >= 4 | (anova_fit[2,8] > .5)){
          if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
          # Breaking for Over factor specification
          message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
          
          return(result_old)
        }
        mirtfit_old <- mirtfit
        try(remove(mirtfit), silent = T)
      }
      
      
      message('Estimating Theta with Multiple Imputation... Please be patient.')
      try(Theta <- fscores(result, full.scores = TRUE, full.scores.SE = F, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP'), silent = T)
      if(exists("Theta")){ # patched 2015/07/21
        
      } else {
        return(result_old)
      }
      
      message('Estimating Model Fit Statistics with Multiple Imputation... Please be patient.')
      try(mirtfit <- findM2(result, Theta=Theta, impute=100, calcNull=TRUE), silent = T)
      
      if(exists("mirtfit")){ # patched 2014/12/05
        try(print(round(mirtfit[1,], digits=4)), silent = T)
        
        #patched 2014/12/2 -- preventing overfatoring
        if(exists("mirtfit_old")){
          if(mirtfit$TLI[1] >= .9 | mirtfit$CFI[1] >= .9) { # i > 1? i > 2? this is heuristic value
            
            if(mirtfit$p[1] >= .5 | mirtfit_old$TLI[1] >= mirtfit$TLI[1] | mirtfit_old$CFI[1] >= mirtfit$CFI[1] | abs(mirtfit$RMSEA[1] - mirtfit_old$RMSEA[1]) <= .001 | (mirtfit$TLI[1] - mirtfit_old$TLI[1]) < .001 | (mirtfit$CFI[1] - mirtfit_old$CFI[1]) < .001) {
              message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
              return(result_old)
            }
          }
          
        }
        
        
        
        ## factor model evaluation ##
        # heywood case check
        if(i > 6){
          if(i==1){
            heywood <- abs(getClass(result)@Fit$F)
            heywood2 <- abs(getClass(result)@Fit$F)
            heywood3 <- abs(getClass(result)@Fit$F)
          } else {
            try(require(GPArotation), silent = T)
            message('Checking Heywood cases by geominQ rotation')
            heywood <- abs(geominQ(getClass(result)@Fit$F, maxit=999999)$loadings)
            message('Checking Heywood cases by bentlerQ rotation')
            heywood2 <- abs(bentlerQ(getClass(result)@Fit$F, maxit=999999)$loadings)
            message('Checking Heywood cases by infomaxQ rotation')
            heywood3 <- abs(infomaxQ(getClass(result)@Fit$F, maxit=999999)$loadings)
          }
          
          communality <- round(result@Fit$h2, 3)
          
          if(sum(communality >= 1) >= 1 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {      
            message("geomin rotation - heywood")
            print(heywood)
            message("bentler rotation - heywood2")
            print(heywood2)
            message("infomax rotation - heywood3")
            print(heywood3)
            message("Communality >= 1"); print(sum(communality >= 1))
            message("\n\n Warning: This model is heywood case(s) occured!\n")
            #           return(result_old)
                       if(i==1) {
                         warning('Solution may not be interprintable')
                       } else {
            result <- result_old # intead of return(result_old) (old method)
                       }
            
            message(" \n ")
            message("\n\n IRT coefficients")
            print(coef(result, rotate = 'geominQ'))
            print(summary(result))
            #print(result, sort=TRUE)
            message("\n\n Use this result. Thanks a lot!")
            return(result)
          }
        }
        if(i <= 6) {
          
          communality <- round(result@Fit$h2, 3)
          
          if(sum(communality >= 1) >= 1) {
            message("\n\n Warning: This model is heywood case(s) occured!\n")
            #           return(result_old)
                       if(i==1) {
                         warning('Solution may not be interprintable')
                       } else {
            result <- result_old # intead of return(result_old) (old method)
                       }
            
            message(" \n ")
            message("\n\n IRT coefficients")
            print(coef(result, rotate = 'geominQ'))
            print(summary(result))
            #print(result, sort=TRUE)
            message("\n\n Use this result. Thanks a lot!")
            return(result)
          }
        }

        
        # Model fit check
        
        # treating abnormal model fit
        if(mirtfit$p[1] >= .5 | mirtfit$TLI[1] >= 1 | mirtfit$CFI[1] >= 1 | mirtfit$TLI[1] < 0 | mirtfit$CFI[1] < 0 | mirtfit$RMSEA_5[1] < 0) {
          
          message("\n\n Warning: This model is something was wrong!")
          if(i==1) {
            
          } else {
            result <- result_old
          }
          
          message(" \n ")
          message("\n\n IRT coefficients")
          print(coef(result, rotate = 'geominQ'))
          print(summary(result))
          #print(result, sort=TRUE)
          message("\n\n Use this result. Thanks a lot!")
          return(result)
        }
        
        # if reached good model fit, try to extract two more factors
        else if(mirtfit$TLI[1] >= .9 && mirtfit$CFI[1] >= .9 && mirtfit$RMSEA_5[1] < RMSEA_Criterion) {
          #else if(getClass(result)@TLI >= .9 && getClass(result)@CFI >= .9 && getClass(result)@RMSEA < .05) {
          message("\n\n Maximum Factor number archieved: ", paste0(i))
          message("\n\n Statistically Perfect! This solution is final solution. It will be save in your value. But, I'll try extract two more factors. Please check can be explain. \n")
          
          # extract one more factor
          i = i + 1
          message('Trying to extract ', paste0(i), ' factor solution')
          if(!exists(as.character(substitute(irtmodel)))) {
            itemtype <- irtmodel
            message('Item type: ', paste0(itemtype))
            
            message('\n')
          } else {
            itemtype <- NULL
          }
          result_onemore <- mirt::mirt(fa_covdata, i, TOL = 1e-3, accelerate = 'squarem', SE = T, itemtype=itemtype, calcNull=T, method = "MHRM", rotate = 'geominQ', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999), control = list(maxit=100000))
          
          #result <- mirt::mirt(fa_covdata, i, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
          #message(" \n ")
          #message("\n\n IRT coefficients")
          #print(coef(result))
          #print(summary(result_onemore))
          
          
          message("\n\n Model Fit Indices")
          #message(" CFI"); print(getClass(result)@CFI)
          #message(" TLI"); print(getClass(result)@TLI)
          #message(" RMSEA"); print(getClass(result)@RMSEA)
          anova_fit <- anova(result, result_onemore)
          print(anova_fit)
          if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
            # Breaking for Over factor specification -- strict
            message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
            
            message("\n\n")
            return(result)
          }
          message('Estimating Theta with Multiple Imputation... Please be patient.')
          
          Theta_onemore <- fscores(result_onemore, full.scores = TRUE, full.scores.SE = F, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
          message("\n\n")
          
          message('Estimating Model Fit Statistics with Multiple Imputation... Please be patient.')
          
          try(mirtfit_onemore <- findM2(result_onemore, Theta=Theta_onemore, impute=100, calcNull=TRUE), silent = T)
          
          #
          if(exists("mirtfit_onemore")){ # patched 2014/12/05
            try(print(round(mirtfit_onemore[1,], digits=4)), silent = T)
            
            #patched 2014/12/2
            #if(i > 1) {
            
            #  if(mirtfit$TLI[1] >= mirtfit_onemore$TLI[1] | mirtfit$CFI[1] >= mirtfit_onemore$CFI[1] | (mirtfit_onemore$TLI[1] - mirtfit$TLI[1]) <= .001 | (mirtfit$CFI_onemore[1] - mirtfit$CFI[1]) <= .001) {
            #    message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
            #    return(result)
            #  }
            #} 
          } else {
            return(result)
          }
          #
          
          if(i==1){
            heywood <- abs(getClass(result_onemore)@Fit$F)
          } else {
            try(require(GPArotation), silent = T)
            message('Checking Heywood cases by geominQ rotation')
            heywood <- abs(geominQ(getClass(result_onemore)@Fit$F, maxit=999999)$loadings)
            message('Checking Heywood cases by bentlerQ rotation')
            heywood2 <- abs(bentlerQ(getClass(result_onemore)@Fit$F, maxit=999999)$loadings)
            message('Checking Heywood cases by infomaxQ rotation')
            heywood3 <- abs(infomaxQ(getClass(result_onemore)@Fit$F, maxit=999999)$loadings)
          }
          
          communality <- round(result@Fit$h2, 3)
          if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
            #if(sum(communality >= .99) > 0) {
            message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
            #rm(result_onemore, Theta_onemore, mirtfit_onemore)
            return(result)
          }
          
          # Model fit check
          else if(mirtfit_onemore$p[1] >= .5 | mirtfit_onemore$TLI[1] >= 1 | mirtfit_onemore$CFI[1] >= 1 | mirtfit$TLI[1] >= mirtfit_onemore$TLI[1] | mirtfit$CFI[1] >= mirtfit_onemore$CFI[1] | mirtfit_onemore$RMSEA_5[1] < 0  | abs(mirtfit_onemore$RMSEA[1] - mirtfit$RMSEA[1]) <= .001 | (mirtfit_onemore$TLI[1] - mirtfit$TLI[1]) <= .01 | (mirtfit_onemore$CFI[1] - mirtfit$CFI[1]) <= .01) {
            message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
            #rm(result_onemore, Theta_onemore, mirtfit_onemore)
            return(result)
          }
          
          # extract two more factor
          i = i + 1
          message('Trying to extract ', paste0(i), ' factor solution')
          if(!exists(as.character(substitute(irtmodel)))) {
            itemtype <- irtmodel
            message('Item type: ', paste0(itemtype))
            message('\n')
          } else {
            itemtype <- NULL
          }
          result_twomore <- mirt::mirt(fa_covdata, i, TOL = 1e-3, accelerate = 'squarem', SE = T, itemtype=itemtype, calcNull=T, method = "MHRM", rotate = 'geominQ', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999), control = list(maxit=100000))
          #result <- mirt::mirt(fa_covdata, i, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
          #message(" \n ")
          #message("\n\n IRT coefficients")
          #print(coef(result))
          #print(summary(result_twomore))
          
          
          message("\n\n Model Fit Indices")
          #message(" CFI"); print(getClass(result)@CFI)
          #message(" TLI"); print(getClass(result)@TLI)
          #message(" RMSEA"); print(getClass(result)@RMSEA)
          anova_fit <- anova(result_onemore, result_twomore)
          print(anova_fit)
          if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
            # Breaking for Over factor specification -- strict
            message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
            #return(result)
            return(result_onemore)
          }
          message("\n\n")
          message('Estimating Theta with Multiple Imputation... Please be patient.')
          
          Theta_twomore <- fscores(result_twomore, full.scores = TRUE, full.scores.SE = F, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
          
          message('Estimating Model Fit Statistics with Multiple Imputation... Please be patient.')
          
          try(mirtfit_twomore <- findM2(result_twomore, Theta=Theta_twomore, impute=100, calcNull=TRUE), silent = T)
          
          #
          if(exists("mirtfit_twomore")){ # patched 2014/12/05
            try(print(round(mirtfit_twomore[1,], digits=4)), silent = T)
            
            #patched 2014/12/2
            #if(i > 1) {
            
            #  if(mirtfit_onemore$TLI[1] >= mirtfit_twomore$TLI[1] | mirtfit_onemore$CFI[1] >= mirtfit_twomore$CFI[1] | (mirtfit_twomore$TLI[1] - mirtfit_onemore$TLI[1]) <= .001 | (mirtfit$CFI_twomore[1] - mirtfit$CFI[1]) <= .001) {
            #    message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
            #    return(result_onemore)
            #  }
            #} 
          } else {
            return(result_onemore)
          }
          #
          
          if(i==1){
            heywood <- abs(getClass(result_twomore)@Fit$F)
          } else {
            try(require(GPArotation), silent = T)
            message('Checking Heywood cases by geominQ rotation')
            heywood <- abs(geominQ(getClass(result_twomore)@Fit$F, maxit=999999)$loadings)
            message('Checking Heywood cases by bentlerQ rotation')
            heywood2 <- abs(bentlerQ(getClass(result_twomore)@Fit$F, maxit=999999)$loadings)
            message('Checking Heywood cases by infomaxQ rotation')
            heywood3 <- abs(infomaxQ(getClass(result_twomore)@Fit$F, maxit=999999)$loadings)
          }
          
          communality <- round(result@Fit$h2, 3)
          if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
            #if(sum(communality >= .99) > 0) {
            message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
            #rm(result_twomore, Theta_twomore, mirtfit_twomore)
            return(result_onemore)
          }
          
          # Model fit check
          else if(mirtfit_twomore$p[1] >= .5 | mirtfit_twomore$TLI[1] >= 1 | mirtfit_twomore$CFI[1] >= 1 | mirtfit_onemore$TLI[1] >= mirtfit_twomore$TLI[1] | mirtfit_onemore$CFI[1] >= mirtfit_twomore$CFI[1] | mirtfit_twomore$RMSEA_5[1] < 0 | abs(mirtfit_twomore$RMSEA[1] - mirtfit_onemore$RMSEA[1]) <= .001 | (mirtfit_twomore$TLI[1] - mirtfit_onemore$TLI[1]) <= .01 | (mirtfit_twomore$CFI[1] - mirtfit_onemore$CFI[1]) <= .01) {
            message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
            #rm(result_twomore, Theta_twomore, mirtfit_twomore)
            return(result_onemore)
          }
          
          # extract three more factor
          i = i + 1
          message('Trying to extract ', paste0(i), ' factor solution')
          if(!exists(as.character(substitute(irtmodel)))) {
            itemtype <- irtmodel
            message('Item type: ', paste0(itemtype))
            
            message('\n')
          } else {
            itemtype <- NULL
          }
          result_threemore <- mirt::mirt(fa_covdata, i, TOL = 1e-3, accelerate = 'squarem', SE = T, itemtype=itemtype, calcNull=T, method = "MHRM", rotate = 'geominQ', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999), control = list(maxit=100000))
          #result <- mirt::mirt(fa_covdata, i, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
          #message(" \n ")
          #message("\n\n IRT coefficients")
          #print(coef(result))
          #print(summary(result_threemore))
          
          
          message("\n\n Model Fit Indices")
          #message(" CFI"); print(getClass(result)@CFI)
          #message(" TLI"); print(getClass(result)@TLI)
          #message(" RMSEA"); print(getClass(result)@RMSEA)
          anova_fit <- anova(result_twomore, result_threemore)
          print(anova_fit)
          if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
            # Breaking for Over factor specification -- strict
            message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
            #return(result)
            return(result_twomore)
          }
          message("\n\n")
          
          message('Estimating Theta with Multiple Imputation... Please be patient.')
          
          Theta_threemore <- fscores(result_threemore, full.scores = TRUE, full.scores.SE = F, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
          
          message('Estimating Model Fit Statistics with Multiple Imputation... Please be patient.')
          
          mirtfit_threemore <- findM2(result_threemore, Theta=Theta_threemore, impute=100, calcNull=TRUE)
          print(round(mirtfit_threemore[1,], 4))
          
          if(i==1){
            heywood <- abs(getClass(result_threemore)@Fit$F)
          } else {
            try(require(GPArotation), silent = T)
            message('Checking Heywood cases by geominQ rotation')
            heywood <- abs(geominQ(getClass(result_threemore)@Fit$F, maxit=999999)$loadings)
            message('Checking Heywood cases by bentlerQ rotation')
            heywood2 <- abs(bentlerQ(getClass(result_threemore)@Fit$F, maxit=999999)$loadings)
            message('Checking Heywood cases by infomaxQ rotation')
            heywood3 <- abs(infomaxQ(getClass(result_threemore)@Fit$F, maxit=999999)$loadings)
          }
          
          communality <- round(result@Fit$h2, 3)
          if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
            #if(sum(communality >= .99) > 0) {
            message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
            #rm(result_threemore, Theta_threemore, mirtfit_threemore)
            return(result_twomore)
          }
          
          # Model fit check
          else if(mirtfit_threemore$p[1] >= .5 | mirtfit_threemore$TLI[1] >= 1 | mirtfit_threemore$CFI[1] >= 1 | mirtfit_twomore$TLI[1] >= mirtfit_threemore$TLI[1] | mirtfit_twomore$CFI[1] >= mirtfit_threemore$CFI[1] | mirtfit_threemore$RMSEA_5[1] < 0  | abs(mirtfit_threemore$RMSEA[1] - mirtfit_twomore$RMSEA[1]) <= .001 | (mirtfit_threemore$TLI[1] - mirtfit_twomore$TLI[1]) <= .01 | (mirtfit_threemore$CFI[1] - mirtfit_twomore$CFI[1]) <= .01) {
            message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
            #rm(result_threemore, Theta_threemore, mirtfit_threemore)
            return(result_twomore)
          }
          
          
          # extract four more factor
          i = i + 1
          message('Trying to extract ', paste0(i), ' factor solution')
          if(!exists(as.character(substitute(irtmodel)))) {
            itemtype <- irtmodel
            message('Item type: ', paste0(itemtype))
            
            message('\n')
          } else {
            itemtype <- NULL
          }
          result_fourmore <- mirt::mirt(fa_covdata, i, TOL = 1e-3, accelerate = 'squarem', SE = T, itemtype=itemtype, calcNull=T, method = "MHRM", rotate = 'geominQ', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999), control = list(maxit=100000))
          #result <- mirt::mirt(fa_covdata, i, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
          #message(" \n ")
          #message("\n\n IRT coefficients")
          #print(coef(result))
          #print(summary(result_fourmore))
          
          
          message("\n\n Model Fit Indices")
          #message(" CFI"); print(getClass(result)@CFI)
          #message(" TLI"); print(getClass(result)@TLI)
          #message(" RMSEA"); print(getClass(result)@RMSEA)
          anova_fit <- anova(result_threemore, result_fourmore)
          print(anova_fit)
          if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
            # Breaking for Over factor specification -- strict
            message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
            #return(result)
            return(result_threemore)
          }
          message("\n\n")
          
          message('Estimating Theta with Multiple Imputation... Please be patient.')
          
          Theta_fourmore <- fscores(result_fourmore, full.scores = TRUE, full.scores.SE = F, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
          
          message('Estimating Model Fit Statistics with Multiple Imputation... Please be patient.')
          
          mirtfit_fourmore <- findM2(result_fourmore, Theta=Theta_fourmore, impute=100, calcNull=TRUE)
          print(round(mirtfit_fourmore[1,], 4))
          
          if(i==1){
            heywood <- abs(getClass(result_fourmore)@Fit$F)
          } else {
            try(require(GPArotation), silent = T)
            message('Checking Heywood cases by geominQ rotation')
            heywood <- abs(geominQ(getClass(result_fourmore)@Fit$F, maxit=999999)$loadings)
            message('Checking Heywood cases by bentlerQ rotation')
            heywood2 <- abs(bentlerQ(getClass(result_fourmore)@Fit$F, maxit=999999)$loadings)
            message('Checking Heywood cases by infomaxQ rotation')
            heywood3 <- abs(infomaxQ(getClass(result_fourmore)@Fit$F, maxit=999999)$loadings)
          }
          
          communality <- round(result@Fit$h2, 3)
          if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
            #if(sum(communality >= .99) > 0) {
            message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
            #rm(result_fourmore, Theta_fourmore, mirtfit_fourmore)
            return(result_threemore)
          }
          
          # Model fit check
          else if(mirtfit_fourmore$p[1] >= .5 | mirtfit_fourmore$TLI[1] >= 1 | mirtfit_fourmore$CFI[1] >= 1 | mirtfit_threemore$TLI[1] >= mirtfit_fourmore$TLI[1] | mirtfit_threemore$CFI[1] >= mirtfit_fourmore$CFI[1] | mirtfit_fourmore$RMSEA_5[1] < 0 | abs(mirtfit_fourmore$RMSEA[1] - mirtfit_threemore$RMSEA[1]) <= .001 | (mirtfit_fourmore$TLI[1] - mirtfit_threemore$TLI[1]) <= .01 | (mirtfit_fourmore$CFI[1] - mirtfit_threemore$CFI[1]) <= .01) {
            message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
            #rm(result_fourmore, Theta_fourmore, mirtfit_fourmore)
            return(result_threemore)
          }
          
          
          # extract five more factor
          i = i + 1
          message('Trying to extract ', paste0(i), ' factor solution')
          if(!exists(as.character(substitute(irtmodel)))) {
            itemtype <- irtmodel
            message('Item type: ', paste0(itemtype))
            
            message('\n')
          } else {
            itemtype <- NULL
          }
          result_fivemore <- mirt::mirt(fa_covdata, i, TOL = 1e-3, accelerate = 'squarem', SE = T, itemtype=itemtype, calcNull=T, method = "MHRM", rotate = 'geominQ', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999), control = list(maxit=100000))
          #result <- mirt::mirt(fa_covdata, i, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
          #message(" \n ")
          #message("\n\n IRT coefficients")
          #print(coef(result))
          #print(summary(result_fivemore))
          
          
          message("\n\n Model Fit Indices")
          #message(" CFI"); print(getClass(result)@CFI)
          #message(" TLI"); print(getClass(result)@TLI)
          #message(" RMSEA"); print(getClass(result)@RMSEA)
          anova_fit <- anova(result_fourmore, result_fivemore)
          print(anova_fit)
          if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
            # Breaking for Over factor specification -- strict
            message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
            #return(result)
            return(result_fourmore)
          }
          message("\n\n")
          
          message('Estimating Theta with Multiple Imputation... Please be patient.')
          
          Theta_fivemore <- fscores(result_fivemore, full.scores = TRUE, full.scores.SE = F, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
          
          message('Estimating Model Fit Statistics with Multiple Imputation... Please be patient.')
          
          mirtfit_fivemore <- findM2(result_fivemore, Theta=Theta_fivemore, impute=100, calcNull=TRUE)
          print(round(mirtfit_fivemore[1,], 4))
          
          if(i==1){
            heywood <- abs(getClass(result_fivemore)@Fit$F)
          } else {
            try(require(GPArotation), silent = T)
            message('Checking Heywood cases by geominQ rotation')
            heywood <- abs(geominQ(getClass(result_fivemore)@Fit$F, maxit=999999)$loadings)
            message('Checking Heywood cases by bentlerQ rotation')
            heywood2 <- abs(bentlerQ(getClass(result_fivemore)@Fit$F, maxit=999999)$loadings)
            message('Checking Heywood cases by infomaxQ rotation')
            heywood3 <- abs(infomaxQ(getClass(result_fivemore)@Fit$F, maxit=999999)$loadings)
          }
          
          communality <- round(result@Fit$h2, 3)
          if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
            #if(sum(communality >= .99) > 0) {
            message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
            #rm(result_fivemore, Theta_fivemore, mirtfit_fivemore)
            return(result_fourmore)
          }
          
          # Model fit check
          else if(mirtfit_fivemore$p[1] >= .5 | mirtfit_fivemore$TLI[1] >= 1 | mirtfit_fivemore$CFI[1] >= 1 | mirtfit_fourmore$TLI[1] >= mirtfit_fivemore$TLI[1] | mirtfit_fourmore$CFI[1] >= mirtfit_fivemore$CFI[1] | mirtfit_fivemore$RMSEA_5[1] < 0 | abs(mirtfit_fivemore$RMSEA[1] - mirtfit_fourmore$RMSEA[1]) <= .001 | (mirtfit_fivemore$TLI[1] - mirtfit_fourmore$TLI[1]) <= .01 | (mirtfit_fivemore$CFI[1] - mirtfit_fourmore$CFI[1]) <= .01) {
            message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
            #rm(result_fivemore, Theta_fivemore, mirtfit_fivemore)
            return(result_fourmore)
          }
          
          else {
            return(result_fivemore) # if can not detect error
          }
          
          # if failed, return original result
          return(result)
        } 
        
        else if(mirtfit$TLI[1] >= .9 | mirtfit$CFI[1] >= .9 | mirtfit$RMSEA_5[1] < RMSEA_Criterion) {
          #else if(getClass(result)@TLI >= .9 | getClass(result)@CFI >= .9 | getClass(result)@RMSEA < .08) {
          message("\n\n Maximum Factor number archieved: ", paste0(i))
          message("\nYou can use this result when you get poor factor solution where next estimation.")
        } else {
          message("\n\n I'll extract more factors. Current factor number(s): ", paste0(i))}
        
        #return(mirtfit)
        #print(residuals(result))
        #par(mfrow=c(2,2))
        message(" \n ") 
        
        
        
        
        message("\n\n\n\n\n")
        
      } else {
        if(i == 1){
          message('Quitting...')
          return(result)
        } else {
          message('Quitting...')
          return(result_old)
        }
      }
      
      
    }
  }
  
  # MIRT (Full-Information Factor Analysis) -- MultipleGroup
  else if(estimator == "multigroup" | estimator == "multiplegroup") {
    
    # Ignore RMSEA evaluation if small items where evaluation of model fits
    if(ncol(fa_covdata) < 10){
      RMSEA_Criterion <- .05
      #fa_covdata <- k.rescale(fa_covdata, 2)
    } else {
      RMSEA_Criterion <- .05
    }
    
    rq <- 0
    
    for(i in request_factors:100) {
      message("\n\n\n Full-Information Item Factor analysis \n")
      #message("Always need to select model : model = c(2PL, 3PL, 4PL, graded)")
      message("\n\nFactor number(s): ", paste0(i))
      message("\nEstimating...\n");
      
      rq <- rq + 1
      
      if(rq > 1) {
        result_old <- result
        remove(result)
      }
      
      if(!exists(as.character(substitute(irtmodel)))) {
        itemtype <- irtmodel
        message('Item type: ', paste0(itemtype))
        
        message('\n')
      } else {
        itemtype <- NULL
      }
      
      if(i == 1) {
        result <-  multipleGroup(fa_covdata, i, group=group, rotate = 'geominQ', itemtype=irtmodel, method = "MHRM", invariance = c('free_means', 'free_var', 'free_cov', 'slopes', 'intercepts'), technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
        
      } else {
        result <-  multipleGroup(fa_covdata, i, group=group, rotate = 'geominQ', itemtype=irtmodel, method = "MHRM", invariance = c('free_means', 'free_var', 'free_cov', 'slopes', 'intercepts'), technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
        
      }
      #result <-  multipleGroup(fa_covdata, i, group=group, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
      #message(" \n ")
      #message("\n\n IRT coefficients")
      #print(coef(result))
      #print(summary(result))
      
      # Print Item Trace Line
      if(i == 1){
        try(print(plot(result)), silent = T)
        try(print(plot(result, type = 'SE')), silent = T)
        try(print(plot(result, type = 'RE')), silent = T)
        try(print(plot(result, type = 'score')), silent = T)
        try(print(plot(result, type = 'trace')), silent = T)}
      
      message("\n\n Model Fit Indices")
      #message(" CFI"); print(getClass(result)@CFI)
      #message(" TLI"); print(getClass(result)@TLI)
      #message(" RMSEA"); print(getClass(result)@RMSEA)
      
      # calculate and evaluation model via AIC, BIC, AICc, SABIC
      if(rq > 1) {
        anova_fit <- anova(result_old, result)
        print(anova_fit)
        message("\n\n")
        if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
          # Breaking for Over factor specification
          message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
          
          return(result_old)
        }
      }
      Theta <- fscores(result, full.scores=TRUE, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
      
      message('Estimating Model Fit Statistics with Multiple Imputation... Please be patient.')
      
      mirtfit <- findM2(result, Theta=Theta, impute=100, calcNull=TRUE)
      print(round(mirtfit[1,], digits=4))
      
      
      
      
      ## factor model evaluation ##
      
      
      # Model fit check
      # treating abnormal model fit
      if(mirtfit$TLI[1] >= 1 | mirtfit$CFI[1] >= 1 | mirtfit$TLI[1] < 0 | mirtfit$CFI[1] < 0 | mirtfit$RMSEA_5[1] < 0) {
        
        message("\n\n Warning: This model is something was wrong!")
        if(i==1) {
          
        } else {
          i = i-1
        }
        if(i == 1) {
          result <- multipleGroup(fa_covdata, i, group=group, rotate = 'geominQ', itemtype=irtmodel, method = "MHRM", invariance = c('free_means', 'free_var', 'free_cov', 'slopes', 'intercepts'), technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
          
        } else {
          result <- multipleGroup(fa_covdata, i, group=group, rotate = 'geominQ', itemtype=irtmodel, method = "MHRM", invariance = c('free_means', 'free_var', 'free_cov', 'slopes', 'intercepts'), technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
          
        }
        message(" \n ")
        message("\n\n IRT coefficients")
        print(coef(result, rotate = 'geominQ'))
        print(summary(result))
        #print(result, sort=TRUE)
        message("\n\n Use this result. Thanks a lot!")
        return(result)
      }
      
      # if reached good model fit, try to extract two more factors
      else if(mirtfit$TLI[1] >= .9 && mirtfit$CFI[1] >= .9 && mirtfit$RMSEA_5[1] < RMSEA_Criterion) {
        #else if(getClass(result)@TLI >= .9 && getClass(result)@CFI >= .9 && getClass(result)@RMSEA < .05) {
        message("\n\n Maximum Factor number archieved: ", paste0(i))
        message("\n\n Statistically Perfect! This solution is final solution. It will be save in your value. But, I'll try extract two more factors. Please check can be explain. \n")
        
        # extract one more factor
        i = i + 1
        message('Trying to extract ', paste0(i), ' factor solution')
        result_onemore <- multipleGroup(fa_covdata, i, group=group, rotate = 'geominQ', method = "MHRM", invariance = c('free_means', 'free_var', 'free_cov', 'slopes', 'intercepts'), technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
        #result <-  multipleGroup(fa_covdata, i, group=group, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
        #message(" \n ")
        #message("\n\n IRT coefficients")
        #print(coef(result))
        #print(summary(result_onemore))
        
        
        message("\n\n Model Fit Indices")
        #message(" CFI"); print(getClass(result)@CFI)
        #message(" TLI"); print(getClass(result)@TLI)
        #message(" RMSEA"); print(getClass(result)@RMSEA)
        anova_fit <- anova(result, result_onemore)
        print(anova_fit)
        if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
          # Breaking for Over factor specification -- strict
          message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
          
          message("\n\n")
          return(result)
        }
        Theta_onemore <- fscores(result_onemore, full.scores=TRUE, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
        message("\n\n")
        mirtfit_onemore <- findM2(result_onemore, Theta=Theta_onemore, impute=100, calcNull=TRUE)
        print(round(mirtfit_onemore[1,], 4))
        
        
        if(i==1){
          heywood <- abs(getClass(result_onemore)@Fit$F)
        } else {
          try(require(GPArotation), silent = T)
          message('Checking Heywood cases by geominQ rotation')
          heywood <- abs(geominQ(getClass(result_onemore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by bentlerQ rotation')
          heywood2 <- abs(bentlerQ(getClass(result_onemore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by infomaxQ rotation')
          heywood3 <- abs(infomaxQ(getClass(result_onemore)@Fit$F, maxit=999999)$loadings)
        }
        
        communality <- round(result@Fit$h2, 3)
        if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
          #if(sum(communality >= .99) > 0) {
          message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
          #rm(result_onemore, Theta_onemore, mirtfit_onemore)
          return(result)
        }
        
        # Model fit check
        else if(mirtfit_onemore$p[1] >= .5 | mirtfit_onemore$TLI[1] >= 1 | mirtfit_onemore$CFI[1] >= 1 | mirtfit$TLI[1] >= mirtfit_onemore$TLI[1] | mirtfit$CFI[1] >= mirtfit_onemore$CFI[1] | mirtfit_onemore$RMSEA_5[1] < 0) {
          message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
          #rm(result_onemore, Theta_onemore, mirtfit_onemore)
          return(result)
        }
        
        # extract two more factor
        i = i + 1
        message('Trying to extract ', paste0(i), ' factor solution')
        result_twomore <- multipleGroup(fa_covdata, i, group=group, rotate = 'geominQ', itemtype=irtmodel, method = "MHRM", invariance = c('free_means', 'free_var', 'free_cov', 'slopes', 'intercepts'), technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
        #result <-  multipleGroup(fa_covdata, i, group=group, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
        #message(" \n ")
        #message("\n\n IRT coefficients")
        #print(coef(result))
        #print(summary(result_twomore))
        
        
        message("\n\n Model Fit Indices")
        #message(" CFI"); print(getClass(result)@CFI)
        #message(" TLI"); print(getClass(result)@TLI)
        #message(" RMSEA"); print(getClass(result)@RMSEA)
        anova_fit <- anova(result_onemore, result_twomore)
        print(anova_fit)
        if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
          # Breaking for Over factor specification -- strict
          message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
          #return(result)
          return(result_onemore)
        }
        message("\n\n")
        Theta_twomore <- fscores(result_twomore, full.scores=TRUE, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
        mirtfit_twomore <- findM2(result_twomore, Theta=Theta_twomore, impute=100, calcNull=TRUE)
        print(round(mirtfit_twomore[1,], 4))
        
        if(i==1){
          heywood <- abs(getClass(result_twomore)@Fit$F)
        } else {
          try(require(GPArotation), silent = T)
          message('Checking Heywood cases by geominQ rotation')
          heywood <- abs(geominQ(getClass(result_twomore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by bentlerQ rotation')
          heywood2 <- abs(bentlerQ(getClass(result_twomore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by infomaxQ rotation')
          heywood3 <- abs(infomaxQ(getClass(result_twomore)@Fit$F, maxit=999999)$loadings)
        }
        
        communality <- round(result@Fit$h2, 3)
        if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
          #if(sum(communality >= .99) > 0) {
          message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
          #rm(result_twomore, Theta_twomore, mirtfit_twomore)
          return(result_onemore)
        }
        
        # Model fit check
        else if(mirtfit_twomore$p[1] >= .5 | mirtfit_twomore$TLI[1] >= 1 | mirtfit_twomore$CFI[1] >= 1 | mirtfit_onemore$TLI[1] >= mirtfit_twomore$TLI[1] | mirtfit_onemore$CFI[1] >= mirtfit_twomore$CFI[1] | mirtfit_twomore$RMSEA_5[1] < 0) {
          message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
          #rm(result_twomore, Theta_twomore, mirtfit_twomore)
          return(result_onemore)
        }
        else {
          return(result_twomore) # if can not detect error
        }
        
        # if failed, return original result
        return(result)
      } 
      
      else if(mirtfit$TLI[1] >= .9 | mirtfit$CFI[1] >= .9 | mirtfit$RMSEA_5[1] < RMSEA_Criterion) {
        #else if(getClass(result)@TLI >= .9 | getClass(result)@CFI >= .9 | getClass(result)@RMSEA < .08) {
        message("\n\n Maximum Factor number archieved: ", paste0(i))
        message("\nYou can use this result when you get poor factor solution where next estimation.")
      } else {
        message("\n\n I'll extract more factors. Current factor number(s): ", paste0(i))}
      
      #return(mirtfit)
      #print(residuals(result))
      #par(mfrow=c(2,2))
      message(" \n ") 
      
      
      
      
      message("\n\n\n\n\n")
      
    }
  }
  
  # MIRT (Full-Information Factor Analysis using Multilevel Theta information)
  else if(estimator == "multilevel") {
    
    # Ignore RMSEA evaluation if small items where evaluation of model fits
    if(ncol(fa_covdata) < 13){
      RMSEA_Criterion <- .05
      #fa_covdata <- k.rescale(fa_covdata, 2)
    } else {
      RMSEA_Criterion <- .05
    }
    
    rq <- 0
    
    for(i in request_factors:100) {
      message("\n\n\n Full-Information Item Factor analysis using Multilevel Theta \n")
      #message("Always need to select model : model = c(2PL, 3PL, 4PL, graded)")
      message("\n\nFactor number(s): ", paste0(i))
      message("\nEstimating...\n");
      
      rq <- rq + 1
      
      if(rq > 1) {
        result_old <- result
        remove(result)
      }
      
      if(!exists(as.character(substitute(irtmodel)))) {
        itemtype <- irtmodel
        message('Item type: ', paste0(itemtype))
        
        message('\n')
      } else {
        itemtype <- NULL
      }
      
      # estimating random theta
      covdat <- data.frame(group = group)
      
      ranmod <- mixedmirt(fa_covdata, covdat, model=i, random = list(~ 1|group, ~ 1|items), itemtype = itemtype, technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
      summary(ranmod)
      raneff <- randef(ranmod)
      randTheta <- raneff$Theta
      
      
      result <- mirt::mirt(fa_covdata, i, itemtype=itemtype, calcNull=T, rotate = 'geominQ', method = 'QMCEM', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999, customTheta = randTheta), maxit=100000)
      #
      #message(" \n ")
      #message("\n\n IRT coefficients")
      #print(coef(result))
      #print(summary(result))
      
      # Print Item Trace Line
      if(i == 1){
        #print(plot(result, type = 'infotrace'))
        try(print(plot(result, type = 'infoSE')))
        try(print(plot(result, type = 'infotrace', facet_items = TRUE)))
        try(print(plot(result, type = 'trace')))}
      
      message("\n\n Model Fit Indices")
      #message(" CFI"); print(getClass(result)@CFI)
      #message(" TLI"); print(getClass(result)@TLI)
      #message(" RMSEA"); print(getClass(result)@RMSEA)
      
      # calculate and evaluation model via AIC, BIC, AICc, SABIC
      if(rq > 1) {
        anova_fit <- anova(result_old, result)
        print(anova_fit)
        message("\n\n")
        if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
          # Breaking for Over factor specification
          message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
          
          return(result_old)
        }
      }
      message('Estimating Theta with Multiple Imputation... Please be patient.')
      Theta <- fscores(result, full.scores=TRUE, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
      message('Estimating Model Fit Statistics with Multiple Imputation... Please be patient.')
      mirtfit <- findM2(result, Theta=Theta, impute=100, calcNull=TRUE)
      print(round(mirtfit[1,], digits=4))
      
      
      
      
      ## factor model evaluation ##
      # heywood case check
      
      if(i==1){
        heywood <- abs(getClass(result)@Fit$F)
        heywood2 <- abs(getClass(result)@Fit$F)
        heywood3 <- abs(getClass(result)@Fit$F)
      } else {
        try(require(GPArotation), silent = T)
        message('Checking Heywood cases by geominQ rotation')
        heywood <- abs(geominQ(getClass(result)@Fit$F, maxit=999999)$loadings)
        message('Checking Heywood cases by bentlerQ rotation')
        heywood2 <- abs(bentlerQ(getClass(result)@Fit$F, maxit=999999)$loadings)
        message('Checking Heywood cases by infomaxQ rotation')
        heywood3 <- abs(infomaxQ(getClass(result)@Fit$F, maxit=999999)$loadings)
      }
      
      communality <- round(result@Fit$h2, 3)
      
      if(sum(communality >= 1) >= 1 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {      
        message("geomin rotation - heywood")
        print(heywood)
        message("bentler rotation - heywood2")
        print(heywood2)
        message("infomax rotation - heywood3")
        print(heywood3)
        message("Communality >= 1"); print(sum(communality >= 1))
        message("\n\n Warning: This model is heywood case(s) occured! I'm fixing now!\n")
        if(i==1) {
          
        } else {
#           i = i-1
        }
        
        if(!exists(as.character(substitute(irtmodel)))) {
          itemtype <- irtmodel
          message('Item type: ', paste0(itemtype))
          
          message('\n')
        } else {
          itemtype <- NULL
        }
        
        # estimating random theta
        covdat <- data.frame(group = group)
        
        ranmod <- mixedmirt(fa_covdata, covdat, model=i, random = list(~ 1|group, ~ 1|items), itemtype = itemtype, technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
        summary(ranmod)
        raneff <- randef(ranmod)
        randTheta <- raneff$Theta
        
        result <- mirt::mirt(fa_covdata, i, itemtype=itemtype, calcNull=T, rotate = 'geominQ', method = 'QMCEM', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999, customTheta = randTheta), maxit=100000)
        
        message(" \n ")
        message("\n\n IRT coefficients")
        print(coef(result))
        print(summary(result))
        #print(result, sort=TRUE)
        message("\n\n Use this result. Thanks a lot!")
        return(result)
      }
      
      # Model fit check
      # treating abnormal model fit
      else if(mirtfit$TLI[1] >= 1 | mirtfit$CFI[1] >= 1 | mirtfit$TLI[1] < 0 | mirtfit$CFI[1] < 0 | mirtfit$RMSEA_5[1] < 0) {
        
        message("\n\n Warning: This model is something was wrong!")
        if(i==1) {
          
        } else {
          i = i-1
        }
        if(!exists(as.character(substitute(irtmodel)))) {
          itemtype <- irtmodel
          message('Item type: ', paste0(itemtype))
          
          message('\n')
        } else {
          itemtype <- NULL
        }
        
        # estimating random theta
        covdat <- data.frame(group = group)
        
        ranmod <- mixedmirt(fa_covdata, covdat, model=i, random = list(~ 1|group, ~ 1|items), itemtype = itemtype, technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
        summary(ranmod)
        raneff <- randef(ranmod)
        randTheta <- raneff$Theta
        
        result <- mirt::mirt(fa_covdata, i, itemtype=itemtype, rotate = 'geominQ', method = 'QMCEM', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999, customTheta = randTheta), maxit=100000)
        
        message(" \n ")
        message("\n\n IRT coefficients")
        print(coef(result))
        print(summary(result))
        #print(result, sort=TRUE)
        message("\n\n Use this result. Thanks a lot!")
        return(result)
      }
      
      # if reached good model fit, try to extract two more factors
      else if(mirtfit$TLI[1] >= .9 && mirtfit$CFI[1] >= .9 && mirtfit$RMSEA_5[1] < RMSEA_Criterion) {
        #else if(getClass(result)@TLI >= .9 && getClass(result)@CFI >= .9 && getClass(result)@RMSEA < .05) {
        message("\n\n Maximum Factor number archieved: ", paste0(i))
        message("\n\n Statistically Perfect! This solution is final solution. It will be save in your value. But, I'll try extract two more factors. Please check can be explain. \n")
        
        # extract one more factor
        i = i + 1
        message('Trying to extract ', paste0(i), ' factor solution')
        if(!exists(as.character(substitute(irtmodel)))) {
          itemtype <- irtmodel
          message('Item type: ', paste0(itemtype))
          
          message('\n')
        } else {
          itemtype <- NULL
        }
        
        # estimating random theta
        covdat <- data.frame(group = group)
        
        ranmod <- mixedmirt(fa_covdata, covdat, model=i, random = list(~ 1|group, ~ 1|items), itemtype = itemtype, technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
        summary(ranmod)
        raneff <- randef(ranmod)
        randTheta <- raneff$Theta
        
        
        result_onemore <- mirt::mirt(fa_covdata, i, itemtype=itemtype, rotate = 'geominQ', method = 'QMCEM', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999, customTheta = randTheta), maxit=100000)
        
        #result <- mirt::mirt(fa_covdata, i, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
        #message(" \n ")
        #message("\n\n IRT coefficients")
        #print(coef(result))
        #print(summary(result_onemore))
        
        
        message("\n\n Model Fit Indices")
        #message(" CFI"); print(getClass(result)@CFI)
        #message(" TLI"); print(getClass(result)@TLI)
        #message(" RMSEA"); print(getClass(result)@RMSEA)
        anova_fit <- anova(result, result_onemore)
        print(anova_fit)
        if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
          # Breaking for Over factor specification -- strict
          message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
          
          message("\n\n")
          return(result)
        }
        Theta_onemore <- fscores(result_onemore, full.scores=TRUE, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
        message("\n\n")
        mirtfit_onemore <- findM2(result_onemore, Theta=Theta_onemore, impute=100, calcNull=TRUE)
        print(round(mirtfit_onemore[1,], 4))
        
        
        if(i==1){
          heywood <- abs(getClass(result_onemore)@Fit$F)
        } else {
          try(require(GPArotation), silent = T)
          message('Checking Heywood cases by geominQ rotation')
          heywood <- abs(geominQ(getClass(result_onemore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by bentlerQ rotation')
          heywood2 <- abs(bentlerQ(getClass(result_onemore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by infomaxQ rotation')
          heywood3 <- abs(infomaxQ(getClass(result_onemore)@Fit$F, maxit=999999)$loadings)
        }
        
        communality <- round(result@Fit$h2, 3)
        if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
          #if(sum(communality >= .99) > 0) {
          message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
          #rm(result_onemore, Theta_onemore, mirtfit_onemore)
          return(result)
        }
        
        # Model fit check
        else if(mirtfit_onemore$TLI[1] >= 1 | mirtfit_onemore$CFI[1] >= 1 | mirtfit$TLI[1] >= mirtfit_onemore$TLI[1] | mirtfit$CFI[1] >= mirtfit_onemore$CFI[1] | mirtfit_onemore$RMSEA_5[1] < 0) {
          message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
          #rm(result_onemore, Theta_onemore, mirtfit_onemore)
          return(result)
        }
        
        # extract two more factor
        i = i + 1
        message('Trying to extract ', paste0(i), ' factor solution')
        if(!exists(as.character(substitute(irtmodel)))) {
          itemtype <- irtmodel
          message('Item type: ', paste0(itemtype))
          
          message('\n')
        } else {
          itemtype <- NULL
        }
        
        # estimating random theta
        covdat <- data.frame(group = group)
        
        ranmod <- mixedmirt(fa_covdata, covdat, model=i, random = list(~ 1|group, ~ 1|items), itemtype = itemtype, technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
        summary(ranmod)
        raneff <- randef(ranmod)
        randTheta <- raneff$Theta
        
        result_twomore <- mirt::mirt(fa_covdata, i, itemtype=itemtype, rotate = 'geominQ', method = 'QMCEM', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999, customTheta = randTheta), maxit=100000)
        #result <- mirt::mirt(fa_covdata, i, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
        #message(" \n ")
        #message("\n\n IRT coefficients")
        #print(coef(result))
        #print(summary(result_twomore))
        
        
        message("\n\n Model Fit Indices")
        #message(" CFI"); print(getClass(result)@CFI)
        #message(" TLI"); print(getClass(result)@TLI)
        #message(" RMSEA"); print(getClass(result)@RMSEA)
        anova_fit <- anova(result_onemore, result_twomore)
        print(anova_fit)
        if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
          # Breaking for Over factor specification -- strict
          message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
          #return(result)
          return(result_onemore)
        }
        message("\n\n")
        Theta_twomore <- fscores(result_twomore, full.scores=TRUE, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
        mirtfit_twomore <- findM2(result_twomore, Theta=Theta_twomore, impute=100, calcNull=TRUE)
        print(round(mirtfit_twomore[1,], 4))
        
        if(i==1){
          heywood <- abs(getClass(result_twomore)@Fit$F)
        } else {
          try(require(GPArotation), silent = T)
          message('Checking Heywood cases by geominQ rotation')
          heywood <- abs(geominQ(getClass(result_twomore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by bentlerQ rotation')
          heywood2 <- abs(bentlerQ(getClass(result_twomore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by infomaxQ rotation')
          heywood3 <- abs(infomaxQ(getClass(result_twomore)@Fit$F, maxit=999999)$loadings)
        }
        
        communality <- round(result@Fit$h2, 3)
        if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
          #if(sum(communality >= .99) > 0) {
          message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
          #rm(result_twomore, Theta_twomore, mirtfit_twomore)
          return(result_onemore)
        }
        
        # Model fit check
        else if(mirtfit_twomore$TLI[1] >= 1 | mirtfit_twomore$CFI[1] >= 1 | mirtfit_onemore$TLI[1] >= mirtfit_twomore$TLI[1] | mirtfit_onemore$CFI[1] >= mirtfit_twomore$CFI[1] | mirtfit_twomore$RMSEA_5[1] < 0) {
          message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
          #rm(result_twomore, Theta_twomore, mirtfit_twomore)
          return(result_onemore)
        }
        
        
        # extract two more factor
        i = i + 1
        message('Trying to extract ', paste0(i), ' factor solution')
        if(!exists(as.character(substitute(irtmodel)))) {
          itemtype <- irtmodel
          message('Item type: ', paste0(itemtype))
          
          message('\n')
        } else {
          itemtype <- NULL
        }
        
        # estimating random theta
        covdat <- data.frame(group = group)
        
        ranmod <- mixedmirt(fa_covdata, covdat, model=i, random = list(~ 1|group, ~ 1|items), itemtype = itemtype, technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
        summary(ranmod)
        raneff <- randef(ranmod)
        randTheta <- raneff$Theta
        
        result_threemore <- mirt::mirt(fa_covdata, i, itemtype=itemtype, rotate = 'geominQ', method = 'QMCEM', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999, customTheta = randTheta), maxit=100000)
        #result <- mirt::mirt(fa_covdata, i, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
        #message(" \n ")
        #message("\n\n IRT coefficients")
        #print(coef(result))
        #print(summary(result_threemore))
        
        
        message("\n\n Model Fit Indices")
        #message(" CFI"); print(getClass(result)@CFI)
        #message(" TLI"); print(getClass(result)@TLI)
        #message(" RMSEA"); print(getClass(result)@RMSEA)
        anova_fit <- anova(result_twomore, result_threemore)
        print(anova_fit)
        if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
          # Breaking for Over factor specification -- strict
          message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
          #return(result)
          return(result_twomore)
        }
        message("\n\n")
        Theta_threemore <- fscores(result_threemore, full.scores=TRUE, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
        mirtfit_threemore <- findM2(result_threemore, Theta=Theta_threemore, impute=100, calcNull=TRUE)
        print(round(mirtfit_threemore[1,], 4))
        
        if(i==1){
          heywood <- abs(getClass(result_threemore)@Fit$F)
        } else {
          try(require(GPArotation), silent = T)
          message('Checking Heywood cases by geominQ rotation')
          heywood <- abs(geominQ(getClass(result_threemore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by bentlerQ rotation')
          heywood2 <- abs(bentlerQ(getClass(result_threemore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by infomaxQ rotation')
          heywood3 <- abs(infomaxQ(getClass(result_threemore)@Fit$F, maxit=999999)$loadings)
        }
        
        communality <- round(result@Fit$h2, 3)
        if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
          #if(sum(communality >= .99) > 0) {
          message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
          #rm(result_threemore, Theta_threemore, mirtfit_threemore)
          return(result_twomore)
        }
        
        # Model fit check
        else if(mirtfit_threemore$TLI[1] >= 1 | mirtfit_threemore$CFI[1] >= 1 | mirtfit_twomore$TLI[1] >= mirtfit_threemore$TLI[1] | mirtfit_twomore$CFI[1] >= mirtfit_threemore$CFI[1] | mirtfit_threemore$RMSEA_5[1] < 0) {
          message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
          #rm(result_threemore, Theta_threemore, mirtfit_threemore)
          return(result_twomore)
        }
        
        
        # extract two more factor
        i = i + 1
        message('Trying to extract ', paste0(i), ' factor solution')
        if(!exists(as.character(substitute(irtmodel)))) {
          itemtype <- irtmodel
          message('Item type: ', paste0(itemtype))
          
          message('\n')
        } else {
          itemtype <- NULL
        }
        
        # estimating random theta
        covdat <- data.frame(group = group)
        
        ranmod <- mixedmirt(fa_covdata, covdat, model=i, random = list(~ 1|group, ~ 1|items), itemtype = itemtype, technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
        summary(ranmod)
        raneff <- randef(ranmod)
        randTheta <- raneff$Theta
        
        result_fourmore <- mirt::mirt(fa_covdata, i, itemtype=itemtype, rotate = 'geominQ', method = 'QMCEM', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999, customTheta = randTheta), maxit=100000)
        #result <- mirt::mirt(fa_covdata, i, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
        #message(" \n ")
        #message("\n\n IRT coefficients")
        #print(coef(result))
        #print(summary(result_fourmore))
        
        
        message("\n\n Model Fit Indices")
        #message(" CFI"); print(getClass(result)@CFI)
        #message(" TLI"); print(getClass(result)@TLI)
        #message(" RMSEA"); print(getClass(result)@RMSEA)
        anova_fit <- anova(result_threemore, result_fourmore)
        print(anova_fit)
        if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
          # Breaking for Over factor specification -- strict
          message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
          #return(result)
          return(result_threemore)
        }
        message("\n\n")
        Theta_fourmore <- fscores(result_fourmore, full.scores=TRUE, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
        mirtfit_fourmore <- findM2(result_fourmore, Theta=Theta_fourmore, impute=100, calcNull=TRUE)
        print(round(mirtfit_fourmore[1,], 4))
        
        if(i==1){
          heywood <- abs(getClass(result_fourmore)@Fit$F)
        } else {
          try(require(GPArotation), silent = T)
          message('Checking Heywood cases by geominQ rotation')
          heywood <- abs(geominQ(getClass(result_fourmore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by bentlerQ rotation')
          heywood2 <- abs(bentlerQ(getClass(result_fourmore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by infomaxQ rotation')
          heywood3 <- abs(infomaxQ(getClass(result_fourmore)@Fit$F, maxit=999999)$loadings)
        }
        
        communality <- round(result@Fit$h2, 3)
        if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
          #if(sum(communality >= .99) > 0) {
          message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
          #rm(result_fourmore, Theta_fourmore, mirtfit_fourmore)
          return(result_threemore)
        }
        
        # Model fit check
        else if(mirtfit_fourmore$TLI[1] >= 1 | mirtfit_fourmore$CFI[1] >= 1 | mirtfit_threemore$TLI[1] >= mirtfit_fourmore$TLI[1] | mirtfit_threemore$CFI[1] >= mirtfit_fourmore$CFI[1] | mirtfit_fourmore$RMSEA_5[1] < 0) {
          message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
          #rm(result_fourmore, Theta_fourmore, mirtfit_fourmore)
          return(result_threemore)
        }
        
        # extract two more factor
        i = i + 1
        message('Trying to extract ', paste0(i), ' factor solution')
        if(!exists(as.character(substitute(irtmodel)))) {
          itemtype <- irtmodel
          message('Item type: ', paste0(itemtype))
          
          message('\n')
        } else {
          itemtype <- NULL
        }
        
        # estimating random theta
        covdat <- data.frame(group = group)
        
        ranmod <- mixedmirt(fa_covdata, covdat, model=i, random = list(~ 1|group, ~ 1|items), itemtype = itemtype, technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
        summary(ranmod)
        raneff <- randef(ranmod)
        randTheta <- raneff$Theta
        
        result_fivemore <- mirt::mirt(fa_covdata, i, itemtype=itemtype, rotate = 'geominQ', method = 'QMCEM', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999, customTheta = randTheta), maxit=100000)
        #result <- mirt::mirt(fa_covdata, i, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
        #message(" \n ")
        #message("\n\n IRT coefficients")
        #print(coef(result))
        #print(summary(result_fivemore))
        
        
        message("\n\n Model Fit Indices")
        #message(" CFI"); print(getClass(result)@CFI)
        #message(" TLI"); print(getClass(result)@TLI)
        #message(" RMSEA"); print(getClass(result)@RMSEA)
        anova_fit <- anova(result_fourmore, result_fivemore)
        print(anova_fit)
        if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
          # Breaking for Over factor specification -- strict
          message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
          #return(result)
          return(result_fourmore)
        }
        message("\n\n")
        Theta_fivemore <- fscores(result_fivemore, full.scores=TRUE, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
        mirtfit_fivemore <- findM2(result_fivemore, Theta=Theta_fivemore, impute=100, calcNull=TRUE)
        print(round(mirtfit_fivemore[1,], 4))
        
        if(i==1){
          heywood <- abs(getClass(result_fivemore)@Fit$F)
        } else {
          try(require(GPArotation), silent = T)
          message('Checking Heywood cases by geominQ rotation')
          heywood <- abs(geominQ(getClass(result_fivemore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by bentlerQ rotation')
          heywood2 <- abs(bentlerQ(getClass(result_fivemore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by infomaxQ rotation')
          heywood3 <- abs(infomaxQ(getClass(result_fivemore)@Fit$F, maxit=999999)$loadings)
        }
        
        communality <- round(result@Fit$h2, 3)
        if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
          #if(sum(communality >= .99) > 0) {
          message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
          #rm(result_fivemore, Theta_fivemore, mirtfit_fivemore)
          return(result_fourmore)
        }
        
        # Model fit check
        else if(mirtfit_fivemore$TLI[1] >= 1 | mirtfit_fivemore$CFI[1] >= 1 | mirtfit_fourmore$TLI[1] >= mirtfit_fivemore$TLI[1] | mirtfit_fourmore$CFI[1] >= mirtfit_fivemore$CFI[1] | mirtfit_fivemore$RMSEA_5[1] < 0) {
          message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
          #rm(result_fivemore, Theta_fivemore, mirtfit_fivemore)
          return(result_fourmore)
        }
        
        # extract two more factor
        i = i + 1
        message('Trying to extract ', paste0(i), ' factor solution')
        if(!exists(as.character(substitute(irtmodel)))) {
          itemtype <- irtmodel
          message('Item type: ', paste0(itemtype))
          
          message('\n')
        } else {
          itemtype <- NULL
        }
        
        # estimating random theta
        covdat <- data.frame(group = group)
        
        ranmod <- mixedmirt(fa_covdata, covdat, model=i, random = list(~ 1|group, ~ 1|items), itemtype = itemtype, technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
        summary(ranmod)
        raneff <- randef(ranmod)
        randTheta <- raneff$Theta
        
        result_sixmore <- mirt::mirt(fa_covdata, i, itemtype=itemtype, rotate = 'geominQ', method = 'QMCEM', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999, customTheta = randTheta), maxit=100000)
        #result <- mirt::mirt(fa_covdata, i, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
        #message(" \n ")
        #message("\n\n IRT coefficients")
        #print(coef(result))
        #print(summary(result_sixmore))
        
        
        message("\n\n Model Fit Indices")
        #message(" CFI"); print(getClass(result)@CFI)
        #message(" TLI"); print(getClass(result)@TLI)
        #message(" RMSEA"); print(getClass(result)@RMSEA)
        anova_fit <- anova(result_fivemore, result_sixmore)
        print(anova_fit)
        if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
          # Breaking for Over factor specification -- strict
          message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
          #return(result)
          return(result_fivemore)
        }
        message("\n\n")
        Theta_sixmore <- fscores(result_sixmore, full.scores=TRUE, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
        mirtfit_sixmore <- findM2(result_sixmore, Theta=Theta_sixmore, impute=100, calcNull=TRUE)
        print(round(mirtfit_sixmore[1,], 4))
        
        if(i==1){
          heywood <- abs(getClass(result_sixmore)@Fit$F)
        } else {
          try(require(GPArotation), silent = T)
          message('Checking Heywood cases by geominQ rotation')
          heywood <- abs(geominQ(getClass(result_sixmore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by bentlerQ rotation')
          heywood2 <- abs(bentlerQ(getClass(result_sixmore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by infomaxQ rotation')
          heywood3 <- abs(infomaxQ(getClass(result_sixmore)@Fit$F, maxit=999999)$loadings)
        }
        
        communality <- round(result@Fit$h2, 3)
        if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
          #if(sum(communality >= .99) > 0) {
          message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
          #rm(result_sixmore, Theta_sixmore, mirtfit_sixmore)
          return(result_fivemore)
        }
        
        # Model fit check
        else if(mirtfit_sixmore$TLI[1] >= 1 | mirtfit_sixmore$CFI[1] >= 1 | mirtfit_fivemore$TLI[1] >= mirtfit_sixmore$TLI[1] | mirtfit_fivemore$CFI[1] >= mirtfit_sixmore$CFI[1] | mirtfit_sixmore$RMSEA_5[1] < 0) {
          message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
          #rm(result_sixmore, Theta_sixmore, mirtfit_sixmore)
          return(result_fivemore)
        }
        
        # extract two more factor
        i = i + 1
        message('Trying to extract ', paste0(i), ' factor solution')
        if(!exists(as.character(substitute(irtmodel)))) {
          itemtype <- irtmodel
          message('Item type: ', paste0(itemtype))
          
          message('\n')
        } else {
          itemtype <- NULL
        }
        
        # estimating random theta
        covdat <- data.frame(group = group)
        
        ranmod <- mixedmirt(fa_covdata, covdat, model=i, random = list(~ 1|group, ~ 1|items), itemtype = itemtype, technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999))
        summary(ranmod)
        raneff <- randef(ranmod)
        randTheta <- raneff$Theta
        
        result_sevenmore <- mirt::mirt(fa_covdata, i, itemtype=itemtype, calcNull=T, rotate = 'geominQ', method = 'QMCEM', technical=list(NCYCLES=500000, symmetric_SEM = F, MAXQUAD=9999999999999, customTheta = randTheta), maxit=100000)
        #result <- mirt::mirt(fa_covdata, i, rotate = 'geominQ', calcNull = TRUE, technical=list(NCYCLES=500000, symmetric_SEM = F,  MAXQUAD=9999999999999), itemtype=irtmodel)
        #message(" \n ")
        #message("\n\n IRT coefficients")
        #print(coef(result))
        #print(summary(result_sevenmore))
        
        
        message("\n\n Model Fit Indices")
        #message(" CFI"); print(getClass(result)@CFI)
        #message(" TLI"); print(getClass(result)@TLI)
        #message(" RMSEA"); print(getClass(result)@RMSEA)
        anova_fit <- anova(result_sixmore, result_sevenmore)
        print(anova_fit)
        if((anova_fit[1,3] < anova_fit[2,3]) >= 1 | (anova_fit[2,8] > .5)){
          # Breaking for Over factor specification -- strict
          message("\nError Occured: SABIC is bigger than previous factor model. Returning previous model result.")
          #return(result)
          return(result_sixmore)
        }
        message("\n\n")
        Theta_sevenmore <- fscores(result_sevenmore, full.scores=TRUE, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
        mirtfit_sevenmore <- findM2(result_sevenmore, Theta=Theta_sevenmore, impute=100, calcNull=TRUE)
        print(round(mirtfit_sevenmore[1,], 4))
        
        if(i==1){
          heywood <- abs(getClass(result_sevenmore)@Fit$F)
        } else {
          try(require(GPArotation), silent = T)
          message('Checking Heywood cases by geominQ rotation')
          heywood <- abs(geominQ(getClass(result_sevenmore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by bentlerQ rotation')
          heywood2 <- abs(bentlerQ(getClass(result_sevenmore)@Fit$F, maxit=999999)$loadings)
          message('Checking Heywood cases by infomaxQ rotation')
          heywood3 <- abs(infomaxQ(getClass(result_sevenmore)@Fit$F, maxit=999999)$loadings)
        }
        
        communality <- round(result@Fit$h2, 3)
        if(sum(communality >= 1) > 0 | sum(heywood >= 1) >= 1 | sum(heywood2 >= 1) >= 1 | sum(heywood3 >= 1) >= 1) {
          #if(sum(communality >= .99) > 0) {
          message("\n\n Warning: This model is heywood case(s) occured! It can not be extract more factor. Please use Final solution.\n")
          #rm(result_sevenmore, Theta_sevenmore, mirtfit_sevenmore)
          return(result_sixmore)
        }
        
        # Model fit check
        else if(mirtfit_sevenmore$TLI[1] >= 1 | mirtfit_sevenmore$CFI[1] >= 1 | mirtfit_sixmore$TLI[1] >= mirtfit_sevenmore$TLI[1] | mirtfit_sixmore$CFI[1] >= mirtfit_sevenmore$CFI[1] | mirtfit_sevenmore$RMSEA_5[1] < 0) {
          message("\n\n Warning: This model is something was wrong! It can not be extract more factor. Please use Final solution.")
          #rm(result_sevenmore, Theta_sevenmore, mirtfit_sevenmore)
          return(result_sixmore)
        }
        
        else {
          return(result_sevenmore) # if can not detect error
        }
        
        # if failed, return original result
        return(result)
      } 
      
      else if(mirtfit$TLI[1] >= .9 | mirtfit$CFI[1] >= .9 | mirtfit$RMSEA_5[1] < RMSEA_Criterion) {
        #else if(getClass(result)@TLI >= .9 | getClass(result)@CFI >= .9 | getClass(result)@RMSEA < .08) {
        message("\n\n Maximum Factor number archieved: ", paste0(i))
        message("\nYou can use this result when you get poor factor solution where next estimation.")
      } else {
        message("\n\n I'll extract more factors. Current factor number(s): ", paste0(i))}
      
      #return(mirtfit)
      #print(residuals(result))
      #par(mfrow=c(2,2))
      message(" \n ") 
      
      
      
      
      message("\n\n\n\n\n")
      
    }
  }
  
  # traditional ML & WLS & OLS
  else if(estimator == "minres" | estimator == "gls" | estimator == "GLS" |estimator == "minchi" | estimator == "MLE" | estimator == "WLS" | estimator == "wls" | estimator == "mle" | estimator == "ML" | estimator == "ml") {
    if(estimator == "MLE") { 
      estimator <- "ml"
    }
    
    #fa_covdata <- polychoric(fa_covdata)$rho
    for(i in request_factors:100) {
      message("\nOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n\nFactor number(s): ", paste0(i))
      
      result <- fa(fa_covdata,i,rotate="geominQ",cor='poly',fm=estimator, min.err = 1e-9, max.iter = 1e+6, alpha=.1, p=.05, scores=T, maxit=1e+6)
      print(result, sort=TRUE)
      
      message("\n\n\n\n\n")
      
      # heywood case check
      heywood <- result$loadings[1:ncol(fa_covdata)]
      communality <- result$communality
      if(sum(round(communality, 2) >= 1) > 0) {
        #if(sum(heywood > 1.01) > 0 | sum(communality >= 1) > 0) {
        message("\n\n Warning: This model is heywood case(s) occured! I'm fixing now!\n")
        i = i-1
        result <- fa(fa_covdata,i,rotate="geominQ",cor='poly',fm=estimator, min.err = 1e-9, max.iter = 1e+6,alpha=.1,p=.05,scores=T, maxit=1e+6)
        print(result, sort=TRUE)
        message("\n\n Use this result. Thanks a lot!")
        return(result)
      }
      
      # Model fit check
      if(result$TLI[1] >= .999 | result$TLI[1] == "-Inf" | !sum(is.na(result$RMSEA[2])==1)==0) {
        message("\n\n Warning: This model is something was wrong!")
        i = i-1
        result <- fa(fa_covdata,i,rotate="geominQ",cor='poly',fm=estimator, min.err = 1e-9, max.iter = 1e+6,alpha=.1,p=.05,scores=T, maxit=1e+6)
        print(result, sort=TRUE)
        message("\n\n Use this result. Thanks a lot!")
        return(result)
      }
      else if(result$TLI[1] >= .9 && result$RMSEA[2] <= .05) {
        message("\n\n Maximum Factor number archieved: ", paste0(i))
        message("\n\n Statistically Perfect! This solution is final solution.")
        return(result)
        break()
      } 
      
      else{
        message("\n\n I'll extract more factors. Current factor number(s): ", paste0(i))}
      
    }
  }
  
  # Robust optimization (ADF with Optimizer)
  else if(estimator == "ADF" | estimator == "adf"){
    for(i in 2:100) {
      
      
      man <- make_manifest(fa_covdata, shrink = FALSE, how = "ranks", seed = 12345, bootstrap = 10000, method = "ADF")
      res <- make_restrictions(man, factors = i, model = "EFA", discrepancy = "default")
      
      message("\n\nFactor number(s): ", paste0(i))
      
      message("\n\nThis is factor(s) pre-solution.\nI'll be repeating this task until factors what you requested.")
      efa <- Factanal(manifest = man, restrictions = res, analytic = TRUE)
      
      #print(show(efa)); print(summary(efa))
      
      ## Rotating with Geomin Criteria
      message("\n\nFactor Rotating.........")
      
      #efa.rotated <- Rotate(efa, criteria = list("geomin"), methodArgs = list(matrix = "PP", nfc_threshold = 0.25, delta = .01), normalize = "cureton-mulaik", NelderMead = FALSE)
      efa.rotated <- Rotate(efa, criteria = list("geomin"), methodArgs = list(matrix = "PP", nfc_threshold = 0.25, delta = .01), NelderMead = T)
      print(efa.rotated)
      print(summary(efa.rotated))
      summary(efa.rotated)
      
      message("\n\nFit measures.........")
      
      #fit.measures <- model_comparison(efa.rotated, correction = c("swain", "bartlett", "none"), conf.level = .9, nsim = 1001)
      fit.measures <- model_comparison(efa.rotated, correction = "none", conf.level = .9, nsim = 1001)
      message("\nRMSEA.........")
      rmsea.output <- print(fit.measures$close_fit$RMSEA$estimate)
      message("\nRMSEA 90% CI.........")
      rmsea.ci <- print(fit.measures$close_fit$RMSEA$conf.int)
      message("\nNFI.........")
      nfi.output <- print(fit.measures$fit_indices$NFI)
      message("\nTLI.........")
      tli.output <- print(fit.measures$fit_indices$TLI)
      message("\nCFI.........")
      cfi.output <- print(fit.measures$fit_indices$CFI)
      
      # Model fit check
      if(fit.measures$fit_indices$TLI >= .9 && fit.measures$close_fit$RMSEA$conf.int[1] <= .05) {
        message("\n\n Maximum Factor number archieved: ", paste0(i))
        message("\n\n Statistically Perfect! This solution is final solution.")
        return(efa.rotated)
        break()
      } 
      else if(fit.measures$fit_indices$TLI >= 1 | fit.measures$fit_indices$CFI >= 1 |fit.measures$fit_indices$TLI < 0 | fit.measures$fit_indices$CFI < 0) {
        message("\n\n Warning: This model is something was wrong! Please Use previous result.")
        break()
      }
      else if(fit.measures$fit_indices$TLI >= .9 | fit.measures$close_fit$RMSEA$conf.int[1] < .08) {
        message("\n\n Maximum Factor number archieved: ", paste0(i))
        message("\nYou can use this result when you get poor factor solution where next estimation.")
      } 
      else{
        message("\n\n I'll extract more factors. Current factor number(s): ", paste0(i))}
      
    }}
  
  ## reliability tools
  # ordinal alpha
  else if(estimator == "oa") {
    polycorr <- polychoric(fa_covdata)
    message("ordinal alpha coefficient\nGadermann, A. M., Guhn, M., & Zumbo, B. D. (2012). Estimating ordinal reliability for likert-type and ordinal item response data: A conceptual, empirical, and practical guide. Practical Assessment, Research &Evaluation, 17, 1-13\n")
    alpha(polycorr$rho)}
  
  # cronbach alpha
  else if(estimator == "alpha") {
    
    alpha(fa_covdata)}
  
  # omega
  else if(estimator == "omega") {
    
    omega(fa_covdata)}
  
  # omega ordinal
  else if(estimator == "oo") {
    polycorr <- polychoric(fa_covdata)
    message("ordinal alpha coefficient\nGadermann, A. M., Guhn, M., & Zumbo, B. D. (2012). Estimating ordinal reliability for likert-type and ordinal item response data: A conceptual, empirical, and practical guide. Practical Assessment, Research &Evaluation, 17, 1-13\n")
    omega(polycorr$rho)}
  
  
  else{
    message("\n\nError.........\nPlease select bfa(bayesian factor analysis) or default(Auto selection ADF or MLE) or ADF or MLE, mirt")
    
    readcharacter <- function()
    { 
      n <- readline(prompt="Enter an method: ")
      if(grepl("^[0-9]+$",n))
      {
        return(readcharacter())
      }
      
      method <- as.character(n)
    }
    
    print(readcharacter())
    
    ## define interactive input
    #readinteger <- function()
    #{ 
    #  n <- readline(prompt="Enter an start factor: ")
    #  if(!grepl("^[0-9]+$",n))
    #  {
    #    return(readinteger())
    #  }
    
    #  minimum_factors <- as.integer(n)
    #}
    
    #print(readinteger())
    
    return(k.aefa(dataset = fa_covdata, estimator=method))
    
  }
  
  #detach("package:vegan", unload=TRUE)
  #detach("package:FAiR", unload=TRUE); detach("package:MCMCpack", unload=TRUE); detach("package:psych", unload=TRUE); detach("package:mirt", unload=TRUE); detach("package:latticeExtra", unload=TRUE)
  
}

#######################
# data cleaning tools #
#######################

# imputation
k.imputation <- function(original, ...) {
  if(!require(mirt)) try(install.packages("mirt", dependencies = TRUE), silent=TRUE); require(rrcovNA)
  
  if(!sum(is.na(original)==1)==0) {
    message('Data imputation...\n')
    factor_structure <- k.aefa(original, estimator = 'mirt', ...)
    message('calculate factor scores')
    factor_score <- fscores(object = factor_structure, full.scores = TRUE, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
    message('impute missing values')
    imputated <- imputeMissing(factor_structure, factor_score)
    try(row.names(imputated) <- row.names(original), silent = T)
    
    head(imputated)
    return(data.frame(imputated))
    
  } else {
    return(original)
  }
  
}

k.imputation.multi <- function(original, cluster, itemtype = ..., ...) {
  if(!require(mirt)) try(install.packages("mirt", dependencies = TRUE), silent=TRUE); require(rrcovNA)
  
  if(!sum(is.na(original)==1)==0) {
    
    factor_structure <- k.aefa(original, estimator = 'multilevel', group = as.factor(cluster), ...)
    factor_score <- fscores(factor_structure, full.scores = TRUE, full.scores.SE = F, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
    imputated <- imputeMissing(factor_structure, factor_score)
    row.names(imputated) <- row.names(original)
    
    head(imputated)
    return(data.frame(imputated))
    
  } else {
    return(original)
  }
  
}

k.imputation.multiplegroup <- function(original, cluster, itemtype = ..., ...) {
  if(!require(mirt)) try(install.packages("mirt", dependencies = TRUE), silent=TRUE); require(rrcovNA)
  
  if(!sum(is.na(original)==1)==0) {
    
    factor_structure <- k.aefa(original, estimator = 'multiplegroup', group = as.factor(cluster), ...)
    factor_score <- fscores(factor_structure, full.scores = TRUE, full.scores.SE = F, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP')
    imputated <- imputeMissing(factor_structure, factor_score)
    row.names(imputated) <- row.names(original)
    
    head(imputated)
    return(data.frame(imputated))
    
  } else {
    return(original)
  }
  
}


# reverse coding
k.recode <- function(dname, target, scale) {
  
  i = 0
  x <- vector() ; y <- vector()
  x <- attributes(dname)$names
  y <- grep(target,x)
  
  for(i in 1:length(y)){
    
    dname[,y[i]] <- (scale + 1 - dname[,y[i]])
    
  }
  
  return(dname)
  
}

k.dichotomous <- function(dframe, start, end) {
  for(i in start:end) {
    dframe[,i] <- dframe[,i] >= median(dframe[,i],na.rm=T)
    dframe[,i] <- mapvalues(dframe[,i], c(TRUE, FALSE), c(1, 0))
    dframe[,i] <- as.integer(dframe[,i])
  }
  dframe <- data.frame(dframe)
  return(dframe)
}

# faking bad & aberrant response detection
fastHMM <- function(dat = ..., ...){
  temp_col <- colnames(dat)
  temp_col2 <- paste0(temp_col, "~1")
  resp <- list()
  model <- list()
  for(i in 1:length(temp_col2)){
    resp[[i]] <- as.formula(temp_col2[i])
    model[[i]] <- multinomial(link = 'identity')
  }
  
  for(i in 1:100){
    if(i == 1){
      try(mod <- depmix(response = resp, data = data.frame(dat), nstates = i, family = model), silent = T)
      try(fm <- fit(mod), silent = T)
    } else {
      fm_old <- fm
      rm(mod)
      rm(fm)
      
      try(mod <- depmix(response = resp, data = data.frame(dat), nstates = i, family = model), silent = T)
      try(fm <- fit(mod), silent = T)
      
      if(BIC(fm_old) < BIC(fm)){
        return(fm_old)
      }
      
    }
  }
}

k.faking <- function(dname = ..., formula = NULL, covdata = NULL, IRTonly = F, ...) { # for aberrant & faking response detection
  dataset <- dname
  
  if(IRTonly == F){
    for(j in 1:1000){
      try(HMM_result <- fastHMM(dat = dname))
      if(exists('HMM_result')){
        
        # HMM model check
        if(HMM_result@homogeneous == TRUE){
          message('\nThis data is homogeneous data. Please be careful for check faking response.')
        }
        
        stateDat <- data.frame(count(HMM_result@posterior, 'state'))
        judgementHMMnormal <- vector(length = length(HMM_result@posterior$state))
        judgementHMMnormal[which(HMM_result@posterior$state == stateDat$state[which(max(stateDat$freq) == stateDat$freq)])] <- TRUE
        
        
        # person-fit test in IRT
        if(sum(is.na(dataset)) == 0){
          dataset.mirt <- fastFIFA(x = dataset, covdata = covdata, formula = formula, ...)
          dataset.response <- personfit(dataset.mirt, method='MAP', QMC = T)
 
          
        } else {
                  try(dataset.mirt <- mirt::mirt(data = dataset, model = 1, itemtype = 'gpcm', covdata = covdata, formula = formula, SE = T, SE.type = 'complete', technical = list(SEtol = 1e-10), ...))
                  if(dataset.mirt@OptimInfo$converged == FALSE){
                    try(dataset.mirt <- mirt::mirt(data = dataset, model = 1, covdata = covdata, formula = formula, SE = T, SE.type = 'complete', technical = list(SEtol = 1e-10), ...))
                  }
                  if(is.na(dataset.mirt@OptimInfo$secondordertest)){
                    stop('fail to estimate standard error')
                  }
                  try(dataset <- imputeMissing(dataset.mirt, QMC = T, Theta = fscores(dataset.mirt, full.scores = T, MI = 100, QMC = T, method = 'MAP'), MI = 100))
                  try(dataset.mirt <- fastFIFA(x = dataset, covdata = covdata, formula = formula, ...))
                  try(dataset.response <- personfit(dataset.mirt, method='MAP', QMC = TRUE))
        }
        
        
        print(hist(dataset.response$Zh))
        dataset.response$Zh <- dataset.response$Zh > -2 #(if abnormal, dataset.response$Zh < -2 is right! : See Hyeongjun Kim (2015) @ SNU)
        IRTnormal <- data.frame(dataset.response$Zh)
        
        # reflect results
        if(HMM_result@homogeneous == TRUE){
          message('\nThis results may not refer detection of fake response; just a aberrant response. Like a contents non-relavant random response.')
          #output <- data.frame(dataset.response$Zh)
          
        }
        
        if(sum(is.na(dataset)) == 0){
          output <- data.frame(rowSums(data.frame(judgementHMMnormal, IRTnormal)) == 2)
        } else {
          output <- data.frame(judgementHMMnormal)
        }
        
        
        
        
        names(output) <- "normal"
        row.names(output) <- row.names(dataset)
        result <- data.frame(dataset, output)
        return(result)
      } else {
        
      }
    }
  } else {
    for(j in 1:1000){
        # person-fit test in IRT
        if(sum(is.na(dataset)) == 0){
          dataset.mirt <- fastFIFA(x = dataset, covdata = covdata, formula = formula, ...)
          dataset.response <- personfit(dataset.mirt, method='MAP', QMC = T)
          

          
        } else {
          try(dataset.mirt <- mirt::mirt(data = dataset, model = 1, itemtype = 'gpcm', covdata = covdata, formula = formula, SE = T, SE.type = 'complete', technical = list(SEtol = 1e-10), ...))
          if(dataset.mirt@OptimInfo$converged == FALSE){
            try(dataset.mirt <- mirt::mirt(data = dataset, model = 1, covdata = covdata, formula = formula, SE = T, SE.type = 'complete', technical = list(SEtol = 1e-10), ...))
          }
          if(is.na(dataset.mirt@OptimInfo$secondordertest)){
            stop('fail to estimate standard error')
          }
          try(dataset <- imputeMissing(dataset.mirt, QMC = T, Theta = fscores(dataset.mirt, full.scores = T, MI = 100, QMC = T, method = 'MAP'), MI = 100))
          try(dataset.mirt <- fastFIFA(x = dataset, covdata = covdata, formula = formula, ...))
          try(dataset.response <- personfit(dataset.mirt, method='MAP', QMC = TRUE))

        }
      
      print(hist(dataset.response$Zh))
      dataset.response$Zh <- dataset.response$Zh > -2 #(if abnormal, dataset.response$Zh < -2 is right! : See Hyeongjun Kim (2015) @ SNU)
      IRTnormal <- data.frame(dataset.response$Zh)
        
      if(sum(is.na(dataset)) == 0){
        output <- data.frame(IRTnormal)
      } else {
        stop('Please use HMM')
      }
          
        
        names(output) <- "normal"
        row.names(output) <- row.names(dataset)
        result <- data.frame(dataset, output)
        return(result)
      }
        
      }
}


aberrantZero <- function(data = ..., covdata = ..., formula = ..., ...){
  
  tempData <- data
  colnames(tempData) <- paste0("testVars", 1:ncol(tempData))
  for(i in 1:100000){
    
    if(i == 1){
      mod <- k.faking(dname = tempData, ...) # TRUE search
    } else {
      
      mod_old <- mod
      rm(mod)
      
      message('current number of samples: ', (nrow(mod_old[mod_old$normal == T,1:ncol(mod_old)-1]))) # count TRUE samples
      print(apply(mod_old[mod_old$normal == T,1:ncol(mod_old)-1], 2, table)) # count TRUE sample patterns
      mod <- k.faking(dname = mod_old[mod_old$normal == T,1:ncol(mod_old)-1], ...) # recalculate using TRUE SAMPLE ONLY
      
      if(nrow(mod) == nrow(mod_old)) { # if no difference between mod and mod_old
        # message('final normal cases: ', paste0(rownames(mod_old)))
        return(rownames(mod_old))
      }
    }

  }
}

# multilevel
k.multilevel <- function(variable = ..., group = ..., scale = ..., operation = ..., x = ..., y = ...) {
  
  ranvar_calculate <- (scale^2-1)/12
  ranvarT_calculate <- (scale^2+2*scale-2)/24
  
  #########################
  # rwg family            #
  #########################
  
  #rwg (single measurement variable)
  if(operation == 'rwg.single') {
    RWG.result<-rwg(variable,group,ranvar=ranvar_calculate) # var, group, rectangular function - (A^2-1)/12
    RWG.T.result<-rwg(variable,group,ranvar=ranvarT_calculate) # var, group, rectangular function - (A^2+2A-2)/24
    
    
    #print(summary(RWG.RELIG))
    message('James, Demaree & Wolf (1984) rwg for single item measures\n')
    
    message('sorted rwg(U)')
    print(sort(RWG.result[,2],decreasing=F))
    
    message('\nsorted rwg(T)')
    print(sort(RWG.T.result[,2],decreasing=F))
    
    message('\nmean rwg(U)')
    print(mean(RWG.result[,2]))
    
    message('\nmean rwg(T)')
    print(mean(RWG.T.result[,2]))
    
    message('\nmedian rwg(U)')
    print(median(RWG.result[,2]))
    
    message('\nmedian rwg(T)')
    print(median(RWG.T.result[,2]))
    
    message('\nmean group size')
    print(mean(RWG.result[,3]))
    #max(sort(RWG.RELIG[,2]))
    #mean(sort(RWG.RELIG[,2]))
    #min(sort(RWG.RELIG[,2]))
    hist(sort(RWG.result[,2]), xlab='rwg', main='Histogram of rwg')
    
    message('\nANOVA Table')
    ICC.calculate <- aov(variable~as.factor(group)); print(summary(ICC.calculate))
    
    message('\nEta-squared')
    print(etaSquared(ICC.calculate))
    
    message('\nICC(1)')
    print(ICC1(ICC.calculate))
    
    message('\nICC(2)')    
    print(ICC2(ICC.calculate))
    
    return(RWG.result)
  }
  else if(operation == 'rwg.multiple') {
    RWG.result<-rwg.j(variable,group,ranvar=ranvar_calculate) # var, group, rectangular function - (A^2-1)/12
    RWG.T.result<-rwg.j(variable,group,ranvar=ranvarT_calculate) # var, group, rectangular function - (A^2+2A-2)/24
    
    
    #print(summary(RWG.RELIG))
    message('James, Demaree & Wolf (1984) rwg for multi-item measures\n')
    message('sorted rwg(U)')
    print(sort(RWG.result[,2],decreasing=T))
    
    message('\nsorted rwg(T)')
    print(sort(RWG.T.result[,2],decreasing=T))
    
    message('\nmean rwg(U)')
    print(mean(RWG.result[,2]))
    
    message('\nmean rwg(T)')
    print(mean(RWG.T.result[,2]))
    
    message('\nmedian rwg(U)')
    print(median(RWG.result[,2]))
    
    message('\nmedian rwg(T)')
    print(median(RWG.T.result[,2]))
    
    hist(sort(RWG.result[,2]), xlab='rwg', main='Histogram of rwg(U)')
    
    
    message('\nANOVA Table (composite value)')
    variable2 <- rowMeans(variable)
    ICC.calculate <- aov(variable2~as.factor(group)); print(summary(ICC.calculate))
    
    message('\nEta-squared')
    print(etaSquared(ICC.calculate))
    
    message('\nICC(1)')
    print(ICC1(ICC.calculate))
    
    message('\nICC(2)')    
    print(ICC2(ICC.calculate))
    
    message('\nMulti ICC(1) & ICC(2)')
    rwg.j_result <- mult.icc(x=variable, grpid=group)
    print(rwg.j_result)
    
    message('\nMean group size')
    print(mean(RWG.result[,3]))
    message('\nStandard deviation of group size')
    print(sd(RWG.result[,3]))
    
    return(RWG.result)
  }
  else if(operation == 'rwg.lindell') {
    RWG.result<-rwg.j.lindell(variable,group,ranvar=ranvar_calculate) # var, group, rectangular function - (A^2-1)/12
    RWG.T.result<-rwg.j.lindell(variable,group,ranvar=ranvarT_calculate) # var, group, rectangular function - (A^2-1)/12
    
    #print(summary(RWG.RELIG))
    message('Lindell, & Brandt(1997; 1999) r* wg(j) for multi-item measures\n')
    message('sorted rwg(U)')
    print(sort(RWG.result[,2],decreasing=F))
    
    message('\nsorted rwg(T)')
    print(sort(RWG.T.result[,2],decreasing=F))
    
    message('\nmean rwg(U)')
    print(mean(RWG.result[,2]))
    
    message('\nmean rwg(T)')
    print(mean(RWG.T.result[,2]))
    
    message('\nmedian rwg(U)')
    print(median(RWG.result[,2]))
    
    message('\nmedian rwg(T)')
    print(median(RWG.T.result[,2]))
    hist(sort(RWG.result[,2]), xlab='rwg', main='Histogram of rwg')
    
    message('\nANOVA Table (composite value)')
    variable2 <- rowMeans(variable)
    ICC.calculate <- aov(variable2~as.factor(group)); print(summary(ICC.calculate))
    
    message('\nEta-squared')
    print(etaSquared(ICC.calculate))
    
    message('\nICC(1)')
    print(ICC1(ICC.calculate))
    
    message('\nICC(2)')    
    print(ICC2(ICC.calculate))
    
    message('\nMulti ICC(1) & ICC(2)')
    print(mult.icc(x=variable, grpid=group))
    
    return(RWG.result)
  }
  
  #########################
  # awg family            #
  #########################
  
  else if(operation == 'awg'){
    AWG.result<-awg(variable,group) # var, group
    
    #print(summary(RWG.RELIG))
    message('Brown and Hauenstein (2005) awg for multi-item measures\n')
    message('sorted awg(s)')
    print(sort(AWG.result[,2],decreasing=F))
    message('\nmean awg')
    print(mean(AWG.result[,2]))
    message('\nmean group size')
    print(mean(AWG.result[,4]))
    #max(sort(RWG.RELIG[,2]))
    #mean(sort(RWG.RELIG[,2]))
    #min(sort(RWG.RELIG[,2]))
    hist(sort(AWG.result[,2]), xlab='awg', main='Histogram of awg')
    return(AWG.result)
  }
  
  else if(operation == 'ad'){
    AD.result<-ad.m(variable,group) # var, group
    
    #print(summary(RWG.RELIG))
    message('Burke, Finkelstein and Dusig (1999) Average Deviation (AD) Agreement for multi-item measures\n')
    message('sorted ad')
    print(sort(AD.result[,2],decreasing=F))
    message('\nmean ad')
    print(mean(AD.result[,2]))
    message('\nmean group size')
    print(mean(AD.result[,3]))
    #max(sort(RWG.RELIG[,2]))
    #mean(sort(RWG.RELIG[,2]))
    #min(sort(RWG.RELIG[,2]))
    hist(sort(AD.result[,2]), xlab='ad', main='Histogram of ad')
    return(AD.result)
  }
  
  else if(operation == 'waba'){
    WABA.result<-waba(x, y, group) # var, group
    
    message('Dansereau, Alutto & Yammarino (1984) WABA II\n')
    message('Covariance Theorem')
    print(WABA.result$Cov.Theorem)
    message('\n n size of observation')
    print(WABA.result$n.obs)
    message('\n number of groups')
    print(WABA.result$n.grps)
    
    return(WABA.result)
  }
  
  
}

## meta-analysis Hunter & Schmidt Estimator
k.meta <- function(dframe = ..., # data frame
                   calc = ..., # Calculation mode
                   study = ..., n = ..., Rxy = ..., Rxx = ..., Ryy = ..., u = ..., moderator = ..., # correlation based meta-analysis
                   measure = ..., treatment_positive = ..., treatment_negative = ..., controlled_positive = ..., controlled_negative = ..., # treatment effect
                   reported_r = ..., selected_SD = ..., whole_SD = ..., # correcting range_restriction
                   ...) {
  
  if(calc == "treatment"){
    # example: meta.test <- k.meta(dframe=dat.bcg, calc="treatment", measure="RR", treatment_positive=tpos, treatment_negative=tneg, controlled_positive=cpos, controlled_negative=cneg)
    
    attach(dframe)
    # for 2x2 experiment
    message('Effect size meta analysis for treatment')
    message(' \n ')
    message('If you need help, please ?escalc and ?rma\n')
    
    dat <- escalc(measure=measure, ai=treatment_positive, bi=treatment_negative, ci=controlled_positive, di=controlled_negative, data=dframe)
    
    #detach(dframe)
    #detach(dframe)
    calculate <- rma(yi = dat$yi, vi = dat$vi, method="HS")
    
    print(calculate)
    message("r coefficient")
    print(z2r(calculate$zval))
    
    detach(dframe)
    return(calculate)
    
  }
  else if(calc == "d"){
    
  }
  else if(calc == "r" | calc == "correlation") {
    # example: k.meta(dframe=ABHt32, calc="r", study=ABHt32$study, n=ABHt32$n, Rxy=ABHt32$Rxy, Rxx=ABHt32$Rxx, Ryy=ABHt32$Ryy, u=NA, moderator=NA)
    
    data <- data.frame(study, Rxy, n, Rxx, Ryy, u, moderator)
    colnames(data) <- c("study", "Rxy", "n", "Rxx", "Ryy", "u", "moderator") 
    head(data)
    
    N <- sum(data$n)    
    k <- nrow(data)
    
    calculate <- MetaTable(data)
    
    result <- data.frame(N, k, calculate)
    
    print(result[1:5])
    print(result[6:10])
    print(result[11:15])
    FileDrawer(data)
    FunnelPlot(data)
    return(result)
    
  }
  else if(calc == "range_restriction") {
    # example: k.meta(calc='range_restriction', reported_r=0.3, selected_SD=6, whole_SD=10)
    # You can use this result where 
    # See T. Yoo. & D. Kim. (2003). A Meta-Analysis for Validity Generalization: Integrating Correlation Coefficients. Digital Business Studies, 9, 61-80.
    
    u <- selected_SD / whole_SD
    c <- sqrt((1-u^2) * reported_r^2 + u^2)
    
    message("u")
    print(u)
    message("range restrict corrected rho")
    print(c)
    return(c)
  }
  else {
    
  }
  
}



## K-means++
k.kmpp <- function(X, k) {
  n <- nrow(X)
  C <- numeric(k)
  C[1] <- sample(1:n, 1)
  
  for (i in 2:k) {
    dm <- distmat(X, X[C, ])
    pr <- apply(dm, 1, min); pr[C] <- 0
    C[i] <- sample(1:n, 1, prob = pr)
  }
  
  kmeans(X, X[C, ])
}

## Polynominal Regression
k.polynomial <- function(dependent = ..., independent = ..., moderator = ..., independent_name = ..., moderator_name = ..., dependent_name = ..., rotation_axis = ..., view_axis = ...) {
  #if(centering == TRUE){
  #  Z <- dependent
  #  X <- scale(independent, center=T, scale=T)
  #  M <- scale(moderator, center=T, scale=T)
  #  X <- X[,1]
  #  M <- M[,1]
  #} else {
  Z <- dependent
  X <- independent
  M <- moderator
  #}
  
  #poly.regression <- lm(Z ~ X + M + I(X^2) + I(M^2) + X * M)
  poly.regression <- lm(Z ~ poly(X, M, degree=2))
  message('Polynomial Regression')
  message('
          Z = intercept + X + X^2 + M + X*M + M^2')
  message('
          ---------------------------------------
          in Terms    in Coefficients
          ---------------------------------------
          intercept:  (Intercept)
          X:          poly(X, M, degree = 2)1.0
          X^2:        poly(X, M, degree = 2)2.0
          M:          poly(X, M, degree = 2)0.1
          X*M:        poly(X, M, degree = 2)1.1
          M^2:        poly(X, M, degree = 2)0.2
          ---------------------------------------')
  print(summary(poly.regression))
  #message('VIF')
  #poly.vif <- vif(poly.regression)
  #print(poly.vif)
  
  #if(sum(poly.vif>10) >= 1) {
  #  par(mfrow=c(2,2))
  #  plot(poly.regression)
  
  #  message('\nCorrelation chart')
  #  print(cor(poly.regression$model))
  #  message(' ')
  
  #  stop("STOP: VIF were too high. Please reconsider about this model.")
  #} else {
  par(mfrow=c(1,1))
  try(persp(poly.regression, X ~ M, col = rainbow(50), zlab='Z', contours = "col", theta = rotation_axis, phi = view_axis)) # xlab = independent_name, ylab = moderator_name, zlab = dependent_name, 
  return(poly.regression)
  #}
  
}

#par(mfrow=c(2,2))

k.rescale <- function(original, scalerange) {
  rescaled <- as.data.frame( lapply(original, cut, scalerange, labels=FALSE) )
  return(rescaled)
}

k.rescale.na <- function(original, scalerange, na_to_zero = ..., ...) {
  if(na_to_zero == TRUE) {
    original[is.na(original <- data.frame(original))] <- 0
    rescaled <- as.data.frame( lapply(original, cut, scalerange, labels=FALSE) )
    return(rescaled)
  } else {
    rescaled <- as.data.frame( lapply(original, cut, scalerange, labels=FALSE) )
    return(rescaled)
  }
}

splitdf <- function(dataframe, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/2))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
}

k.split <- function(dataframe){
  split <- splitdf(dataframe, seed=12345)
  str(split)
  lapply(split,nrow)
  lapply(split,head)
  
  training <- split$trainset
  testing <- split$testset
  
  training$sample <- "developmental"
  testing$sample <- "validation"
  
  split <- rbind(training, testing)
  split$sample <- as.factor(split$sample)
  str(split$sample)
  
  return(split)
}



#################################
# Analytical tools for          #
# Likert style (ordinal) scale  #
#################################

# Multiple Imputation
# (require m > 100)
likertMI <- function(model = ..., data = ..., m = 100, fun = 'sem', estimator = 'MML', parameterization = ..., chi = 'lmrr', ...) {
  
  #########################
  # Multiple              #
  # Imputation for        #
  # likert style measured #
  # variables             #
  # in SEM context        #
  #########################
  # Seongho Bae           #
  # seongho@kw.ac.kr      #
  # Nov 6th 2014          #
  # Under GNU / GPL       #
  #########################
  
  # checking packages for multiple impuation in SEM context
  if(!require('semTools')) {
    try(install.packages('semTools', dependencies = TRUE), silent=TRUE)
  }
  
  if(!require('Amelia')) {
    try(install.packages('Amelia', dependencies = TRUE), silent=TRUE)
  }
  
  if(!require('lavaan')) {
    try(install.packages('lavaan', dependencies = TRUE), silent=TRUE)
  }
  
  # loading packages for multiple imputation in SEM context
  require('semTools')
  require('Amelia')
  require('lavaan')
  
  fit <- sem(model=model, data=data.frame(data)) # extract variable names in model syntax
  message("sample size (listwise): ", paste0(nrow(data.frame(fit@Data@X))))
  
  fit_MI <- runMI(model=model, data=data.frame(data[,attributes(fit)$Model@dimNames[[1]][[1]]]), m=m, fun = fun, ordered=names(data[,attributes(fit)$Model@dimNames[[1]][[1]]]), miArgs=list(ords = attributes(fit)$Model@dimNames[[1]][[1]]), estimator = estimator, chi = chi, ...)#, control=list(optim.method="L-BFGS-B"), ...)
  
  cat(summary(fit_MI, standardize=T))
  print(inspect(fit_MI, 'fit'))
  Sys.sleep(60)
  
  fit_data <- data.frame(fit_MI@Data@X)
  colnames(fit_data) <- attributes(fit_MI)$Model@dimNames[[1]][[1]]
  message("sample size (imputated): ", paste0(nrow(data.frame(fit_MI@Data@X))))
  
  if(fun == 'cfa' | fun == 'CFA') {
    fit_MI_theta <- cfa(model=model, data=fit_data, ordered=names(fit_data), estimator = estimator, ...) #, control=list(optim.method="L-BFGS-B")
  } else if(fun == 'sem' | fun == 'SEM') {
    fit_MI_theta <- sem(model=model, data=fit_data, ordered=names(fit_data), estimator = estimator, ...) #, control=list(optim.method="L-BFGS-B")
  } else if(fun == 'growth' | fun == 'GROWTH') {
    fit_MI_theta <- growth(model=model, data=fit_data, ordered=names(fit_data), estimator = estimator, ...) #, control=list(optim.method="L-BFGS-B")
  }
  
  #   return(fit_MI)
  return(fit_MI_theta)
  
  # fit_cfaMI <- likertMI(model=model.sem, data=rawdata, m=100, estimator='WLSMV', fun='cfa', verbose=T)
  # fit_semMI <- likertMI(model=model.sem, data=rawdata, m=100, estimator='WLSMV', fun='sem', verbose=T)
}


# Full-automatically information item factor analysis until investigating simple factor structure
likertFA <- function(data = ..., ...) {
  
  fa_covdata <- data.frame(data)
  # fa_covdata <- k.imputation(fa_covdata, ...) # data imputation
  
  message('\nStage 1 of calcuation: Discriminant coefficient based Evaluation')
  ncol_stage1 <- ncol(fa_covdata)
  message('\nCurrent number of Items: ', paste0(ncol_stage1))
  result <- k.aefa(fa_covdata, estimator='mirt', ...)
  
  #   test <- data.frame(coef(result))
  #   b <- abs(test[1, grep("a[0-9]$", colnames(test))]) >= .5 # F. B. Baker (2001; p. 159; .5 <= a < 2)
  #   
  #   test <- test[1, colnames(b)] # replace
  #   
  #   c <- vector() # saving names
  #   d <- 0 # name counter
  #   require(stringr)
  #   
  #   for(i in 1:length(colnames(b))) {
  #     if((abs(test[i]) >= .5) == TRUE){
  #       d <- d+1
  #       c[d] <- str_sub(colnames(test[i]), end = -4)
  #     } else {
  #       
  #     }
  #   }
  
  find_a_lower <- which(data.frame(MDISC(result)) >= .5)
  #   find_a_upper <- which(data.frame(MDISC(result)) < 2.00) # F. B. Baker (2001; p. 159; .5 <= a < 2)
  #   find_a <- c(find_a_lower)#, find_a_upper)
  
  fa_covdata_temp <- data.frame(fa_covdata[,find_a_lower])
  
  #   #test <- vector()
  #   test <- try(mod2values(result), silent = T)
  #   if(exists("test")){
  #     test <- data.frame(test[grep("a[0-9]+", test$name),])
  #     test <- subset(test, test$value >= .4)
  #   } else {
  #     test <- data.frame()
  #   }
  #   
  #   
  #   if(nrow(test)==0){
  #     message('Passing........')
  #   } else {
  #     test$item <- as.character(test$item)
  #     test$item <- as.factor(test$item)
  #     include <- as.character(test$item)
  #     myvars <- names(fa_covdata) %in% c(include)
  #     # fa_covdata_temp <- fa_covdata[myvars]
  #     fa_covdata_temp <- subset(fa_covdata, select=myvars)
  fa_covdata <- fa_covdata_temp
  #   }
  
  # evaluation of IRT itemfit
  message('\nStage 2 of calcuation: Item Fit Based Evaluation')
  
  
  ncol_stage2 <- ncol(fa_covdata)
  if(ncol_stage2 == 0){
    stop('All items are deleted')
  } else {
    message('\nCurrent number of Items: ', paste0(ncol_stage2))
  }
  
  
  # screening items (stage 1)
  if(ncol_stage1 == ncol_stage2){
    # evaluating model
    #result <- k.aefa(fa_covdata, estimator='mirt', ...)
    #print(summary(result))
  } else {
    # evaluating model
    result <- k.aefa(fa_covdata, estimator='mirt', ...)
    #print(summary(result))
  }
  
  message('estimating Theta')
  try(test2_fscores <- fscores(object = result, rotate = 'geominQ', full.scores = T, plausible.draws = 100, QMC = TRUE, method = 'MAP', MI = 100))
  message('estimating itemfit')
  try(test2 <- itemfit(result, method = 'MAP', Theta = test2_fscores, impute = 100, QMC = TRUE), silent = T)
  
  if(exists("test2")) { # patch for mirt function bug...
    test2$cal <- test2$S_X2/test2$df.S_X2
    test2 <- subset(test2, test2$cal >= 3)
  } else {
    test2 <- data.frame()
    #test2 <- 0
  }
  
  if(nrow(test2)==0){
    message('Passing........')
  } else {
    exclude <- as.character(test2$item)
    exclude <- as.factor(exclude)
    exclude <- as.character(exclude)
    
    myvars <- names(fa_covdata) %in% c(exclude)
    fa_covdata <- fa_covdata[!myvars]
  }
  
  
  message('\nStage 3 of calcuation: Factor Loading Based Evaluation')
  
  for(i in 1:10000){    
    
    ncol_stage3 <- ncol(fa_covdata)
    if(ncol_stage3 == 0){
      stop('All items are deleted')
    } else {
      message('\nCurrent number of Items: ', paste0(ncol_stage3))
    }
    
    if(ncol_stage1 == ncol_stage3 | ncol_stage2 == ncol_stage3){
      # evaluating model
      #result <- k.aefa(fa_covdata, estimator='mirt', ...)
      #print(summary(result))
    } else {
      # evaluating model
      result <- k.aefa(fa_covdata, estimator='mirt', ...)
      #print(summary(result))
    }
    
    
    if(ncol(result@Fit$F)==1){
      rotF_geomin <- data.frame(result@Fit$F)
      h2 <- data.frame(result@Fit$h2)
    } else {
      # getting factor loadings
      message('Rotating Factor Solution now')
      rotF_geomin <- geominQ(result@Fit$F, maxit = 100000)
      rotF_geomin <- data.frame(rotF_geomin$loadings)
      h2 <- data.frame(result@Fit$h2)
    }
    
    
    if(i > 1){
      exclude_rownum1 <- NULL
      exclude_rownum1_low <- NULL
      exclude_rownum2 <- NULL
      exclude_rownum2_low <- NULL
      exclude_rownum3 <- NULL
      exclude_rownum3_low <- NULL
      exclude <- NULL
      exclude_rownum_low <- NULL
      myvars <- NULL
      
      if(length(exclude_rownum1) != 0){
        try(rm(exclude_rownum1, exclude_rownum2, exclude_rownum2_low, exclude_rownum3, exclude_rownum3_low, exclude, exclude_rownum_low, myvars), silent = T)
        
      } else {
        try(rm(exclude_rownum1, exclude_rownum1_low, exclude_rownum2, exclude_rownum2_low, exclude_rownum3, exclude_rownum3_low, exclude, exclude_rownum_low, myvars), silent = T)
        
      }
    }
    
    # evaluation of factor structure (stage 2) -- evaluating cross loadings and very small loadngs
    exclude_rownum1 <- which(rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) > 0.4) >=2) # cross loadings based delection
    if(length(exclude_rownum1) != 0){
      exclude_rownum1_low <- as.numeric(names(which(min(rowSums(abs(rotF_geomin[rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) > 0.4) >=2,]))) == rowSums(abs(rotF_geomin[rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) > 0.4) >=2,])))))
    }    
    exclude_rownum2 <- which(h2 < .5) # communality based delection
    exclude_rownum2_low <- which(h2 == min(h2))
    exclude_rownum3 <- which(rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) < 0.4) == ncol(rotF_geomin)) # fail to load any factors
    exclude_rownum3_low <- which(rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)])) == min(rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]))))
    
    #     exclude_rownum <- as.numeric(paste(c(exclude_rownum1, exclude_rownum2, exclude_rownum3))) # list of doing delection
    
    if(length(exclude_rownum3) > 0){
      message('[WARN] No loadings to any factor(s) occured!')
      exclude_rownum_low <- as.numeric(paste(c(exclude_rownum3_low))) # list of doing delection
      #       message(paste0(exclude_rownum_low))
    } else if(length(exclude_rownum1) > 0){
      message('[WARN] Cross Loadings occured!')
      exclude_rownum_low <- as.numeric(paste(c(exclude_rownum1_low))) # list of doing delection
    } else if(length(exclude_rownum2) > 0){
      message('[WARN] Communality problem occured!')
      exclude_rownum_low <- as.numeric(paste(c(exclude_rownum2_low))) # list of doing delection
    } else {
      message('[Done!]')
      exclude_rownum_low <- NULL
      #       exclude_rownum_low <- as.numeric(paste(c(exclude_rownum1_low, exclude_rownum2_low, exclude_rownum3_low))) # list of doing delection
    }
    
    #     exclude_rownum <- as.numeric(paste(c(exclude_rownum1, exclude_rownum3))) # list of doing delection
    
    # the start of variable delection
    if(length(exclude_rownum_low) > 0) { # if number of delection is non-zero
      exclude <- vector() # make vectors
      #       j <- 0 # set to zero exclude list counter
      
      #       for(i in c(exclude_rownum)){ # saving list of delection variable
      #       j <- j + 1
      #print(j)
      #message('Trying to extract ', paste0(i), ' factor solution') -- for debug code
      #         print(names(fa_covdata)[i])
      exclude <- names(data.frame(result@Data$data))[c(exclude_rownum_low)]
      message(paste0(exclude))
      #         exclude[j] <- names(fa_covdata)[c(exclude_rownum_low)]
      #       }
      
      myvars <- names(data.frame(result@Data$data)) %in% c(exclude)
      fa_covdata <- data.frame(result@Data$data)[!myvars]
      
    } else { # (length(exclude_rownum == 0))
      try(print(findM2(result, calcNull = T, Theta = fscores(result, full.scores = T, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP'), impute = 100)), silent = T)
      #cat('number of iteration : paste0(a))
      return(result)
    }
    
    
  }
}

# Full-automatically information item factor analysis until investigating simple factor structure
betaFA <- function(data = ..., ...) {
  
  fa_covdata <- data.frame(data)
  # fa_covdata <- k.imputation(fa_covdata, ...) # data imputation
  
  message('\nStage 1 of calcuation: Discriminant coefficient based Evaluation')
  ncol_stage1 <- ncol(fa_covdata)
  message('\nCurrent number of Items: ', paste0(ncol_stage1))
  result <- k.aefa(fa_covdata, estimator='mirt', ...)
  
#   test <- data.frame(coef(result))
#   b <- abs(test[1, grep("a[0-9]$", colnames(test))]) >= .5 # F. B. Baker (2001; p. 159; .5 <= a < 2)
#   
#   test <- test[1, colnames(b)] # replace
#   
#   c <- vector() # saving names
#   d <- 0 # name counter
#   require(stringr)
#   
#   for(i in 1:length(colnames(b))) {
#     if((abs(test[i]) >= .5) == TRUE){
#       d <- d+1
#       c[d] <- str_sub(colnames(test[i]), end = -4)
#     } else {
#       
#     }
#   }
  
  find_a_lower <- which(data.frame(MDISC(result)) >= .5)
#   find_a_upper <- which(data.frame(MDISC(result)) < 2.00) # F. B. Baker (2001; p. 159; .5 <= a < 2)
#   find_a <- c(find_a_lower)#, find_a_upper)
  
    fa_covdata_temp <- data.frame(fa_covdata[,find_a_lower])
  
  #   #test <- vector()
  #   test <- try(mod2values(result), silent = T)
  #   if(exists("test")){
  #     test <- data.frame(test[grep("a[0-9]+", test$name),])
  #     test <- subset(test, test$value >= .4)
  #   } else {
  #     test <- data.frame()
  #   }
  #   
  #   
  #   if(nrow(test)==0){
  #     message('Passing........')
  #   } else {
  #     test$item <- as.character(test$item)
  #     test$item <- as.factor(test$item)
  #     include <- as.character(test$item)
  #     myvars <- names(fa_covdata) %in% c(include)
  #     # fa_covdata_temp <- fa_covdata[myvars]
  #     fa_covdata_temp <- subset(fa_covdata, select=myvars)
  fa_covdata <- fa_covdata_temp
  #   }
  
  # evaluation of IRT itemfit
  message('\nStage 2 of calcuation: Item Fit Based Evaluation')
  
  
  ncol_stage2 <- ncol(fa_covdata)
  if(ncol_stage2 == 0){
    stop('All items are deleted')
  } else {
    message('\nCurrent number of Items: ', paste0(ncol_stage2))
  }
  
  
  # screening items (stage 1)
  if(ncol_stage1 == ncol_stage2){
    # evaluating model
    #result <- k.aefa(fa_covdata, estimator='mirt', ...)
    #print(summary(result))
  } else {
    # evaluating model
    result <- k.aefa(fa_covdata, estimator='mirt', ...)
    #print(summary(result))
  }
  
  message('estimating Theta')
  try(test2_fscores <- fscores(object = result, rotate = 'geominQ', full.scores = T, plausible.draws = 100, QMC = TRUE, method = 'MAP', MI = 100))
  message('estimating itemfit')
  try(test2 <- itemfit(result, method = 'MAP', Theta = test2_fscores, impute = 100, QMC = TRUE), silent = T)
  
  if(exists("test2")) { # patch for mirt function bug...
    test2$cal <- test2$S_X2/test2$df.S_X2
    test2 <- subset(test2, test2$cal >= 3)
  } else {
    test2 <- data.frame()
    #test2 <- 0
  }
  
  if(nrow(test2)==0){
    message('Passing........')
  } else {
    exclude <- as.character(test2$item)
    exclude <- as.factor(exclude)
    exclude <- as.character(exclude)
    
    myvars <- names(fa_covdata) %in% c(exclude)
    fa_covdata <- fa_covdata[!myvars]
  }
  
  
  message('\nStage 3 of calcuation: Factor Loading Based Evaluation')
  
  for(i in 1:10000){    
    
    ncol_stage3 <- ncol(fa_covdata)
    if(ncol_stage3 == 0){
      stop('All items are deleted')
    } else {
      message('\nCurrent number of Items: ', paste0(ncol_stage3))
    }
    
    if(ncol_stage1 == ncol_stage3 | ncol_stage2 == ncol_stage3){
      # evaluating model
      #result <- k.aefa(fa_covdata, estimator='mirt', ...)
      #print(summary(result))
    } else {
      # evaluating model
      result <- k.aefa(fa_covdata, estimator='mirt', ...)
      #print(summary(result))
    }
    
    
    if(ncol(result@Fit$F)==1){
      rotF_geomin <- data.frame(result@Fit$F)
      h2 <- data.frame(result@Fit$h2)
    } else {
      # getting factor loadings
      message('Rotating Factor Solution now')
      rotF_geomin <- geominQ(result@Fit$F, maxit = 100000)
      rotF_geomin <- data.frame(rotF_geomin$loadings)
      h2 <- data.frame(result@Fit$h2)
    }
    
    
    if(i > 1){
      exclude_rownum1 <- NULL
      exclude_rownum1_low <- NULL
      exclude_rownum2 <- NULL
      exclude_rownum2_low <- NULL
      exclude_rownum3 <- NULL
      exclude_rownum3_low <- NULL
      exclude <- NULL
      exclude_rownum_low <- NULL
      myvars <- NULL
      
      if(length(exclude_rownum1) != 0){
        try(rm(exclude_rownum1, exclude_rownum2, exclude_rownum2_low, exclude_rownum3, exclude_rownum3_low, exclude, exclude_rownum_low, myvars), silent = T)
        
      } else {
        try(rm(exclude_rownum1, exclude_rownum1_low, exclude_rownum2, exclude_rownum2_low, exclude_rownum3, exclude_rownum3_low, exclude, exclude_rownum_low, myvars), silent = T)
        
      }
    }
    
    # evaluation of factor structure (stage 2) -- evaluating cross loadings and very small loadngs
    exclude_rownum1 <- which(rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) > 0.4) >=2) # cross loadings based delection
    if(length(exclude_rownum1) != 0){
      exclude_rownum1_low <- as.numeric(names(which(min(rowSums(abs(rotF_geomin[rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) > 0.4) >=2,]))) == rowSums(abs(rotF_geomin[rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) > 0.4) >=2,])))))
    }    
    exclude_rownum2 <- which(h2 < (.5)^2) # communality based delection
    exclude_rownum2_low <- which(h2 == min(h2))
    exclude_rownum3 <- which(rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) < 0.4) == ncol(rotF_geomin)) # fail to load any factors
    exclude_rownum3_low <- which(rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)])) == min(rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]))))
    
#     exclude_rownum <- as.numeric(paste(c(exclude_rownum1, exclude_rownum2, exclude_rownum3))) # list of doing delection
    
    if(length(exclude_rownum3) > 0){
      message('[WARN] No loadings to any factor(s) occured!')
      exclude_rownum_low <- as.numeric(paste(c(exclude_rownum3_low))) # list of doing delection
      #       message(paste0(exclude_rownum_low))
    } else if(length(exclude_rownum1) > 0){
      message('[WARN] Cross Loadings occured!')
      exclude_rownum_low <- as.numeric(paste(c(exclude_rownum1_low))) # list of doing delection
    } else if(length(exclude_rownum2) > 0){
      message('[WARN] Communality problem occured!')
      exclude_rownum_low <- as.numeric(paste(c(exclude_rownum2_low))) # list of doing delection
    } else {
      message('[Done!]')
      exclude_rownum_low <- NULL
#       exclude_rownum_low <- as.numeric(paste(c(exclude_rownum1_low, exclude_rownum2_low, exclude_rownum3_low))) # list of doing delection
    }
    
    #     exclude_rownum <- as.numeric(paste(c(exclude_rownum1, exclude_rownum3))) # list of doing delection
    
    # the start of variable delection
    if(length(exclude_rownum_low) > 0) { # if number of delection is non-zero
      exclude <- vector() # make vectors
#       j <- 0 # set to zero exclude list counter
      
      #       for(i in c(exclude_rownum)){ # saving list of delection variable
#       j <- j + 1
      #print(j)
      #message('Trying to extract ', paste0(i), ' factor solution') -- for debug code
      #         print(names(fa_covdata)[i])
      exclude <- names(data.frame(result@Data$data))[c(exclude_rownum_low)]
      message(paste0(exclude))
      #         exclude[j] <- names(fa_covdata)[c(exclude_rownum_low)]
      #       }
      
      myvars <- names(data.frame(result@Data$data)) %in% c(exclude)
      fa_covdata <- data.frame(result@Data$data)[!myvars]
      
    } else { # (length(exclude_rownum == 0))
      # try(print(findM2(result, calcNull = T, Theta = fscores(result, full.scores = T, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP'), impute = 100)), silent = T)
      #cat('number of iteration : paste0(a))
      return(result)
    }
    
    
  }
}


# Full-automatically information item factor analysis when investigating simple factor structure automatically failed
redo_FA <- function(model = ..., ...) {
  
  
  message('Redo FA (trying to fix exploratory factor model automatically when you get strange solution using k.aefa() native function)\n')
  
  fa_covdata <- data.frame(model@Data$data)
  init <- 0
  
  for(i in 1:10000){    
    
    ncol_stage3 <- ncol(fa_covdata)
    init <- init + 1
    
    message('\niteration: ', paste0(init), '\n')
    
    if(ncol_stage3 == 0){
      stop('All items are deleted')
    } else {
      message('\nCurrent number of Items: ', paste0(ncol_stage3))
    }
    
    if(init == 1){
      result <- model
    } else {
      result <- k.aefa(fa_covdata, estimator='mirt', ...)
    }
    
    if(ncol(result@Fit$F)==1){
      rotF_geomin <- data.frame(result@Fit$F)
      h2 <- data.frame(result@Fit$h2)
    } else {
      # getting factor loadings
      message('Rotating Factor Solution now')
      rotF_geomin <- geominQ(result@Fit$F, maxit = 100000)
      rotF_geomin <- data.frame(rotF_geomin$loadings)
      h2 <- data.frame(result@Fit$h2)
    }
    
    
    
    
    # evaluation of factor structure (stage 2) -- evaluating cross loadings and very small loadngs
    exclude_rownum1 <- which(rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) > 0.4) >=2) # cross loadings based delection
    if(length(exclude_rownum1) != 0){
      exclude_rownum1_low <- as.numeric(names(which(min(rowSums(abs(rotF_geomin[rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) > 0.4) >=2,]))) == rowSums(abs(rotF_geomin[rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) > 0.4) >=2,])))))
    }    
    exclude_rownum2 <- which(h2 < (.25)) # communality based delection
    exclude_rownum2_low <- which(h2 == range(h2)[1])
    exclude_rownum3 <- which(rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) < 0.4) == ncol(rotF_geomin)) # fail to load any factors
    exclude_rownum3_low <- which(abs(rotF_geomin[1:ncol(rotF_geomin)]) == range(abs(rotF_geomin[1:ncol(rotF_geomin)]))[1])
    
    exclude_rownum <- as.numeric(paste(c(exclude_rownum1, exclude_rownum2, exclude_rownum3))) # list of doing delection
    if(length(exclude_rownum1) == 0){
      exclude_rownum_low <- as.numeric(paste(c(exclude_rownum2_low, exclude_rownum3_low))) # list of doing delection
    } else {
      exclude_rownum_low <- as.numeric(paste(c(exclude_rownum1_low, exclude_rownum2_low, exclude_rownum3_low))) # list of doing delection
    }
    
    #     exclude_rownum <- as.numeric(paste(c(exclude_rownum1, exclude_rownum3))) # list of doing delection
    
    # the start of variable delection
    if(length(exclude_rownum) > 0) { # if number of delection is non-zero
      exclude <- vector() # make vectors
      j <- 0 # set to zero exclude list counter
      
      #       for(i in c(exclude_rownum)){ # saving list of delection variable
      j <- j + 1
      #print(j)
      #message('Trying to extract ', paste0(i), ' factor solution') -- for debug code
      #         print(names(fa_covdata)[i])
      exclude <- names(fa_covdata)[c(exclude_rownum_low)]
      #         exclude[j] <- names(fa_covdata)[c(exclude_rownum_low)]
      #       }
      
      myvars <- names(fa_covdata) %in% c(exclude)
      fa_covdata <- fa_covdata[!myvars]
      
    } else { # (length(exclude_rownum == 0))
      try(print(findM2(result, calcNull = T, Theta = fscores(result, full.scores = T, plausible.draws = 100, rotate = "geominQ", MI = 100, QMC = TRUE, method = 'MAP'), impute = 100)), silent = T)
      #cat('number of iteration : paste0(a))
      return(result)
    }
  }
}

redo_FA_TAM <- function(model = ..., ...) {
  
  
  message('Redo FA (trying to fix exploratory factor model automatically when you get strange solution using k.aefa() native function)\n')
  
  fa_covdata <- data.frame(model$resp)
  init <- 0
  
  for(i in 1:10000){    
    
    ncol_stage3 <- ncol(fa_covdata)
    init <- init + 1
    
    message('\niteration: ', paste0(init), '\n')
    
    if(ncol_stage3 == 0){
      stop('All items are deleted')
    } else {
      message('\nCurrent number of Items: ', paste0(ncol_stage3))
    }
    
    if(init == 1){
      result <- model
    } else {
      result <- k.aefa(fa_covdata, estimator='TAM', ...)
    }
    
    if(ncol(result$B.stand)==1){
      rotF_geomin <- data.frame(result$B.stand)
      h2 <- data.frame(as.vector(rowSums(result$B.stand^2)))
    } else {
      # getting factor loadings
      message('Rotating Factor Solution now')
      rotF_geomin <- geominQ(result$B.stand, maxit = 100000)
      rotF_geomin <- data.frame(rotF_geomin$loadings)
      h2 <- data.frame(as.vector(rowSums(result$B.stand^2)))
    }
    
    
    
    
    # evaluation of factor structure (stage 2) -- evaluating cross loadings and very small loadngs
    exclude_rownum1 <- which(rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) > 0.4) >=2) # cross loadings based delection
    if(length(exclude_rownum1) != 0){
      exclude_rownum1_low <- as.numeric(names(which(min(rowSums(abs(rotF_geomin[rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) > 0.4) >=2,]))) == rowSums(abs(rotF_geomin[rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) > 0.4) >=2,])))))
    }    
    exclude_rownum2 <- which(h2 < (.25)) # communality based delection
    exclude_rownum2_low <- which(h2 == range(h2)[1])
    exclude_rownum3 <- which(rowSums(abs(rotF_geomin[1:ncol(rotF_geomin)]) < 0.4) == ncol(rotF_geomin)) # fail to load any factors
    exclude_rownum3_low <- which(abs(rotF_geomin[1:ncol(rotF_geomin)]) == range(abs(rotF_geomin[1:ncol(rotF_geomin)]))[1])
    
    exclude_rownum <- as.numeric(paste(c(exclude_rownum1, exclude_rownum2, exclude_rownum3))) # list of doing delection
    if(length(exclude_rownum1) == 0){
      exclude_rownum_low <- as.numeric(paste(c(exclude_rownum2_low, exclude_rownum3_low))) # list of doing delection
    } else {
      exclude_rownum_low <- as.numeric(paste(c(exclude_rownum1_low, exclude_rownum2_low, exclude_rownum3_low))) # list of doing delection
    }
    
    #     exclude_rownum <- as.numeric(paste(c(exclude_rownum1, exclude_rownum3))) # list of doing delection
    
    # the start of variable delection
    if(length(exclude_rownum) > 0) { # if number of delection is non-zero
      exclude <- vector() # make vectors
      j <- 0 # set to zero exclude list counter
      
      #       for(i in c(exclude_rownum)){ # saving list of delection variable
      j <- j + 1
      #print(j)
      #message('Trying to extract ', paste0(i), ' factor solution') -- for debug code
      #         print(names(fa_covdata)[i])
      exclude <- names(fa_covdata)[c(exclude_rownum_low)]
      #         exclude[j] <- names(fa_covdata)[c(exclude_rownum_low)]
      #       }
      
      myvars <- names(fa_covdata) %in% c(exclude)
      fa_covdata <- fa_covdata[!myvars]
      
    } else { # (length(exclude_rownum == 0))
      
      return(result)
    }
  }
}

# transform lavaan syntax to
# Exploratory Full-information item factor analysis
lavaan2aefa <- function(model = ..., data = ..., ...) {
  fit <- sem(model=model, data=data) # extract variable names in model syntax
  aefa_result <- k.aefa(dataset = data[,attributes(fit)$Model@dimNames[[1]][[1]]], estimator = 'mirt', ...)
  return(aefa_result)
}

# transform lavaan syntax to
# Exploratory Full-automatically information item factor analysis for analytical trouble shooting
lavaan2likertFA <- function(model = ..., data = ..., ...) {
  fit <- sem(model=model, data=data) # extract variable names in model syntax
  likertfa_result <- likertFA(data=data[,attributes(fit)$Model@dimNames[[1]][[1]]], ...)
  return(likertfa_result)
}

k.lca <- function(data) {
    message(paste0("Kwangwoon Automated Exploratory Factor Analysis [k.aefa] 3 -- under GNU GPL 2 license.\nk.lca: automated latent class analysis\n"))
  for(i in 1:1000) {
    message(paste0("trying to extract ", i, " classes"))
    
    if(i>1){
      mod_old <- mod_
      log_old <- mod_@Fit$logLik
      aic_old <- mod_@Fit$AIC
      bic_old <- mod_@Fit$BIC
      aicc_old <- mod_@Fit$AICc
      sabic_old <- mod_@Fit$SABIC
      dic_old <- mod_@Fit$DIC
      
      rm(mod_)
    }
    
    if(max(data, na.rm = T) - min(data, na.rm = T) == 1){
      itemtypeINPUT <- 'lca'
    } else {
      itemtypeINPUT <- 'nlca'
    }
    
    try(mod_ <- mdirt(data, i, itemtype = itemtypeINPUT, nruns = 100, GenRandomPars = F, return_max = T, QMC = T, technical = list(NCYCLES = 500)), silent = T)
#     if(exists('mod_') == F){ # failover
#       try(mod_ <- mdirt(data, i, itemtype = itemtypeINPUT, nruns = 100, GenRandomPars = T, return_max = T, QMC = T, optimizer = 'solnp', technical = list(NCYCLES = 500)), silent = T)
#     }
    if(i==1){
      
    } else {
      if(dic_old < mod_@Fit$DIC | mod_@OptimInfo$converged == F | exists('mod_') == F ){
        
        if(mod_@Model$nfact == 1){ # failover
          if(itemtypeINPUT == 'nlca'){
            itemtypeINPUT <- 'lca'
          } else {
            itemtypeINPUT <- 'nlca'
          }
          
          try(mod_ <- mdirt(data, i, itemtype = itemtypeINPUT, nruns = 100, GenRandomPars = T, return_max = T, QMC = T, technical = list(NCYCLES = 1000)))
          
          if(mod_@Model$nfact != 1){
            mod_old <- mod_
          }
        }
        
        return(mod_old)
      }
    }
  }
}

findCluster_iter <- function(data){
  for(i in 1:100){
    mod <- k.lca(data)
    if(as.integer(mod@Model$nfact) >= as.integer(round(nrow(mod@Data$data)/200, 0))){
      return(mod)
    } else {
      rm(mod)
    }
  }
}

findCluster <- function(data){

  mod <- findCluster_iter(data)
    
  print(plot(mod, facet_items = FALSE))
  print(plot(mod))
  
  fs <- fscores(mod)
  
  class_prob <- data.frame(apply(fs, 1, function(x) sample(1:mod@Model$nfact, 1, prob=x)))
  colnames(class_prob) <- "Class"
  return(class_prob)
}

k.sampling <- function(data, n){
  set.seed(1234)
  data <- data[sample(nrow(data), n), ]
  return(data)
}

# k.select <- function(model, items){
#   select <- data.frame(model@Data$data[,rownames(subset(model@Fit$F, model@Fit$F >= min(sort(model@Fit$F, decreasing = T)[1:items])))])
#   return(select)
#   
# 
# }

k.select <- function(model, items){
  select <- model@Data$data[,c(subset(colnames(model@Data$data), model@Fit$h2 >= min(sort(model@Fit$h2, decreasing = T)[1:items])))]
  return(select)
}

# finding problemistic data
k.fixdata <- function(data, start, end, bioend){
  
  data <- data[!duplicated(data[start:bioend]),] # remove duplicate cases automatically
  
  dropVector <- vector()
  j <- 0
  for(i in 1:nrow(data)){
    if(std(as.numeric(data[i,start:end])) == 0){
      j <- j + 1
      dropVector[j] <- i
      print(dropVector)
    } else {
      
    }
  }
  return(data[c(-dropVector),])
}

# surveyFA addon
fastFIFA <- function(x, covdata = NULL, formula = NULL, SE = F, SE.type = "crossprod", skipNominal = T, forceGRSM = F, assumingFake = F, masterThesis = F, forceRasch = F, unstable = F, forceMHRM = F, ...){
  for(i in 1:100){
    if (i == 1){
      message('\nfactor number: ', paste0(i))
    } else {
      message('\nfactor numbers: ', paste0(i))
    }
    
    # optimizer config
    if(length(covdata) == 0){ # if no covariate variables
      if(forceMHRM == T | forceGRSM == T | assumingFake == T | masterThesis == T){
        estimationMETHOD <- 'MHRM'
        optimINPUT <- NULL
        optimCTRL  <- NULL  #list(control = list(trace = F))
      } else if(i < 2){
        if(unstable == T){
          estimationMETHOD <- 'QMCEM'
          optimINPUT <- NULL
          optimCTRL  <- NULL  #list(control = list(trace = F))
        } else {
          estimationMETHOD <- 'EM'
          optimINPUT <- NULL
          optimCTRL  <- NULL  #list(control = list(trace = F))
        }
      } else {
        if(unstable == T){
          estimationMETHOD <- 'QMCEM'
          optimINPUT <- NULL
          optimCTRL  <- NULL  #list(control = list(trace = F))
        } else {
          estimationMETHOD <- 'MHRM'
          optimINPUT <- NULL
          optimCTRL <- NULL
        }
      }
      covdataINPUT <- NULL
      formulaINPUT <- NULL
      
      message('estimation method: ', paste0(estimationMETHOD))
      
    } else {
      
      covdataINPUT <- covdata
      formulaINPUT <- formula
      
      if(forceMHRM == T | forceGRSM == T | assumingFake == T | masterThesis == T){
        message('MHRM currently not supported with latent regressors')
        estimationMETHOD <- 'QMCEM'
        optimINPUT <- 'nlminb'
        optimCTRL  <- NULL  #list(control = list(trace = F))
      } else if(i < 2){
        if(unstable == T){
          estimationMETHOD <- 'QMCEM'
          optimINPUT <- NULL
          optimCTRL <- NULL
        } else{
          estimationMETHOD <- 'EM'
          optimINPUT <- 'nlminb'
          optimCTRL <- NULL
        }
      } else {
        estimationMETHOD <- 'QMCEM'
        optimINPUT <- 'nlminb' # NULL
        optimCTRL <- NULL
      }
      
      message('estimation method: ', paste0(estimationMETHOD))
      message('latent regression formula: ', paste0(formulaINPUT))
    }
    
    # forcing SE estimation activate
    if((sum(is.na(x)) != 0) && SE.type == 'crossprod'){
      SE <- T
      if(length(covdata) == 0){
        if(estimationMETHOD == 'MHRM'){ # Richadson (BL) isn't support MHRM estimation method
          SE.type <- 'MHRM'
          } else {
          SE.type <- 'Richardson'
          }
      } else {
        SE.type <- 'complete'
      }
    }
    
    # standard error estimation config
    if((SE == T && SE.type == 'SEM') == T){
      accelerateINPUT <- 'none'
      TOLINPUT <- NULL
      SEtolINPUT <- 1e-4
      symmetric_SEMINPUT <- TRUE
      
      message('TOL: ', 'default', ' / SEtol: ', SEtolINPUT, ' / SE.type: ', SE.type, ' / Accelerator: ', accelerateINPUT, ' / Symmetric SEM: ', symmetric_SEMINPUT)
    } else if((SE == T && estimationMETHOD == 'MHRM') == T){
      accelerateINPUT <- 'squarem'
      TOLINPUT <- NULL
      SEtolINPUT <- 1e-9
      symmetric_SEMINPUT <- TRUE
      
      message('TOL: ', 'default', ' / SEtol: ', SEtolINPUT, ' / SE.type: ', SE.type, ' / Accelerator: ', accelerateINPUT, ' / Symmetric SEM: ', symmetric_SEMINPUT)
    } else {
      accelerateINPUT <- 'squarem'
      TOLINPUT <- NULL
      SEtolINPUT <- 1e-10
      symmetric_SEMINPUT <- FALSE
      
      if(SE == T){
        message('Standard Error Estimation: On')
        message('TOL: ', 'default', ' / SEtol: ', SEtolINPUT, ' / SE.type: ', SE.type, ' / Accelerator: ', accelerateINPUT, ' / Symmetric SEM: ', symmetric_SEMINPUT)
      }
    }
    
    if(max(x, na.rm = T) - min(x, na.rm = T) == 1){ # dichotomous items
      
      message('\nMIRT model: ideal point')
      try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'ideal', method = estimationMETHOD, accelerate = accelerateINPUT, calcNull = T, technical = list(symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT, removeEmptyRows = T), TOL = TOLINPUT, covdata = covdataINPUT, formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,  SE.type = SE.type, ... = ...), silent = T)
      
      if(exists('modTEMP') == F | modTEMP@OptimInfo$converged != 1){
        
        message('\nMIRT model: Noncompensatory 2PL')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = '2PL', method = estimationMETHOD, accelerate = accelerateINPUT, calcNull = T, technical = list(symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT, removeEmptyRows = T), TOL = TOLINPUT, covdata = covdataINPUT, formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,  SE.type = SE.type, ... = ...), silent = T)
      }
      
      if(exists('modTEMP') == F | modTEMP@OptimInfo$converged != 1){
        
        message('\nMIRT model: Partially compensatory 2PL')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'partcomp', method = estimationMETHOD, accelerate = accelerateINPUT, calcNull = T, technical = list(symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT, removeEmptyRows = T), TOL = TOLINPUT, covdata = covdataINPUT, formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,  SE.type = SE.type, ... = ...), silent = T)
      }
      
      if(exists('modTEMP') == F){
        stop('Fail to find Factor solutions')
      }
      
    } else { # polytomous items
      
      # forceRasch
      if(forceRasch == T){
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'Rasch', method = estimationMETHOD, accelerate = accelerateINPUT, calcNull = T, technical = list(symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT, removeEmptyRows = T), TOL = TOLINPUT, covdata = covdataINPUT, formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE, SE.type = SE.type, ...), silent = F)
        try(return(modTEMP))
      }
      
      # for detecting fake responses
      if(SE == T){ # if grsm, Standard error estimation was unsuccessful.
        
      } else {
        #if((max(describe(x)$range) - min(describe(x)$range)) == 0 | forceGRSM == T | assumingFake == T | masterThesis == T){
        if(forceGRSM == T | assumingFake == T | masterThesis == T){
          x <- data.frame(x)
          
          if(length(which(describe(x)$max > median(describe(x)$max))) != 0){ # preventing weird input (e.g: 6 in 5 point scale)
            for(ii in which(describe(x)$max > median(describe(x)$max))){
              x[,ii] <- mapvalues(x[,ii], describe(x[,which(describe(x)$max > median(describe(x)$max))])$max, NA)
            }
          }
          
          x <- x[,which(describe(x)$range == max(describe(x)$range))]
          k <- vector()
          for(iii in 1:length(x)){
            for(j in range(na.omit(x[iii]))[1]:range(na.omit(x[iii]))[2]){
              if((sum(na.omit(x[iii]) == j) == 0) == TRUE){
                k[length(k)+1] <- colnames(x[iii])              
              }
            }
          }
          
          x <- (x[,!colnames(x) %in% k])
          for(iiii in 1:ncol(x)){
            x[,iiii] <- as.integer(x[,iiii])
          }
          
          message('\nMIRT model: graded rating scale')
          if(i == 1){
            
            try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'grsmIRT', method = estimationMETHOD, accelerate = accelerateINPUT, calcNull = T, technical = list(symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT, removeEmptyRows = T), TOL = TOLINPUT, covdata = covdataINPUT, formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE, SE.type = SE.type, ...), silent = F)
            
          } else {
            
            try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'grsm', method = estimationMETHOD, accelerate = accelerateINPUT, calcNull = T, technical = list(symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT, removeEmptyRows = T), TOL = TOLINPUT, covdata = covdataINPUT, formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE, SE.type = SE.type, ...), silent = F)
            
          }
          
          if(modTEMP@OptimInfo$converged == 1){
            skipNominal <- T
          }
        }
      }

      
      # nominal model
      if(skipNominal == F){
        message('\nMIRT model: nominal')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'nominal', method = estimationMETHOD, accelerate = accelerateINPUT, calcNull = T, technical = list(symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT, removeEmptyRows = T), TOL = TOLINPUT, covdata = covdataINPUT, formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE, SE.type = SE.type, ...), silent = F)
      } else {
        
        if(exists('modTEMP') == F){# | (max(describe(x)$range) - min(describe(x)$range)) != 0){
          message('\nMIRT model: Generalized partial credit')
          try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'gpcm', method = estimationMETHOD, accelerate = accelerateINPUT, calcNull = T, technical = list(symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT, removeEmptyRows = T), TOL = TOLINPUT, covdata = covdataINPUT, formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE, SE.type = SE.type, ...), silent = F)
        }
        
      }
      
      # generalized partial credit model (non-sequential)
      if(exists('modTEMP') == F | ((modTEMP@OptimInfo$converged != 1) && sum(modTEMP@Model$itemtype == 'grsm') == ncol(modTEMP@Data$data)) | ((modTEMP@OptimInfo$converged != 1) && sum(modTEMP@Model$itemtype == 'grsmIRT') == ncol(modTEMP@Data$data)) | (modTEMP@OptimInfo$converged != 1 && skipNominal == F)){
        message('\nMIRT model: Generalized partial credit')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'gpcm', method = estimationMETHOD, accelerate = accelerateINPUT, calcNull = T, technical = list(symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT, removeEmptyRows = T), TOL = TOLINPUT, covdata = covdataINPUT, formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE, SE.type = SE.type, ...), silent = F)
      }
      
      # graded response model (sequential)
      if(exists('modTEMP') == F | modTEMP@OptimInfo$converged != 1){
        message('\nMIRT model: Graded response')
        try(modTEMP <- mirt::mirt(data = x, model = i, method = estimationMETHOD, accelerate = accelerateINPUT, calcNull = T, technical = list(symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT, removeEmptyRows = T), TOL = TOLINPUT, covdata = covdataINPUT, formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE, SE.type = SE.type, ...), silent = T)
      }
      
      # finally, if can not converge
      if(exists('modTEMP') == F){
        stop('Fail to find Factor solutions')
      }
    }
    
    if(i == 1){ # ICC printing
      try(print(plot(modTEMP, type = 'infoSE')))
      try(print(plot(modTEMP, type = 'infotrace', facet_items = TRUE)))
      try(print(plot(modTEMP, type = 'trace')))
    }
    
    if(i == 1 && modTEMP@OptimInfo$converged != 1){
      stop('No convergence')
    }
    
    
    if(i > 1){
      if(ncol(modTEMP@Fit$F) == 1){
        rotMat <- modTEMP@Fit$F
      } else {
        rotMat <- geominQ(modTEMP@Fit$F, maxit = 1e+5)$loadings
      }
      
      if(modTEMP@Fit$DIC > modOLD@Fit$DIC | modTEMP@OptimInfo$converged != 1 | sum(round(rotMat, 2) >= 1.00) != 0 | sum(round(modTEMP@Fit$h2, 2) >= .99) != 0){ # modTEMP@Fit$AICc > modOLD@Fit$AICc | 
        message('optimal factor numbers: ', paste0(i-1))
        return(modOLD)
      } else if(sum(colSums(round(abs(rotMat), 2) > .4) < 2) != 0) {
        message('optimal factor numbers: ', paste0(i-1))
        return(modOLD)
      }
      
    }
    if(exists('modTEMP') == T){
      modOLD <- modTEMP
      try(rm(modTEMP))
    } else {
      stop('Fail to convergence')
    }
  }
}

surveyFA <- function(data = ..., covdata = NULL, formula = NULL, SE = F, SE.type = "crossprod", skipNominal = T, forceGRSM = F, assumingFake = F, masterThesis = F, forceRasch = F, unstable = F, forceMHRM = F, printFactorStructureRealtime = F, ...) {
  message('---------------------------------------------------------')
  message(' k.aefa: kwangwoon automated exploratory factor analysis ')
  message('---------------------------------------------------------\n')
  
  message('Calculating Initial Factor model')
  iteration_num <- 1
  message('Iteration: ', iteration_num, '\n')
  surveyFixMod <- fastFIFA(x = data, covdata = covdata, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, ...)
  
  itemFitDone <- FALSE
  while (!itemFitDone) {
    
    message('\nChecking item local independence assumption')
    iteration_num <- iteration_num + 1
    message('Iteration: ', iteration_num, '\n')
    
    if(sum(is.na(surveyFixMod@Data$data)) == 0){
      surveyFixMod_itemFit <- itemfit(x = surveyFixMod, Zh = T,
                                      method = 'MAP',
                                      QMC = T,
                                      fscores(surveyFixMod, method = 'MAP',
                                              QMC = T))
    } else {
      
      mirtCluster()
      surveyFixMod_itemFit <- itemfit(x = surveyFixMod, Zh = T,
                                      impute = 100,
                                      method = 'MAP',
                                      QMC = T,
                                      fscores(surveyFixMod, method = 'MAP',
                                              QMC = T, impute = 100))
      
      mirtCluster(remove = T)
      
    }
    if(sum(is.na(surveyFixMod_itemFit$p.S_X2)) == 0 && length(which(surveyFixMod_itemFit$S_X2[1:surveyFixMod@Data$nitems]/surveyFixMod_itemFit$df.S_X2[1:surveyFixMod@Data$nitems] > 3)) != 0){
      surveyFixMod <- fastFIFA(surveyFixMod@Data$data[,-which(max(surveyFixMod_itemFit$S_X2[1:surveyFixMod@Data$nitems]/surveyFixMod_itemFit$df.S_X2[1:surveyFixMod@Data$nitems]) == surveyFixMod_itemFit$S_X2[1:surveyFixMod@Data$nitems]/surveyFixMod_itemFit$df.S_X2[1:surveyFixMod@Data$nitems])], covdata = covdata, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, ...)
    } else if(length(which(surveyFixMod_itemFit$Zh[1:surveyFixMod@Data$nitems] < -2)) != 0){
      surveyFixMod <- fastFIFA(surveyFixMod@Data$data[,-which(min(surveyFixMod_itemFit$Zh[1:surveyFixMod@Data$nitems]) == surveyFixMod_itemFit$Zh[1:surveyFixMod@Data$nitems])], covdata = covdata, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, ...)
    } else {
      itemFitDone <- TRUE
    }
    
    if(printFactorStructureRealtime == T){
      message('\Realtime Factor Structure after iteration')
      print(round(GPArotation::geominQ(surveyFixMod@Fit$F, maxit = 10000)$loadings, 2))
    }
  }
  
  message('\nChecking aberrant responses')
  iteration_num <- iteration_num + 1
  message('Iteration: ', iteration_num, '\n')
  
  noAberrant <- k.faking(surveyFixMod@Data$data, IRTonly = T, covdata = covdata, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, ...)
  if(length(covdata) == 0){ # anyway, covdata is NULL
    surveyFixMod <- fastFIFA(surveyFixMod@Data$data[which(noAberrant$normal==TRUE),], covdata = covdata, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, ...)
  } else {
    covdata_workout <- covdata
    surveyFixMod <- fastFIFA(surveyFixMod@Data$data[which(noAberrant$normal==TRUE),], covdata = covdata_workout[which(noAberrant$normal==TRUE),], formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, ...)
  }
  
  if(printFactorStructureRealtime == T){
    message('\Realtime Factor Structure after iteration')
    print(round(GPArotation::geominQ(surveyFixMod@Fit$F, maxit = 10000)$loadings, 2))
  }
  
  # autofix
  fixFactorStructure_Done <- FALSE
  surveyFixMod_Workout <- surveyFixMod
  while (!fixFactorStructure_Done) {
    message('\nFixing Factor Model')
    iteration_num <- iteration_num + 1
    message('Iteration: ', iteration_num, '\n')
    
    LowCommunalities <- surveyFixMod_Workout@Fit$h2[which(min(surveyFixMod_Workout@Fit$h2) == surveyFixMod_Workout@Fit$h2)]
    
    if(ncol(surveyFixMod_Workout@Fit$F) == 1){
      Fmatrix <- surveyFixMod_Workout@Fit$F
    } else {
      Fmatrix <- GPArotation::geominQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings
    }
    
    NoLoadings <- surveyFixMod_Workout@Fit$h2[which(rowSums(abs(round(Fmatrix, 2)) < .4) == ncol(surveyFixMod_Workout@Fit$F))]
    
    
    # h2 have to >= .3
    if(LowCommunalities < .3^2){
      
      surveyFixMod_New <- fastFIFA(surveyFixMod_Workout@Data$data[,-which(min(surveyFixMod_Workout@Fit$h2) == surveyFixMod_Workout@Fit$h2)], covdata = attr(surveyFixMod_Workout@ParObjects$lrPars, 'df'), formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, ...)
      surveyFixMod_Workout <- surveyFixMod_New
      
      if(printFactorStructureRealtime == T){
        message('\Realtime Factor Structure after iteration')
        print(round(GPArotation::geominQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings, 2))
      }
      
    } else if(length(NoLoadings) != 0){
      if(as.logical(length(names(which(NoLoadings == min(NoLoadings))) != 0))){
        surveyFixMod_New <- fastFIFA(surveyFixMod_Workout@Data$data[,!colnames(surveyFixMod_Workout@Data$data) %in% names(which(NoLoadings == min(NoLoadings)))], covdata = attr(surveyFixMod_Workout@ParObjects$lrPars, 'df'), formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, ...)
        surveyFixMod_Workout <- surveyFixMod_New
        
        if(printFactorStructureRealtime == T){
          message('\Realtime Factor Structure after iteration')
          print(round(GPArotation::geominQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings, 2))
        }
        
      }
    } else {
      fixFactorStructure_Done <- TRUE
    }
    
  }
  
  if(printFactorStructureRealtime == T){
    message('\Final Factor Structure')
    print(round(GPArotation::geominQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings, 2))
  }
  return(surveyFixMod_Workout)
  
}
