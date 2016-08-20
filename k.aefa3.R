# Kwangwoon Automated Exploratory Factor Analysis (K.AEFA)
# Seongho Bae (seongho@kw.ac.kr)


##############
# aefa frontend #
##############

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
try(library(depmixS4), silent = T)
try(library(Rsolnp), silent = T)
try(library(Cairo), silent = T)
try(library(cairoDevice), silent = T)
try(library(stringr), silent = T)
try(library(SQUAREM), silent = T)
try(library(psychometric), silent = T)
try(library(psych), silent = T)
try(library(FAiR), silent = T)
try(library(bfa), silent = T)
try(library(mirt), silent = T)
try(library(latticeExtra), silent = T)
try(library(pracma), silent = T)
try(library(multilevel), silent = T)
try(library(nlme), silent = T)
try(library(lsr), silent = T)
try(library(rsm), silent = T)
try(library(car), silent = T)
try(library(TAM), silent = T)
try(library(GPArotation), silent = T)
try(library(lavaan), silent = T)
try(library(semTools), silent = T)



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

## fixTYPO for likert scaling
fixTYPO <- function(cleaningData){
  if((length(which(median(psych::describe(cleaningData)$rnage) != psych::describe(cleaningData)$range)) != 0) | length(which(median(psych::describe(cleaningData)$max) != psych::describe(cleaningData)$max)) != 0){
    for(ii in which(describe(cleaningData)$max > median(describe(cleaningData)$max))){
      cleaningData[,ii] <- mapvalues(cleaningData[,ii], describe(cleaningData[,which(describe(cleaningData)$max > median(describe(cleaningData)$max))])$max, NA)
    }
    for(ii in which(describe(cleaningData)$min < median(describe(cleaningData)$min))){
      cleaningData[,ii] <- mapvalues(cleaningData[,ii], describe(cleaningData[,which(describe(cleaningData)$min < median(describe(cleaningData)$min))])$min, NA)
    }
  }
  return(cleaningData)
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

k.faking <- function(data = ..., covdata = NULL, formula = NULL, SE = F, SE.type = "crossprod", skipNominal = T, forceGRSM = F, assumingFake = F, masterThesis = F, forceRasch = F, unstable = F, forceMHRM = F, printFactorStructureRealtime = F, itemkeys = NULL, survey.weights = NULL, IRTonly = F, ...) { # for aberrant & faking response detection
  dataset <- data
  dname <- data
  
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
          dataset.mirt <- fastFIFA(x = as.data.frame(data), covdata = as.data.frame(covdata), formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, itemkeys = itemkeys, survey.weights = survey.weights, ...)
          dataset.response <- personfit(dataset.mirt, method='MAP', QMC = T)
          
          
        } else {
          dataset.mirt <- fastFIFA(x = as.data.frame(data), covdata = as.data.frame(covdata), formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, itemkeys = itemkeys, survey.weights = survey.weights, ...)
          dataset_temp <- imputeMissing(x = dataset.mirt, Theta = fscores(dataset.mirt, method = 'MAP', QMC = T), QMC = T, impute = 100)
          dataset.mirt2 <- fastFIFA(x = as.data.frame(dataset_temp), covdata = as.data.frame(covdata), formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, itemkeys = itemkeys, survey.weights = survey.weights, ...)
          dataset.response <- personfit(dataset.mirt2, method='MAP', QMC = T)
        }
        
        
        print(hist(dataset.response$Zh))
        dataset.response$Zh <- dataset.response$Zh > -2.58 #(if abnormal, dataset.response$Zh < -2 is right! : See Hyeongjun Kim (2015) @ SNU)
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
        dataset.mirt <- fastFIFA(x = as.data.frame(data), covdata = as.data.frame(covdata), formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, itemkeys = itemkeys, survey.weights = survey.weights, ...)
        dataset.response <- personfit(dataset.mirt, method='MAP', QMC = T)
        
        
      } else {
        message('imputing missing values')
        dataset.mirt <- fastFIFA(x = as.data.frame(data), covdata = as.data.frame(covdata), formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, itemkeys = itemkeys, survey.weights = survey.weights, ...)
        dataset_temp <- imputeMissing(x = dataset.mirt, Theta = fscores(dataset.mirt, method = 'MAP', QMC = T), QMC = T, impute = 100)
        message('re-estimating parameters')
        dataset.mirt2 <- fastFIFA(x = as.data.frame(dataset_temp), covdata = as.data.frame(covdata), formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, itemkeys = itemkeys, survey.weights = survey.weights, ...)
        dataset.response <- personfit(dataset.mirt2, method='MAP', QMC = T)
      }
      
      print(hist(dataset.response$Zh))
      dataset.response$Zh <- dataset.response$Zh > -2 #(if abnormal, dataset.response$Zh < -2 is right! : See Hyeongjun Kim (2015) @ SNU)
      IRTnormal <- data.frame(dataset.response$Zh)
      
      # if(sum(is.na(dataset)) == 0){
      output <- data.frame(IRTnormal)
      # } else {
      # stop('Please use HMM')
      # }
      
      
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
      mod <- k.faking(data = tempData, ...) # TRUE search
    } else {
      
      mod_old <- mod
      rm(mod)
      
      message('current number of samples: ', (nrow(mod_old[mod_old$normal == T,1:ncol(mod_old)-1]))) # count TRUE samples
      print(apply(mod_old[mod_old$normal == T,1:ncol(mod_old)-1], 2, table)) # count TRUE sample patterns
      mod <- k.faking(data = mod_old[mod_old$normal == T,1:ncol(mod_old)-1], ...) # recalculate using TRUE SAMPLE ONLY
      
      if(nrow(mod) == nrow(mod_old)) { # if no difference between mod and mod_old
        # message('final normal cases: ', paste0(rownames(mod_old)))
        return(rownames(mod_old))
      }
    }
    
  }
}

# CTT multilevel
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
  library('semTools')
  library('Amelia')
  library('lavaan')
  
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
fastFIFA <- function(x, covdata = NULL, formula = NULL, SE = F, SE.type = "crossprod", skipNominal = T,
                     forceGRSM = F, assumingFake = F, masterThesis = F, forceRasch = F, unstable = F,
                     forceMHRM = F, forceNormalEM = F, itemkeys = NULL, survey.weights = NULL, allowMixedResponse = T,
                     forceUIRT = F, skipIdealPoint = F, MHRM_SE_draws = 1e+4, forceNRM = F, ...){
  
  for(i in 1:100){
    try(invisible(gc()), silent = T) # garbage cleaning
    
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
        optimCTRL  <- NULL
        empiricalhist <- FALSE
        NCYCLES <- 4000
      } else if(length(survey.weights) != 0) {
        estimationMETHOD <- 'QMCEM'
        optimINPUT <- NULL
        optimCTRL  <- NULL
        empiricalhist <- FALSE
        NCYCLES <- NULL
      } else if(i < 2){
        if(unstable == T){
          estimationMETHOD <- 'QMCEM'
          optimINPUT <- NULL
          optimCTRL  <- NULL
          empiricalhist <- FALSE
          NCYCLES <- NULL
        } else if (forceNormalEM == T) {
          estimationMETHOD <- 'EM'
          optimINPUT <- NULL
          optimCTRL  <- NULL
          empiricalhist <- FALSE
          NCYCLES <- 1e+5
        } else {
          estimationMETHOD <- 'EM'
          optimINPUT <- NULL
          optimCTRL  <- NULL
          empiricalhist <- TRUE
          NCYCLES <- 1e+5
        }
      } else {
        if(unstable == T){
          estimationMETHOD <- 'QMCEM'
          optimINPUT <- NULL
          optimCTRL  <- NULL
          empiricalhist <- FALSE
          NCYCLES <- NULL
        } else {
          estimationMETHOD <- 'MHRM'
          optimINPUT <- NULL
          optimCTRL <- NULL
          empiricalhist <- FALSE
          NCYCLES <- 4000
        }
      }
      covdataINPUT <- NULL
      formulaINPUT <- NULL
      
      message('estimation method: ', paste0(estimationMETHOD))
      
    } else { # with covariate
      
      covdataINPUT <- covdata
      formulaINPUT <- formula
      
      if(forceMHRM == T | forceGRSM == T | assumingFake == T | masterThesis == T){
        message('MHRM currently not supported with latent regressors')
        estimationMETHOD <- 'QMCEM'
        optimINPUT <- NULL
        optimCTRL  <- NULL
        empiricalhist <- FALSE
        NCYCLES <- 1e+4
        
      } else if(length(survey.weights) != 0) {
        estimationMETHOD <- 'QMCEM'
        optimINPUT <- NULL
        optimCTRL  <- NULL
        empiricalhist <- FALSE
        NCYCLES <- NULL
      } else if(i < 2){
        if(unstable == T){
          estimationMETHOD <- 'QMCEM'
          optimINPUT <- NULL
          optimCTRL <- NULL
          empiricalhist <- FALSE
          NCYCLES <- NULL
        } else if (forceNormalEM == T) {
          estimationMETHOD <- 'EM'
          optimINPUT <- NULL
          optimCTRL  <- NULL
          empiricalhist <- FALSE
          NCYCLES <- 1e+5
        } else {
          estimationMETHOD <- 'EM'
          optimINPUT <- NULL
          optimCTRL <- NULL
          empiricalhist <- TRUE
          NCYCLES <- 1e+5
        }
      } else {
        estimationMETHOD <- 'QMCEM'
        optimINPUT <- NULL # NULL
        optimCTRL <- NULL
        empiricalhist <- FALSE
        NCYCLES <- NULL
      }
      
      message('estimation method: ', paste0(estimationMETHOD))
      message('latent regression formula: ', paste0(formulaINPUT))
    }
    if(estimationMETHOD == 'EM'){
      message('Empirical Histogram for find Prior distribution: ', empiricalhist)
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
        if(estimationMETHOD == 'MHRM'){ # Richadson (BL) isn't support MHRM estimation method
          SE.type <- 'MHRM'
        } else {
          SE.type <- 'complete'
        }
      }
    }
    
    # standard error estimation config
    if((SE == T && SE.type == 'SEM') == T){
      accelerateINPUT <- 'none'
      TOLINPUT <- NULL
      SEtolINPUT <- 1e-5
      symmetric_SEMINPUT <- FALSE
      
      message('TOL: ', 'default', ' / SEtol: ', SEtolINPUT, ' / SE.type: ', SE.type, ' / Accelerator: ',
              accelerateINPUT, ' / Symmetric SEM: ', symmetric_SEMINPUT)
    } else if((SE == T && estimationMETHOD == 'MHRM') == T){
      accelerateINPUT <- 'squarem'
      TOLINPUT <- NULL
      SEtolINPUT <- 1e-3
      symmetric_SEMINPUT <- TRUE
      SE.type <- 'MHRM'
      
      message('TOL: ', 'default', ' / SEtol: ', SEtolINPUT, ' / SE.type: ', SE.type, ' / Accelerator: ',
              accelerateINPUT, ' / Symmetric SEM: ', symmetric_SEMINPUT)
    } else {
      accelerateINPUT <- 'squarem'
      TOLINPUT <- NULL
      SEtolINPUT <- 1e-10
      symmetric_SEMINPUT <- FALSE
      
      if(SE == T){
        message('Standard Error Estimation: On')
        message('TOL: ', 'default', ' / SEtol: ', SEtolINPUT, ' / SE.type: ', SE.type, ' / Accelerator: ',
                accelerateINPUT, ' / Symmetric SEM: ', symmetric_SEMINPUT)
      }
    }
    
    # removeEmptyRows config
    if(length(covdata) != 0){
      removeEmptyRowsConf <- FALSE
    } else {
      removeEmptyRowsConf <- TRUE
    }
    
    if(max(x, na.rm = T) - min(x, na.rm = T) == 1){ # dichotomous items
      
      # forceRasch (dichotomous)
      if(forceRasch == T){
        message('\nMIRT model: Rasch')
        
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'Rasch', method = estimationMETHOD,
                                  accelerate = accelerateINPUT, calcNull = T,
                                  technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                   removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                  formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ...), silent = F)
        try(return(modTEMP))
      }
      
      if(nrow(x) >= 2000){
        
        if(nrow(x) >= 5000){
          message('\nMIRT model: Noncompensatory 4PL')
          try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = '4PL', method = estimationMETHOD,
                                    accelerate = accelerateINPUT, calcNull = T,
                                    technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                     removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                    formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                    SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ... = ...), silent = T)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
        }
        
        
        if(exists('modTEMP') == F){
          
          message('\nMIRT model: Noncompensatory 3PL')
          try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = '3PL', method = estimationMETHOD,
                                    accelerate = accelerateINPUT, calcNull = T,
                                    technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                     removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                    formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                    SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ... = ...), silent = T)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
        }
        
        if(exists('modTEMP') == F){
          
          message('\nMIRT model: Noncompensatory 3PL with upper asymptote estimated')
          try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = '3PLu', method = estimationMETHOD,
                                    accelerate = accelerateINPUT, calcNull = T,
                                    technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                     removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                    formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                    SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ... = ...), silent = T)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
        }
        
        if(exists('modTEMP') == F){
          
          message('\nMIRT model: Partially compensatory 3PL')
          try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'PC3PL', method = estimationMETHOD,
                                    accelerate = accelerateINPUT, calcNull = T,
                                    technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                     removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                    formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                    SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ... = ...), silent = T)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
        }
      }
      
      if(exists('modTEMP') == F && skipIdealPoint == F){
        
        message('\nMIRT model: ideal point')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'ideal', method = estimationMETHOD,
                                  accelerate = accelerateINPUT, calcNull = T,
                                  technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                   removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                  formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ... = ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      if(exists('modTEMP') == F){
        
        message('\nMIRT model: Noncompensatory 2PL')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = '2PL', method = estimationMETHOD,
                                  accelerate = accelerateINPUT, calcNull = T,
                                  technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                   removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                  formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ... = ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      if(exists('modTEMP') == F){
        
        message('\nMIRT model: Partially compensatory 2PL')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'PC2PL', method = estimationMETHOD,
                                  accelerate = accelerateINPUT, calcNull = T,
                                  technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                   removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                  formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ... = ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      if(exists('modTEMP') == F && i == 1){
        
        message('\nMIRT model: spline response')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'spline', method = estimationMETHOD,
                                  accelerate = accelerateINPUT, calcNull = T,
                                  technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                   removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                  formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ... = ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      if(exists('modTEMP') == F){
        if(i == 1){
          stop('Fail to find Factor solutions: Model didn\'t converge.')
        } else {
          return(modOLD)
        }
      }
      
    } else if(length(itemkeys) != 0){ # 2-4PLNRM
      if(nrow(x) >= 2000){
        message('\nMIRT model: Noncompensatory 4PL Nominal response')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = '4PLNRM', method = estimationMETHOD,
                                  key = itemkeys, accelerate = accelerateINPUT, calcNull = T,
                                  technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                   removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                  formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ... = ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
        
        if(exists('modTEMP') == F){
          
          message('\nMIRT model: Noncompensatory 3PL Nominal response')
          try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = '3PLNRM', method = estimationMETHOD,
                                    key = itemkeys, accelerate = accelerateINPUT, calcNull = T,
                                    technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                     removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                    formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                    SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ... = ...), silent = T)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
        }
        
        
        if(exists('modTEMP') == F){
          
          message('\nMIRT model: Noncompensatory 3PL Nominal response with upper asymptote estimated')
          try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = '3PLuNRM', method = estimationMETHOD,
                                    key = itemkeys, accelerate = accelerateINPUT, calcNull = T,
                                    technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                     removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                    formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                    SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ... = ...), silent = T)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
        }
      }
      
      if(exists('modTEMP') == F){
        
        message('\nMIRT model: Noncompensatory 2PL Nominal response')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = '2PLNRM', method = estimationMETHOD,
                                  key = itemkeys, accelerate = accelerateINPUT, calcNull = T,
                                  technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                   removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                  formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ... = ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      if(exists('modTEMP') == F){
        message('\nMIRT model: Nominal response without keys')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'nominal', method = estimationMETHOD,
                                  accelerate = accelerateINPUT, calcNull = T,
                                  technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                   removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                  formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, key = NULL, ...), silent = F)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      if(exists('modTEMP') == F){
        if(i == 1){
          stop('Fail to find Factor solutions: Model didn\'t converge.')
        } else {
          return(modOLD)
        }
      }
      
    } else if((sum(psych::describe(x)$range == 1) != 0) && allowMixedResponse == T) { # mixed format (Construct Responses + Multiple Choices)
      if(skipNominal == F){
        if(skipIdealPoint == F){
          
          itemtype_mixed <- vector()
          for(i in 1:ncol(x)){
            if(psych::describe(x[,i])$range == 1){
              itemtype_mixed[i] <- 'ideal'
              dichotomous_type <- 'ideal point'
            } else {
              itemtype_mixed[i] <- 'nominal'
            }
          }
          message('\nMIRT model: nominal response + ', paste0(dichotomous_type))
          try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = itemtype_mixed, method = estimationMETHOD,
                                    accelerate = accelerateINPUT, calcNull = T,
                                    technical = list(MAXQUAD = 20000000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                     removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                    formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                    SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ...), silent = F)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
        }
        
        if(exists('modTEMP') == F){
          itemtype_mixed <- vector()
          for(i in 1:ncol(x)){
            if(psych::describe(x[,i])$range == 1){
              if(nrow(x) >= 5000){
                itemtype_mixed[i] <- '4PL'
                dichotomous_type <- '4PL'
              } else if(nrow(x) >= 2000){
                itemtype_mixed[i] <- '3PL'
                dichotomous_type <- '3PL'
              } else {
                itemtype_mixed[i] <- '2PL'
                dichotomous_type <- '2PL'
              }
            } else {
              itemtype_mixed[i] <- 'nominal'
            }
          }
          message('\nMIRT model: nominal response + ', paste0(dichotomous_type), ' (automatically set by sample size)')
          try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = itemtype_mixed, method = estimationMETHOD,
                                    accelerate = accelerateINPUT, calcNull = T,
                                    technical = list(MAXQUAD = 20000000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                     removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                    formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                    SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ...), silent = F)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
        }
        
      }
      
      # generalized partial credit model (non-sequential)
      if(exists('modTEMP') == F && forceNRM == F){
        itemtype_mixed <- vector()
        for(i in 1:ncol(x)){
          if(psych::describe(x[,i])$range == 1){
            if(nrow(x) >= 5000){
              itemtype_mixed[i] <- '4PL'
            } else if(nrow(x) >= 2000){
              itemtype_mixed[i] <- '3PL'
            } else {
              itemtype_mixed[i] <- 'ideal'
            }
          } else {
            itemtype_mixed[i] <- 'gpcm'
          }
        }
        message('\nMIRT model: Generalized partial credit + ideal or 3-4PL')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = itemtype_mixed, method = estimationMETHOD,
                                  accelerate = accelerateINPUT, calcNull = T,
                                  technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                   removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                  formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ...), silent = F)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      # graded response model (sequential)
      if(exists('modTEMP') == F && forceNRM == F){
        itemtype_mixed <- vector()
        for(i in 1:ncol(x)){
          if(psych::describe(x[,i])$range == 1){
            if(nrow(x) >= 5000){
              itemtype_mixed[i] <- '4PL'
            } else if(nrow(x) >= 2000){
              itemtype_mixed[i] <- '3PL'
            } else {
              itemtype_mixed[i] <- 'ideal'
            }
          } else {
            itemtype_mixed[i] <- 'graded'
          }
        }
        message('\nMIRT model: Graded response + ideal or 3-4PL')
        try(modTEMP <- mirt::mirt(data = x, model = i, method = estimationMETHOD, itemtype = itemtype_mixed, accelerate = accelerateINPUT,
                                  calcNull = T, technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT,
                                                                 SEtol = SEtolINPUT, removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES),
                                  TOL = TOLINPUT, covdata = covdataINPUT, formula = formulaINPUT,
                                  optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        } else {
          stop('retry to add forceMHRM = TRUE argument')
        }
      }
      
      # finally, if can not converge
      if(exists('modTEMP') == F){
        if(i == 1){
          stop('Fail to find Factor solutions: Model didn\'t converge.')
        } else {
          return(modOLD)
        }
      }
    } else { # polytomous items
      
      # forceRasch (PCM)
      if(forceRasch == T){
        message('\nMIRT model: Partial Credit')
        
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'Rasch', method = estimationMETHOD,
                                  accelerate = accelerateINPUT, calcNull = T,
                                  technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                   removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                  formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ...), silent = F)
        try(return(modTEMP))
      }
      
      # forceGRSM (polytomous)
      # for detecting fake responses
      # if(SE == T){ # if grsm, Standard error estimation was unsuccessful.
      
      # } else {
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
          
          try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'grsmIRT', method = estimationMETHOD,
                                    accelerate = accelerateINPUT, calcNull = T,
                                    technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                     removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                    formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                    SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ...), silent = F)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
          
          if(exists('modTEMP') == F){
            if(i == 1){
              stop('Fail to find Factor solutions: Model didn\'t converge.')
            } else {
              return(modOLD)
            }
          }
          
        } else {
          
          try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'grsm', method = estimationMETHOD,
                                    accelerate = accelerateINPUT, calcNull = T,
                                    technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                     removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                    formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                    SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ...), silent = F)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
          
          if(exists('modTEMP') == F){
            if(i == 1){
              stop('Fail to find Factor solutions: Model didn\'t converge.')
            } else {
              return(modOLD)
            }
          }
          
        }
        
        if(modTEMP@OptimInfo$converged == 1){
          skipNominal <- T
        }
      }
      # }
      
      # polytomous
      # nominal model
      if(skipNominal == F){
        message('\nMIRT model: nominal response')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'nominal', method = estimationMETHOD,
                                  accelerate = accelerateINPUT, calcNull = T,
                                  technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                   removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                  formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ...), silent = F)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      # generalized partial credit model (non-sequential)
      if(exists('modTEMP') == F && forceNRM == F){
        message('\nMIRT model: Generalized partial credit')
        try(modTEMP <- mirt::mirt(data = x, model = i, itemtype = 'gpcm', method = estimationMETHOD,
                                  accelerate = accelerateINPUT, calcNull = T,
                                  technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                   removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES), TOL = TOLINPUT, covdata = covdataINPUT,
                                  formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ...), silent = F)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      # graded response model (sequential)
      if(exists('modTEMP') == F && forceNRM == F){
        message('\nMIRT model: Graded response')
        try(modTEMP <- mirt::mirt(data = x, model = i, method = estimationMETHOD, accelerate = accelerateINPUT,
                                  calcNull = T, technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT,
                                                                 SEtol = SEtolINPUT, removeEmptyRows = removeEmptyRowsConf, NCYCLES = NCYCLES),
                                  TOL = TOLINPUT, covdata = covdataINPUT, formula = formulaINPUT,
                                  optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                  SE.type = SE.type, survey.weights = survey.weights, empiricalhist = empiricalhist, ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      # finally, if can not converge
      if(exists('modTEMP') == F){
        if(i == 1){
          stop('Fail to find Factor solutions: Model didn\'t converge.')
        } else {
          return(modOLD)
        }
      }
    }
    
    if(i == 1){ # ICC printing
      try(print(plot(modTEMP, type = 'infoSE')))
      try(print(plot(modTEMP, type = 'infotrace', facet_items = TRUE)))
      try(print(plot(modTEMP, type = 'trace')))
    }
    
    if(forceUIRT == T){
      return(modTEMP)
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
      
      if(modTEMP@Fit$DIC > modOLD@Fit$DIC | modTEMP@OptimInfo$converged != 1 | sum(round(modTEMP@Fit$h2, 3) >= .999) != 0){ # modTEMP@Fit$AICc > modOLD@Fit$AICc | 
        message('optimal factor numbers: ', paste0(i-1))
        return(modOLD)
      } #else if(sum(colSums(round(abs(rotMat), 2) > .4) < 2) != 0) {
        #message('optimal factor numbers: ', paste0(i-1))
        #return(modOLD)
      #}
      
    }
    if(exists('modTEMP') == T){
      modOLD <- modTEMP
      try(rm(modTEMP))
    } else {
      stop('Fail to convergence')
    }
  }
}


surveyFA <- function(data = ..., covdata = NULL, formula = NULL, SE = F,
                     SE.type = "crossprod", skipNominal = T, forceGRSM = F,
                     assumingFake = F, masterThesis = F, forceRasch = F,
                     unstable = F, forceNormalEM = F, forceMHRM = F,
                     printFactorStructureRealtime = F, itemkeys = NULL,
                     survey.weights = NULL, allowMixedResponse = T, autofix = F,
                     forceUIRT = F, skipIdealPoint = F, MHRM_SE_draws = 1e+4, bifactorSolution = F, skipS_X2 = F, forceNRM = F, needGlobalOptimal = T, ...) {
  message('---------------------------------------------------------')
  message(' k.aefa: kwangwoon automated exploratory factor analysis ')
  message('---------------------------------------------------------\n')
  
  if(bifactorSolution) {
    rotateCriteria <- 'bifactorQ'
  } else {
    rotateCriteria <- 'geominQ'
  }
  
  message('Calculating Initial Factor model')
  iteration_num <- 1
  message('Iteration: ', iteration_num, '\n')
  surveyFixMod <- fastFIFA(x = as.data.frame(data), covdata = as.data.frame(covdata),
                           formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal,
                           forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis,
                           forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM,
                           itemkeys = itemkeys, survey.weights = survey.weights, allowMixedResponse = allowMixedResponse,
                           autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
  if(needGlobalOptimal == T && forceUIRT == F){
    surveyFixMod <- deepFA(surveyFixMod)
  }
  itemFitDone <- FALSE
  while (!itemFitDone) {
    surveyFixModRAW <- data.frame(mirt::extract.mirt(surveyFixMod, 'data'))
    surveyFixModCOV <- data.frame(attr(surveyFixMod@ParObjects$lrPars, "df"))
    
    if(ncol(surveyFixModRAW) > 3){
      
      
      message('\nChecking item local independence assumption')
      iteration_num <- iteration_num + 1
      message('Iteration: ', iteration_num, '\n')
      
      if(sum(is.na(surveyFixMod@Data$data)) == 0){
        try(surveyFixMod_itemFitTest <- itemfit(x = surveyFixMod, S_X2 = T, Zh = T, infit = T,
                                                method = 'MAP',
                                                QMC = T,
                                                Theta = fscores(surveyFixMod, method = 'MAP',
                                                                QMC = T, maxit = 1e+5, rotate = rotateCriteria), rotate = rotateCriteria, maxit = 1e+5))
      } else {
        mirtCluster()
        try(surveyFixMod_itemFitTest <- itemfit(x = surveyFixMod, S_X2 = T, Zh = T, infit = T,
                                                impute = 100,
                                                method = 'MAP',
                                                QMC = T,
                                                Theta = fscores(surveyFixMod, method = 'MAP',
                                                                QMC = T, maxit = 1e+5, rotate = rotateCriteria), rotate = rotateCriteria, maxit = 1e+5))
        mirtCluster(remove = T)
      }
      if(!exists('surveyFixMod_itemFitTest')){
        S_X2 <- FALSE
        Zh <- TRUE
        infit <- TRUE
        if (sum(surveyFixMod@Model$itemtype == 'Rasch') != 0) {
          activateInfitOnly <- TRUE
          activateZhOnly <- FALSE
        } else {
          activateInfitOnly <- FALSE
          activateZhOnly <- TRUE
        }
      } else { # for normal conditions (e.g. gpcm, 2-4PL)
        S_X2 <- TRUE
        Zh <- TRUE
        infit <- TRUE
        activateInfitOnly <- FALSE
        activateZhOnly <- FALSE
      }
      
      # item fit calculation for detection of weird item(s)
      if(sum(is.na(surveyFixMod@Data$data)) == 0){
        if (exists('surveyFixMod_itemFitTest')) {
          surveyFixMod_itemFit <- surveyFixMod_itemFitTest # re use upon calculation
          rm(surveyFixMod_itemFitTest)
        } else {
          
          surveyFixMod_itemFit <- itemfit(x = surveyFixMod, S_X2 = S_X2, Zh = Zh, infit = infit,
                                          method = 'MAP',
                                          QMC = T,
                                          Theta = fscores(surveyFixMod, method = 'MAP',
                                                          QMC = T, maxit = 1e+5, rotate = rotateCriteria), maxit = 1e+5, rotate = rotateCriteria)
        }
        
      } else {
        if (exists('surveyFixMod_itemFitTest')) {
          surveyFixMod_itemFit <- surveyFixMod_itemFitTest # re use upon calculation
          rm(surveyFixMod_itemFitTest)
        } else {
          mirtCluster()
          surveyFixMod_itemFit <- itemfit(x = surveyFixMod, S_X2 = S_X2, Zh = Zh, infit = infit,
                                          impute = 100,
                                          method = 'MAP',
                                          QMC = T,
                                          Theta = fscores(surveyFixMod, method = 'MAP',
                                                          QMC = T, maxit = 1e+5, rotate = rotateCriteria), maxit = 1e+5, rotate = rotateCriteria)
          
          mirtCluster(remove = T)
        }
        
        
      }
      
      print(surveyFixMod_itemFit)
      
      # item evaluation
      if(activateInfitOnly == T){ # if can't calculate S-X2 fit when itemtype = 'Rasch'
        if(length(c(
          union(which(max((surveyFixMod_itemFit$infit)) >= 1.5),
                which(max((surveyFixMod_itemFit$outfit)) >= 1.5)))) > 0){
          
          message('\nRasch infit & outfit (.5 ~ 1.5): beta version')
          surveyFixMod <- fastFIFA(surveyFixModRAW[,-c(
            union(which(max((surveyFixMod_itemFit$infit)) == (surveyFixMod_itemFit$infit)),
                  which(max((surveyFixMod_itemFit$outfit)) == (surveyFixMod_itemFit$outfit))))],
            itemkeys = itemkeys[,-c(
              union(which(max((surveyFixMod_itemFit$infit)) == (surveyFixMod_itemFit$infit)),
                    which(max((surveyFixMod_itemFit$outfit)) == (surveyFixMod_itemFit$outfit))))], covdata = surveyFixModCOV, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, survey.weights = survey.weights, allowMixedResponse = allowMixedResponse, autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
          
        } else {
          itemFitDone <- TRUE
        }
      } else if(activateZhOnly == FALSE){ # normal IRT condition
        if(sum(is.na(surveyFixMod_itemFit$df.S_X2)) != 0){
          message('\nremoving items df is NA')
          
          surveyFixMod <- fastFIFA(surveyFixModRAW[,-which(is.na(surveyFixMod_itemFit$df.S_X2) == TRUE)], itemkeys = itemkeys[-which(is.na(surveyFixMod_itemFit$df.S_X2) == TRUE)], covdata = surveyFixModCOV, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, survey.weights = survey.weights, allowMixedResponse = allowMixedResponse, autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
          if(needGlobalOptimal == T && forceUIRT == F){
            surveyFixMod <- deepFA(surveyFixMod)
          }
        } else if(sum(na.omit(surveyFixMod_itemFit$df.S_X2) == 0) != 0){
          message('\nremoving items df is 0')
          
          surveyFixMod <- fastFIFA(surveyFixModRAW[,-which(surveyFixMod_itemFit$df.S_X2 == 0)], itemkeys = itemkeys[-which(surveyFixMod_itemFit$df.S_X2 == 0)], covdata = surveyFixModCOV, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, survey.weights = survey.weights, allowMixedResponse = allowMixedResponse, autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
          if(needGlobalOptimal == T && forceUIRT == F){
            surveyFixMod <- deepFA(surveyFixMod)
          }
        } else if(length(which(surveyFixMod_itemFit$Zh[1:surveyFixMod@Data$nitems] < -1.96)) != 0){ # Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness measurement with polychotomous item response models and standardized indices. British Journal of Mathematical and Statistical Psychology, 38(1), 67-86.
          message('\nDrasgow, F., Levine, M. V., & Williams, E. A. (1985)')
          surveyFixMod <- fastFIFA(surveyFixModRAW[,-which(min(surveyFixMod_itemFit$Zh[1:surveyFixMod@Data$nitems]) == surveyFixMod_itemFit$Zh[1:surveyFixMod@Data$nitems])], itemkeys = itemkeys[-which(min(surveyFixMod_itemFit$Zh[1:surveyFixMod@Data$nitems]) == surveyFixMod_itemFit$Zh[1:surveyFixMod@Data$nitems])], covdata = surveyFixModCOV, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, survey.weights = survey.weights, allowMixedResponse = allowMixedResponse, autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
          if(needGlobalOptimal == T && forceUIRT == F){
            surveyFixMod <- deepFA(surveyFixMod)
          }
        } else if(nrow(data) <= 5000 && sum(is.na(surveyFixMod_itemFit$p.S_X2[1:surveyFixMod@Data$nitems])) == 0 && length(which(surveyFixMod_itemFit$S_X2[1:surveyFixMod@Data$nitems]/surveyFixMod_itemFit$df.S_X2[1:surveyFixMod@Data$nitems] >= 3)) != 0 && skipS_X2 == F && surveyFixMod@Model$nfact == 1){ # Drasgow, F., Levine, M. V., Tsien, S., Williams, B., & Mead, A. D. (1995). Fitting polytomous item response theory models to multiple-choice tests. Applied Psychological Measurement, 19(2), 143-166.
          message('\nDrasgow, F., Levine, M. V., Tsien, S., Williams, B., & Mead, A. D. (1995)')
          surveyFixMod <- fastFIFA(surveyFixModRAW[,-which(max(surveyFixMod_itemFit$S_X2[1:surveyFixMod@Data$nitems]/surveyFixMod_itemFit$df.S_X2[1:surveyFixMod@Data$nitems]) == surveyFixMod_itemFit$S_X2[1:surveyFixMod@Data$nitems]/surveyFixMod_itemFit$df.S_X2[1:surveyFixMod@Data$nitems])], itemkeys = itemkeys[-which(max(surveyFixMod_itemFit$S_X2[1:surveyFixMod@Data$nitems]/surveyFixMod_itemFit$df.S_X2[1:surveyFixMod@Data$nitems]) == surveyFixMod_itemFit$S_X2[1:surveyFixMod@Data$nitems]/surveyFixMod_itemFit$df.S_X2[1:surveyFixMod@Data$nitems])], covdata = surveyFixModCOV, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, survey.weights = survey.weights, allowMixedResponse = allowMixedResponse, autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
          if(needGlobalOptimal == T && forceUIRT == F){
            surveyFixMod <- deepFA(surveyFixMod)
          }
        } else if(nrow(data) <= 5000 && sum(is.na(surveyFixMod_itemFit$p.S_X2[1:surveyFixMod@Data$nitems])) == 0 && length(which(surveyFixMod_itemFit$p.S_X2[1:surveyFixMod@Data$nitems] < .05)) != 0 && skipS_X2 == F && surveyFixMod@Model$nfact == 1){ # Kang, T., & Chen, T. T. (2008). Performance of the Generalized S‐X2 Item Fit Index for Polytomous IRT Models. Journal of Educational Measurement, 45(4), 391-406.; Reise, S. P. (1990). A comparison of item- and person-fit methods of assessing model-data fit in IRT. Applied Psychological Measurement, 14, 127-137.
          message('\nKang, T., & Chen, T. T. (2008); Reise, S. P. (1990)')
          surveyFixMod <- fastFIFA(surveyFixModRAW[,-which(max(surveyFixMod_itemFit$S_X2[1:surveyFixMod@Data$nitems]/surveyFixMod_itemFit$df.S_X2[1:surveyFixMod@Data$nitems]) == surveyFixMod_itemFit$S_X2[1:surveyFixMod@Data$nitems]/surveyFixMod_itemFit$df.S_X2[1:surveyFixMod@Data$nitems])], itemkeys = itemkeys[-which(max(surveyFixMod_itemFit$S_X2[1:surveyFixMod@Data$nitems]/surveyFixMod_itemFit$df.S_X2[1:surveyFixMod@Data$nitems]) == surveyFixMod_itemFit$S_X2[1:surveyFixMod@Data$nitems]/surveyFixMod_itemFit$df.S_X2[1:surveyFixMod@Data$nitems])], covdata = surveyFixModCOV, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, survey.weights = survey.weights, allowMixedResponse = allowMixedResponse, autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
          if(needGlobalOptimal == T && forceUIRT == F){
            surveyFixMod <- deepFA(surveyFixMod)
          }
        } else if (forceRasch == T) {
          if(length(c(
            union(which(max((surveyFixMod_itemFit$infit)) >= 1.5),
                  which(max((surveyFixMod_itemFit$outfit)) >= 1.5)))) > 0){
            
            message('\nRasch infit & outfit (.5 ~ 1.5): beta version')
            surveyFixMod <- fastFIFA(surveyFixModRAW[,-c(
              union(which(max((surveyFixMod_itemFit$infit)) == (surveyFixMod_itemFit$infit)),
                    which(max((surveyFixMod_itemFit$outfit)) == (surveyFixMod_itemFit$outfit))))],
              itemkeys = itemkeys[,-c(
                union(which(max((surveyFixMod_itemFit$infit)) == (surveyFixMod_itemFit$infit)),
                      which(max((surveyFixMod_itemFit$outfit)) == (surveyFixMod_itemFit$outfit))))], covdata = surveyFixModCOV, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, survey.weights = survey.weights, allowMixedResponse = allowMixedResponse, autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
            
          } else {
            itemFitDone <- TRUE
          }
        } else {
          itemFitDone <- TRUE
        }
      } else if(length(which(surveyFixMod_itemFit$Zh[1:surveyFixMod@Data$nitems] < -1.96)) != 0){ # if can't calculate S-X2 fit when itemtype = 'ideal'; Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness measurement with polychotomous item response models and standardized indices. British Journal of Mathematical and Statistical Psychology, 38(1), 67-86.
        message('\nDrasgow, F., Levine, M. V., & Williams, E. A. (1985)')
        surveyFixMod <- fastFIFA(surveyFixModRAW[,-which(min(surveyFixMod_itemFit$Zh[1:surveyFixMod@Data$nitems]) == surveyFixMod_itemFit$Zh[1:surveyFixMod@Data$nitems])], itemkeys = itemkeys[-which(min(surveyFixMod_itemFit$Zh[1:surveyFixMod@Data$nitems]) == surveyFixMod_itemFit$Zh[1:surveyFixMod@Data$nitems])], covdata = surveyFixModCOV, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, survey.weights = survey.weights, allowMixedResponse = allowMixedResponse, autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
        if(needGlobalOptimal == T && forceUIRT == F){
          surveyFixMod <- deepFA(surveyFixMod)
        }
      } else {
        itemFitDone <- TRUE
      }
    } else {
      itemFitDone <- TRUE
    }
    
    
    if(printFactorStructureRealtime == T){
      message('\nRealtime Factor Structure after iteration')
      if(ncol(surveyFixMod@Fit$F)==1){
        print(round(surveyFixMod@Fit$F, 2))
      } else {
        if(bifactorSolution){
          print(round(GPArotation::bifactorQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings, 2))
        } else {
          print(round(GPArotation::geominQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings, 2))
        }
      }
    }
  } # the end of while loop
  
  if(autofix == T){
    
    message('\nChecking aberrant responses')
    iteration_num <- iteration_num + 1
    message('Iteration: ', iteration_num, '\n')
    
    surveyFixModRAW <- data.frame(mirt::extract.mirt(surveyFixMod, 'data'))
    surveyFixModCOV <- data.frame(attr(surveyFixMod@ParObjects$lrPars, "df"))
    
    noAberrant <- k.faking(surveyFixModRAW, IRTonly = T, itemkeys = itemkeys, covdata = surveyFixModCOV, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, survey.weights = survey.weights, allowMixedResponse = allowMixedResponse, autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
    if(length(covdata) == 0){ # anyway, covdata is NULL
      surveyFixMod <- fastFIFA(surveyFixModRAW[which(noAberrant$normal==TRUE),], itemkeys = itemkeys, covdata = surveyFixModCOV, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, survey.weights = survey.weights[which(noAberrant$normal==TRUE)], allowMixedResponse = allowMixedResponse, autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
    } else {
      covdata_workout <- surveyFixModCOV
      surveyFixMod <- fastFIFA(surveyFixModRAW[which(noAberrant$normal==TRUE),], itemkeys = itemkeys, covdata = covdata_workout[which(noAberrant$normal==TRUE),], formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, survey.weights = survey.weights[which(noAberrant$normal==TRUE)], allowMixedResponse = allowMixedResponse, autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
    }
    
    if(printFactorStructureRealtime == T){
      message('\nRealtime Factor Structure after iteration')
      if(ncol(surveyFixMod@Fit$F)==1){
        print(round(surveyFixMod@Fit$F, 2))
      } else {
        if(bifactorSolution){
          print(round(GPArotation::bifactorQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings, 2))
        } else {
          print(round(GPArotation::geominQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings, 2))
        }
      }
    }
    
    # autofix config
    fixFactorStructure_Done <- FALSE
    surveyFixMod_Workout <- surveyFixMod
    
    while (!fixFactorStructure_Done) { # start of while loop
      surveyFixModRAW <- data.frame(mirt::extract.mirt(surveyFixMod_Workout, 'data')) # update data.frame
      surveyFixModCOV <- data.frame(attr(surveyFixMod_Workout@ParObjects$lrPars, "df"))
      
      message('\nFixing Factor Model')
      iteration_num <- iteration_num + 1
      message('Iteration: ', iteration_num, '\n')
      
      tempG <- mirt::extract.mirt(surveyFixMod_Workout, 'itemtype')
      if(sum(tempG != 'Rasch') != 0){
        LowCommunalities <- surveyFixMod_Workout@Fit$h2[which(min(surveyFixMod_Workout@Fit$h2) == surveyFixMod_Workout@Fit$h2)]
      }
      
      
      if(ncol(surveyFixMod_Workout@Fit$F) == 1){
        Fmatrix <- surveyFixMod_Workout@Fit$F
      } else {
        if(bifactorSolution){
          Fmatrix <- GPArotation::bifactorQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings
        } else {
          Fmatrix <- GPArotation::geominQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings
        }
      }
      
      NoLoadings <- surveyFixMod_Workout@Fit$h2[which(rowSums(abs(round(Fmatrix, 2)) < .4) == ncol(surveyFixMod_Workout@Fit$F))]
      
      
      # h2 have to >= .3
      if(sum(tempG != 'Rasch') == 0){
        fixFactorStructure_Done <- TRUE
      } else if(LowCommunalities < .3^2){
        
        surveyFixMod_New <- fastFIFA(surveyFixModRAW[,-which(min(surveyFixMod_Workout@Fit$h2) == surveyFixMod_Workout@Fit$h2)], itemkeys = itemkeys[-which(min(surveyFixMod_Workout@Fit$h2) == surveyFixMod_Workout@Fit$h2)], covdata = surveyFixModCOV, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, survey.weights = survey.weights, allowMixedResponse = allowMixedResponse, autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
        surveyFixMod_Workout <- surveyFixMod_New
        
        if(printFactorStructureRealtime == T){
          message('\nRealtime Factor Structure after iteration')
          if(ncol(surveyFixMod_Workout@Fit$F)==1){
            print(round(surveyFixMod_Workout@Fit$F, 2))
          } else {
            if(bifactorSolution){
              print(round(GPArotation::bifactorQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings, 2))
            } else {
              print(round(GPArotation::geominQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings, 2))
            }
          }
        }
        
      } else if(length(NoLoadings) != 0){ # noloadings
        if(as.logical(length(names(which(NoLoadings == min(NoLoadings))) != 0))){
          surveyFixMod_New <- fastFIFA(surveyFixModRAW[,!colnames(surveyFixModRAW) %in% names(which(NoLoadings == min(NoLoadings)))], itemkeys = itemkeys[,!colnames(surveyFixModRAW) %in% names(which(NoLoadings == min(NoLoadings)))], covdata = surveyFixModCOV, formula = formula, SE = SE, SE.type = SE.type, skipNominal = skipNominal, forceGRSM = forceGRSM, assumingFake = assumingFake, masterThesis = masterThesis, forceRasch = forceRasch, unstable = unstable, forceMHRM = forceMHRM, survey.weights = survey.weights, allowMixedResponse = allowMixedResponse, autofix = autofix, forceUIRT = forceUIRT, skipIdealPoint = skipIdealPoint, forceNRM = forceNRM, forceNormalEM = forceNormalEM, ...)
          surveyFixMod_Workout <- surveyFixMod_New
          
          if(printFactorStructureRealtime == T){
            message('\nRealtime Factor Structure after iteration')
            if(ncol(surveyFixMod_Workout@Fit$F)==1){
              print(round(surveyFixMod_Workout@Fit$F, 2))
            } else {
              if(bifactorSolution){
                print(round(GPArotation::bifactorQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings, 2))
              } else {
                print(round(GPArotation::geominQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings, 2))
              }
            }
          }
          
        }
      } else {
        fixFactorStructure_Done <- TRUE
      }
      
    } # the end of while loop
    
    if(printFactorStructureRealtime == T){
      message('\nFinal Factor Structure')
      if(ncol(surveyFixMod_Workout@Fit$F)==1){
        print(round(surveyFixMod_Workout@Fit$F, 2))
      } else {
        if(bifactorSolution){
          print(round(GPArotation::bifactorQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings, 2))
        } else {
          print(round(GPArotation::geominQ(surveyFixMod_Workout@Fit$F, maxit = 10000)$loadings, 2))
        }
      }
    }
    return(surveyFixMod_Workout)
  } else {
    return(surveyFixMod)
  }
  
}

fastMultipleGroup <- function(x, covdata = NULL, formula = NULL, SE = F, SE.type = "crossprod", skipNominal = T,
                              forceGRSM = F, assumingFake = F, masterThesis = F, forceRasch = F, unstable = F,
                              forceMHRM = F, itemkeys = NULL, survey.weights = NULL, group = ...,
                              invariance = c('free_means', 'free_var', colnames(x)), ...){
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
      } else if(length(survey.weights) != 0) {
        estimationMETHOD <- 'QMCEM'
        optimINPUT <- 'nlminb'
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
      
    } else { # with covariate
      
      covdataINPUT <- covdata
      formulaINPUT <- formula
      
      if(forceMHRM == T | forceGRSM == T | assumingFake == T | masterThesis == T){
        message('MHRM currently not supported with latent regressors')
        estimationMETHOD <- 'QMCEM'
        optimINPUT <- 'nlminb'
        optimCTRL  <- NULL  #list(control = list(trace = F))
      } else if(length(survey.weights) != 0) {
        estimationMETHOD <- 'QMCEM'
        optimINPUT <- 'nlminb'
        optimCTRL  <- NULL  #list(control = list(trace = F))
      } else if(i < 2){
        if(unstable == T){
          estimationMETHOD <- 'QMCEM'
          optimINPUT <- NULL
          optimCTRL <- NULL
        } else {
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
        if(estimationMETHOD == 'MHRM'){ # Richadson (BL) isn't support MHRM estimation method
          SE.type <- 'MHRM'
        } else {
          SE.type <- 'complete'
        }
      }
    }
    
    # standard error estimation config
    if((SE == T && SE.type == 'SEM') == T){
      accelerateINPUT <- 'none'
      TOLINPUT <- NULL
      SEtolINPUT <- 1e-4
      symmetric_SEMINPUT <- TRUE
      
      message('TOL: ', 'default', ' / SEtol: ', SEtolINPUT, ' / SE.type: ', SE.type, ' / Accelerator: ',
              accelerateINPUT, ' / Symmetric SEM: ', symmetric_SEMINPUT)
    } else if((SE == T && estimationMETHOD == 'MHRM') == T){
      accelerateINPUT <- 'squarem'
      TOLINPUT <- NULL
      SEtolINPUT <- 1e-3
      symmetric_SEMINPUT <- TRUE
      SE.type <- 'MHRM'
      
      message('TOL: ', 'default', ' / SEtol: ', SEtolINPUT, ' / SE.type: ', SE.type, ' / Accelerator: ',
              accelerateINPUT, ' / Symmetric SEM: ', symmetric_SEMINPUT)
    } else {
      accelerateINPUT <- 'squarem'
      TOLINPUT <- NULL
      SEtolINPUT <- 1e-10
      symmetric_SEMINPUT <- FALSE
      
      if(SE == T){
        message('Standard Error Estimation: On')
        message('TOL: ', 'default', ' / SEtol: ', SEtolINPUT, ' / SE.type: ', SE.type, ' / Accelerator: ',
                accelerateINPUT, ' / Symmetric SEM: ', symmetric_SEMINPUT)
      }
    }
    
    # removeEmptyRows config
    if(length(covdata) != 0){
      removeEmptyRowsConf <- FALSE
    } else {
      removeEmptyRowsConf <- TRUE
    }
    
    if(max(x, na.rm = T) - min(x, na.rm = T) == 1){ # dichotomous items
      
      if(nrow(x) >= 5000){
        message('\nMIRT model: Noncompensatory 4PL')
        try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = '4PL', method = estimationMETHOD,
                                           accelerate = accelerateINPUT, calcNull = T,
                                           technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                            removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                           formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                           SE.type = SE.type, survey.weights = survey.weights, group = group, invariance = invariance, ... = ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
        
        if(exists('modTEMP') == F){
          
          message('\nMIRT model: Noncompensatory 3PL')
          try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = '3PL', method = estimationMETHOD,
                                             accelerate = accelerateINPUT, calcNull = T,
                                             technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                              removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                             formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                             SE.type = SE.type, survey.weights = survey.weights, group = group, invariance = invariance, ... = ...), silent = T)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
        }
        
        if(exists('modTEMP') == F){
          
          message('\nMIRT model: Noncompensatory 3PL with lower or upper asymptote estimated')
          try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = '3PLu', method = estimationMETHOD,
                                             accelerate = accelerateINPUT, calcNull = T,
                                             technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                              removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                             formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                             SE.type = SE.type, survey.weights = survey.weights, group = group, invariance = invariance, ... = ...), silent = T)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
        }
        
        if(exists('modTEMP') == F){
          
          message('\nMIRT model: Partially compensatory 3PL')
          try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = 'PC3PL', method = estimationMETHOD,
                                             accelerate = accelerateINPUT, calcNull = T,
                                             technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                              removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                             formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                             SE.type = SE.type, survey.weights = survey.weights, group = group, invariance = invariance, ... = ...), silent = T)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
        }
      }
      
      if(exists('modTEMP') == F){
        
        message('\nMIRT model: ideal point')
        try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = 'ideal', method = estimationMETHOD,
                                           accelerate = accelerateINPUT, calcNull = T,
                                           technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                            removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                           formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                           SE.type = SE.type, survey.weights = survey.weights, group = group, invariance = invariance, ... = ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      if(exists('modTEMP') == F){
        
        message('\nMIRT model: Noncompensatory 2PL')
        try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = '2PL', method = estimationMETHOD,
                                           accelerate = accelerateINPUT, calcNull = T,
                                           technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                            removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                           formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                           SE.type = SE.type, survey.weights = survey.weights, group = group, invariance = invariance, ... = ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      if(exists('modTEMP') == F | modTEMP@OptimInfo$converged != 1){
        
        message('\nMIRT model: Partially compensatory 2PL')
        try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = 'PC2PL', method = estimationMETHOD,
                                           accelerate = accelerateINPUT, calcNull = T,
                                           technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                            removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                           formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                           SE.type = SE.type, survey.weights = survey.weights, group = group, invariance = invariance, ... = ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      if(exists('modTEMP') == F){
        if(i == 1){
          stop('Fail to find Factor solutions: Model didn\'t converge.')
        } else {
          return(modOLD)
        }
      }
      
    } else if(length(itemkeys) != 0){ # 2-4PLNRM when given item keys (nested logit)
      if(nrow(x) >= 5000){
        message('\nMIRT model: Noncompensatory 4PL Nominal response')
        try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = '4PLNRM', method = estimationMETHOD,
                                           key = itemkeys, accelerate = accelerateINPUT, calcNull = T,
                                           technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                            removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                           formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                           SE.type = SE.type, survey.weights = survey.weights, group = group, invariance = invariance, ... = ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
        
        if(exists('modTEMP') == F){
          
          message('\nMIRT model: Noncompensatory 3PL Nominal response')
          try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = '3PLNRM', method = estimationMETHOD,
                                             key = itemkeys, accelerate = accelerateINPUT, calcNull = T,
                                             technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                              removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                             formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                             SE.type = SE.type, survey.weights = survey.weights, group = group, invariance = invariance, ... = ...), silent = T)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
        }
        
        
        if(exists('modTEMP') == F){
          
          message('\nMIRT model: Noncompensatory 3PL Nominal response with lower or upper asymptote estimated')
          try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = '3PLuNRM', method = estimationMETHOD,
                                             key = itemkeys, accelerate = accelerateINPUT, calcNull = T,
                                             technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                              removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                             formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                             SE.type = SE.type, survey.weights = survey.weights, group = group, invariance = invariance, ... = ...), silent = T)
          if(exists('modTEMP')){
            if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
          }
        }
      }
      
      if(exists('modTEMP') == F){
        
        message('\nMIRT model: Noncompensatory 2PL Nominal response')
        try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = '2PLNRM', method = estimationMETHOD,
                                           key = itemkeys, accelerate = accelerateINPUT, calcNull = T,
                                           technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                            removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                           formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                           SE.type = SE.type, survey.weights = survey.weights, group = group, invariance = invariance, ... = ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      if(exists('modTEMP') == F){
        message('\nMIRT model: Nominal response without keys')
        try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = 'nominal', method = estimationMETHOD,
                                           accelerate = accelerateINPUT, calcNull = T,
                                           technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                            removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                           formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                           SE.type = SE.type, survey.weights = survey.weights, key = NULL, ...), silent = F)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      if(exists('modTEMP') == F){
        if(i == 1){
          stop('Fail to find Factor solutions: Model didn\'t converge.')
        } else {
          return(modOLD)
        }
      }
      
    } else { # polytomous items
      
      # forceRasch (PCM)
      if(forceRasch == T){
        message('\nMIRT model: Partial Credit')
        
        try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = 'Rasch', method = estimationMETHOD,
                                           accelerate = accelerateINPUT, calcNull = T,
                                           technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                            removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                           formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                           SE.type = SE.type, survey.weights = survey.weights, group = group, invariance = invariance, ... = ...), silent = F)
        try(return(modTEMP))
      }
      
      # forceGRSM (polytomous)
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
            
            try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = 'grsmIRT', method = estimationMETHOD,
                                               accelerate = accelerateINPUT, calcNull = T,
                                               technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                                removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                               formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                               SE.type = SE.type, survey.weights = survey.weights,
                                               group = group, invariance = invariance, ... = ...), silent = F)
            if(exists('modTEMP')){
              if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
            }
            
            if(exists('modTEMP') == F){
              if(i == 1){
                stop('Fail to find Factor solutions: Model didn\'t converge.')
              } else {
                return(modOLD)
              }
            }
            
          } else {
            
            try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = 'grsm', method = estimationMETHOD,
                                               accelerate = accelerateINPUT, calcNull = T,
                                               technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                                removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                               formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                               SE.type = SE.type, survey.weights = survey.weights,
                                               group = group, invariance = invariance, ... = ...), silent = F)
            if(exists('modTEMP')){
              if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
            }
            
            if(exists('modTEMP') == F){
              if(i == 1){
                stop('Fail to find Factor solutions: Model didn\'t converge.')
              } else {
                return(modOLD)
              }
            }
            
          }
          
          if(modTEMP@OptimInfo$converged == 1){
            skipNominal <- T
          }
        }
      }
      
      # polytomous
      # nominal model
      if(skipNominal == F){
        message('\nMIRT model: nominal response')
        try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = 'nominal', method = estimationMETHOD,
                                           accelerate = accelerateINPUT, calcNull = T,
                                           technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                            removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                           formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                           SE.type = SE.type, survey.weights = survey.weights,
                                           group = group, invariance = invariance, ... = ...), silent = F)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      # generalized partial credit model (non-sequential)
      if(exists('modTEMP') == F){
        message('\nMIRT model: Generalized partial credit')
        try(modTEMP <- mirt::multipleGroup(data = x, model = i, itemtype = 'gpcm', method = estimationMETHOD,
                                           accelerate = accelerateINPUT, calcNull = T,
                                           technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT, SEtol = SEtolINPUT,
                                                            removeEmptyRows = removeEmptyRowsConf), TOL = TOLINPUT, covdata = covdataINPUT,
                                           formula = formulaINPUT, optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                           SE.type = SE.type, survey.weights = survey.weights,
                                           group = group, invariance = invariance, ... = ...), silent = F)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      # graded response model (sequential)
      if(exists('modTEMP') == F){
        message('\nMIRT model: Graded response')
        try(modTEMP <- mirt::multipleGroup(data = x, model = i, method = estimationMETHOD, accelerate = accelerateINPUT,
                                           calcNull = T, technical = list(MAXQUAD = 2000000, MHRM_SE_draws = MHRM_SE_draws, symmetric_SEM = symmetric_SEMINPUT,
                                                                          SEtol = SEtolINPUT, removeEmptyRows = removeEmptyRowsConf),
                                           TOL = TOLINPUT, covdata = covdataINPUT, formula = formulaINPUT,
                                           optimizer = optimINPUT, solnp_args = optimCTRL, SE = SE,
                                           SE.type = SE.type, survey.weights = survey.weights,
                                           group = group, invariance = invariance, ... = ...), silent = T)
        if(exists('modTEMP')){
          if(modTEMP@OptimInfo$converged != 1){rm(modTEMP)}
        }
      }
      
      # finally, if can not converge
      if(exists('modTEMP') == F){
        if(i == 1){
          stop('Fail to find Factor solutions: Model didn\'t converge.')
        } else {
          return(modOLD)
        }
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

numericMI <- function(model = ..., data = ..., m = 100, fun = 'sem', estimator = 'MLMV', parameterization = 'delta', chi = 'mplus', ...) {
  
  #########################
  # Multiple              #
  # Imputation for        #
  # variables             #
  # in SEM context        #
  #########################
  # Seongho Bae           #
  # seongho@kw.ac.kr      #
  # May 6th 2016          #
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
  library('semTools')
  library('Amelia')
  library('lavaan')
  
  fit <- sem(model=model, data=data.frame(data)) # extract variable names in model syntax
  message("sample size (listwise): ", paste0(nrow(data.frame(fit@Data@X))))
  
  fit_MI <- runMI(model=model, data=data.frame(data[,attributes(fit)$Model@dimNames[[1]][[1]]]), m=m, fun = fun, estimator = estimator, chi = chi, ...)#, control=list(optim.method="L-BFGS-B"), ...)
  
  cat(summary(fit_MI, standardize=T))
  print(inspect(fit_MI, 'fit'))
  
  return(fit_MI)
  message('inspect(fit, "impute")')
  
}
cmvFA <- function(x, MHRM = F){
  initialModel <- surveyFA(x, autofix = F, forceMHRM = MHRM, bifactorSolution = T)
  workModel <- initialModel
  STOP <- FALSE
  
  # rotation
  while (!STOP) {
    try(invisible(gc()), silent = T)
    if(ncol(workModel@Fit$F) > 1){
      if(workModel@Model$nfact == 2){
        rotSumMat <- GPArotation::bifactorQ(workModel@Fit$F, maxit = 1e+6)$loadings[,2]
      } else {
        rotSumMat <- data.frame(GPArotation::bifactorQ(workModel@Fit$F, maxit = 1e+6)$loadings[,2:workModel@Model$nfact])
      }
      
      # evaluation
      evaluationMat <- abs(rotSumMat) < .1
      print(rotSumMat)
      print(evaluationMat)
      
      if(is.data.frame(rotSumMat)){
        evaluationMat <- rowSums(evaluationMat)
        rotSumMat <- rowSums(abs(rotSumMat))
      } else {
        evaluationMat <- as.numeric(evaluationMat)
        rotSumMat <- (abs(rotSumMat))
      }
      print(rotSumMat)
      print(evaluationMat)
      print(sum(evaluationMat))
      print(ncol(workModel@Data$data))
      
      if(sum(evaluationMat) == ncol(workModel@Data$data)){ # number of inappropreate items should equal with number of items in model
        message('Done')
        return(workModel)
      } else {
        message(paste0(colnames(workModel@Data$data[,which(abs(rotSumMat) == min(abs(rotSumMat)))])))
        workModel <- fastFIFA(workModel@Data$data[,-which(abs(rotSumMat) == min(abs(rotSumMat)))], forceMHRM = MHRM)
      }
      
    } else { # if one dimensional model
      rotSumMat <- workModel@Fit$F
      evaluationMat <- abs(rotSumMat) < .1
      print(rotSumMat)
      print(evaluationMat)
      evaluationMat <- as.numeric(evaluationMat)
      rotSumMat <- (abs(rotSumMat))
      if(sum(evaluationMat) == ncol(workModel@Data$data) | sum(evaluationMat) == 0 | ncol(workModel@Data$data) <= 5){ # number of inappropreate items should equal with number of items in model
        message('Done')
        return(workModel)
      } else {
        message(paste0(colnames(workModel@Data$data[,which(abs(rotSumMat) == min(abs(rotSumMat)))])))
        workModel <- fastFIFA(workModel@Data$data[,-which(abs(rotSumMat) == min(abs(rotSumMat)))], forceMHRM = MHRM)
      }
      
      # return(workModel)
    }
  }
  
}

bifactorFA <- function(data = ..., skipS_X2 = F, forceMHRM = F, covdata = NULL, formula = NULL, skipNominal = F, allowMixedResponse = T, itemkeys = NULL) {
  mod <- surveyFA(data = data, bifactorSolution = T, skipS_X2 = skipS_X2, forceMHRM = forceMHRM, autofix = F, covdata = covdata, formula = formula, skipNominal = skipNominal, allowMixedResponse = allowMixedResponse, itemkeys = itemkeys)
  STOP <- FALSE
  while (!STOP) {
    if(ncol(mod@Fit$F) == 1){
      rotMAT <- abs(mod@Fit$F)
    } else {
      rotMAT <- abs(GPArotation::bifactorQ(mod@Fit$F, maxit = 1e+5)$loadings)[,1]
    }
    
    print(rotMAT)
    
    if(sum(rotMAT < .999) != ncol(mod@Data$data)){
      mod <- surveyFA(data = mod@Data$data[,-which(rotMAT == max(rotMAT))], bifactorSolution = T, skipS_X2 = skipS_X2, forceMHRM = forceMHRM, autofix = F, covdata = covdata, formula = formula, skipNominal = skipNominal, allowMixedResponse = allowMixedResponse, itemkeys = itemkeys[-which(rotMAT == max(rotMAT))])
    } else if(sum(rotMAT > .1) != ncol(mod@Data$data)){
      mod <- surveyFA(data = mod@Data$data[,-which(rotMAT == min(rotMAT))], bifactorSolution = T, skipS_X2 = skipS_X2, forceMHRM = forceMHRM, autofix = F, covdata = covdata, formula = formula, skipNominal = skipNominal, allowMixedResponse = allowMixedResponse, itemkeys = itemkeys[-which(rotMAT == max(rotMAT))])
    } else {
      return(mod)
    }
    
  }
}


findMLCA <- function(data = ..., startN = 1, empiricalhist = F, group = NULL, empiricaloptimal = T){
  # try(invisible(gc()), silent = T)
  DICindices <- vector()
  j <- 0
  nfact <- vector()
  
  if(sum(psych::describe(data)$range == 0) == 0){
    workData <- data
  } else {
    workData <- data[,-which(psych::describe(data)$range == 0)]
  }
  message('starting find global optimal of latent class')
  for(i in startN:ncol(workData)){
    try(invisible(tempModel <- mirt::mdirt(data = workData, model = i, empiricalhist = empiricalhist, group = group, nruns = 100)), silent = F)
    if(exists('tempModel')){
      if(tempModel@OptimInfo$converged){
        message('class: ', i, ' / DIC: ', tempModel@Fit$DIC)
        j <- j + 1
        DICindices[j] <- tempModel@Fit$DIC
        nfact[j] <- i
        rm(tempModel)
      }
    }
  }
  bestModel <- which(min(DICindices) == DICindices)
  message('find global optimal: ', nfact[bestModel])
  
  if(empiricaloptimal){
    
    testFS <- fscores(mirt::mdirt(data = workData, model = nfact[bestModel], empiricalhist = empiricalhist, group = group, nruns = 100), QMC = T)
    
    membership <- vector()
    for(i in 1:nrow(testFS)){
      membership[i] <- which(testFS[i,] == max(testFS[i,]))
    }
    
    message('empirical optimal: ', max(membership))
    
    return(mirt::mdirt(data = workData, model = max(membership), empiricalhist = empiricalhist, group = group, nruns = 100))
  } else {
    return(mirt::mdirt(data = workData, model = nfact[bestModel], empiricalhist = empiricalhist, group = group, nruns = 100))
  }
}

doMLCA <- function(data = ..., startN = 1, empiricalhist = F, group = NULL){
  mirtCluster()
  if(is.data.frame(data) | is.matrix(data)){
    workModel <- findMLCA(data = data, startN = startN, empiricalhist = T, empiricaloptimal = F, group = group)
  } else {
    workModel <- data
  }
  initData <- workModel@Data$data
  workData <- workModel@Data$data
  
  STOP <- FALSE
  while(!STOP){
    
    # item fit evaluation
    workModelFit <- mirt::itemfit(workModel, impute = 100, QMC = T)
    FitSize <- workModelFit$S_X2/workModelFit$df.S_X2
    
    print(cbind(workModelFit, FitSize))
    
    if(sum(na.omit(workModelFit$S_X2[1:ncol(workData)] == "NaN")) != 0){
      workModel <- findMLCA(workData[,-which(workModelFit$S_X2[1:ncol(workData)] == "NaN")], empiricalhist = F, empiricaloptimal = T)
    } else if(sum(workModelFit$S_X2[1:ncol(workData)]/workModelFit$df.S_X2[1:ncol(workData)] > 10) != 0){
      workModel <- findMLCA(workData[,-which(max(workModelFit$S_X2[1:ncol(workData)]/workModelFit$df.S_X2[1:ncol(workData)]) == workModelFit$S_X2[1:ncol(workData)]/workModelFit$df.S_X2[1:ncol(workData)])], empiricalhist = F, empiricaloptimal = T)
    } else {
      STOP <- TRUE
    }
    rm(workData)
    workData <- initData[,colnames(workModel@Data$data)] # update work data
  }
  
  if(sum(is.na(initData)) != 0){
    workData <- mirt::imputeMissing(workModel, fscores(workModel, QMC = T))
    workModel <- findMLCA(workData, empiricalhist = F, empiricaloptimal = T)
  }
  
  print(plot(workModel, facet_items = FALSE))
  print(plot(workModel))
  
  fs <- fscores(workModel, QMC = T)
  
  class_prob <- data.frame(apply(fs, 1, function(x) sample(1:workModel@Model$nfact, 1, prob=x)))
  colnames(class_prob) <- "Class"
  mirtCluster(remove = T)
  
  return(class_prob)
}

deepFA <- function(mirtModel){ # for search more factors with prevent local optimal
  DICindices <- vector()
  DICindices[1] <- mirtModel@Fit$DIC
  
  
  j <- 1
  
  if(mirtModel@Options$method == 'EM' && length(attr(mirtModel@ParObjects$lrPars, 'formula')[[1]]) == 0) {
    method <- 'MHRM'
  } else if(mirtModel@Options$method == 'EM' && length(attr(mirtModel@ParObjects$lrPars, 'formula')[[1]]) != 0) {
    method <- 'QMCEM'
  } else {
    method <- mirtModel@Options$method
  }
  
  message('searching global optimal...')
  start <- mirtModel@Model$nfact + 1
  end <- mirtModel@Model$nfact + 10 # see http://www.tandfonline.com/doi/abs/10.1080/00273171.2012.710386
  
  nfact <- vector()
  nfact[1] <- mirtModel@Model$nfact
  for(i in start:end){
    try(invisible(gc()), silent = T)
    
    try(invisible(tempModel <- mirt::mirt(data = mirtModel@Data$data, model = i, itemtype = mirtModel@Model$itemtype, SE = mirtModel@Options$SE, SE.type = mirtModel@Options$SE.type, covdata = attr(mirtModel@ParObjects$lrPars, 'df'), formula = attr(mirtModel@ParObjects$lrPars, 'formula')[[1]], method = method, optimizer = mirtModel@Options$Moptim, accelerate = mirtModel@Options$accelerate, verbose = F, technical = list(NCYCLES = mirtModel@Options$technical$NCYCLES, MAXQUAD = mirtModel@Options$technical$MAXQUAD, SEtol = mirtModel@Options$technical$SEtol, symmetric_SEM = mirtModel@Options$technical$symmetric_SEM, removeEmptyRows = mirtModel@Options$technical$removeEmptyRows, MHRM_SE_draws = mirtModel@Options$technical$MHRM_SE_draws))), silent = T)
    if(exists('tempModel')){
      if(tempModel@OptimInfo$converged){
        message(i, ' factors were converged; DIC: ', tempModel@Fit$DIC)
        j <- j + 1
        DICindices[j] <- tempModel@Fit$DIC
        nfact[j] <- i
        rm(tempModel)
      }
    }
  }
  bestModel <- which(min(DICindices) == DICindices)
  message('find global optimal: ', nfact[bestModel])
  if(bestModel == 1){
    return(mirtModel)
  } else {
    return(mirt::mirt(data = mirtModel@Data$data, model = nfact[bestModel], itemtype = mirtModel@Model$itemtype, SE = mirtModel@Options$SE, SE.type = mirtModel@Options$SE.type, covdata = attr(mirtModel@ParObjects$lrPars, 'df'), formula = attr(mirtModel@ParObjects$lrPars, 'formula')[[1]], method = method, optimizer = mirtModel@Options$Moptim, accelerate = mirtModel@Options$accelerate, verbose = F, technical = list(NCYCLES = mirtModel@Options$technical$NCYCLES, MAXQUAD = mirtModel@Options$technical$MAXQUAD, SEtol = mirtModel@Options$technical$SEtol, symmetric_SEM = mirtModel@Options$technical$symmetric_SEM, removeEmptyRows = mirtModel@Options$technical$removeEmptyRows, MHRM_SE_draws = mirtModel@Options$technical$MHRM_SE_draws)))
  }
}
