source('https://raw.githubusercontent.com/seonghobae/k.aefa/master/k.aefa3.R')

autoFIPC <- function(newformXData = ..., oldformYData = ..., newformCommonItemNames = ..., oldformCommonItemNames = ..., itemtype = '3PL', newformBILOGprior = NULL, oldformBILOGprior = NULL, tryFitwholeNewItems = T, tryFitwholeOldItems = T, checkIPD = T, tryEM = F, freeMEAN = T, forceNormalZeroOne = F, parameterOverwrite = F, ...){
  
  # print credits
  message('automated Fixed Item Parameter Calibration: aFIPC 0.2')
  message('Seongho Bae (seongho@kw.ac.kr)\n')
  
  # checking configure
  if(length(newformCommonItemNames) != length(oldformCommonItemNames)){
    stop('Common Items are not equal')
  }
  
  if(length(newformCommonItemNames) == 0 | length(oldformCommonItemNames) == 0){
    stop('Please provide common item names')
  }
  
  # checking common items
  message('Checking correspond common item names')
  to <- rep("<<<", length(newformCommonItemNames))
  print(data.frame(cbind(newformCommonItemNames, to, oldformCommonItemNames)))
  correspondItems <- data.frame(cbind(newformCommonItemNames, oldformCommonItemNames))
  
  
  checkCorrect <- function()
  { 
    n <- readline(prompt="Is it correct? (1: Yes 2: No) : ")
    if(!grepl("^[0-9]+$",n))
    {
      return(checkCorrect())
    }
    
    return(as.integer(n))
  }
  confirm <- checkCorrect()
  if(confirm != 1){
    stop('Please write down pairs correctly')
  }
  
  # estimate models for calibration
  if(!is.data.frame(oldformYData) && !is.matrix(oldformYData)){ # if Data is mirt model
    oldFormModel <- oldformYData
    oldformYDataK <- data.frame(oldFormModel@Data$data)
  } else { # if Data is data.frame
    oldformYDataK <- oldformYData
    if(itemtype == '3PL' && length(oldformBILOGprior) == 0){
      checkoldformBILOGprior <- function()
      { 
        n <- readline(prompt="Do you want to use default BILOG-MG priors for oldform Data? (1: Yes 2: No) : ")
        if(!grepl("^[0-9]+$",n))
        {
          return(readinteger())
        }
        
        return(as.integer(n))
      }
      oldformBILOGprior <- checkoldformBILOGprior()
      if(oldformBILOGprior == 1){
        oldformBILOGprior <- TRUE
      } else {
        oldformBILOGprior <- FALSE
      }
    }
    
    message('\nestimating oldForm (Y) parameters')
    if(itemtype == '3PL' && oldformBILOGprior == T){
      message('with traditional MMLE/EM approach')
      oldFormModelSyntax <- mirt::mirt.model(paste0('F1 = 1-',ncol(oldformYData),'\n',
                                                    'PRIOR = (1-',ncol(oldformYData),', a1, lnorm, 1, 1.6487), ', '(1-',ncol(oldformYData),', g, norm, .22, .08)'))
      try(oldFormModel <- mirt::mirt(data = oldformYData, model = oldFormModelSyntax, itemtype = itemtype, SE = T, accelerate = 'squarem'))
    } else { # try to search priors automatically. if it fail, try to bayesian approaches
      message('with estimate prior distribution using an empirical histogram approach. please be patient.')
      try(oldFormModel <- mirt::mirt(data = oldformYData, model = 1, itemtype = itemtype, SE = T, accelerate = 'squarem', empiricalhist = T, technical = list(NCYCLES = 1e+5), GenRandomPars = F))
      
    }
    
    if(tryFitwholeOldItems == T){
      if(!oldFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
        message('Estimation failed. estimating new parameters with no prior distribution using quasi-Monte Carlo EM estimation. please be patient.')
        
        try(rm(oldFormModel))
        try(oldFormModel <- mirt::mirt(data = oldformYDataK, 1, itemtype = itemtype, SE = T, method = 'QMCEM', accelerate = 'squarem', technical = list(NCYCLES = 1e+5), GenRandomPars = F))
      }
      
      if(!oldFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
        message('Estimation failed. estimating new parameters with no prior distribution using  Cai\'s (2010) Metropolis-Hastings Robbins-Monro (MHRM) algorithm. please be patient.')
        
        try(rm(oldFormModel))
        while (!exists('oldFormModel')) {
          try(oldFormModel <- mirt::mirt(data = oldformYDataK, 1, itemtype = itemtype, SE = T, method = 'MHRM', accelerate = 'squarem', technical = list(NCYCLES = 1e+5, MHRM_SE_draws = 200000), GenRandomPars = F))
        }
      }
    }
    
    if(!oldFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
      message('Estimation failed. trying to remove weird items by itemfit statistics')
      try(rm(oldFormModel))
      
      oldFormModel <- surveyFA(oldformYData, autofix = F, SE = T, forceUIRT = T)
    }
    
    if(!oldFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
      message('Estimation failed. trying to remove weird items by itemfit statistics by normal MMLE/EM')
      try(rm(oldFormModel))
      
      oldFormModel <- surveyFA(oldformYData, autofix = F, SE = T, forceUIRT = T, forceNormalEM = T)
    }
    
    if(!oldFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
      message('Estimation failed. trying to remove weird items by itemfit statistics by MMLE/QMCEM')
      try(rm(oldFormModel))
      
      oldFormModel <- surveyFA(oldformYData, autofix = F, SE = T, forceUIRT = T, unstable = T)
    }
    
    if(!oldFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
      message('Estimation failed. trying to remove weird items by itemfit statistics by MMLE/MHRM')
      try(rm(oldFormModel))
      
      oldFormModel <- surveyFA(oldformYData, autofix = F, SE = T, forceUIRT = T, forceMHRM = T)
    }
    
    if(!oldFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
      stop('Estimation failed. Please check test quality.')
    }
    
  }
  
  if(!is.data.frame(newformXData) && !is.matrix(newformXData)){ # if Data is mirt model
    newFormModel <- newformXData
    newformXDataK <- data.frame(newFormModel@Data$data)
  } else {
    newformXDataK <- newformXData
    if(itemtype == '3PL' && length(newformBILOGprior) == 0){
      checknewformBILOGprior <- function()
      { 
        n <- readline(prompt="Do you want to use default BILOG-MG priors for newform Data? (1: Yes 2: No) : ")
        if(!grepl("^[0-9]+$",n))
        {
          return(checknewformBILOGprior())
        }
        
        return(as.integer(n))
      }
      newformBILOGprior <- checknewformBILOGprior()
      if(newformBILOGprior == 1){
        newformBILOGprior <- TRUE
      } else {
        newformBILOGprior <- FALSE
      }
    }
    
    
    message('\nestimating newForm (X) parameters')
    if(itemtype == '3PL' && newformBILOGprior == T){
      message('with traditional MMLE/EM approach')
      newFormModelSyntax <- mirt::mirt.model(paste0('F1 = 1-',ncol(newformXData),'\n',
                                                    'PRIOR = (1-',ncol(newformXData),', a1, lnorm, 1, 1.6487), ', '(1-',ncol(newformXData),', g, norm, .22, .08)'))
      try(newFormModel <- mirt::mirt(data = newformXData, model = newFormModelSyntax, itemtype = itemtype, SE = T, accelerate = 'squarem'))
    } else {
      message('with estimate prior distribution using an empirical histogram approach. please be patient.')
      try(newFormModel <- mirt::mirt(data = newformXDataK, 1, itemtype = itemtype, SE = T, empiricalhist = T, accelerate = 'squarem', technical = list(NCYCLES = 1e+5), GenRandomPars = F))
    }
    
    if (tryFitwholeNewItems) {
      
      if(!newFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
        message('Estimation failed. estimating new parameters with no prior distribution using quasi-Monte Carlo EM estimation. please be patient.')
        
        try(rm(newFormModel))
        try(newFormModel <- mirt::mirt(data = newformXDataK, 1, itemtype = itemtype, SE = T, method = 'QMCEM', accelerate = 'squarem', technical = list(NCYCLES = 1e+5), GenRandomPars = F))
      }
      
      if(!newFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
        message('Estimation failed. estimating new parameters with no prior distribution using  Cai\'s (2010) Metropolis-Hastings Robbins-Monro (MHRM) algorithm. please be patient.')
        
        try(rm(newFormModel))
        while (!exists('newFormModel')) {
          try(newFormModel <- mirt::mirt(data = newformXDataK, 1, itemtype = itemtype, SE = T, method = 'MHRM', accelerate = 'squarem', technical = list(NCYCLES = 1e+5, MHRM_SE_draws = 200000), GenRandomPars = F))
        }
      }
      
    }
    
    if(!newFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
      message('Estimation failed. trying to remove weird items by itemfit statistics')
      try(rm(newFormModel))
      
      newFormModel <- surveyFA(newformXData, autofix = F, SE = T, forceUIRT = T)
    }
    
    if(!newFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
      message('Estimation failed. trying to remove weird items by itemfit statistics again by normal MMLE/EM')
      try(rm(newFormModel))
      
      newFormModel <- surveyFA(newformXData, autofix = F, SE = T, forceUIRT = T, forceNormalEM = T)
    }
    
    if(!newFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
      message('Estimation failed. trying to remove weird items by itemfit statistics again by MMLE/QMCEM')
      try(rm(newFormModel))
      
      newFormModel <- surveyFA(newformXData, autofix = F, SE = T, forceUIRT = T, unstable = T)
    }
    
    if(!newFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
      message('Estimation failed. trying to remove weird items by itemfit statistics again by MMLE/MHRM')
      try(rm(newFormModel))
      
      newFormModel <- surveyFA(newformXData, autofix = F, SE = T, forceUIRT = T, forceMHRM = T)
    }
    
    if(!newFormModel@OptimInfo$secondordertest && !itemtype == 'ideal'){
      stop('Estimation failed. Please check test quality.')
    }
  }
  
  
  # do FIPC
  NewScaleParms <- mirt::mod2values(newFormModel)
  OldScaleParms <- mirt::mod2values(oldFormModel)
  
  if(!parameterOverwrite){
    NewScaleParms[, "est"] <- TRUE
    OldScaleParms[, "est"] <- TRUE
  }
  
  NewScaleParms[which(NewScaleParms$item == paste0('GROUP')), "est"] <- FALSE
  OldScaleParms[which(OldScaleParms$item == paste0('GROUP')), "est"] <- FALSE
  
  NewScaleParms[which(NewScaleParms$name == "COV_11"), "est"] <- TRUE
  OldScaleParms[which(OldScaleParms$name == "COV_11"), "est"] <- TRUE
  
  if(itemtype == 'Rasch'){
    NewScaleParms[which(NewScaleParms$name == "a1"), "est"] <- FALSE
    OldScaleParms[which(OldScaleParms$name == "a1"), "est"] <- FALSE
  }
  
  #IPD
  if(checkIPD == T){
    # config
    IPDgroup <- as.factor(c(rep('oldForm', nrow(oldformYDataK)), rep('newForm', nrow(newformXDataK))))
    IPDItemCount <- 0
    IPDItemNamesOldForm <- vector()
    IPDItemNamesNewForm <- vector()
    
    # IPD target item checking
    for(i in 1:length(oldformCommonItemNames)){
      if((length(grep(paste0('^',newformCommonItemNames[i],'$'), colnames(newformXDataK[colnames(newFormModel@Data$data)]))) == 1) == TRUE && (length(grep(paste0('^',oldformCommonItemNames[i],'$'), colnames(oldformYDataK[colnames(oldFormModel@Data$data)]))) == 1) == TRUE){
        IPDItemCount <- IPDItemCount + 1
        IPDItemNamesOldForm[IPDItemCount] <- names(oldformYDataK[oldformCommonItemNames[i]])
        IPDItemNamesNewForm[IPDItemCount] <- names(newformXDataK[newformCommonItemNames[i]])
        
      } else {
        
      }
    }
    
    # IPD Data generation
    IPDItemList <- data.frame(rbind(IPDItemNamesOldForm, IPDItemNamesNewForm))
    
    IPDData <- data.frame(matrix(nrow = length(IPDgroup), ncol = IPDItemCount))
    colnames(IPDData) <- paste0('X', 1:IPDItemCount)
    print(IPDItemNamesOldForm)
    print(IPDItemNamesNewForm)
    IPDData[1:nrow(oldformYDataK),] <- oldformYDataK[,IPDItemNamesOldForm]
    IPDData[nrow(oldformYDataK)+1:nrow(newformXDataK),] <- newformXDataK[,IPDItemNamesNewForm]
    
    # IPD estimation
    IPDParmNames <- OldScaleParms$name
    IPDParmNames <- IPDParmNames[!duplicated(IPDParmNames)]
    IPDParmNames <- IPDParmNames[-c(grep("^MEAN", IPDParmNames), grep("^COV", IPDParmNames), grep("^ak", IPDParmNames), grep("^d0$", IPDParmNames))]
    IPDParmNames <- as.character(IPDParmNames)
    
    mirt::mirtCluster()
    message('Discovering IPD')
    if(itemtype == 'nominal' | tryEM == T){
      modIPD_MG <- multipleGroup(IPDData, model = 1, group = IPDgroup,
                                 itemtype = itemtype, method = 'EM', invariance = names(IPDData), empiricalhist = T, technical = list(NCYCLES = 1e+5, removeEmptyRows=TRUE))
      try(modIPD_DIF <- DIF(modIPD_MG, IPDParmNames, scheme = 'drop_sequential', method = 'EM', empiricalhist = T, technical = list(NCYCLES = 1e+5)))
    } else {
      modIPD_MG <- multipleGroup(IPDData, model = 1, group = IPDgroup,
                                 itemtype = itemtype, method = 'MHRM', invariance = names(IPDData), technical = list(NCYCLES = 1e+5, removeEmptyRows=TRUE))
      try(modIPD_DIF <- DIF(modIPD_MG, IPDParmNames, scheme = 'drop_sequential', method = 'MHRM', technical = list(NCYCLES = 1e+5)))
    }
    mirt::mirtCluster(remove = T)
    
    if(exists('modIPD_DIF')){
      
      modIPD_IPDItem <- names(modIPD_DIF)
      CommonItemList_NOIPD <- colnames(IPDData)[!colnames(IPDData) %in% modIPD_IPDItem]
      print(modIPD_DIF)
      print(CommonItemList_NOIPD)
      
      ActualoldFormCommonItem <- vector(length = length(CommonItemList_NOIPD))
      ActualnewFormCommonItem <- vector(length = length(CommonItemList_NOIPD))
      for(i in 1:length(CommonItemList_NOIPD)){
        ActualoldFormCommonItem[i] <- as.character(IPDItemList[CommonItemList_NOIPD][1,i])
        ActualnewFormCommonItem[i] <- as.character(IPDItemList[CommonItemList_NOIPD][2,i])
      }
      
      message('ActualoldFormCommonItem: ', ActualoldFormCommonItem)
      message('ActualnewFormCommonItem: ', ActualnewFormCommonItem)
      
      oldformCommonItemNames <- ActualoldFormCommonItem
      newformCommonItemNames <- ActualnewFormCommonItem
      message('oldformCommonItemNames: ', ActualoldFormCommonItem)
      message('newformCommonItemNames: ', ActualnewFormCommonItem)
    } else {
      message('No IPD detected')
    }
  }
  
  for(i in 1:length(oldformCommonItemNames)){
    
    if((length(grep(paste0('^',newformCommonItemNames[i],'$'), colnames(newformXDataK[colnames(newFormModel@Data$data)]))) == 1) == TRUE && (length(grep(paste0('^',oldformCommonItemNames[i],'$'), colnames(oldformYDataK[colnames(oldFormModel@Data$data)]))) == 1) == TRUE && (length(levels(as.factor(newFormModel@Data$data[,grep(paste0('^',newformCommonItemNames[i],'$'), colnames(newformXDataK[colnames(newFormModel@Data$data)]))]))) == length(levels(as.factor(oldFormModel@Data$data[,grep(paste0('^',oldformCommonItemNames[i],'$'), colnames(oldformYDataK[colnames(oldFormModel@Data$data)]))]))))){
      message('applying ', paste0(newformCommonItemNames[i]), ' <<< ', paste0(oldformCommonItemNames[i]), ' as common item use')
      
      message('   Newform Parms: ', paste0(NewScaleParms[which(NewScaleParms$item == paste0(newformCommonItemNames[i])), "value"], ' '))
      message('   Oldform Parms: ', paste0(OldScaleParms[which(OldScaleParms$item == paste0(oldformCommonItemNames[i])), "value"], ' '))
      
      NewScaleParms[which(NewScaleParms$item == paste0(newformCommonItemNames[i])), "value"] <- OldScaleParms[which(OldScaleParms$item == paste0(oldformCommonItemNames[i])), "value"]
      message('   Linkedform Parms: ', paste0(NewScaleParms[which(NewScaleParms$item == paste0(newformCommonItemNames[i])), "value"], ' '), '\n')
      
      NewScaleParms[which(NewScaleParms$item == paste0(newformCommonItemNames[i])), "est"] <- FALSE
      
    } else {
      message('skipping ', paste0(newformCommonItemNames[i]), ' <<< ', paste0(oldformCommonItemNames[i]), ' as common item use')
    }
  }
  
  if(length(attr(newFormModel@ParObjects$lrPars, 'parnum')) != 0 && length(attr(oldFormModel@ParObjects$lrPars, 'parnum')) != 0){
    NewScaleParms[which(NewScaleParms$item == paste0('BETA')), "value"] <- OldScaleParms[which(OldScaleParms$item == paste0('BETA')), "value"]
    NewScaleParms[which(NewScaleParms$item == paste0('BETA')), "est"] <- FALSE
    
    message('applying BETA parameter as linking')
    
    message('   Linkedform Parms: ', paste0(NewScaleParms[which(NewScaleParms$item == paste0('BETA')), "value"], ' '), '\n')
    betaFormula <- attr(newFormModel@ParObjects$lrPars, 'formula')[[1]]
    betaCOVdata <- attr(newFormModel@ParObjects$lrPars, 'df')
    betaSE <- FALSE
    betaEmpiricalhist <- FALSE
  } else {
    betaFormula <- NULL
    betaCOVdata <- NULL
    betaSE <- TRUE
    betaEmpiricalhist <- TRUE
    
  }
  
  message('\nestimating Linked Form Eq(X) parameters')
  if(forceNormalZeroOne){
    freeMEAN <- F
    
    NewScaleParms[which(NewScaleParms$name == "COV_11"), "est"] <- FALSE
    OldScaleParms[which(OldScaleParms$name == "COV_11"), "est"] <- FALSE
    NewScaleParms[which(NewScaleParms$name == "MEAN_11"), "est"] <- FALSE
    OldScaleParms[which(OldScaleParms$name == "MEAN_11"), "est"] <- FALSE
    NewScaleParms[which(NewScaleParms$name == "COV_11"), "value"] <- 1
    OldScaleParms[which(OldScaleParms$name == "MEAN_11"), "value"] <- 0
    
    
  }
  if(freeMEAN == T){
    LinkedModelSyntax <- mirt::mirt.model(paste0('F1 = 1-',ncol(newformXDataK[colnames(newFormModel@Data$data)]),'\n',
                                                 'MEAN = F1'))
    
    NewScaleParms[which(NewScaleParms$name == "MEAN_1"), "est"] <- TRUE
    OldScaleParms[which(OldScaleParms$name == "MEAN_1"), "est"] <- TRUE
    
  } else {
    LinkedModelSyntax <- mirt::mirt.model(paste0('F1 = 1-',ncol(newformXDataK[colnames(newFormModel@Data$data)]),'\n' ))
  }
  
  print(NewScaleParms)
  
  
  if(itemtype == 'nominal' | tryEM == T){
    if(betaEmpiricalhist){
      message('with MMLE/EM + empirical histogram approach. please be patient.')
      
    } else {
      message('with MMLE/EM approach. please be patient.')
      
    }
    if(sum(NewScaleParms$est) == 0){
      # LinkedModel <- oldFormModel
      
      LinkedModel <- mirt::mirt(data = newformXDataK[colnames(newFormModel@Data$data)], LinkedModelSyntax, itemtype = newFormModel@Model$itemtype, method = 'EM', SE = F, accelerate = 'squarem', empiricalhist = betaEmpiricalhist, technical = list(NCYCLES = 1e+6, SEtol = 1e-4, MHRM_SE_draws = 1e+5), pars = NewScaleParms, GenRandomPars = F, covdata = betaCOVdata, formula = betaFormula)
      
      
    } else {
      LinkedModel <- mirt::mirt(data = newformXDataK[colnames(newFormModel@Data$data)], LinkedModelSyntax, itemtype = newFormModel@Model$itemtype, method = 'EM', SE = betaSE, accelerate = 'squarem', empiricalhist = betaEmpiricalhist, technical = list(NCYCLES = 1e+6, SEtol = 1e-4, MHRM_SE_draws = 1e+5), pars = NewScaleParms, GenRandomPars = F, covdata = betaCOVdata, formula = betaFormula)
      
    }
    
  } else {
    message('with Cai\'s (2010) Metropolis-Hastings Robbins-Monro (MHRM) approach. please be patient.')
    
    if(sum(NewScaleParms$est) == 0){
      # LinkedModel <- oldFormModel
      LinkedModel <- mirt::mirt(data = newformXDataK[colnames(newFormModel@Data$data)], LinkedModelSyntax, itemtype = newFormModel@Model$itemtype, method = 'MHRM', SE = F, accelerate = 'squarem', TOL = .0005, technical = list(NCYCLES = 1e+6, SEtol = 1e-4, MHRM_SE_draws = 1e+5), pars = NewScaleParms, GenRandomPars = F, covdata = betaCOVdata, formula = betaFormula)
      
    } else {
      LinkedModel <- mirt::mirt(data = newformXDataK[colnames(newFormModel@Data$data)], LinkedModelSyntax, itemtype = newFormModel@Model$itemtype, method = 'MHRM', SE = betaSE, accelerate = 'squarem', TOL = .0005, technical = list(NCYCLES = 1e+6, SEtol = 1e-4, MHRM_SE_draws = 1e+5), pars = NewScaleParms, GenRandomPars = F, covdata = betaCOVdata, formula = betaFormula)
      
    }
    
  }
  
  # if(!LinkedModel@OptimInfo$secondordertest){
  #   message('Estimation failed. estimating new parameters with no prior distribution using quasi-Monte Carlo EM estimation. please be patient.')
  #   
  #   rm(LinkedModel)
  #   try(LinkedModel <- mirt::mirt(data = newformXDataK[colnames(newFormModel@Data$data)], LinkedModelSyntax, itemtype = newFormModel@Model$itemtype, SE = T, method = 'QMCEM', accelerate = 'squarem', technical = list(NCYCLES = 1e+5), pars = NewScaleParms, GenRandomPars = F))
  # }
  # 
  # if(!LinkedModel@OptimInfo$secondordertest){
  #   message('Estimation failed. estimating new parameters with no prior distribution using  Cai\'s (2010) Metropolis-Hastings Robbins-Monro (MHRM) algorithm. please be patient.')
  #   
  #   rm(LinkedModel)
  #   while (!exists('LinkedModel')) {
  #     try(LinkedModel <- mirt::mirt(data = newformXDataK[colnames(newFormModel@Data$data)], LinkedModelSyntax, itemtype = newFormModel@Model$itemtype, SE = T, method = 'MHRM', accelerate = 'squarem', technical = list(NCYCLES = 1e+5, MHRM_SE_draws = 200000), pars = NewScaleParms, GenRandomPars = T))
  #   }
  # }
  
  # if(!LinkedModel@OptimInfo$secondordertest){
  #   stop('Estimation failed. Please check test quality.')
  # }
  
  # calculate expected score
  ExpectedScoreOldform <- mirt::expected.test(x = oldFormModel, Theta = fscores(oldFormModel, method = 'MAP'))
  ExpectedScoreLinkedform <- mirt::expected.test(x = LinkedModel, Theta = fscores(LinkedModel, method = 'MAP'))
  ExpectedScoreNewform <- mirt::expected.test(x = newFormModel, Theta = fscores(newFormModel, method = 'MAP'))
  
  # calculate theta
  ThetaOldform <- fscores(oldFormModel, method = 'MAP')
  ThetaLinkedform <- fscores(LinkedModel, method = 'MAP')
  ThetaNewform <- fscores(newFormModel, method = 'MAP')
  
  # save results as object
  modelReturn <- new.env()
  modelReturn$oldFormModel <- oldFormModel
  modelReturn$newFormModel <- newFormModel
  modelReturn$LinkedModel <- LinkedModel
  modelReturn$ExpectedScoreOldform <- ExpectedScoreOldform
  modelReturn$ExpectedScoreLinkedform <- ExpectedScoreLinkedform
  modelReturn$ExpectedScoreNewform <- ExpectedScoreNewform
  modelReturn$ThetaOldform <- ThetaOldform
  modelReturn$ThetaNewform <- ThetaNewform
  modelReturn$ThetaLinkedform <- ThetaLinkedform
  if(checkIPD){
    modelReturn$IPDData <- data.frame(IPDData, IPDgroup)
    if(exists('CommonItemList_NOIPD')){
      modelReturn$IPDCommonItemList <- IPDItemList[CommonItemList_NOIPD]
    }
  }
  
  return(as.list(modelReturn))
}
