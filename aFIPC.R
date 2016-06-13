source('https://raw.githubusercontent.com/seonghobae/k.aefa/master/k.aefa3.R')

autoFIPC <- function(newformXData = ..., oldformYData = ..., newformCommonItemNames = ..., oldformCommonItemNames = ..., itemtype = '3PL', newformBILOGprior = NULL, oldformBILOGprior = NULL, tryFitwholeNewItems = T, tryFitwholeOldItems = T, ...){
  
  # print credits
  message('automated Fixed Item Parameter Calibration: aFIPC 0.1')
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
    oldFormModel <- mirt::mirt(data = oldformYData, model = oldFormModelSyntax, itemtype = itemtype, SE = T, SE.type = 'complete', accelerate = 'squarem')
  } else { # try to search priors automatically. if it fail, try to bayesian approaches
    message('with estimate prior distribution using an empirical histogram approach. please be patient.')
    oldFormModel <- mirt::mirt(data = oldformYData, model = 1, itemtype = itemtype, SE = T, SE.type = 'complete', accelerate = 'squarem', empiricalhist = T, technical = list(NCYCLES = 1e+5), GenRandomPars = F)
    
  }
  
  if(tryFitwholeOldItems == T){
    if(!oldFormModel@OptimInfo$secondordertest){
      message('Estimation failed. estimating new parameters with no prior distribution using quasi-Monte Carlo EM estimation. please be patient.')
      
      try(rm(oldFormModel))
      try(oldFormModel <- mirt::mirt(data = oldformYData, 1, itemtype = itemtype, SE = T, SE.type = 'complete', method = 'QMCEM', accelerate = 'squarem', technical = list(NCYCLES = 1e+5), GenRandomPars = F))
    }
    
    if(!oldFormModel@OptimInfo$secondordertest){
      message('Estimation failed. estimating new parameters with no prior distribution using  Cai\'s (2010) Metropolis-Hastings Robbins-Monro (MHRM) algorithm. please be patient.')
      
      try(rm(oldFormModel))
      while (!exists('oldFormModel')) {
        try(oldFormModel <- mirt::mirt(data = oldformYData, 1, itemtype = itemtype, SE = T, method = 'MHRM', SE.type = 'MHRM', accelerate = 'squarem', technical = list(NCYCLES = 1e+5, MHRM_SE_draws = 200000), GenRandomPars = F))
      }
    }
  }
  
  if(!oldFormModel@OptimInfo$secondordertest){
    message('Estimation failed. trying to remove weird items by itemfit statistics')
    try(rm(oldFormModel))
    
    oldFormModel <- surveyFA(oldformYData, autofix = F, SE = T)
  }
  
  if(!oldFormModel@OptimInfo$secondordertest){
    stop('Estimation failed. Please check test quality.')
  }
  
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
    newFormModel <- mirt::mirt(data = newformXData, model = newFormModelSyntax, itemtype = itemtype, SE = T, SE.type = 'complete', accelerate = 'squarem')
  } else {
    message('with estimate prior distribution using an empirical histogram approach. please be patient.')
    newFormModel <- mirt::mirt(data = newformXData, 1, itemtype = itemtype, SE = T, SE.type = 'complete', empiricalhist = T, accelerate = 'squarem', technical = list(NCYCLES = 1e+5), GenRandomPars = F)
  }
  
  if (tryFitwholeNewItems) {
    
    if(!newFormModel@OptimInfo$secondordertest){
      message('Estimation failed. estimating new parameters with no prior distribution using quasi-Monte Carlo EM estimation. please be patient.')
      
      try(rm(newFormModel))
      try(newFormModel <- mirt::mirt(data = newformXData, 1, itemtype = itemtype, SE = T, method = 'QMCEM', SE.type = 'complete', accelerate = 'squarem', technical = list(NCYCLES = 1e+5), GenRandomPars = F))
    }
    
    if(!newFormModel@OptimInfo$secondordertest){
      message('Estimation failed. estimating new parameters with no prior distribution using  Cai\'s (2010) Metropolis-Hastings Robbins-Monro (MHRM) algorithm. please be patient.')
      
      try(rm(newFormModel))
      while (!exists('newFormModel')) {
        try(newFormModel <- mirt::mirt(data = newformXData, 1, itemtype = itemtype, SE = T, method = 'MHRM', SE.type = 'MHRM', accelerate = 'squarem', technical = list(NCYCLES = 1e+5, MHRM_SE_draws = 200000), GenRandomPars = F))
      }
    }
    
  }
  
  if(!newFormModel@OptimInfo$secondordertest){
    message('Estimation failed. trying to remove weird items by itemfit statistics')
    try(rm(newFormModel))
    
    newFormModel <- surveyFA(newformXData, autofix = F, SE = T)
  }
  
  if(!newFormModel@OptimInfo$secondordertest){
    stop('Estimation failed. Please check test quality.')
  }
  
  # do FIPC
  NewScaleParms <- mod2values(newFormModel)
  OldScaleParms <- mod2values(oldFormModel)
  
  for(i in 1:length(oldformCommonItemNames)){
    if((length(grep(paste0('^',newformCommonItemNames[i],'$'), colnames(newformXData[colnames(newFormModel@Data$data)]))) == 1) == TRUE && (length(grep(paste0('^',oldformCommonItemNames[i],'$'), colnames(oldformYData[colnames(oldFormModel@Data$data)]))) == 1) == TRUE){
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
  
  message('\nestimating Linked Form Eq(X) parameters')
  message('with estimate prior distribution using an Cai\'s (2010) Metropolis-Hastings Robbins-Monro (MHRM) approach. please be patient.')
  LinkedModelSyntax <- mirt::mirt.model(paste0('F1 = 1-',ncol(newformXData[colnames(newFormModel@Data$data)]),'\n',
                                               'MEAN = F1'))
  LinkedModel <- mirt::mirt(data = newformXData[colnames(newFormModel@Data$data)], LinkedModelSyntax, itemtype = newFormModel@Model$itemtype, SE = T, SE.type = 'MHRM', method = 'MHRM', accelerate = 'squarem', technical = list(NCYCLES = 1e+5), pars = NewScaleParms, GenRandomPars = F, technical = list(MHRM_SE_draws = 1e+5))
  
  # if(!LinkedModel@OptimInfo$secondordertest){
  #   message('Estimation failed. estimating new parameters with no prior distribution using quasi-Monte Carlo EM estimation. please be patient.')
  #   
  #   rm(LinkedModel)
  #   try(LinkedModel <- mirt::mirt(data = newformXData[colnames(newFormModel@Data$data)], LinkedModelSyntax, itemtype = newFormModel@Model$itemtype, SE = T, method = 'QMCEM', SE.type = 'complete', accelerate = 'squarem', technical = list(NCYCLES = 1e+5), pars = NewScaleParms, GenRandomPars = F))
  # }
  # 
  # if(!LinkedModel@OptimInfo$secondordertest){
  #   message('Estimation failed. estimating new parameters with no prior distribution using  Cai\'s (2010) Metropolis-Hastings Robbins-Monro (MHRM) algorithm. please be patient.')
  #   
  #   rm(LinkedModel)
  #   while (!exists('LinkedModel')) {
  #     try(LinkedModel <- mirt::mirt(data = newformXData[colnames(newFormModel@Data$data)], LinkedModelSyntax, itemtype = newFormModel@Model$itemtype, SE = T, method = 'MHRM', SE.type = 'MHRM', accelerate = 'squarem', technical = list(NCYCLES = 1e+5, MHRM_SE_draws = 200000), pars = NewScaleParms, GenRandomPars = T))
  #   }
  # }
  
  # if(!LinkedModel@OptimInfo$secondordertest){
  #   stop('Estimation failed. Please check test quality.')
  # }
  
  # calculate expected score
  ExpectedScoreOldform <- mirt::expected.test(x = oldFormModel, Theta = fscores(oldFormModel, method = 'MAP'))
  ExpectedScoreLinkedform <- mirt::expected.test(x = LinkedModel, Theta = fscores(LinkedModel, method = 'MAP'))
  ExpectedScoreNewform <- mirt::expected.test(x = newFormModel, Theta = fscores(newFormModel, method = 'MAP'))
  
  
  # save results as object
  modelReturn <- new.env()
  modelReturn$oldFormModel <- oldFormModel
  modelReturn$newFormModel <- newFormModel
  modelReturn$LinkedModel <- LinkedModel
  modelReturn$ExpectedScoreOldform <- ExpectedScoreOldform
  modelReturn$ExpectedScoreLinkedform <- ExpectedScoreLinkedform
  modelReturn$ExpectedScoreNewform <- ExpectedScoreNewform
  
  
  return(as.list(modelReturn))
}
