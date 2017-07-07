# k.aefa
Kwangwoon Automated Exploratory Factor Analysis (k.aefa): automatically find the optimal number of factors with different exploratory factoring methods since September 2013. In 2016, it turns nearby 3!

## What can do using this source code?
Just find optimal factor numbers in the context of exploratory factor analysis. With some functions, they have some useful toys for psychologists like aberrant data point detection and deletion with fully automated by statistical criteria and FULLY AUTOMATION of Exploratory Factor Analysis what include checking the violence of item interdependence assumption in IRT and SEM, automatically variable deletion after automatically model estimations.

### Supported factor analytical strategies
Currently, Morden True-score Theory based model approaches; Compensatory Full-information Exploratory IFA (from mirt() in mirt). Exploratory Bifactor model, Exploratory n-dimensional model are available to use.

## Why did this project?
This project built and updated for personal conveniences during survey data analysis since September 2013. However, some imported packages require the GPL 2 or GPL 3 licence. So decided opens this source code even skeleton code status; not a library. Will improve this just a source code to the library.

## How can use it?
```R
  # load source code
  library(RCurl)
  destfile = "k.aefa3.R"
  x = getBinaryURL("https://raw.githubusercontent.com/seonghobae/k.aefa/master/k.aefa3.R",
                    followlocation = TRUE, ssl.verifypeer = FALSE)
  writeBin(x, destfile, useBytes = TRUE)
  source(paste("k.aefa3.R", sep = ""))
  rm(destfile, x)
  
  # find optimal factor numbers
  mod1 <- fastFIFA(your_data_frame)
  
  # doing fully automated exploratory factor analysis
  mod2 <- surveyFA(your_data_frame)
  
  # find optimal factor numbers with person covariates (latent regression of fixed effects)
  mod3 <- fastFIFA(your_data_frame, covdata = your_demographic_data_frame,
                  formula = ~1 + your + variable + names + in + demographic + data + frame)
  
  # doing fully automated exploratory factor analysis with person covariates (latent regression of fixed effects)
  mod4 <- surveyFA(your_data_frame, covdata = your_demographic_data_frame,
                  formula = ~1 + your + variable + names + in + demographic + data + frame)
                  
  # doing fully automated exploratory factor analysis with survey weights
  mod5 <- surveyFA(your_data_frame, survey.weights = your_survey_weights_where_get_from_Finite_population)
```
