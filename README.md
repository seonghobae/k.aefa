# k.aefa
Kwangwoon Automated Exploratory Factor Anlaysis (k.aefa): automatically find optimal number of factors with various exploratory factoring methods since September, 2013. in 2016, it turns nearby 3!

## What can I do using this source code?
Just find optimal factor numbers in context of exploratory factor analysis. With some functions, they have some useful toys for psychologists like aberrant data point detection and delection with fully automated by statistical criteria, and *FULLY AUTOMATION* of Exploratory Factor Analysis what include check violence of item interdependence assumption in IRT and SEM, automatically variable delection after automatically model estimations.

### Supported factor analytical strategies
#### Classical True-score Theory by model based approach (Traditional Structural Equation Model)
It supports traditional Classical True-score Theory based model approaches by dependes correlation only; ML, minres, wls, gls (from ```fa()``` in psych and ```efaUnrotate()``` in semTools), ADF (from ```FAiR```).
#### Morden True-socre Theory by model based approach (Multidimensional IRT or Categorical Structural Equation Model)
The other one is Morden True-score Theory based model approaches; (Non-Compensatory) Full-information Exploratory IFA (from ```mirt()``` in mirt and ```tam.fa()``` in TAM), Limited-information exploratory IFA using WLSMV estimation in Mplus7 (from ```efaUnrotate()``` in semTools). 

## Why are you made this project?
This project made and updated for personal convinence during survey data analysis using gathering Likert rating scale style survey forms since September 2013. but some imported packages requires the GPL 2 or GPL 3 licence. So I decided open this source code even skeleton code status; not a library. I'll improve this just a source code to library.

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
