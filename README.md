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
