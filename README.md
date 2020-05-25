# Joint model of count and ordered variables
 R codes for random parameters negative binomial fractional split ordered logit (probit) regression model

This code has been produced as a part of my doctoral research. Please cite the following article if you use this code in any kind:

Afghari, A.P., 2019. Detecting motor vehicle crash blackspots based on their underlying behavioural, engineering, and spatial causes (Doctoral dissertation, University of Queensland).

This R code corresponds with the recently developed advanced joint econometric model of crash count and crash severity. The model is specified with Random Parameters to accout for unobserved heterogeneity in data. In the existing code, the model only has an observed correlation term between the count model and the ordered model. However, the same logic may be followed to define an unobserved correlation term between the two models.

Although this code has been written for crash count and severity data, it may be used for any dependent variables that are of count and ordered nature (e.g. number of patients at hospitals and severity of disease).

Prior to running the code, make sure to have "maxLik" and "randtoolbox" packages installed.

