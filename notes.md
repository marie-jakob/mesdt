# Notes 



### Input Syntax

+ Use R formula stuff?


### GLMM Estimation in R (lme4) and Julia (MixedModels)

Julia: https://juliastats.org/MixedModels.jl/v4.13/optimization/#Generalized-Linear-Mixed-Effects-Models

lme4: https://www.rdocumentation.org/packages/lme4/versions/1.1-32/topics/glmer

+ fast = true option in Julia corresponds to nAQG = 0 in lme4!

+ "somewhat less accurate but way faster"

+ nAQG > 1 also possible but mostly infeasible

„The distinction between the "fast" and "slow" algorithms in the MixedModels package (nAGQ=0 or nAGQ=1 in lme4) is whether the fixed-effects parameters are optimized in PIRLS or in the nonlinear optimizer. In a call to the pirls! function the first argument is a GeneralizedLinearMixedModel, which is modified during the function call.“

+ could be used for selecting the random-effects structure, power analyses or in general simulations

+ final fit with fast = false or nAQG = 1

### Determining the Random-Effects Structure

+ Implement two strategies: 
	
	+ Barr et al. (2013): Maximal Model -> Start with maximal model, reduce sequentially and use the first model that fits the data

	+ Matuschek et al. (2017): Balancing power and Type I error rate -> determine maximal fitting model and reduce random-effects structure sequentially until there is a significantly worse model fit (determined by LRTs)

+ ::buildmer package could work for this -> Julia?


### Dealing with Predictors

+ Automatically use effect coding for categorical predictors?

+ automatically center continous predictors? 

+ -> TODO: check how afex does it

### SDT Parameter Estimates

+ for categorical predictors: transform GLM parameters (population-level and individual) to SDT

+ for continuous predictors: return SDT estimates for mean values of predictors?

### Confidence Intervals

##### Fixed Effects

+ -> Compute standard errors based on the Hessian -> returned by lme4 through vcov(model_fit)

+ transformation to SEs for SDT parameters through addition of Varianzadditionssatz -> account for all entries of the Fisher information matrix (covariances!)

+ TODO: only for the fixed parameters or also for the individual parameters? 


##### Random Effects (Individual Parameters)

+ no CIs for individual parameters? probably




