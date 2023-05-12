# Notes 



### Input Syntax

+ Use R formula stuff? -> afex: uses 3 different syntax variants for the ANOVA (might also be an option)

+ basic variant: specify formula for mu, lambda and a random term

+ input is used to construct two model matrices (for lambda and mu) and a random term which are then given to lme4::glmer()

	+ allows to estimate interactions without the presence of main effects in the model (which is necessary for predictors affecting sensitivity and not response bias and the LRTs)



### GLMM Estimation with different backends


__Julia__: https://juliastats.org/MixedModels.jl/v4.13/optimization/#Generalized-Linear-Mixed-Effects-Models

__lme4__: https://www.rdocumentation.org/packages/lme4/versions/1.1-32/topics/glmer

+ fast = true option in Julia corresponds to nAQG = 0 in lme4!

+ "somewhat less accurate but way faster"

+ nAQG > 1 also possible but mostly infeasible

„The distinction between the "fast" and "slow" algorithms in the MixedModels package (nAGQ=0 or nAGQ=1 in lme4) is whether the fixed-effects parameters are optimized in PIRLS or in the nonlinear optimizer. In a call to the pirls! function the first argument is a GeneralizedLinearMixedModel, which is modified during the function call.“

+ could be used for selecting the random-effects structure, power analyses or in general simulations

+ final fit with fast = false or nAQG = 1


__glmmTMB__ package: 

+ New R package using automatic differentiation to estimate GLMMs

+ may be faster than lme4 (but has less functionality at the moment)

+ https://glmmtmb.github.io/glmmTMB__glmmTMB__ package: 

+ New R package using automatic differentiation to estimate GLMMs

+ may be faster than lme4 (but has less functionality at the moment)

+ https://glmmtmb.github.io/glmmTMB/

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

+ different variants for CIs: Wald CIs, Profile CIs, Parametric Bootstrapped CIs -> all possible in lme4

+ Wald CIs: based on the inverse Hessian

	+ very easy to compute from lme4 output (see paper)

	+ relatively strong assumption of a quadratic likelihood surface

+ Likelihood Profile CIs: ???

	+ more time-consuming but less severe assumptions (something with the asymptotic distribution of the deviance)

+ Parametric Bootstrap CIs: take a lot of time but avoid all asymptotic assumptions 

	+ not really feasible in practice 


+ -> Compute standard errors based on the Hessian -> returned by lme4 through vcov(model_fit)

	+ might not be accurate because assumptionsa are not met --> good enough?
	+ lme4 also offers profile CIs but those take forever to compute and it might be difficult to transform these to the SDT parameter space

+ transformation to SEs for SDT parameters through addition of Varianzadditionssatz -> account for all entries of the Fisher information matrix (covariances!)

+ TODO: only for the fixed parameters or also for the individual parameters? 


##### Random Effects (Individual Parameters)

+ no CIs for individual parameters? probably




