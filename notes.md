# Notes 


### Input Syntax

+ Use R formula stuff? -> afex: uses 3 different syntax variants for the ANOVA (might also be an option)

+ basic variant: specify formulas for mu and lambda including fixed and random parts

+ input is used to construct two model matrices (for lambda and mu) for fixed and for random effects respectively which are then given to lme4::glmer()
	+ ("regular" lm() model matrices are used for the random terms, not the real lme4 random effects model matrices as these are pretty complicated and seem easy to break)

	+ allows to estimate interactions without the presence of main effects in the model (which is necessary for predictors affecting sensitivity and not response bias and the LRTs)
	
+ "||" notation in input allowed if correlations are suppressed for both mu and lambda
	+ works for factors as well (even though it does not work in lme4) because data are entered numerically via the model matrices (problem in lme4 is that factors with > 2 levels are not split into separate terms, afex allows this via expand_re)
	+ "||" does weird stuff with matrices -> columns of a given model matrix are always correlated in the model (e.g., `~ ... + (m_matrix || ID)` -> `~ ... + (m_matrix[, 1] + m_matrix[, 2] + ... || ID)`


### GLMM Estimation with different backends

#### Estimation
+ (from Julia MixedModels documentation):
+ not the exact deviance, but the Laplace approximation to the deviance is optimized
+ Laplace approximation:
  + "the Laplace approximation entails finding a Gaussian approximation to a continuous probability density."
  + used in GLMM to approximate the integral given through random effects
+ optimization w.r.t. the fixed effects coefficients $\beta$ and the covariance parameters $\theta$
+ 


#### Backends

__Julia__: https://juliastats.org/MixedModels.jl/v4.13/optimization/#Generalized-Linear-Mixed-Effects-Models

+ integration of Julia MixedModels into R: https://github.com/mikabr/jglmm

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

+ afex: automatically sets effect coding and gives a warning if continuous predictors are not centered

### SDT Parameter Estimates

+ standard regression output but separated for sensitivity and response bias

+ SDT estimates per category: similar to emmeans package (or in combination with the package)

+ for categorical predictors: transform GLM parameters (population-level and individual) to SDT

+ for continuous predictors: return SDT estimates for mean values of predictors?

### Confidence Intervals

##### Fixed Effects

+ different variants for CIs: Wald CIs, Profile CIs, Parametric Bootstrapped CIs -> all possible in lme4

+ Wald CIs: based on the inverse Hessian (observed Fisher information)

	+ very easy to compute from lme4 output (see paper) -> variance addition of the regression parameters, including the nondiagonal elements for the covariance terms

	+ relatively strong assumption of a quadratic likelihood surface -> adequacy can be diagnosed via likelihood profiles (see below)

+ Likelihood Profile CIs: 

	+ gauge shape of the likelihood for a specific parameter of interest to compute CIs on that basis -> possible based on one fit of the model

	+ more time-consuming but less severe assumptions (something with the asymptotic distribution of the deviance) -> quadratic shape of log likelihood surface is not required

	+ not trivial for SDT since the SDT parameters are not part of the likelihood function -> shape of the likelihood function "in multiple directions" needed (analogous to nondiagonal elements of the Hessian)

	+ profile zeta plots (see lme4 paper) can be used to diagnose to what extent the assumption of a quadratic likelihood surface holds -> if this holds, Wald CIs are fine

+ Parametric Bootstrap CIs: take a lot of time but avoid all asymptotic assumptions 

	+ not really feasible in practice 

+ TODO: only for the fixed parameters or also for the individual parameters? 


##### Random Effects (Individual Parameters)

+ no CIs for individual parameters? probably


### Prediction & Simulation

+ predict() function would be nice -> could run into problems with the way the model matrices are generated. Adding the variables from the model matrices to the data might help

+ simulating from known parameter values -> basis for future simulation-based power analysis functionality

## Resources


+ lme4 estimation:

	+ http://lme4.r-forge.r-project.org/slides/2011-03-16-Amsterdam/

	+ https://www.alexejgossmann.com/Dissect_lmer_part1/

+ profile likelihood stuff

	+ https://stats.stackexchange.com/questions/27976/how-can-i-find-the-standard-deviation-of-the-sample-standard-deviation-from-a-no/27981#27981

	+ https://stats.stackexchange.com/questions/28671/what-is-the-exact-definition-of-profile-likelihood



