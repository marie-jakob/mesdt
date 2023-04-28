# mlsdt

Fit Multi-Level Signal Detection Theory Models Using Maximum Likelihood Estimation. 


### GLMM Estimation in R (lme4) and Julia (MixedModels)

Julia: https://juliastats.org/MixedModels.jl/v4.13/optimization/#Generalized-Linear-Mixed-Effects-Models

lme4: https://www.rdocumentation.org/packages/lme4/versions/1.1-32/topics/glmer

+ fast = true option in Julia corresponds to nAQG = 0 in lme4!

+ "somewhat less accurate but way faster"

+ nAQG > 1 also possible but mostly infeasible

„The distinction between the "fast" and "slow" algorithms in the MixedModels package (nAGQ=0 or nAGQ=1 in lme4) is whether the fixed-effects parameters are optimized in PIRLS or in the nonlinear optimizer. In a call to the pirls! function the first argument is a GeneralizedLinearMixedModel, which is modified during the function call.“

+ could be used for selecting the random-effects structure, power analyses or in general simulations

+ final fit with fast = false or nAQG = 1
