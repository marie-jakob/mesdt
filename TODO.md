# TODOs


+ support for "||" notation -> check
+ support for mix of correlated and uncorrelated random effects
+ type 2 LRTs -> check
+ parallelization
+ pretty output
+ documentation
+ optionally test parameters instead of whole factors
+ parametric bootstrap
+ reduced random-effects structure
  + strategy 1: "keep it maximal"
  + strategy 2: Matuschek et al. 

+ optional other backends:
  + Julia MixedModels
  + glmmTMB  -> check

+ emmeans-like estimates + SEs

+ Plots

+ Check user input

+ Tests for customized SDT user output

### Open Questions

+ How to credit Henrik Singmann (a lot of things are similar to afex)
+ How to handle contrast coding? Automatically use sum contrasts and notify the user? 
+ How to handle overparametrized models (when lme4 drops columns)? (fixed-effect model matrix is rank deficient so dropping XX columns / coefficients)
+ How to handle reverse-coded trial_type variable? 
