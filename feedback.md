# Feedback MathPsych

### Why 3 backends?
+ Idea: different trade-offs of ease of installation and estimation speed
+ still a valid point -> I'd have to test how much the estimates differ between the different backends
+ Maybe only use lme4 if this is fast enough?


### What about convergence issues?
+ Why / how does the approach then still make sense? 
+ -> Convergence issues often indicate that the random-effects structure is too complex to be supported by the data -> variance and correlation estimates close to 0
+ -> Reduction should not alter the results 
+ Bayesian models would still be overparametrized (but constrained by the priors)

### Why not Bayesian?
+ -> find good arguments for this for the manuscript later
+ speed of estimation (and for model comparisons)
+ statistical inference through LRTs and parametric bootstrapping -> debate about appropriate strategies for Bayesian model selection
+ Bayesian model comparison might also work here with Savage-Dickey ratios (of the fixed-effects estimates)

