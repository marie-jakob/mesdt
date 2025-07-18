# Notes 

### TODOs

Functionality / Technical stuff:
(+ gumbel-max + tests)
+ Aggregation
+ tests continuous predictors
+ tests bootstrap with many samples / compare with afex
+ gumbel-max

Post-Processing / Pretty Stuff
+ apa_print_mlsdt() method
  + can have a dependency for papaja
+ pretty output
  + correct names for the tests
+ df for random effects
+ progress bars PB
+ reformulas

Documentations / Conventions / Naming / Practicality
+ (Tests for customized SDT user output)
+ Vignette extreme value SDT
+ set sum contrasts
+ center continuous predictors


Maybe: 
+ optionally test parameters instead of whole factors
+ support for mix of correlated and uncorrelated random effects
  + That should also be easier when the translation of the input syntax is changed
+ Plots
+ test only a subset of random effects
+ tests random effects with random slopes with more than two levels

Probably not: 
+ reduced random-effects structure
  + strategy 1: "keep it maximal"
  + strategy 2: Matuschek et al. 

Done: 
+ Model translation
+ Fitting
+ support for "||" notation
+ crossed random effects
+ LRTs: type 3 and 2
+ parametric bootstrap (lme4 & glmmTMB)
+ test random effects for significance (type 3)
+ test only a selected subset of fixed effects
+ optional other backends:
  + glmmTMB
+ parallelization
+ mesdt model class
+ renaming stuff
+ simulate() method
+ emmeans method
+ input checks fit_mesdt()
+ change mu / lambda to c / d' -> or just use "sensitivity" and "response bias"  -> discriminability
+ different distributions: gumbel-min link function (gumbel max?)
  + -> adapt to do the correct link stuff (reverse-coding response variable or something
  + -> see DeCarlo 1998 paper
+ skip tests if packages not installed
  + summary() (+ print.summary()) method
    + -> print random effects
  + print() method
  + print() etc. method for the tests
+ README
+ datasets: detecting bias project?
  + or use simulated data for a toy question?
+ Fix random seed thing for bootstrap
+ documentation
+ examples

### Questions Co-Authors

+ naming conventions for SDT parameters - what to use here? 
+ df for random effects tests
+ leave random effects out for now?
+ How to handle contrast coding? Automatically use sum contrasts and notify the user? 
+ default for compute_tests() test argument



### Open Questions

+ How to credit Henrik Singmann (a lot of things are similar to afex + he helped)
+ How to handle overparametrized models (when lme4 drops columns)? (fixed-effect model matrix is rank deficient so dropping XX columns / coefficients)
+ How to handle reverse-coded trial_type variable? -> warning to user


### Open Questions - Raphi & Constantin

+ Remove correlations for all random effects grouping factors at once or sequentially or depending on there being correlations = 1 or = -1? 
+ Crossed random effects: do two separate selections and combine (and then do another selection)?
+ ANOVA stuff


### Input Syntax

+ Use R formula stuff? -> afex: uses 3 different syntax variants for the ANOVA (might also be an option)

+ basic variant: specify formulas for mu and lambda including fixed and random parts

+ input is used to construct two model matrices (for lambda and mu) for fixed and for random effects respectively which are then given to lme4::glmer()
	+ ("regular" lm() model matrices are used for the random terms, not the real lme4 random effects model matrices as these are pretty complicated and seem easy to break)

	+ allows to estimate interactions without the presence of main effects in the model (which is necessary for predictors affecting sensitivity and not response bias and the LRTs)
	
+ "||" notation in input allowed if correlations are suppressed for both mu and lambda
	+ works for factors as well (even though it does not work in lme4) because data are entered numerically via the model matrices (problem in lme4 is that factors with > 2 levels are not split into separate terms, afex allows this via expand_re)
	+ "||" does weird stuff with matrices -> columns of a given model matrix are always correlated in the model (e.g., `~ ... + (m_matrix || ID)` -> `~ ... + (m_matrix[, 1] + m_matrix[, 2] + ... || ID)`)


### GLMM Estimation with different backends

#### Estimation
+ (from Julia MixedModels documentation):
+ not the exact deviance, but the Laplace approximation to the deviance is optimized
+ Laplace approximation:
  + "the Laplace approximation entails finding a Gaussian approximation to a continuous probability density."
  + used in GLMM to approximate the integral given through random effects
+ optimization w.r.t. the fixed effects coefficients $\beta$ and the covariance parameters $\theta$



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

+ may be faster than lme4 (but has less functionality at the moment) -> definitely faster by a lot

+ https://glmmtmb.github.io/glmmTMB__glmmTMB__ package: 


__glmmADMB__ package:

+ also uses automatic differentiation

+ older than glmmTMB and does not seem to be actively supported anymore

### Determining the Random-Effects Structure

+ Implement two strategies: 
	
	+ Barr et al. (2013): Maximal Model -> Start with maximal model, reduce sequentially and use the first model that fits the data -> adapted to SDT models (do this "separately" for Mu and Lambda -> see below)

	+ Matuschek et al. (2017): Balancing power and Type I error rate -> determine maximal fitting model and reduce random-effects structure sequentially until there is a significantly worse model fit (determined by LRTs)
	+ -> simplify this by (1) finding the maximal model (in a backward selection), (2) test all remaining random effects for significance and (3) remove all non-significant predictors (optionally (4) do this as long as there are nonsignificant random effects in the model)

+ ::buildmer package is probably not that easy to include -> do this manually

+ Constantin: adhere to marginality in the random-effects selection (wegen Reviewer)
	+ optimal: optional forward & backward selection implementieren + testing strategies (e.g., LRT, AIC)
	+ -> Empfehlung: forward selection for complex models, backward selection for less complex models

##### "Keep it maximal"

+ User gives the maximum random-effects structure justified by the design

+ estimate the full model with nAGQ = 0 and terms are removed until the model converges to a non-singular fit
	+ for removal, adhere to principle of marginality _within_ sensitivity and response bias
	+ first, remove correlations
	+ then: check variance of slopes for the highest-order predictors on sensitivity and response bias and remove the slope with the lowest estimated variance from the model
	+ if the model still does not converge when only random intercepts for sensitivity and response bias are in the model, warn the user accordingly and suggest a different modeling approach (e.g., fitting on aggregated data) -> Still fit the model
+ Store the formula and fit the converging model again with nAGQ = 1
+ if the model does not converge, increase the iterations iteratively until `max_iter` is exceeded
+ if the model does not converge by then -> reduce random-effects structure further
	+ again: if the random intercept for either sensitivity or response bias would be removed in the next step, the user is warned about this
+ Return: 
	+ the final estimated model
	+ Table with sensitivity and response bias estimates
	+ the final model formula
	+ some string in nice formatting explaining the model?


##### lme4
+ 
+ check if model is singular: `isSingular()`
+ check if model converged: https://search.r-project.org/CRAN/refmans/performance/html/check_convergence.html ?
+ or simply check the messages: https://stackoverflow.com/questions/64114841/capturing-convergence-message-from-lme4-package-in-r


### Dealing with Predictors

+ Automatically use effect coding for categorical predictors?

+ automatically center continuous predictors? 

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


### Significance Testing

##### Fixed-Effects Hypothesis Tests

+ the standard test --> "Does my predictor influence the criterion"? 

+ likelihood ratio tests or parametric bootstrapping -> parametric bootstrapping can take ages

+ entails fitting a submodel for every to-be-tested predictor

+ Type III SS: fuck marginality, only remove the to-be-tested fixed effect from the model

+ Type II-like SS for this case: adhere to marginality _within_ sensitivity and response bias, respectively
	+ i.e., a submodel testing predictors for response bias still includes everything related to sensitivity, but includes only same- or lower-order effects of the predictor in question on the response bias
	+ Models to fit: 
 	+ For every order in the model:
    	+ Fit model with all effects with <= order (i.e., a model with all main effects, a model with all main effects + two-way interactions etc.)
    	+ Reference model for the respective order
    	+ For every effect of that order:
      	+ Fit a model with all effects with <= order _without_ the to-be-tested effect
      	+ Compare that model with the reference model

+ TODO: balanced vs. strict null stuff -> I don't think this "strict null" stuff makes a lot of sense -> Instead there will be the possibility to test for variability in the fixed effect separately


##### Random-Effects Tests

+ test variability in fixed effects -> LRTs or parametric bootstrapping for random slopes
+ needed for the Bates et al. (2015) / Matuschek et al. (2017) strategy (at least the one based on the LRT, the other strategies simply compare all models via AIC)

+ Stoel et al. (2006): LRTs for parameters with boundary constraints are tested with mixtures of Chi^2 distributions

+ Scheipl et al. (2008): propose a number of tests based on REML estimation, which is (to my knowledge) not possible for GLMM

+ LRTs: dfs are difficult since we are at the edge of the parameter space (i.e., var = 0) -> sampling distribution is a mixture of chi^2 distributions
	+ Davis-Stober (2009) references Chant, 1974; Chernoff, 1954; Robertson, Wright, & Dykstra, 1988; Shapiro, 1985; Wilks, 1938
	+ formula for df = 1
	+ from Thorstens foundations lectures: If the models differ by one variance and _k_ covariance parameters, the reference distribution is 0.5 * chi^2(1 + _k_) + 0.5 * chi^2(1 + _k_)
	+ 50:50 mixture

+ without adhering to marginality: Type III-like (remove exactly the to-be-tested random effect from the model)
+ does this make sense for type II? -> yes, but similar to the fixed effects, new full models need to be estimated
+ crossed random effects: for both type II and type III, all random effects for a different grouping factor are left in the model

+ implemented at the moment: type III like (since type II is a ton of work)
	+ for correlated random effects, every term is removed together with its corresponding covariance terms (-> dfs > 1)

##### Parametric Bootstrapping

Procedure:
+ fit the full and reduced model
+ simulate data from the reduced model and fit both the full and the reduced model to the simulated datasets
+ difference in -2LL between the full and reduced model on the simulated datasets provides a reference distribution for distribution for likelihood ratio values under H0
+ p value: percentage of simulated likelihood ratio values that are larger than the observed likelihood ratio value

+ afex does this using pbkrtest (Halekoh \& Hoojsgaard, 2014)
	+ built on lme4
	+ implements pb for glmm as well (but focus on LMM)
+ manual implementation for glmmTMB


### Prediction & Simulation

+ predict() function would be nice -> could run into problems with the way the model matrices are generated. Adding the variables from the model matrices to the data might help

+ simulating from known parameter values -> basis for future simulation-based power analysis functionality


### Simulation Study

Goals: 
+ Parameter Recovery
+ Power + Error Rate for Two-Step Procedure &  hierarchical approach (backward selection + parsimonious strategy + one-way ANOVA evtl.?)
	+ LRTs + Wald tests


+ two-step approach: two different corrections for 0 and 1 cells (Hautus, 1995)
	+ log-linear rule: add 0.5 to each cell (row and column totals + 1)
	+ --> less biased estimates (that always underestimate population d')
	+ 1/(2N) rule: 0 is replaced by 1 / (2N), 1 by 1 - 1/(2N)
	+ in general, Hautus (1995) showed that the correction methods lead to biased estimates (especially for small trial numbers)
	+ -> nice argument for the paper
	+ hierarchical approach deals with this through shrinkage <3

##### Design

+ Number of trials: (= number of signal trials in one cell)
	+ Axt et al. (2020): N = 200, 16 trials; N = 160, 16 trials; N = 300, 
	+ Detecting Bias studies: N = 250, minimum of 24 trials (32 trials when excluding the balanced condition)
	+ 


##### Barr et al. (2013) Setup

+ compared various LMMs and traditional methods (min F and F_1 times F_2) regarding Type I error rate and power
	+ similar to what I want to do
+ Design: item manipulation (between vs. within) x number of items (12 vs. 24) x fixed effect (present vs. absent)
+ manipulated between vs. within item manipulations, different numbers of items, different random-effect sizes and covariances
+ a single two-level treatment factor
+ 24 subjects, 12 vs. 24 items
+ Sampling model: random intercepts + random slopes slopes for items and subjects, respectively, correlated
+ residuals drawn from a normal distribution
+ population parameters for each dataset were drawn from distributions, except for the fixed effect 
	+ for tests of power: set to .8, yielding maximal power of about .5
	+ for tests of type I error rate: set to 0
+ different test statistics (t vs. chi^2 vs. MCMC)
+ 100000 repetitions
+ Analyses / Models:
	+ maximal model
	+ random intercepts only
	+ random slopes only
	+ maximal without correlations
	+ 6 model selection strategies
		+ random intercepts were always kept in the model
		+ Forward vs. backward
		+ by-item slope added / removed first vs. by-participant slope added / removed first vs. best path strategy (slope that leads to the best model is added / removed first)
+ Non-convergence: 
	+ non-converging models were progressively simplified until convergence was reached
	+ slope associated with the smaller variance was dropped (justification in the online appendix)
	+ reduced models were included in the metrics for the more complex models (e.g., if a maximal model is reduced to a random-intercept only model, it still "counts" as a maximal model)
	+ when the random-intercepts models did not converge, they were discarded (0.002 % of cases)
+ They don't always mention how they dealt with the correlations? 
+ Performance metrics:
	+ Type I error rate, power and Power', power rate corrected for anticonservativity
	+ alpha = .05 (results similar for alpha = .01 and alpha = .1)
	+ analyzed performance metrics as a function of "critical variances" (i.e., variances that )

+ Additional "worst case" analysis for the maximal model strategy: set random slope variances to 0 and estimate maximal model


##### Matuschek et al. (2017)

+ Sampling model: random intercepts + random slopes for items and subjects, respectively, correlated (.6)
+ 20 items, 50 subjects, 2 conditions
+ intercept fixed at 2000, effect of condition fixed at 25 vs. 0 (presence vs. absence of an effect)
+ SD of random intercepts (items and subjects) fixed to 100
+ SD of random slopes varied between 0 and 120
+ residual SD ~ N(0, 300^2)
+ Models:
	+ 1. Maximal model
	+ 2. Maximal model without correlations
	+ 3. Random intercepts + random slope for subjects
	+ 4. Random intercepts + random slope for item
	+ 5. Random intercepts
+ 2 Simulation scenarios: 
	+ 1. effect-related (co-)variances are set to 0 -> maximal cost of the maximal model specification
	+ 2. maximal model matches the data-generating process but variance components are below the threshold of detecability (i.e., small variances, small sample size)
+ Non-convergence: if any model did not converge, the the sample was redrawn
+ Model selection: 
	+ LRT criterion with alpha = .2
	+ backwards selection (without a "best-path" strategy) (1 -> 2 -> 3 -> 4 -> 5)

##### Boehm et al. (2018)

+ Comparison of different analysis strategies for hierarchically structured (response time)data
	+ two-step approach (my terminonology)
	+ fitting hierarchical models and conducting t tests on the resulting parameters
	+ fully hierarchical approach
+ Simulation design: 
	+ number of simulated trials per participant (2 vs. 5 vs. 15 vs. 30)
	+ number of participants (2, 5, 15, 30)
	+ effect size (0, 0.1, 0.5, 1)
+ 64 combinations, 200 repetitions


##### Setup

+ intercept, one binary predictor for response bias and sensitivity, crossed random effects for items and subjects with random intercepts, slopes and correlations
+ presence vs. absence of an effect

+ Models:
	+ two-step analysis procedure
	+ maximal
	+ parsimonious


### Simulation Study 2 - Crossed Random Effects

+ don't manipulate effects on mu and lambda separately anymore 
+ 


### Simulation Study 3 - unequal cell sizes and missing data


## Examples for the Two-Step Approach
+ Batailler et al. (2021)
+ Axt et al. (2018)
+ Nahon et al. (2024)
+ Fort & Shulman (2024)



## Resources


+ lme4 estimation:

	+ http://lme4.r-forge.r-project.org/slides/2011-03-16-Amsterdam/

	+ https://www.alexejgossmann.com/Dissect_lmer_part1/

+ profile likelihood stuff

	+ https://stats.stackexchange.com/questions/27976/how-can-i-find-the-standard-deviation-of-the-sample-standard-deviation-from-a-no/27981#27981

	+ https://stats.stackexchange.com/questions/28671/what-is-the-exact-definition-of-profile-likelihood



