# mesdt

Estimate and test mixed-effects signal detection theory models with maximum likelihood estimation. This package leverages the equivalence between a subclass of SDT and a subclass of generalized linear models first shown by DeCarlo (1998) to estimate these models using software for generalized linear mixed models, that is, the `lme4` or `glmmTMB` package. 
Mixed-effects SDT models can be specified on the level of SDT parameters discriminability and response bias using the typical R formula syntax:

` discriminability ~ X + (X | ID)`

` bias ~ X + (X | ID)` 

`mesdt` than translates the model formulas to a GLMM, using either `lme4` or `glmmTMB` to estimate the model (based on the user-specified backend), and transforms the parameters back to SDT space, such that all post-processing (e.g., estimating marginal means) can take place on the level of SDT parameters as well. Hypothesis tests related to fixed and random effects can be conducted with Wald tests, likelihood ratio tests (type II and type III), and tests based on parametric bootstrapping (type II and type III). 

## Installation

also install lme4 or glmmTMB

## Example

### Step 0: Prepare the Data

+ trial_type_var must be coded as either 1 and -1 or 



### Step 1: Specify and Fit a Mixed-Effects SDT Model

### Step 2: Compute Hypothesis Tests for Selected Parameters

### Step 3: Estimate marginal means

### Step 4: Report results with `papaja` methods

## References

glmmTMB package
lme4 package
DeCarlo 1998
