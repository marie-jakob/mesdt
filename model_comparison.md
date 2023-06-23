# Model Comparison

### Type II SS

+ Models to fit: 
  + For every order in the model:
    + Fit model with all effects with <= order (i.e., a model with all main effects, a model with all main effects + two-way interactions etc.)
    + Reference model for the respective order
    + For every effect of that order:
      + Fit a model with all effects with <= order _without_ the to-be-tested effect
      + Compare that model with the reference model
