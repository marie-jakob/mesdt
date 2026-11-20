# mesdt (development version)

# mesdt 0.2.0

## Bug fixes
* Numeric trial types coded 0/1 are now converted correctly
* Type II tests: correct df and p-values for interactions involving factors with more than two levels

## New features
* `dv_signal_level` and `trial_type_signal_level` to name the signal level of factors explicitly.
* `fit_mesdt()` reports which factor level is treated as signal trials / signal responses,
  and `summary()` shows the coding.
* `fit_mesdt()` stops with an informative error if the model is rank-deficient, 
  naming the redundant terms.
* The warning about the random-effects structure is now only given if the structure
  cannot be fitted as specified (previously a message was shown whenever uncorrelated
  terms were present).
