library(mesdt)
# install.packages("tibble")
library(tibble)
set_backend("lme4")

debi3subset <- as_tibble(mesdt::debi3subset)
contrasts(debi3subset$committee) <- contr.sum(2)
contrasts(debi3subset$emp_gender) <- contr.sum(2)


debi3subset <- debi3subset[debi3subset$committee == "granted", ]

start <- Sys.time()
mod <- fit_mesdt(
  discriminability =~ committee * emp_gender * participant_gender + (committee * emp_gender || id),
  response_bias =~ committee * emp_gender * participant_gender + (committee * emp_gender || id),
  data = debi3subset,
  trial_type = "status",
  dv = "assessment",
  correlate_sdt_params = FALSE
)
end <- Sys.time()



# lme 4: prints
# fixed-effect model matrix is rank deficient so dropping 8 columns / coefficients
attr(mod$fit_obj@pp$X, "msgRankdrop")
# "fixed-effect model matrix is rank deficient so dropping 8 columns / coefficients"
# When nothing is dropped: NULL

# glmmTMB: 
attr(getME(mod$fit_obj, "X"), "col.dropped")
#mm[["lambda"]]committee1                     mm[["lambda"]]committee1:emp_gender1 
#2                                                        5 
#mm[["lambda"]]committee1:participant_gender1 mm[["lambda"]]committee1:emp_gender1:participant_gender1 
#6                                                        8 
#mm[["mu"]]committee1                         mm[["mu"]]committee1:emp_gender1 
#10                                                       13 
#mm[["mu"]]committee1:participant_gender1     mm[["mu"]]committee1:emp_gender1:participant_gender1 
#14                                                       16 
# glmmTMB: When nothing is dropped: NULL


# glm()

mod <- fit_mesdt(
  discriminability =~ committee * emp_gender * participant_gender,
  response_bias =~ committee * emp_gender * participant_gender,
  data = debi3subset,
  trial_type = "status",
  dv = "assessment",
  correlate_sdt_params = FALSE
)


