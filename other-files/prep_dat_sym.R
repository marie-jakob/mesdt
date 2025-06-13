debi3 <- readRDS("other-files/dat_sym.rds")

names(debi3)

debi3 %>%
  dplyr::select(id, assessment, status, objective, committee, emp_gender, file_name,
         participant_gender, age) -> debi3
debi3$status <- ifelse(debi3$status == "signal", 1, -1)
debi3$assessment <- factor(debi3$assessment)
contrasts(debi3$assessment) <- contr.sum(2)
contrasts(debi3$committee) <- contr.sum(2)
contrasts(debi3$emp_gender) <- contr.sum(2)
debi3$participant_gender <- factor(ifelse(debi3$participant_gender == "w", "f", "m"))
contrasts(debi3$participant_gender) <- contr.sum(2)
debi3$id <- factor(debi3$id)
debi3$file_name <- factor(debi3$file_name)

use_data(debi3, overwrite = TRUE)

