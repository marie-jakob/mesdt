#' Data from "On Bias in Detecting Bias: A Signal Detection Analysis of Attributions to Gender Discrimination (Experiment 3)
#'
#' Data reported in Jakob, Shechter, Calanchini, & Klauer (in prep.)
#' investigating attributions to gender discrimination with
#' a signal detection theory approach.
#' Participants had to judge 256 fictional pay raise decisions as biased
#' or unbiased. The cases involved male and female employees (`emp_gender`), who
#' were either granted or denied a pay raise. A quarter of all cases was biased,
#' that is, the actual committee decision (`committee`) did not match the
#' objectively fair decision (`objective`). Across all biased and unbiased cases,
#' male and female employees were denied and granted a pay raise in equal
#' proportions; thus, there was no overall gender bias in the committee
#' decisions. Still, we found that participants' response bias was aligned with
#' cultural stereotypes about gender bias in pay raise decisions, that is, they
#' were more liberal to judge a decision as biased if a female employee was
#' denied or a male employee was granted a pay raise.
#'
#' @format
#' A data frame with 124672 rows and 9 columns.
#' \describe{
#'   \item{id}{Factor indicating the participant ID}
#'   \item{assessment}{Factor indicating participants' judgments whether the
#'   presented pay raise decision was fair or likely biased}
#'   \item{status}{Factor indicating whether the trial was a signal (biased
#'   decision) or a noise (unbiased decision) trial}
#'   \item{objective}{Factor indicating the correct / unbiased decision}
#'   \item{committee}{Factor indicating whether the trial involved a decision
#'   where a pay raise was granted or denied}
#'   \item{emp_gender}{Factor indicating the perceived gender of the employee
#'   in a given trial}
#'   \item{file_name}{Factor indicating the picture shown in a given trial}
#'   \item{participant_gender}{Factor indicating the gender of the participant
#'   that gave the response in the trial}
#'   \item{age}{Numerical value indicating the age of the participant that gave
#'   the response in the trial}
#'   ...
#' }
#' @source TODO (once we have updated the preprint)
"debi3"
