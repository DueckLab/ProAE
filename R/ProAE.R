#' ProAE.
#'
#' This package is a collection of tools for analyzing and plotting PRO-CTCAE data
#'
#' @name ProAE
#' @docType package
NULL

#' PRO-CTCAE variable / label crosswalk
#'
#' A crosswalk / look-up table of expected variable names for associated PRO-CTCAE symptom items.
#'
#' \itemize{
#'   \item name. Expected variable name - item number/letter corresponds to the NCI-PRO-CTCAE (English) Item Library Version 1.0
#'   \item short_label. Item label including the item symptom group and frequency, severity, interference, or presence component
#' }
#'
#' @name PROCTCAE_table
#' @docType data
#' @author Blake Langlais
#' @references \url{https://healthcaredelivery.cancer.gov/pro-ctcae/instrument-pro.html}
#' @keywords data
#' @format A data frame with 124 rows and 2 variables
NULL

#' PRO-CTCAE data reflecting acute drug toxicity
#'
#' Simulated example data where the drug group experiences acute toxicity followed by symptom abatement over the course of treatment.
#'
#' \itemize{
#'   \item id. Subject identifier
#'   \item Cycle. Time variable denoting visits/cycles (1-10)
#'   \item arm. Treatment groups include drug and placebo
#'   \item PROCTCAE_78A_SC. PRO-CTCAE frequency item for nosebleeds
#'   \item PROCTCAE_78B_SC. PRO-CTCAE severity item for nosebleeds
#' }
#'
#' @name tox_acute
#' @docType data
#' @author Blake Langlais
#' @keywords data
#' @format A data frame with 1400 rows and 5 variables
NULL

#' PRO-CTCAE data reflecting chronic drug toxicity
#'
#' Simulated example data where the drug group experiences chronic toxicity over the course of treatment.
#'
#' \itemize{
#'   \item id. Subject identifier
#'   \item Cycle. Time variable denoting visits/cycles (1-10)
#'   \item arm. Treatment groups include drug and placebo
#'   \item PROCTCAE_78A_SC. PRO-CTCAE frequency item for nosebleeds
#'   \item PROCTCAE_78B_SC. PRO-CTCAE severity item for nosebleeds
#' }
#'
#' @name tox_chronic
#' @docType data
#' @author Blake Langlais
#' @keywords data
#' @format A data frame with 1400 rows and 5 variables
NULL

#' PRO-CTCAE data reflecting cumulative drug toxicity
#'
#' Simulated example data where drug toxicity is cumulative over the course of treatment.
#'
#' \itemize{
#'   \item id. Subject identifier
#'   \item Cycle. Time variable denoting visits/cycles (1-10)
#'   \item arm. Treatment groups include drug and placebo
#'   \item PROCTCAE_78A_SC. PRO-CTCAE frequency item for nosebleeds
#'   \item PROCTCAE_78B_SC. PRO-CTCAE severity item for nosebleeds
#' }
#'
#' @name tox_cumulative
#' @docType data
#' @author Blake Langlais
#' @keywords data
#' @format A data frame with 1400 rows and 5 variables
NULL

#' PRO-CTCAE data reflecting cyclical drug toxicity
#'
#' Simulated example data where the drug group experiences cyclicly toxicity post-treatment administration
#'
#' \itemize{
#'   \item id. Subject identifier
#'   \item Cycle. Time variable denoting visits/cycles (1-10)
#'   \item arm. Treatment groups include drug and placebo0
#'   \item PROCTCAE_78A_SC. PRO-CTCAE frequency item for nosebleeds
#'   \item PROCTCAE_78B_SC. PRO-CTCAE severity item for nosebleeds
#' }
#'
#' @name tox_cyclic
#' @docType data
#' @author Blake Langlais
#' @keywords data
#' @format A data frame with 1400 rows and 5 variables
NULL
#' PRO-CTCAE data reflecting late incipient drug toxicity
#'
#' Simulated example data where the drug group experiences late incipient toxicity towards the end of the treatment period.
#'
#' \itemize{
#'   \item id. Subject identifier
#'   \item Cycle. Time variable denoting visits/cycles (1-10)
#'   \item arm. Treatment groups include drug and placebo
#'   \item PROCTCAE_78A_SC. PRO-CTCAE frequency item for nosebleeds
#'   \item PROCTCAE_78B_SC. PRO-CTCAE severity item for nosebleeds
#' }
#'
#' @name tox_late
#' @docType data
#' @author Blake Langlais
#' @keywords data
#' @format A data frame with 1400 rows and 5 variables
NULL


