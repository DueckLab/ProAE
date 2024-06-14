# The `ProAE` Package <a href="https://duecklab.github.io/proae.html"><img src="man/images/logo.png" alt="ProAE logo" style="float:right;height:232.25px" align="right" height="232.25"></a>

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ProAE)](https://CRAN.R-project.org/package=ProAE)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/ProAE)](https://CRAN.R-project.org/package=ProAE)
[![Downloads](http://cranlogs.r-pkg.org/badges/ProAE)](https://CRAN.R-project.org/package=ProAE)
<!-- badges: end -->

## Overview

`library(proae)` is a collection of functions to facilitate standardized analysis and graphical procedures when using the National Cancer Institute’s Patient-Reported Outcomes version of the Common Terminology Criteria for Adverse Events (PRO-CTCAE).

To use the functions, data needs to be long (multiple observations per patient (e.g one line per cycle) and one column per PRO-CTCAE item). Additionally, PRO-CTCAE items need to follow specific naming structure. PRO-CTCAE variable names are made up of FOUR components: 1)’PROCTCAE’, 2) number [1,2,3, ..., i, ..., 80], 3) ’A’, ’B’, or ’C’ component of the i-th PRO-CTCAE field, 4) and ’SCL’ (if severity, interference, or frequency) or ’IND’ (if yes/no variable). Each component must be delimited with an underscore (_).

`ProAE::PROCTCAE_table` provides a crosswalk/look-up table of expected variable names for associated PRO-CTCAE symptom items. 

PRO-CTCAE Reference: https://healthcaredelivery.cancer.gov/pro-ctcae/instrument-pro.html

## The `toxScores()` Function

`toxScores()` is a function that easiliy reformats PRO-CTCAE items that were collected as text responses to numeric, applies zero-imputation procedures, and constructs PRO-CTCAE composite grades. This function accepts 1 or up to all 124 PRO-CTCAE items. This function returns a data frame with the respective numerical re-coding and computations. This function should be used prior to other functions in this package.

Composite grading reference: https://pubmed.ncbi.nlm.nih.gov/33258687/

## The `toxFigures` Function

`toxFigures()` is a function that creates publication quality longitudinal bar charts of symptomatic AE profiles. Frequency distributions of PRO-CTCAE scores can be shown over the course of the trial and can be stratified by treatment arm. An R list object is returned with the available PRO-CTCAE items indexed as elements of the list beside a reference table. See toxFigures vignette for more detail about how to customize these figures extensively using function parameters. 

## The `toxTables()` Function

`toxTables()` is a function that will produce PRO-CTCAE tables similar to those typically created for CTCAE data. These tables show arm-level sample size and frequency distributions for PRO-CTCAE items. P-values for statistical comparisons using chi-squared or Fisher's exact tests can be calculated as well as risk-difference between arms and alpha specified confidence intervals. See the toxTables() vignette for more details about how to use this function. 

## The `toxSummary()` Function

`toxSummary()` is a function that creates patient-level adjusted scores and group-level summary statistics. When requesting patient-level summary measures, the output is a data frame object with one observation per patient and include the `id_var` and PRO-CTCAE items are replaced with the adjusted score.

`summary_measure` options are:

*  `max` - patient-level - the maximum score over the treatment period for each item

*  `max_post_bl` - patient-level - the maximum score of all of the post-baseline timepoints for each item

*  `bl_adjusted` - patient-level - the maximum score post-baseline if the maximum score is more severe than baseline. If the baseline score is equivalent or more severe than all post-baseline scores then zero (0) is used as the adjusted score

* `toxicity_index` - patient-level - the toxicity index for each patient. Reference: https://pubmed.ncbi.nlm.nih.gov/34371138/

* `AUC_worsening' - group-level - For each level of the `arm_var` or the entire sample, it calculates the modified Area Under the Curve (AUC) where PRO-CTCAE scores are worsening from baseline

* `AUC_improvement' - group-level - For each level of the `arm_var` or the entire sample, it calculates the modified Area Under the Curve (AUC) where PRO-CTCAE scores are improvement from baseline

## The `toxAUC()` Function

`toxAUC()` is function that summarizes patient's symptomatic experience by calculating incremental AUC adjusted for baseline symptoms. Arm-level AUCs can be compared statistically using function parameters. See the toxAUC() vignette for more details about how to use this function.


## Other objects

`ProAE` also includes simulated data of common symptomatic AE profiles such as acute, cumulative, cyclical, and late toxicity at various effect sizes for power or exploratory analysis.

* `tox_acute` is a simulated data frame where the drug group experiences acute toxicity followed by symptom abatement over the course of treatment.

* `tox_chronic` is simulated data frame where the drug group experiences chronic toxicity over the course of treatment.

* `tox_cummulative` is a simulated data frame where drug toxicity is cumulative over the course of treatment.

* `tox_late` is a simulated data frame where the drug group experiences late incipient toxicity towards the end of the treatment period.

