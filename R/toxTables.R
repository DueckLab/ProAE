#' Create toxicity tables for individual and composite PRO-CTCAE survey items.
#'
#'
#' 	  Data format should be in 'long' format, where each PRO-CTCAE item is a
#' 	  variable/column.
#'
#' @param dsn A data.frame object with PRO-CTCAE data.
#' @param id_var A character string. Name of ID variable differentiating each
#'   PRO-CTCAE survey/participant entered as a quoted string.
#' @param cycle_var A character string. Name of variable differentiating one
#'   longitudinal/repeated. PRO-CTCAE survey from another, within an individual
#'   ID.
#' @param baseline_val A number indicating the expected baseline cycle/time
#'   point.
#' @param arm_var A character string. Name of arm variable differentiating
#'   treatment groups. Must be character or factor class. Overall frequencies
#'   will be reported if no arm/grouping variable is provided. Defaults to
#'   \code{NA}.
#' @param type A character string. Type of summary measure to be be used.
#'   Options include: \code{"max_post_bl"} = Use subjects' maximum score
#'   post-baseline visit. \code{"bl_adjusted"} = Use subjects' baseline adjusted
#'   score over the study period. \code{"max"} = Use subjects' maximum score
#'   over the study period. The baseline adjusted score is derived by the
#'   following: If the maximum score post-baseline is more severe than the
#'   baseline score, then the use maximum score post-baseline is used as the
#'   adjusted score. Otherwise, if the maximum score post-baseline is the same
#'   or less serve than the baseline score, then zero (0) is used as the
#'   adjusted score. Defaults to \code{"bl_adjusted"}.
#' @param test A character string. Specify the statistical test to apply where
#'   comparing rates among arms. Options include: \code{"c"} = chi square,
#'   \code{"f"} = fisher's exact. Defaults to \code{"c"}.
#' @param cycle_limit A number. Limit the data to be analyzed up to and
#'   including a given cycle number or time point. Defaults to \code{NA}.
#' @param riskdiff Logical. Calculates risk differences between two arms. Valid
#'   if there are only two arms in the data.frame specified. This option will
#'   countermand options called with the \code{test} parameter. Defaults to
#'   \code{FALSE}.
#' @param risk_ci A character string. Specify the confidence interval type
#'   to be constructed for risk differences. Options include: \code{"wald"},
#'   and \code{"agresti-caffo"}. Defaults to \code{"wald"}.
#' @param risk_ci_alpha A number between 0 and 1. Specify the alpha level of
#'   the risk difference confidence intervals. Defaults to \code{0.05}.
#' @return A list object with data.frame elements for individual items and
#'   composite grades.
#' @importFrom magrittr %>%
#' @examples
#' toxTables(dsn=ProAE::tox_acute, id_var="id", cycle_var="Cycle", baseline_val=1)
#' @export

toxTables = function(dsn,
                    id_var,
                    cycle_var,
                    baseline_val,
                    type="bl_adjusted",
                    test="c",
                    riskdiff=FALSE,
                    risk_ci = "wald",
                    risk_ci_alpha = 0.05,
                    arm_var=NA,
                    cycle_limit=NA){

  # ----------------------------------------------------------------
  # -- Checks 1/2
  # ----------------------------------------------------------------

  ## -- Assign binding for data.frame variables used within dplyr functions as global variables

  base_val = NULL
  max_val = NULL
  bl_adjusted = NULL
  max_post_bl = NULL
  max = NULL

  ## -- Required parameters

  if(exists("dsn")){
    if(!is.data.frame(dsn)){
      stop("param dsn must be provided as a data.frame object")
    }
  } else {stop("param dsn not provided")}

  if(exists("id_var")){
    if(!is.character(id_var)){
      stop("param id_var must be provided as a character string")
    } else if (!(id_var %in% colnames(dsn))){
      stop(paste0("param id_var (", id_var, ") not found as a variable in dsn (", deparse(substitute(dsn)), ")"))
    }
  } else {stop("param id_var not provided")}

  if(exists("cycle_var")){
    if(!is.character(cycle_var)){
      stop("param cycle_var must be provided as a character string")
    } else if (!(cycle_var %in% colnames(dsn))){
      stop(paste0("param cycle_var (", cycle_var, ") not found as a variable in dsn (", deparse(substitute(dsn)), ")"))
    }
  } else {stop("param cycle_var not provided")}

  if(exists("baseline_val")){
    if(!(is.numeric(baseline_val) | is.integer(baseline_val) | length(baseline_val)==1)){
      stop("param baseline_val must be provided as a single number, of class numeric or integer")
    }
  } else {stop("param baseline_val not provided")}

  if(!(risk_ci %in% c("wald", "agresti-caffo"))){
    stop("param risk_ci must be one of the fallowing; 'wald' 'agresti-caffo'")
  }

  if(!(is.numeric(risk_ci_alpha))){
    stop("param risk_ci_alpha must be a numeric value between 0 and 1")
  } else if(!(0 < risk_ci_alpha & risk_ci_alpha < 1)){
    stop("param risk_ci_alpha must be a numeric value between 0 and 1")
  }

  ## -- Check for any duplicate individuals within cycles

  if(any(duplicated(dsn[,c(id_var, cycle_var)]))){
    stop(paste0("Duplicate observations were found within id_var and cycle_var combinations (", id_var, " and ", cycle_var, ")"))
  }

  ## -------------------------------------------------------------
  ## --- Reference data sets
  ## -------------------------------------------------------------

  # -- Get existing PROCTCAE variables in dsn

  # - Individual items
  dsn_items0 = toupper(names(dsn)[toupper(names(dsn)) %in% proctcae_vars$name])
  dsn_items = dsn_items0[! dsn_items0 %in% as.character(proctcae_vars$name[proctcae_vars$fmt %in% c("yn_2_fmt", "yn_3_fmt", "yn_4_fmt")])]

  # - Composites
  proctcae_vars_comp0 = proctcae_vars[,-1] %>% dplyr::mutate_if(is.factor, as.character)
  proctcae_vars_comp0 = proctcae_vars_comp0[!proctcae_vars_comp0$name %in% as.character(proctcae_vars_comp0$name[proctcae_vars_comp0$fmt %in% c("yn_2_fmt",
                                                                                                                                                "yn_3_fmt",
                                                                                                                                                "yn_4_fmt")]),]
  proctcae_vars_comp = c()
  proctcae_vars_comp$name = paste0(substr(proctcae_vars_comp0$name, 1, nchar(proctcae_vars_comp0$name)-5), "_COMP")
  proctcae_vars_comp$short_label = sub(proctcae_vars_comp0$short_label, pattern = " [[:alpha:]]*$", replacement = "")
  proctcae_vars_comp = unique(data.frame(proctcae_vars_comp))
  dsn_comps = toupper(names(dsn)[toupper(names(dsn)) %in% proctcae_vars_comp$name])

  # ----------------------------------------------------------------
  # -- Checks 2/2
  # ----------------------------------------------------------------

  ## -- Confirm there are available PRO-CTCAE variables within dsn with expected naming convention
  if(identical(dsn_items, character(0)) & identical(dsn_comps, character(0))){
    stop(paste0("No PRO-CTCAE variables found within dsn (", deparse(substitute(dsn)), ") meeting the expected naming convention"))
  }
  # ----------------------------------------------------------------
  # --- Limit cycles
  # ----------------------------------------------------------------

  if(!is.na(cycle_limit) & is.numeric(cycle_limit)){
    dsn = dsn[dsn[,cycle_var] <= cycle_limit,]
  }

  # ----------------------------------------------------------------
  # --- Individual items
  # ----------------------------------------------------------------
  tab_combined0 = data.frame()

  for(item in dsn_items){

    # --- Construct summary scores
    # -- Maximum post baseline
    dsn_max0 = dsn[dsn[,cycle_var] != baseline_val & !is.na(dsn[,cycle_var]), c(id_var, item)] %>%
      dplyr::group_by(get(id_var)) %>%
      dplyr::mutate(max_post_bl = ifelse(all(is.na(get(item))), NA, max(get(item), na.rm=TRUE)))
    dsn_max = unique(dsn_max0[,c(id_var, "max_post_bl")])

    # -- Baseline adjusted
    # - Subjects with an observed AE grade at baseline
    base_ids = dsn[dsn[,cycle_var] == baseline_val & !is.na(dsn[,cycle_var]), id_var]

    dsn_adj0 = dsn[dsn[,id_var] %in% base_ids, c(id_var, cycle_var,item)] %>%
      dplyr::group_by(get(id_var)) %>%
      dplyr::mutate(base_val = get(item)[get(cycle_var) == baseline_val],
             max_val = ifelse(all(is.na(get(item))), NA, max(get(item), na.rm=TRUE)),
             bl_adjusted = ifelse(base_val >= max_val, 0, max_val))
    dsn_adj = unique(dsn_adj0[,c(id_var, "bl_adjusted", "max_val")])

    dsn1 = merge(x=merge(x=dsn, y=dsn_max, by=id_var, all.x=TRUE),
                 y=dsn_adj, by=id_var, all.x=TRUE)

    dsn1$bl_adjusted[is.na(dsn1$max_post_bl)] = NA
    dsn1 = as.data.frame(dsn1)

    # ----------------------------------------------------------------
    # --- Overall (no arm provided)
    # ----------------------------------------------------------------
    if(is.na(arm_var)){
      dsn2 = dsn1 %>%
        dplyr::select(id_var, bl_adjusted, max_post_bl, max_val) %>%
        dplyr::group_by(get(id_var)) %>%
        dplyr::slice(1)

      dsn2 = as.data.frame(dsn2)

      # -- Choose baseline adjusted score or maximum score post-baseline
      if(tolower(type) == "bl_adjusted"){
        dsn2$present = as.numeric(dsn2$bl_adjusted > 0)
        dsn2$severe = as.numeric(dsn2$bl_adjusted >= 3)
      } else if(tolower(type) == "max_post_bl"){
        dsn2$present = as.numeric(dsn2$max_post_bl > 0)
        dsn2$severe = as.numeric(dsn2$max_post_bl >= 3)
      } else if(tolower(type) == "max"){
        dsn2$present = as.numeric(dsn2$max_val > 0)
        dsn2$severe = as.numeric(dsn2$max_val >= 3)
      }

      # ----------------------------------------------------------------
      # -- Statistics
      # ----------------------------------------------------------------

      dsn2$severe = factor(dsn2$severe, levels = c(0,1))
      dsn2$present = factor(dsn2$present, levels = c(0,1))

      # - Frequency tables
      tab_pres = table(dsn2$present)
      tab_sev = table(dsn2$severe)

      # -- Overall sample size;
      n0 = sum(tab_pres)
      n = data.frame(n0)
      colnames(n) = "overall"

      # -- Present
      col_pres = cbind.data.frame(n=n0, freq=tab_pres[2])
      col_pres$perc = 100*round(col_pres$freq / col_pres$n, 2)
      col_pres$text = paste0(col_pres$freq, " (",col_pres$perc,"%)")
      col_pres$text = gsub("NaN", "", col_pres$text, fixed=TRUE)

      col_pres_out = t(col_pres$text)
      colnames(col_pres_out) = "overall_pres"

      # -- Severe
      col_sev = cbind.data.frame(n=n0, freq=tab_sev[2])
      col_sev$perc = 100*round(col_sev$freq / col_sev$n, 2)
      col_sev$text = paste0(col_sev$freq, " (",col_sev$perc,"%)")
      col_sev$text = gsub("NaN", "", col_sev$text, fixed=TRUE)

      col_sev_out = t(col_sev$text)
      colnames(col_sev_out) = "overall_sev"

      row_out = data.frame(item_lab = as.character(proctcae_vars$short_label[proctcae_vars$name==item]),
                           n,
                           col_pres_out,
                           col_sev_out)

    } else {
      dsn2 = dsn1 %>%
        dplyr::select(id_var, arm_var, bl_adjusted, max_post_bl, max_val) %>%
        dplyr::group_by(get(id_var)) %>%
        dplyr::slice(1)

      dsn2 = as.data.frame(dsn2[,-5])

      # -- Choose baseline adjusted score or maximum score post-baseline

      if(tolower(type) == "bl_adjusted"){
        dsn2$present = as.numeric(dsn2$bl_adjusted > 0)
        dsn2$severe = as.numeric(dsn2$bl_adjusted >= 3)
      } else if(tolower(type) == "max_post_bl"){
        dsn2$present = as.numeric(dsn2$max_post_bl > 0)
        dsn2$severe = as.numeric(dsn2$max_post_bl >= 3)
      } else if(tolower(type) == "max"){
        dsn2$present = as.numeric(dsn2$max_post_bl > 0)
        dsn2$severe = as.numeric(dsn2$max_post_bl >= 3)
      }

      dsn2$severe = factor(dsn2$severe, levels = c(0,1))
      dsn2$present = factor(dsn2$present, levels = c(0,1))

      # ----------------------------------------------------------------
      # -- Statistics
      # ----------------------------------------------------------------

      # - Frequency tables
      tab_pres = table(dsn2[,arm_var], dsn2$present)
      tab_sev = table(dsn2[,arm_var], dsn2$severe)

      # - Risk differences (if called and there are two arms)
      if(riskdiff==TRUE & nrow(table(dsn[,arm_var]))==2){
        # Consider letting user specify alpha level
        alpha = risk_ci_alpha

        p1_pres = tab_pres[1,2]/sum(tab_pres[1,])
        p2_pres = tab_pres[2,2]/sum(tab_pres[2,])
        rdiff_est_pres = p1_pres - p2_pres

        if(risk_ci == "wald"){
          conf_ints = DescTools::BinomDiffCI(x1=tab_pres[1,2],
                                             n1=sum(tab_pres[1,]),
                                             x2=tab_pres[2,2],
                                             n2=sum(tab_pres[2,]),
                                             conf.level = 1 - alpha,
                                             sides = "two.sided",
                                             method = "wald")
          rdiff_lcl_pres = conf_ints[2]
          rdiff_ucl_pres = conf_ints[3]
        } else if(risk_ci == "agresti-caffo"){
          conf_ints = DescTools::BinomDiffCI(x1=tab_pres[1,2],
                                             n1=sum(tab_pres[1,]),
                                             x2=tab_pres[2,2],
                                             n2=sum(tab_pres[2,]),
                                             conf.level = 1 - alpha,
                                             sides = "two.sided",
                                             method = "ac")
          rdiff_lcl_pres = conf_ints[2]
          rdiff_ucl_pres = conf_ints[3]
        }

        rdiff_pres = paste0(100*round(rdiff_est_pres,2),"% (",100*round(rdiff_lcl_pres,2),"%, ",100*round(rdiff_ucl_pres,2),"%)")

        p1_sev = tab_sev[1,2]/sum(tab_sev[1,])
        p2_sev = tab_sev[2,2]/sum(tab_sev[2,])
        rdiff_est_sev = p1_sev - p2_sev

        if(risk_ci == "wald"){
          conf_ints = DescTools::BinomDiffCI(x1=tab_sev[1,2],
                                             n1=sum(tab_sev[1,]),
                                             x2=tab_sev[2,2],
                                             n2=sum(tab_sev[2,]),
                                             conf.level = 1 - alpha,
                                             sides = "two.sided",
                                             method = "wald")
          rdiff_lcl_sev = conf_ints[2]
          rdiff_ucl_sev = conf_ints[3]
        } else if(risk_ci == "agresti-caffo"){
          conf_ints = DescTools::BinomDiffCI(x1=tab_sev[1,2],
                                             n1=sum(tab_sev[1,]),
                                             x2=tab_sev[2,2],
                                             n2=sum(tab_sev[2,]),
                                             conf.level = 1 - alpha,
                                             sides = "two.sided",
                                             method = "ac")
          rdiff_lcl_sev = conf_ints[2]
          rdiff_ucl_sev = conf_ints[3]
        }

        rdiff_sev = paste0(100*round(rdiff_est_sev,2),"% (",100*round(rdiff_lcl_sev,2),"%, ",100*round(rdiff_ucl_sev,2),"%)")

      }else if(riskdiff==TRUE & nrow(table(dsn[,arm_var]))>2){
        test = "C"
        riskdiff = FALSE
      }

      # - Testing
      if(riskdiff==FALSE & tolower(test)=="c"){
        pv_pres = round(stats::chisq.test(tab_pres, correct = FALSE)$p.value, 4)
        pv_sev = round(stats::chisq.test(tab_sev, correct = FALSE)$p.value, 4)
      } else if(riskdiff==FALSE & tolower(test)=="f"){
        pv_pres = round(stats::fisher.test(tab_pres)$p.value, 4)
        pv_sev = round(stats::fisher.test(tab_sev)$p.value, 4)

      }

      # -- Group sample size;
      arm_n0 = rowSums(tab_pres)
      arm_n = t(arm_n0)
      colnames(arm_n) = paste0(gsub(" ", "", rownames(tab_pres), fixed = TRUE),"_n")

      # - If there is an arm with n=0 then remove the pvalues from the output
      for(j in colnames(arm_n)){
        if(arm_n[,j]==0){
          pv_pres=NA
          pv_sev=NA
        }
      }

      # -- Present
      col_pres = cbind.data.frame(n=arm_n0, freq=tab_pres[,2])
      col_pres$perc = 100*round(col_pres$freq / col_pres$n, 2)
      col_pres$text = paste0(col_pres$freq, " (",col_pres$perc,"%)")
      col_pres$text = gsub("NaN", "", col_pres$text, fixed=TRUE)

      col_pres_out = t(col_pres$text)
      colnames(col_pres_out) = paste0(gsub(" ", "", rownames(col_pres), fixed = TRUE),"_pres")

      if(all(col_pres$perc==0)){
        pv_pres=NA
      }

      # -- Severe
      col_sev = cbind.data.frame(n=arm_n0, freq=tab_sev[,2])
      col_sev$perc = 100*round(col_sev$freq / col_sev$n, 2)
      col_sev$text = paste0(col_sev$freq, " (",col_sev$perc,"%)")
      col_sev$text = gsub("NaN", "", col_sev$text, fixed=TRUE)

      col_sev_out = t(col_sev$text)
      colnames(col_sev_out) = paste0(gsub(" ", "", rownames(col_sev), fixed = TRUE),"_sev")

      if(all(col_sev$perc==0)){
        pv_sev=NA
      }

      # -- Combine for row output
      if(riskdiff==TRUE & nrow(table(dsn[,arm_var]))==2){
        row_out = data.frame(item_lab = as.character(proctcae_vars$short_label[proctcae_vars$name==item]),
                             arm_n,
                             col_pres_out,
                             rdiff_pres,
                             col_sev_out,
                             rdiff_sev)

      } else if(riskdiff==FALSE){
        row_out = data.frame(item_lab = as.character(proctcae_vars$short_label[proctcae_vars$name==item]),
                             arm_n,
                             col_pres_out,
                             pv_pres,
                             col_sev_out,
                             pv_sev)
      }
    }

    tab_combined0 = rbind(tab_combined0, row_out)

  }

  # -- sort respective tables by their PROCTCAE_##
  tab_combined = tab_combined0 %>% dplyr::mutate_if(is.factor, as.character)
  tab_combined = tab_combined[order(match(tab_combined$item_lab, proctcae_vars$short_label)),]
  rownames(tab_combined) = NULL

  # ----------------------------------------------------------------
  # --- Composite scores
  # ----------------------------------------------------------------

  tab_combined_comps0 = data.frame()

  for(item in dsn_comps){
    # ----------------------------------------------------------------
    # -- Construct summary scores
    # ----------------------------------------------------------------

    # --- Construct summary scores
    # -- Maximum post baseline
    dsn_max0 = dsn[dsn[,cycle_var] != baseline_val & !is.na(dsn[,cycle_var]), c(id_var, item)] %>%
      dplyr::group_by(get(id_var)) %>%
      dplyr::mutate(max_post_bl = ifelse(all(is.na(get(item))), NA, max(get(item), na.rm=TRUE)))
    dsn_max = unique(dsn_max0[,c(id_var, "max_post_bl")])

    # -- Baseline adjusted
    # - Subjects with an observed AE grade at baseline
    base_ids = dsn[dsn[,cycle_var] == baseline_val & !is.na(dsn[,cycle_var]), id_var]

    dsn_adj0 = dsn[dsn[,id_var] %in% base_ids, c(id_var, cycle_var,item)] %>%
      dplyr::group_by(get(id_var)) %>%
      dplyr::mutate(base_val = get(item)[get(cycle_var) == baseline_val],
             max_val = ifelse(all(is.na(get(item))), NA, max(get(item), na.rm=TRUE)),
             bl_adjusted = ifelse(base_val >= max_val, 0, max_val))
    dsn_adj = unique(dsn_adj0[,c(id_var, "bl_adjusted", "max_val")])

    dsn1 = merge(x=merge(x=dsn, y=dsn_max, by=id_var, all.x=TRUE),
                 y=dsn_adj, by=id_var, all.x=TRUE)

    dsn1$bl_adjusted[is.na(dsn1$max_post_bl)] = NA
    dsn1 = as.data.frame(dsn1)

    if(is.na(arm_var)){
      # ----------------------------------------------------------------
      # --- Overall (no arm provided)
      # ----------------------------------------------------------------

      dsn2 = dsn1 %>%
        dplyr::select(id_var, bl_adjusted, max_post_bl, max_val) %>%
        dplyr::group_by(get(id_var)) %>%
        dplyr::slice(1)

      dsn2 = as.data.frame(dsn2[,-4])

      # -- Choose baseline adjusted score or maximum score post-baseline
      if(tolower(type) == "bl_adjusted"){
        dsn2$present = as.numeric(dsn2$bl_adjusted > 0)
        dsn2$severe = as.numeric(dsn2$bl_adjusted >= 3)
      } else if(tolower(type) == "max_post_bl"){
        dsn2$present = as.numeric(dsn2$max_post_bl > 0)
        dsn2$severe = as.numeric(dsn2$max_post_bl >= 3)
      } else if(tolower(type) == "max"){
        dsn2$present = as.numeric(dsn2$max_post_bl > 0)
        dsn2$severe = as.numeric(dsn2$max_post_bl >= 3)
      }

      # ----------------------------------------------------------------
      # -- Statistics
      # ----------------------------------------------------------------

      dsn2$severe = factor(dsn2$severe, levels = c(0,1))
      dsn2$present = factor(dsn2$present, levels = c(0,1))

      # - Frequency tables
      tab_pres = table(dsn2$present)
      tab_sev = table(dsn2$severe)

      # -- Overall sample size;
      n0 = sum(tab_pres)
      n = data.frame(n0)
      colnames(n) = "overall"

      # -- Present
      col_pres = cbind.data.frame(n=n0, freq=tab_pres[2])
      col_pres$perc = 100*round(col_pres$freq / col_pres$n, 2)
      col_pres$text = paste0(col_pres$freq, " (",col_pres$perc,"%)")
      col_pres$text = gsub("NaN", "", col_pres$text, fixed=TRUE)

      col_pres_out = t(col_pres$text)
      colnames(col_pres_out) = "overall_pres"

      # -- Severe
      col_sev = cbind.data.frame(n=n0, freq=tab_sev[2])
      col_sev$perc = 100*round(col_sev$freq / col_sev$n, 2)
      col_sev$text = paste0(col_sev$freq, " (",col_sev$perc,"%)")
      col_sev$text = gsub("NaN", "", col_sev$text, fixed=TRUE)

      col_sev_out = t(col_sev$text)
      colnames(col_sev_out) = "overall_sev"

      row_out = data.frame(item_lab = as.character(proctcae_vars_comp$short_label[proctcae_vars_comp$name==item]),
                           n,
                           col_pres_out,
                           col_sev_out)

    } else {

      dsn2 = dsn1 %>%
        dplyr::select(id_var, arm_var, bl_adjusted, max_post_bl, max_val) %>%
        dplyr::group_by(get(id_var)) %>%
        dplyr::slice(1)

      dsn2 = as.data.frame(dsn2[,-5])

      # -- Choose baseline adjusted score or maximum score post-baseline

      if(tolower(type) == "bl_adjusted"){
        dsn2$present = as.numeric(dsn2$bl_adjusted > 0)
        dsn2$severe = as.numeric(dsn2$bl_adjusted >= 3)
      } else if(tolower(type) == "max_post_bl"){
        dsn2$present = as.numeric(dsn2$max_post_bl > 0)
        dsn2$severe = as.numeric(dsn2$max_post_bl >= 3)
      } else if(tolower(type) == "max"){
        dsn2$present = as.numeric(dsn2$max_post_bl > 0)
        dsn2$severe = as.numeric(dsn2$max_post_bl >= 3)
      }

      # ----------------------------------------------------------------
      # -- Statistics
      # ----------------------------------------------------------------

      dsn2$severe = factor(dsn2$severe, levels = c(0,1))
      dsn2$present = factor(dsn2$present, levels = c(0,1))

      # - Frequency tables
      tab_pres = table(dsn2[,arm_var], dsn2$present)
      tab_sev = table(dsn2[,arm_var], dsn2$severe)

      # - Risk differences (if called and there are two arms)
      if(riskdiff==TRUE & nrow(table(dsn[,arm_var]))==2){
        # Consider letting user specify alpha level
        alpha = risk_ci_alpha

        p1_pres = tab_pres[1,2]/sum(tab_pres[1,])
        p2_pres = tab_pres[2,2]/sum(tab_pres[2,])
        rdiff_est_pres = p1_pres - p2_pres

        if(risk_ci == "wald"){
          conf_ints = DescTools::BinomDiffCI(x1=tab_pres[1,2],
                                             n1=sum(tab_pres[1,]),
                                             x2=tab_pres[2,2],
                                             n2=sum(tab_pres[2,]),
                                             conf.level = 1 - alpha,
                                             sides = "two.sided",
                                             method = "wald")
          rdiff_lcl_pres = conf_ints[2]
          rdiff_ucl_pres = conf_ints[3]
        } else if(risk_ci == "agresti-caffo"){
          conf_ints = DescTools::BinomDiffCI(x1=tab_pres[1,2],
                                             n1=sum(tab_pres[1,]),
                                             x2=tab_pres[2,2],
                                             n2=sum(tab_pres[2,]),
                                             conf.level = 1 - alpha,
                                             sides = "two.sided",
                                             method = "ac")
          rdiff_lcl_pres = conf_ints[2]
          rdiff_ucl_pres = conf_ints[3]
        }

        rdiff_pres = paste0(100*round(rdiff_est_pres,2),"% (",100*round(rdiff_lcl_pres,2),"%, ",100*round(rdiff_ucl_pres,2),"%)")

        p1_sev = tab_sev[1,2]/sum(tab_sev[1,])
        p2_sev = tab_sev[2,2]/sum(tab_sev[2,])
        rdiff_est_sev = p1_sev - p2_sev

        if(risk_ci == "wald"){
          conf_ints = DescTools::BinomDiffCI(x1=tab_sev[1,2],
                                             n1=sum(tab_sev[1,]),
                                             x2=tab_sev[2,2],
                                             n2=sum(tab_sev[2,]),
                                             conf.level = 1 - alpha,
                                             sides = "two.sided",
                                             method = "wald")
          rdiff_lcl_sev = conf_ints[2]
          rdiff_ucl_sev = conf_ints[3]
        } else if(risk_ci == "agresti-caffo"){
          conf_ints = DescTools::BinomDiffCI(x1=tab_sev[1,2],
                                             n1=sum(tab_sev[1,]),
                                             x2=tab_sev[2,2],
                                             n2=sum(tab_sev[2,]),
                                             conf.level = 1 - alpha,
                                             sides = "two.sided",
                                             method = "ac")
          rdiff_lcl_sev = conf_ints[2]
          rdiff_ucl_sev = conf_ints[3]
        }

        rdiff_sev = paste0(100*round(rdiff_est_sev,2),"% (",100*round(rdiff_lcl_sev,2),"%, ",100*round(rdiff_ucl_sev,2),"%)")

      }else if(riskdiff==TRUE & nrow(table(dsn[,arm_var]))>2){
        test = "C"
        riskdiff = FALSE
      }

      # - Testing
      if(riskdiff==FALSE & tolower(test)=="c"){
        pv_pres = round(stats::chisq.test(tab_pres, correct = FALSE)$p.value, 4)
        pv_sev = round(stats::chisq.test(tab_sev, correct = FALSE)$p.value, 4)
      } else if(riskdiff==FALSE & tolower(test)=="f"){
        pv_pres = round(stats::fisher.test(tab_pres)$p.value, 4)
        pv_sev = round(stats::fisher.test(tab_sev)$p.value, 4)
      }

      # -- Group sample size;
      arm_n0 = rowSums(tab_pres)
      arm_n = t(arm_n0)
      colnames(arm_n) = paste0(gsub(" ", "", rownames(tab_pres), fixed = TRUE),"_n")

      # - If there is an arm with n=0 then remove the pvalues from the output
      for(j in colnames(arm_n)){
        if(arm_n[,j]==0){
          pv_pres=NA
          pv_sev=NA
        }
      }

      # -- Present
      col_pres = cbind.data.frame(n=arm_n0, freq=tab_pres[,2])
      col_pres$perc = 100*round(col_pres$freq / col_pres$n, 2)
      col_pres$text = paste0(col_pres$freq, " (",col_pres$perc,"%)")
      col_pres$text = gsub("NaN", "", col_pres$text, fixed=TRUE)

      col_pres_out = t(col_pres$text)
      colnames(col_pres_out) = paste0(gsub(" ", "", rownames(col_pres), fixed = TRUE),"_pres")

      if(all(col_pres$perc==0)){
        pv_pres=NA
      }

      # -- Severe
      col_sev = cbind.data.frame(n=arm_n0, freq=tab_sev[,2])
      col_sev$perc = 100*round(col_sev$freq / col_sev$n, 2)
      col_sev$text = paste0(col_sev$freq, " (",col_sev$perc,"%)")
      col_sev$text = gsub("NaN", "", col_sev$text, fixed=TRUE)

      col_sev_out = t(col_sev$text)
      colnames(col_sev_out) = paste0(gsub(" ", "", rownames(col_sev), fixed = TRUE),"_sev")

      if(all(col_sev$perc==0)){
        pv_sev=NA
      }

      # -- Combine for row output
      if(riskdiff==TRUE & nrow(table(dsn[,arm_var]))==2){
        row_out = data.frame(item_lab = as.character(proctcae_vars_comp$short_label[proctcae_vars_comp$name==item]),
                             arm_n,
                             col_pres_out,
                             rdiff_pres,
                             col_sev_out,
                             rdiff_sev)

      } else if(riskdiff==FALSE){
        row_out = data.frame(item_lab = as.character(proctcae_vars_comp$short_label[proctcae_vars_comp$name==item]),
                             arm_n,
                             col_pres_out,
                             pv_pres,
                             col_sev_out,
                             pv_sev)
      }
    }

    tab_combined_comps0 = rbind(tab_combined_comps0, row_out)

  }

  # -- sort respective tables by their PROCTCAE_##
  tab_combined_comps = tab_combined_comps0 %>% dplyr::mutate_if(is.factor, as.character)
  tab_combined_comps = tab_combined_comps[order(match(tab_combined_comps$item_lab, proctcae_vars_comp$short_label)),]
  rownames(tab_combined_comps) = NULL

  out = list(tab_combined, tab_combined_comps)
  names(out) = c("individual", "composite")

  # --------------------------------------------------------------------------
  # --- Exit
  # --------------------------------------------------------------------------

  return(out)
}
