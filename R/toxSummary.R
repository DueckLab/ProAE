#' Create patient-level and group-level summary statistics.
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
#' @param summary_measure A character string. Type of summary statistic to be
#'   used. Please consult current literature for appropriate interpretations of
#'   the summary measure selected and suitable analysis procedures for comparing
#'   groups. Options include: \code{"max"} = Use subjects'
#'   maximum score. \code{"max_post_bl"} = Use subjects' maximum score
#'   post-baseline visit. \code{"bl_adjusted"} = Use subjects' baseline adjusted
#'   score over the study period. The baseline adjusted score is derived by the
#'   following: If the maximum score post-baseline is more severe than the
#'   baseline score, then the use maximum score post-baseline is used as the
#'   adjusted score. Otherwise, if the maximum score post-baseline is the same
#'   or less serve than the baseline score, then zero (0) is used as the
#'   adjusted score. \code{"toxicity_idex"} = Construct patient-level toxicity
#'   index. \code{"AUC_worsening"} = Calculate group-level AUC describing
#' @param baseline_val A number indicating the expected baseline cycle/time
#'   point.
#' @param arm_var A character string. Name of arm variable differentiating
#'   treatment arms or other grouping factor. Required for group-level
#'   summary measures.
#' @return A data.frame with only the id and PRO-CTCAE variables being summarized.
#'   Each subject will now only have 1 observation (PRO-CTCAE variables are now the summary measure value).
#' @importFrom magrittr %>%
#' @examples
#' toxSummary(dsn=ProAE::tox_acute,
#' id_var="id",
#' cycle_var="Cycle",
#' baseline_val=1,
#' summary_measure = "max")
#' @export

toxSummary <- function(dsn,
                       id_var,
                       cycle_var,
                       summary_measure,
                       baseline_val=NA,
                       arm_var=NA){

  # ----------------------------------------------------------------
  # --- Checks 1/2
  # ----------------------------------------------------------------

  ## -- Assign binding for objects
  baseline = NULL

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

  if(!(summary_measure %in% c("max", "max_post_bl", "bl_adjusted", "toxicity_index", "AUC_worsening", "AUC_improvement"))){
    stop("param summary_measure must be one of the fallowing; 'max', 'max_post_bl', 'toxicity_idex', 'bl_adjusted', 'AUC_worsening', 'AUC_improvement'")
  }

  ## -- Checks for summary measures requiring baseline data information.
  if(summary_measure %in% c("max_post_bl", "bl_adjusted", "AUC_worsening", "AUC_improvement")){
    if(is.na(baseline_val)){
      stop("param baseline_val must be provided for this measure")
    } else{
      if(!(is.numeric(baseline_val) | is.integer(baseline_val) | length(baseline_val)==1)){
        stop("param baseline_val must be provided for this measure as a single number, of class numeric or integer")
      }
      if(min(dsn[, cycle_var]) != baseline_val){
        stop(paste0("The value of the param baseline_val (", baseline_val, ") is not the smallest ", cycle_var, "."))
      }
    }
  }

  if(summary_measure %in% c("AUC_worsening", "AUC_improvement")){
    if(is.na(arm_var) | !(arm_var %in% colnames(dsn))){
      stop(paste0("A valid arm or grouping variable within dsn must be provided for the summary measure selected"))
    }
  }

  # ----------------------------------------------------------------
  # --- Get existing PRO-CTCAE variables in dsn
  # ----------------------------------------------------------------

  ## -- Individual items
  dsn_items0 = toupper(names(dsn)[toupper(names(dsn)) %in% proctcae_vars$name])
  dsn_items = dsn_items0[! dsn_items0 %in% as.character(proctcae_vars$name[proctcae_vars$fmt %in% c("yn_2_fmt", "yn_3_fmt", "yn_4_fmt")])]

  ## -- Composites
  proctcae_vars_comp0 = proctcae_vars[,-1] %>% dplyr::mutate_if(is.factor, as.character)
  proctcae_vars_comp0 = proctcae_vars_comp0[!proctcae_vars_comp0$name %in% as.character(proctcae_vars_comp0$name[proctcae_vars_comp0$fmt %in% c("yn_2_fmt",
                                                                                                                                                "yn_3_fmt",
                                                                                                                                                "yn_4_fmt")]),]
  proctcae_vars_comp = c()
  proctcae_vars_comp$name = paste0(substr(proctcae_vars_comp0$name, 1, nchar(proctcae_vars_comp0$name)-5), "_COMP")
  proctcae_vars_comp$short_label = sub(proctcae_vars_comp0$short_label, pattern = " [[:alpha:]]*$", replacement = "")
  proctcae_vars_comp = unique(data.frame(proctcae_vars_comp, stringsAsFactors=FALSE))
  dsn_comps = toupper(names(dsn)[toupper(names(dsn)) %in% proctcae_vars_comp$name])

  dsn_items = c(dsn_items, dsn_comps)

  # ----------------------------------------------------------------
  # --- Checks 2/2
  # ----------------------------------------------------------------

  ## -- Confirm there are available PRO-CTCAE variables within dsn with expected naming convention
  if(identical(dsn_items, character(0))){
    stop(paste0("No PRO-CTCAE variables found within dsn (",
                deparse(substitute(dsn)),
                ") meeting the expected naming convention to summarize"))
  }

  ## -- Check for any duplicate individuals within cycles
  if(any(duplicated(dsn[,c(id_var, cycle_var)]))){
    stop(paste0("Duplicate observations were found within id_var and cycle_var combinations (", id_var, " and ", cycle_var, ")"))
  }

  ## -- Check to make sure scores are elements of {0, 1, 2, 3, 4}
  if(sum(!unlist(dsn[, dsn_items], use.names = FALSE)%in%c(0, 1, 2, 3, 4, NA))>0){
    stop("PRO-CTCAE item scores must be an integer between 0 and 4.")
  }

  # ----------------------------------------------------------------
  # --- Individual-level summary measures
  # ----------------------------------------------------------------

  if(summary_measure %in% c("max", "max_post_bl", "bl_adjusted", "toxicity_index")){
    ## -----------------------------
    ## -- Max
    ## -- This is the largest score across all timepoints for any observed data
    ## -----------------------------

    if(summary_measure == "max"){
      dsn_summary <- dsn %>%
        dplyr::group_by(get(id_var)) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(dsn_items), ~ifelse(all(is.na(.)), NA, max(., na.rm = TRUE))), .groups = "drop")
    }

    ## -----------------------------
    ## -- Max post-baseline
    ## -- This is the largest score across all timepoints after the baseline. Must have data after baseline to
    ## -- calculate this summary measure, otherwise return NA
    ## -----------------------------

    if(summary_measure == "max_post_bl") {
      dsn_summary <- dsn %>%
        dplyr::mutate(baseline = dplyr::if_else(get(cycle_var) == baseline_val, TRUE, FALSE)) %>%
        dplyr::group_by(get(id_var)) %>%
        dplyr::filter(!baseline) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(dsn_items), ~ifelse(all(is.na(.)), NA, max(., na.rm = TRUE)), .groups = "drop"))
    }

    ## -----------------------------
    ## -- Baseline adjusted max
    ## -- This is the largest score across all timepoints after the baseline, including only scores that were
    ## -- worse than the patient's baseline score. Must have data at baseline and after baseline to
    ## -- calculate this summary measure, otherwise return NA

    ## -- Function to calculate baseline adjustment summary statistic.
    ## -- scores := a numeric vector for each cycle for a given id
    ## -- baseline := a boolean vector identifying which score is the baseline score
    ## -----------------------------

    bl_adjusted <- function(scores, baseline) {
      if(length(scores) == 1 | sum(baseline) == 0){
        return(NA)
      } else {
        base_val = scores[which(baseline)]
        max_post_bl_val = ifelse(all(is.na(scores[which(!baseline)])), NA, max(scores[which(!baseline)], na.rm = TRUE))
        ifelse(base_val >= max_post_bl_val, 0, max_post_bl_val)
      }
    }

    if(summary_measure == "bl_adjusted"){
      dsn_summary = dsn %>%
        dplyr::mutate(baseline = (get(cycle_var) == baseline_val)) %>%
        dplyr::group_by(get(id_var)) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(dsn_items), ~bl_adjusted(., baseline)), .groups = "drop")
    }

    ## -----------------------------
    ## -- Toxicity index
    ## -----------------------------

    toxicity_index = function(x){
      # -- Remove NA obs and sort descending
      x_tmp = sort(x[!is.na(x)], decreasing = TRUE)
      if(length(x_tmp)==1){
        # only one grade
        ti = x_tmp[1]
      } else if(sum(x_tmp[-1])==0){
        # only one none-zero grade
        ti = x_tmp[1]
      } else {
        # compute
        ti = x_tmp[1]
        for(i in 1:(length(x_tmp)-1)){
          ti = ti + (x_tmp[i+1] / prod((1+x_tmp[1:i])))
          # -- prevent default rounding for large decimal portions
          if(ti-x_tmp[1] >= 0.9999){
            ti = x_tmp[1] + 0.9999
            # send warning to console / log that some estimates of the toxicity index were seen to
            # round out of a logical range, therefore subsequent round was prevented and returned
            # with estimates of [integer].9999
            break
          }
        }
      }
      return(ti)
    }

    if(summary_measure == "toxicity_index"){
      dsn_summary = dsn %>%
        dplyr::group_by(get(id_var)) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(dsn_items), ~toxicity_index(.)), .groups = "drop")
    }

    ## -----------------------------
    ## -----------------------------

    ## -- update id_var name
    colnames(dsn_summary)[1] = id_var

    ## -- add back in any previous id's that were filtered out
    dsn_summary_complete = dsn %>%
      dplyr::select(dplyr::all_of(id_var)) %>%
      dplyr::distinct() %>%
      dplyr::left_join(dsn_summary, by = dplyr::all_of(id_var))
  }


  ## ----------------------------------------------------------------
  ## --- Group-level summary measures
  ## ----------------------------------------------------------------

  if(summary_measure %in% c("AUC_worsening", "AUC_improvement")){

    ## -----------------------------
    ## -- AUC: AE burden worsening/improvement from baseline
    ## -----------------------------

    # -- Function for line-line intersection
    intersects_fun = function(blmean, x, y){
      intersects = c()
      for(i in 1:(length(y)-1)){
        y_ith_pair = y[i:(i+1)]
        x_ith_pair = x[i:(i+1)]
        # -- Get intersection if baseline value falls between the two adjacent points
        if(sum(blmean <= y_ith_pair)==1){
          # x-y swap here to accommodate vertical line intersection
          int_fun = stats::approxfun(y=x_ith_pair, x=y_ith_pair)
          int = int_fun(v = blmean)
          if(!is.na(int)){
            intersects = c(intersects, int)
          }
        }
      }
      x2_0 = c(x, intersects)
      y2 = c(y,rep(rep(blmean, length(intersects))))
      # Reorder by time with additional time at intersections
      x2 = x2_0[order(x2_0)]
      y2 = y2[order(x2_0)]
      return(list(x2,y2))
    }

    # -- Group-level AE means data
    # -- Only subjects with existing baseline values are considered for AUC calculation.
    for(ref_i in dsn_items){
      dsn_slice = dsn[,c(id_var, cycle_var, arm_var, ref_i)]
      dsn_slice_no_NA = dsn_slice[dsn_slice[,cycle_var]==baseline_val & !is.na(dsn_slice[,ref_i]),]
      dsn_slice_ids = unique(dsn_slice_no_NA[,id_var])

      dsn_slice = dsn_slice[dsn_slice[,id_var] %in% dsn_slice_ids,]

      if(which(dsn_items == ref_i) == 1){
        group_auc = stats::aggregate(data = dsn_slice[,c(arm_var, cycle_var, ref_i)],
                                     stats::formula(paste0(". ~", arm_var,"+", cycle_var)),
                                     mean)
      }else if(which(dsn_items == ref_i) > 1){
        group_auc = merge(group_auc,
                          stats::aggregate(data = dsn_slice[,c(arm_var, cycle_var, ref_i)],
                                           stats::formula(paste0(". ~", arm_var,"+", cycle_var)),
                                           mean),
                          by = c(arm_var, cycle_var),
                          all = TRUE)
      }
    }

    group_auc = group_auc[with(group_auc, order(group_auc[,cycle_var], group_auc[,arm_var])), ]

    # -- Create data frame for AUC to be attached
    group_out1 = data.frame(arm = unique(group_auc[,arm_var]))
    names(group_out1) = arm_var
    group_out2 = as.data.frame(matrix(as.numeric(),
                                      nrow = nrow(group_out1),
                                      ncol = length(dsn_items)))
    names(group_out2) = dsn_items
    auc_summary = cbind(group_out1, group_out2)

    for (item in dsn_items){
      for(i in unique(group_auc[,arm_var])){
        bl_val = group_auc[group_auc[,arm_var]==i & !is.na(group_auc[,item]) & group_auc[,cycle_var]==baseline_val,item]
        ints = intersects_fun(bl_val,
                              x = unique(stats::na.omit(group_auc[group_auc[,arm_var]==i & !is.na(group_auc[,item]),cycle_var])),
                              y = group_auc[group_auc[,arm_var]==i & !is.na(group_auc[,item]),item])
        t_i = ints[[1]]
        y_i = ints[[2]] - bl_val

        if(summary_measure == "AUC_worsening"){
          y_i[y_i < 0] = 0
        }

        if(summary_measure == "AUC_improvement"){
          y_i[y_i > 0] = 0
        }
        auc_i = sum( (t_i[2:length(t_i)] - t_i[2:length(t_i)-1]) * (y_i[2:length(t_i)] + y_i[2:length(t_i)-1])/2 )
        auc_summary[auc_summary[,arm_var]==i, item] = auc_i
      }
    }

    ## -----------------------------
    ## -----------------------------

    dsn_summary_complete = auc_summary

  }

  ## -----------------------------
  ## -----------------------------

  return(dsn_summary_complete)

}

