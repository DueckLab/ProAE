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
#'   adjusted score. \code{"toxicity_index"} = Construct patient-level toxicity
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

  proctcae_vars_comp0 = proctcae_vars[,-1] #gets rid of the fmt column
  fact_cols = sapply(proctcae_vars_comp0, is.factor) # Identify all factor columns
  proctcae_vars_comp0[fact_cols] = lapply(proctcae_vars_comp0[fact_cols], as.character)
  # if they are factor, change them to character


  proctcae_vars_comp0 = proctcae_vars_comp0[!proctcae_vars_comp0$name %in%
                                              as.character(proctcae_vars_comp0$name[proctcae_vars_comp0$fmt %in%
                                                                                      c("yn_2_fmt","yn_3_fmt","yn_4_fmt")]),]

  proctcae_vars_comp = c()
  proctcae_vars_comp$name = paste0(substr(proctcae_vars_comp0$name, 1, nchar(proctcae_vars_comp0$name)-5), "_COMP")
  proctcae_vars_comp$short_label = sub(proctcae_vars_comp0$short_label, pattern = " [[:alpha:]]*$", replacement = "")
  proctcae_vars_comp = unique(data.frame(proctcae_vars_comp, stringsAsFactors=FALSE))
  dsn_comps = toupper(names(dsn)[toupper(names(dsn)) %in% proctcae_vars_comp$name])

  dsn_items = c(dsn_items, dsn_comps)
  dsn_items =dsn_items[sort.list(dsn_items)]

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

  if(summary_measure %in% c("max", "max_post_bl", "bl_adjusted","toxicity_index")){


    # -------------APPLIES TO ALL SUMMARY MEASURES---------------------------------------------------
    refset = data.frame(name = dsn_items) #change to refset (now a df)
    out = unique(dsn[id_var]) #df with only id. will be merged onto for each dsn_item column


    for(i in 1:nrow(refset)){ #loop through each dsn_item

      item = as.character(refset$name[i]) #item is current dsn_item


      # -------------TOXICITY INDEX---------------------------------------------------------------

      toxicity_index = function(x){ #function
        # -- Remove NA obs and sort descending
        x_tmp = sort(x[!is.na(x)], decreasing = TRUE)
        if(length(x_tmp)==1){# length would be 1 if there was only 1 num in the vector
          # only one grade
          ti = x_tmp[1] #ti is just the only num in the vector then
        } else if(sum(x_tmp[-1])==0){ #-1 takes away the first number (the max).
          #this tests if there is more than 1 nonzero
          # only one none-zero grade
          ti = x_tmp[1] #therefore ti is just the max (if only nonzero)
        } else {
          # compute
          ti = x_tmp[1]
          for(i in 1:(length(x_tmp)-1)){
            ti = ti + (x_tmp[i+1] / prod((1+x_tmp[1:i]))) #is there a rounding issue?
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

        tox_vector = vector('list',nrow(unique(dsn[id_var]))) #vector to add tox index scores
        id_vector = vector('list',nrow(unique(dsn[id_var]))) #vector to keep track of which id on

        id_list = as.list(unique(dsn[,id_var]))


        for(j in 1:length(id_list)){ #loop through each unique id variable

          id = id_list[j] #id is the id actual value. J is th number to be subset
          x = dsn[dsn[,id_var] == id,c(item)] #numeric of all values in the "item" column for id = j
          t_index = toxicity_index(x) #using the previously defined function: outputs toxicity index
          tox_vector[[j]] = t_index #add tox score to score vector
          id_vector[[j]]  = id #add j to id vector

        }



        tox_df = cbind(as.data.frame(unlist(tox_vector)),as.data.frame(unlist(id_vector))) #bind the two vectors together
        colnames(tox_df) = c(item,id_var) #change col names

        out = merge(x = tox_df,y = out, by = id_var) #every loop through items, add one more column for each item

      }else{



        #----------------MAX VALUE UNAJUSTED-----------------------------------------------------------------------

        if(summary_measure == "max"){

          suppressWarnings({
          max_overall = stats::aggregate(dsn[item],
                                         by = dsn[id_var],
                                         FUN = max, na.rm = T)

          max_overall[,2] = as.integer(max_overall[,2])
          })

          #the final output!
          out = merge(x =max_overall,y = out, by = id_var)


        }else{#------------APPLIES TO MAX POST BASELINE OR BASELINE ADJUSTED---------------#\

          #create df with data only when cycle_var = baseline var
          #end result: one id column, and one item col that shows the baseline value per each id
          base_adj0 = dsn[dsn[,cycle_var] == baseline_val, c(id_var,item)]
          colnames(base_adj0)[2] = "base_val"

          #merge this with regular dsn
          #end result: a new df with all dsn cols + base_val col that has the same value for all
          #rows of the same id
          base_adj1 = merge(x = dsn, y = base_adj0, by = id_var, all.x = T)


          #this essentially creates max post baseline
          suppressWarnings({
          base_adj2 = stats::aggregate(base_adj1[base_adj1[,cycle_var]>baseline_val, item],
                                       by = list(base_adj1[base_adj1[,cycle_var]>baseline_val, id_var]),
                                       FUN = max,na.rm = T)

          base_adj2[,2] = as.integer(base_adj2[,2])
          })
          #base_adj2[,2] = as.integer(base_adj2[,2])
          #update names
          colnames(base_adj2) = c(id_var,"max_post_bl") #max post baseline

          #have to make sure the NA's are handled the right way
          base_adj2$max_pbl_updated =
            ifelse(is.na(base_adj2$max_post_bl)==T | is.infinite(base_adj2$max_post_bl)==T,NA,
                   base_adj2$max_post_bl) #as long as max post bl is not NA, then keep

          #drop max_post_bl variable
          base_adj2 = subset(base_adj2, select = -max_post_bl )


          #----------------Max Post-Baseline-------------------#

          if(summary_measure == "max_post_bl"){

            #create final df to use for max pbl
            max_pbl_final = base_adj2
            colnames(max_pbl_final) = c(id_var,item)

            #merge for each df. Out should be the dsn_summary
            out = merge(x = max_pbl_final,y = out, by = id_var)

            #----------------Adjusted Baseline-------------------#

          }else{
            #create a data frame that has baseline value and max_post_bl value
            base_adj3 = merge(x = base_adj1, y = unique(base_adj2), by = id_var, all.x = T)

            #NA check: if max post bl is na and bl is bigger than max post bl = 0
            base_adj3$adj_bl2=ifelse(!is.na(base_adj3$max_pbl_updated) & base_adj3$base_val>=base_adj3$max_pbl_updated,0,
                                     ifelse(is.na(base_adj3$base_val)==TRUE,NA,
                                            base_adj3$max_pbl_updated))

            #get it back to only one id and change the name so that it is the dsn_item name
            base_adj4 = base_adj3[base_adj3[,cycle_var] == baseline_val, c(id_var,"adj_bl2")]
            colnames(base_adj4)[colnames(base_adj4) == "adj_bl2"] = item

            #merge for each df. Out should be the dsn_summary
            out = merge(x = base_adj4,y = out, by = id_var)
          }

        }

      }
    }

    ## -- rearrange dsn_items columns in the in right order!!

    #new code
    order_new = c(id_var,dsn_items)
    out <- out[,order_new]
    dsn_summary = out

    ## -- add back in any previous id's that were filtered out

    dsn_summary_complete = merge(dsn_summary,unique(dsn[id_var]), by = id_var)
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
