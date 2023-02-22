#' Create longitudinal mean score line plots for PRO-CTCAE data with modified
#' Area Under the Curve (AUC) estimates showing descriptive symptomatic adverse
#' event burden worsening and improvement from baseline.
#'
#'
#' 	  Data format should be in 'long' format, where each PRO-CTCAE item is a
#' 	  variable/column. AUC calculations will only include subjects with non-
#' 	  missing baseline values (within each PRO-CTCAE item).
#'
#' @param dsn A data.frame object with PRO-CTCAE data
#' @param id_var A character string. Name of ID variable differentiating each
#'   unique patient.
#' @param cycle_var A character string. Name of variable differentiating one
#'   longitudinal/repeated. PRO-CTCAE survey from another, within an individual
#'   ID.
#' @param arm_var A character string. Name of arm variable differentiating
#'   treatment arms or other grouping factor. Up to 4 arms can be supported.
#' @param baseline_val A number indicating the expected baseline cycle/time
#'   point.
#' @param arm_var A character string. Name of arm variable differentiating
#'   treatment groups. Must be character or factor class. Overall frequencies
#'   will be reported if no arm/grouping variable is provided. Defaults to
#'   \code{NA}.
#' @param cycle_limit A number. Limit the number of cycles to be use to
#'   calculate the AUC metrics up to and including a given cycle number.
#'   All available cycle time points are used if no cycle number is provided.
#'   Defaults to \code{NA}.
#' @param y_limit A number. Y axis limit for plots. Defaults to \code{4}.
#' @param tab_ymin A number. Y axis coordinate for adjusting the vertical
#'   placement of the AUC table within the figure. Defaults to \code{NA}.
#' @param tab_ymax A number. Y axis coordinate for adjusting the vertical
#'   placement of the AUC table within the figure. Defaults to \code{NA}.
#' @param round_dec A number. Number of decimal places to be shown within
#'   the AUC table. Defaults to \code{2}.
#' @param overwrite_title A character string. Add main title to plots. Defaults
#'   to \code{NA}.
#' @param permute_tests Logical. Calls to calculate p values comparing the
#'   difference in AUC between two arms using a permutation test. Typical two-
#'   sided null hypothesis for a permutation test is applied. That is, assigning
#'   subjects to groups is interchangeable when calculating AUC. Computation
#'   time may be extensive depending on data size, number of PRO-CTCAE items,
#'   and number of permutations called. Consider staring out an open window or
#'   crafting a haiku during this time. Defaults to \code{FALSE}.
#' @param permute_n A number. The number of permutations to be used for
#'   permutation tests. Defaults to \code{2000}.
#' @param bootstrap_ci Logical. Calls to construct alpha-level confidence
#'   intervals for the difference in AUC between arms. Similar considerations
#'   for computation time as \code{permute_tests} are recommended here. Defaults
#'   to \code{FALSE}.
#' @param bootstrap_ci_alpha A number. Specifies the alpha level for bootstrap
#'   confidence intervals. Must be between 0 and 1. Defaults to \code{0.05}.
#' @param bootstrap_n A number. The number of bootstrap iterations to be used
#'   for bootstrap confidence intervals. Defaults to \code{2000}.
#'
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' AUC=toxAUC(dsn = ProAE::tox_acute,
#' id_var = "id",
#' cycle_var = "Cycle",
#' baseline_val = 1)
#' AUC[[1]]
#' }
#' @export


toxAUC = function(dsn,
                  id_var,
                  cycle_var,
                  baseline_val,
                  arm_var=NA,
                  cycle_limit=NA,
                  y_limit=4,
                  tab_ymin=NA,
                  tab_ymax=NA,
                  round_dec = 2,
                  overwrite_title = NA,
                  permute_tests = FALSE,
                  permute_n = 2000,
                  bootstrap_ci = FALSE,
                  bootstrap_ci_alpha = 0.05,
                  bootstrap_n = 2000){


    # ----------------------------------------------------------------
    # --- Checks 1
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

    ## -------------------------------------------------------------
    ## --- Reference data sets
    ## -------------------------------------------------------------

    # -- Get all scale items
    ref0 = proctcae_vars$name
    ref1 = ref0[ref0 %in% grep("_SCL", ref0, value=TRUE)]
    ref2 = rbind(matrix(ref1),
                 matrix(unique(paste0(substr(ref1,1,nchar(ref1)-5),"_COMP"))))

    # -- Get list of valid items within dsn
    ref = toupper(names(dsn)[toupper(names(dsn)) %in% ref2])

    # -- Restrict / limit analysis to a subset of cycles/timepoints
    if(!is.na(cycle_limit)){
      dsn = dsn[dsn[,cycle_var] <= cycle_limit,]
    } else{
      cycle_limit = max(dsn[,cycle_var])
    }

    # ----------------------------------------------------------------
    # -- Defaults
    # ----------------------------------------------------------------

    if(is.na(arm_var)){
      arm_var = "overall_arm___"
      dsn$overall_arm___ = "Overall"
    }

    # ------------------------------------------------------------------------------
    # -- AUC at the group level (area under the mean AE grade curve)
    # ------------------------------------------------------------------------------

    # -- Group-level AE means data
    # -- Only subjects with existing baseline values are considered for AUC calculation.
    for(ref_i in ref){
      dsn_slice = dsn[,c(id_var, cycle_var, arm_var, ref_i)]
      dsn_slice_no_NA = dsn_slice[dsn_slice[,cycle_var]==baseline_val & !is.na(dsn_slice[,ref_i]),]
      dsn_slice_ids = unique(dsn_slice_no_NA[,id_var])

      dsn_slice = dsn_slice[dsn_slice[,id_var] %in% dsn_slice_ids,]

      if(which(ref == ref_i) == 1){
        group_auc = stats::aggregate(data = dsn_slice[,c(arm_var, cycle_var, ref_i)],
                                     stats::formula(paste0(". ~", arm_var,"+", cycle_var)),
                                     mean)
      }else if(which(ref == ref_i) > 1){
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
                                      ncol = length(ref)))
    names(group_out2) = ref

    # --------
    # --- Reference lists for use in plot reference data sets
    # --------

    ref_labs0 = proctcae_vars[which(proctcae_vars$name %in% ref),-1]
    comp_tab0 = ref_labs0
    comp_tab0$name = gsub("([0-9]+).*$", "\\1", comp_tab0$name)
    comp_tab0$short_label = gsub("Frequency", "Composite", comp_tab0$short_label)
    comp_tab0$short_label = gsub("Severity", "Composite", comp_tab0$short_label)
    comp_tab0$short_label = gsub("Interference", "Composite", comp_tab0$short_label)
    comp_tab0$name = paste0(comp_tab0$name, "_COMP")
    comp_tab = unique(comp_tab0)

    ref_labs = rbind(ref_labs0, comp_tab)
    ref_labs = ref_labs[which(ref_labs$name %in% ref),]
    ref_labs$index = seq.int(nrow(ref_labs))


    # --------
    # --- Repeat the baseline adjusted with line intersection improvement
    # --------

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

    # --- WORSENING - Baseline adjusted (floor AUC is zero)
    bladj_group_out_floor_inter = cbind(group_out1, group_out2)
    ribbon_plot = vector(mode='list', length(ref))
    for (item in ref_labs$name){
      j = ref_labs[ref_labs$name ==item,"index"]

      ribbon_plot_arms = data.frame()
      for(i in unique(group_auc[,arm_var])){
        bl_val = group_auc[group_auc[,arm_var]==i & !is.na(group_auc[,item]) & group_auc[,cycle_var]==baseline_val,item]
        ints = intersects_fun(bl_val,
                              x = unique(stats::na.omit(group_auc[group_auc[,arm_var]==i & !is.na(group_auc[,item]),cycle_var])),
                              y = group_auc[group_auc[,arm_var]==i & !is.na(group_auc[,item]),item])
        t_i = ints[[1]]
        y_i = ints[[2]] - bl_val
        y_i[y_i < 0] = 0
        auc_i = sum( (t_i[2:length(t_i)] - t_i[2:length(t_i)-1]) * (y_i[2:length(t_i)] + y_i[2:length(t_i)-1])/2 )
        bladj_group_out_floor_inter[bladj_group_out_floor_inter[,arm_var]==i, item] = auc_i
        ribbon_plot_i = data.frame(arm = rep(i, length(t_i)),
                                   time = t_i,
                                   value = ints[[2]],
                                   bl_val = bl_val)
        names(ribbon_plot_i) = c(arm_var, cycle_var, item, "bl_val")
        ribbon_plot_arms = rbind.data.frame(ribbon_plot_arms, ribbon_plot_i)

      }
      ribbon_plot[[j]] = ribbon_plot_arms[order(ribbon_plot_arms[,cycle_var],ribbon_plot_arms[,arm_var]),]
    }

    # --- IMPROVEMENT - Baseline adjusted (ceiling AUC is zero)
    bladj_group_out_ceiling_inter = cbind(group_out1, group_out2)
    for (item in ref_labs$name){
      for(i in unique(group_auc[,arm_var])){
        bl_val = group_auc[group_auc[,arm_var]==i & !is.na(group_auc[,item]) & group_auc[,cycle_var]==baseline_val,item]
        ints = intersects_fun(bl_val,
                              x = unique(stats::na.omit(group_auc[group_auc[,arm_var]==i & !is.na(group_auc[,item]),cycle_var])),
                              y = group_auc[group_auc[,arm_var]==i & !is.na(group_auc[,item]),item])
        t_i = ints[[1]]
        y_i = ints[[2]] - bl_val
        y_i[y_i > 0] = 0
        auc_i = sum( (t_i[2:length(t_i)] - t_i[2:length(t_i)-1]) * (y_i[2:length(t_i)] + y_i[2:length(t_i)-1])/2 )
        bladj_group_out_ceiling_inter[bladj_group_out_ceiling_inter[,arm_var]==i, item] = auc_i
      }
    }

    # ------------------------------------------------------------------------------
    # --- Allow for permutation tests for the difference in AUC between 2 arms (iff 2 arms)
    # ------------------------------------------------------------------------------

    if(nrow(data.frame(table(dsn[,arm_var]))) == 2 & permute_tests == TRUE){

      perm_dsn = dsn[, c(id_var, arm_var, cycle_var, item)]
      ids = unique(perm_dsn[,c(id_var, arm_var)])

      for(l in 1:permute_n){

        # randomly assign patient ids to group a or group b
        ids$arm_rand = sample(ids[,arm_var], length(ids[,arm_var]), replace = FALSE)
        perm_dat = merge(perm_dsn[,-grep(arm_var, names(perm_dsn))], ids[,c(id_var, "arm_rand")], id_var)
        arm_var2 = "arm_rand"
        # # -- Get mean score at each cycle for each arm
        group_auc_perm = stats::aggregate(data = perm_dat[,c(arm_var2, cycle_var, item)],
                                   stats::formula(paste0(". ~", arm_var2,"+", cycle_var)),
                                   mean)

        # -- Initiates data.frames
        if(l==1){
          group_out1 = data.frame(arm = unique(group_auc_perm[,arm_var2]))
          names(group_out1) = arm_var2
          group_out2 = as.data.frame(matrix(as.numeric(),
                                            nrow = nrow(group_out1),
                                            ncol = 1))
          names(group_out2) = item
          bladj_group_out_floor_inter_auc = cbind(group_out1, group_out2) #### left off here
          bladj_group_out_ceiling_inter_auc = cbind(group_out1, group_out2)
          bladj_group_out_floor_inter_auc_diffs = data.frame()
          bladj_group_out_ceiling_inter_auc_diffs = data.frame()
        }

        # --- WORSENING - Baseline adjusted (floor AUC is zero)
        for(i in unique(group_auc_perm[,arm_var2])){
          bl_val = group_auc_perm[group_auc_perm[,arm_var2]==i & !is.na(group_auc_perm[,item]) & group_auc_perm[,cycle_var]==baseline_val,item]
          ints = intersects_fun(bl_val,
                                x = unique(stats::na.omit(group_auc_perm[group_auc_perm[,arm_var2]==i & !is.na(group_auc_perm[,item]),cycle_var])),
                                y = group_auc_perm[group_auc_perm[,arm_var2]==i & !is.na(group_auc_perm[,item]),item])
          t_i = ints[[1]]
          y_i = ints[[2]] - bl_val
          y_i[y_i < 0] = 0
          auc_i = sum( (t_i[2:length(t_i)] - t_i[2:length(t_i)-1]) * (y_i[2:length(t_i)] + y_i[2:length(t_i)-1])/2 )
          if(l==1){
            bladj_group_out_floor_inter_auc[bladj_group_out_floor_inter_auc[,arm_var2]==i, item] = auc_i
            bladj_group_out_floor_inter_auc$iteration = 1
          } else{
            run_l = data.frame(arm = i,
                               item = auc_i,
                               iteration = l)
            names(run_l)[1:2] = c(arm_var2, item)
            bladj_group_out_floor_inter_auc = rbind.data.frame(bladj_group_out_floor_inter_auc, run_l)
          }
        }
        arm1_auc = bladj_group_out_floor_inter_auc[bladj_group_out_floor_inter_auc$iteration==l,item][1]
        arm2_auc = bladj_group_out_floor_inter_auc[bladj_group_out_floor_inter_auc$iteration==l,item][2]
        bladj_group_out_floor_inter_auc_diffs = rbind.data.frame(bladj_group_out_floor_inter_auc_diffs,
                                                                 data.frame(iteration = l,
                                                                            diff_auc = arm1_auc - arm2_auc))

        # --- IMPROVEMENT - Baseline adjusted (ceiling AUC is zero)
        for(i in unique(group_auc_perm[,arm_var2])){
          bl_val = group_auc_perm[group_auc_perm[,arm_var2]==i & !is.na(group_auc_perm[,item]) & group_auc_perm[,cycle_var]==baseline_val,item]
          ints = intersects_fun(bl_val,
                                x = unique(stats::na.omit(group_auc_perm[group_auc_perm[,arm_var2]==i & !is.na(group_auc_perm[,item]),cycle_var])),
                                y = group_auc_perm[group_auc_perm[,arm_var2]==i & !is.na(group_auc_perm[,item]),item])
          t_i = ints[[1]]
          y_i = ints[[2]] - bl_val
          y_i[y_i > 0] = 0
          auc_i = sum( (t_i[2:length(t_i)] - t_i[2:length(t_i)-1]) * (y_i[2:length(t_i)] + y_i[2:length(t_i)-1])/2 )
          if(l==1){
            bladj_group_out_ceiling_inter_auc[bladj_group_out_ceiling_inter_auc[,arm_var2]==i, item] = auc_i
            bladj_group_out_ceiling_inter_auc$iteration = 1
          } else{
            run_l = data.frame(arm = i,
                               item = auc_i,
                               iteration = l)
            names(run_l)[1:2] = c(arm_var2, item)
            bladj_group_out_ceiling_inter_auc = rbind.data.frame(bladj_group_out_ceiling_inter_auc, run_l)
          }
        }
        arm1_auc = bladj_group_out_ceiling_inter_auc[bladj_group_out_ceiling_inter_auc$iteration==l,item][1]
        arm2_auc = bladj_group_out_ceiling_inter_auc[bladj_group_out_ceiling_inter_auc$iteration==l,item][2]
        bladj_group_out_ceiling_inter_auc_diffs = rbind.data.frame(bladj_group_out_ceiling_inter_auc_diffs,
                                                                   data.frame(iteration = l,
                                                                              diff_auc = arm1_auc - arm2_auc))
      }
    }

    # ------------------------------------------------------------------------------
    # --- Allow for bootstrap for alpha-level confidence intervals for difference in AUC between 2 arms (diff 2 arms)
    # ------------------------------------------------------------------------------

    if(nrow(data.frame(table(dsn[,arm_var]))) == 2 & bootstrap_ci == TRUE){

      unique_ids_arm1 = unique(dsn[dsn[,arm_var] == unique(dsn[,arm_var])[1],id_var])
      unique_ids_arm2 = unique(dsn[dsn[,arm_var] == unique(dsn[,arm_var])[2],id_var])

      for(l in 1:bootstrap_n){

        # -- Create bootstrap sample
        boot_ids_arm1 = sample(x = unique_ids_arm1, size = length(unique_ids_arm1), replace = TRUE)
        boot_ids_arm2 = sample(x = unique_ids_arm2, size = length(unique_ids_arm2), replace = TRUE)
        boot_ids = c(boot_ids_arm1, boot_ids_arm2)
        boot_dat = c()
        for(k in 1:length(boot_ids)){
          id = boot_ids[k]
          boot_dat = rbind.data.frame(boot_dat, dsn[dsn[,id_var] == id,])
        }

        # -- Get mean score at each cycle for each arm
        group_auc_boot = stats::aggregate(data = boot_dat[,c(arm_var, cycle_var, item)],
                                          stats::formula(paste0(". ~", arm_var,"+", cycle_var)),
                                   mean)

        if(l==1){
          group_out1 = data.frame(arm = unique(group_auc_boot[,arm_var]))
          names(group_out1) = arm_var
          group_out2 = as.data.frame(matrix(as.numeric(),
                                            nrow = nrow(group_out1),
                                            ncol = 1))
          names(group_out2) = item
          bladj_group_out_floor_inter_boot = cbind(group_out1, group_out2)
          bladj_group_out_ceiling_inter_boot = cbind(group_out1, group_out2)
        }


        # --- WORSENING - Baseline adjusted (floor AUC is zero)
        for(i in unique(group_auc_boot[,arm_var])){
          bl_val = group_auc_boot[group_auc_boot[,arm_var]==i & !is.na(group_auc_boot[,item]) & group_auc_boot[,cycle_var]==baseline_val,item]
          ints = intersects_fun(bl_val,
                                x = unique(stats::na.omit(group_auc_boot[group_auc_boot[,arm_var]==i & !is.na(group_auc_boot[,item]),cycle_var])),
                                y = group_auc_boot[group_auc_boot[,arm_var]==i & !is.na(group_auc_boot[,item]),item])
          t_i = ints[[1]]
          y_i = ints[[2]] - bl_val
          y_i[y_i < 0] = 0
          auc_i = sum( (t_i[2:length(t_i)] - t_i[2:length(t_i)-1]) * (y_i[2:length(t_i)] + y_i[2:length(t_i)-1])/2 )
          if(l==1){
            bladj_group_out_floor_inter_boot[bladj_group_out_floor_inter_boot[,arm_var]==i, item] = auc_i
            bladj_group_out_floor_inter_boot$iteration = 1
          } else{
            run_l = data.frame(arm = i,
                               item = auc_i,
                               iteration = l)
            names(run_l)[1:2] = c(arm_var, item)
            bladj_group_out_floor_inter_boot = rbind.data.frame(bladj_group_out_floor_inter_boot, run_l)
          }
        }

        # --- IMPROVEMENT - Baseline adjusted (ceiling AUC is zero)
        for(i in unique(group_auc_boot[,arm_var])){
          bl_val = group_auc_boot[group_auc_boot[,arm_var]==i & !is.na(group_auc_boot[,item]) & group_auc_boot[,cycle_var]==baseline_val,item]
          ints = intersects_fun(bl_val,
                                x = unique(stats::na.omit(group_auc_boot[group_auc_boot[,arm_var]==i & !is.na(group_auc_boot[,item]),cycle_var])),
                                y = group_auc_boot[group_auc_boot[,arm_var]==i & !is.na(group_auc_boot[,item]),item])
          t_i = ints[[1]]
          y_i = ints[[2]] - bl_val
          y_i[y_i > 0] = 0
          auc_i = sum( (t_i[2:length(t_i)] - t_i[2:length(t_i)-1]) * (y_i[2:length(t_i)] + y_i[2:length(t_i)-1])/2 )
          if(l==1){
            bladj_group_out_ceiling_inter_boot[bladj_group_out_ceiling_inter_boot[,arm_var]==i, item] = auc_i
            bladj_group_out_ceiling_inter_boot$iteration = 1
          } else{
            run_l = data.frame(arm = i,
                               item = auc_i,
                               iteration = l)
            names(run_l)[1:2] = c(arm_var, item)
            bladj_group_out_ceiling_inter_boot = rbind.data.frame(bladj_group_out_ceiling_inter_boot, run_l)
          }
        }
      }
    }


    # ------------------------------------------------------------------------------
    # --- Line plots with AUC
    # ------------------------------------------------------------------------------

    list_out = vector(mode='list', length(ref))

    for(item in ref_labs$name){
      i = ref_labs[ref_labs$name ==item,"index"]

      plot_auc = group_auc[,c(arm_var, cycle_var, item)]

      bl_dat = plot_auc[plot_auc[,cycle_var] == baseline_val,]

      ## Overall / 1-arm accommodation
      if(nrow(bl_dat)==1){
        bl_arm1 = bl_dat[1,item]
        name_arm1 = bl_dat[1,arm_var]

        plot_auc$bl_val = ifelse(plot_auc[,arm_var] == name_arm1, bl_arm1)

        # -- AUC table
        group_bl_floor_inter = bladj_group_out_floor_inter[,c(arm_var, item)]
        group_bl_ceiling_inter =  bladj_group_out_ceiling_inter[,c(arm_var, item)]

        group_bl_floor_inter = round(group_bl_floor_inter[,2], round_dec)
        group_bl_ceiling_inter = round(group_bl_ceiling_inter[,2], round_dec)

        anno_tab = cbind(data.frame(bl_dat[,1]),
                         group_bl_floor_inter,
                         group_bl_ceiling_inter)

        names(anno_tab) = c("",
                            "Worsening [solid area]",
                            "Improvement [striped area]")

        ribbon_dat = ribbon_plot[[i]]

        ribbon_pattern_dat0 = ribbon_dat[ribbon_dat[,item] <= ribbon_dat$bl_val, ]
        ribbon_pattern_dat0 = ribbon_pattern_dat0[order(ribbon_pattern_dat0[,arm_var], ribbon_pattern_dat0[,cycle_var]),]
        ribbon_pattern_dat1 = ribbon_pattern_dat0[ribbon_pattern_dat0[,item] <= ribbon_pattern_dat0$bl_val & ribbon_pattern_dat0[,arm_var] == name_arm1, ]

        ribbon_pattern_dat = rbind.data.frame(ribbon_pattern_dat1[-1,])

        item_title = ref_labs[ref_labs[,"name"]==item, "short_label"]

        df1 <- data.frame(x1 = baseline_val, x2 = cycle_limit, y1 = bl_arm1, y2 = bl_arm1, arm = name_arm1)

        df = rbind.data.frame(df1)
        names(df)[5] = arm_var

        # -- Adjustments in vertical placement of the AUC table within the figure
        if(is.na(tab_ymin) | is.na(tab_ymax)){
          tab_ymin = 2.75 - (4 - y_limit)
          tab_ymax = 3.95 - (4 - y_limit)
        }

        # -- Force the fill colors to be the same for all plots
        plot_auc[,arm_var] = as.factor(plot_auc[,arm_var])
        ribbon_dat[,arm_var] = as.factor(ribbon_dat[,arm_var])

        # -- Solid cols
        ribbon_dat$arm_col = ifelse(ribbon_dat[,arm_var] == name_arm1, "black", "#E69F00")
        arm_fill_cols = unique(ribbon_dat$arm_col)

        # -- Pattern cols
        ribbon_pattern_dat$arm_col = ifelse(ribbon_pattern_dat[,arm_var] == name_arm1, "black", "#E69F00")
        arm_fill_pattern_cols = unique(ribbon_pattern_dat$arm_col)

        figure_i = ggplot2::ggplot(plot_auc, ggplot2::aes_string(x=cycle_var, y=item, group=arm_var, fill=arm_var)) +
          ggplot2::geom_line(ggplot2::aes_string(color=arm_var, linetype=arm_var)) +
          ggplot2::geom_point(ggplot2::aes_string(color=arm_var)) +

          # -- Worsening
          ggplot2::geom_ribbon(data = ribbon_dat[ribbon_dat[,item] >= ribbon_dat$bl_val & ribbon_dat[,arm_var] == name_arm1, ],
                      ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, fill=arm_var), alpha = 0.4, inherit.aes = FALSE) +

           ggpattern::geom_ribbon_pattern(
            data = ribbon_pattern_dat,
            ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, pattern_fill = arm_var, pattern_color = arm_var),
            inherit.aes = FALSE,
            pattern_alpha = 0.4,
            fill            = NA,
            pattern         = 'stripe',
            pattern_spacing = 0.02,
            pattern_density = 0.5,
            show.legend = FALSE
          ) +

          # --- Add horizontal line indicating baseline value
          ggplot2::geom_segment(ggplot2::aes_string(x = "x1", y = "y1", xend = "x2", yend = "y2", colour = arm_var),
                                linetype="dashed",
                                data = df,
                                inherit.aes = FALSE) +

          ggplot2::scale_x_continuous(name=cycle_var, limits=c(baseline_val, cycle_limit), breaks = baseline_val:cycle_limit) +
          ggplot2::scale_y_continuous(name="Mean Score", limits=c(0, y_limit)) +
          ggplot2::scale_linetype_manual(values=c("solid")) +
          ggplot2::scale_color_manual(values=arm_fill_cols) +
          ggplot2::scale_fill_manual(values=arm_fill_cols) +
          ggpattern::scale_pattern_color_manual(values=arm_fill_pattern_cols) +
          ggpattern::scale_pattern_fill_manual(values=arm_fill_pattern_cols) +
          ggplot2::annotation_custom(gridExtra::tableGrob(t(anno_tab)), xmin=0.375*cycle_limit, xmax=0.8125*cycle_limit, ymin=tab_ymin, ymax=tab_ymin) +
          ggplot2::theme_bw() +
          ggplot2::ggtitle(ifelse(is.na(overwrite_title),
                         paste0(item_title, ": Baseline adjusted AUC"),
                         overwrite_title)) +
          ggplot2::labs(caption = paste0("Dashed horizontal line indicates arm-level baseline symptomatic adverse event level\n",
                                "Positive AUC represents increased symptomatic adverse events from baseline (Worsening)\n",
                                "Negative AUC represents decreased symptomatic adverse events from baseline (Improvement)\n")) +

          ggplot2::theme(legend.position = "none",
                legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                legend.title = ggplot2::element_blank(),
                plot.caption = ggplot2::element_text(hjust = 0, size = 10))

        auc_tab_out = as.data.frame(t(anno_tab)[2:3,], stringsAsFactors = FALSE)
        names(auc_tab_out) = t(anno_tab)[1,]
        rownames(auc_tab_out) = rownames(t(anno_tab)[2:3,])

      } else if(nrow(bl_dat)==2){
        bl_arm1 = bl_dat[1,item]
        name_arm1 = bl_dat[1,arm_var]

        bl_arm2 = bl_dat[2,item]
        name_arm2 = bl_dat[2,arm_var]

        plot_auc$bl_val = ifelse(plot_auc[,arm_var] == name_arm1, bl_arm1, bl_arm2)

        # -- AUC table
        group_bl_floor_inter = bladj_group_out_floor_inter[,c(arm_var, item)]
        group_bl_ceiling_inter =  bladj_group_out_ceiling_inter[,c(arm_var, item)]

        group_bl_floor_inter = round(group_bl_floor_inter[,2], round_dec)
        group_bl_ceiling_inter = round(group_bl_ceiling_inter[,2], round_dec)

        anno_tab = cbind(data.frame(bl_dat[,1]),
                         group_bl_floor_inter,
                         group_bl_ceiling_inter)

        names(anno_tab) = c("",
                            "Worsening [solid area]",
                            "Improvement [striped area]")

        auc_tab_out = as.data.frame(t(anno_tab)[2:3,], stringsAsFactors = FALSE)
        names(auc_tab_out) = t(anno_tab)[1,]
        rownames(auc_tab_out) = rownames(t(anno_tab)[2:3,])

        auc_tab_out2 = auc_tab_out
        anno_tab_dup = anno_tab

        if(permute_tests == TRUE){
          auc_tab_out[] <- lapply(as.data.frame(t(anno_tab)[2:3,], stringsAsFactors = FALSE), function(x) as.numeric(as.character(x)))

          worse_diff = auc_tab_out[1,1] - auc_tab_out[1,2]
          worse_pv = sum(abs(bladj_group_out_floor_inter_auc_diffs$diff_auc) >= abs(worse_diff)) / permute_n

          improv_diff = auc_tab_out[2,1] - auc_tab_out[2,2]
          improv_pv = sum(abs(bladj_group_out_ceiling_inter_auc_diffs$diff_auc) >= abs(improv_diff)) / permute_n

          temp_i = data.frame(auc_tab_out)
          rownames(temp_i) = NULL
          AUC_table = cbind.data.frame(item = sub(pattern = sub('^.* ([[:alnum:]]+)$', ' \\1', ref_labs[ref_labs[,"name"]==item, "short_label"]),
                                                  replacement = "", ref_labs[ref_labs[,"name"]==item, "short_label"]),
                                       type = sub('^.* ([[:alnum:]]+)$', '\\1', ref_labs[ref_labs[,"name"]==item, "short_label"]),
                                       AUC_type = sub(" .*$", "", rownames(auc_tab_out)),
                                       temp_i,
                                       P = c(worse_pv, improv_pv))


          AUC_table_display = AUC_table
          AUC_table_display$p_fmt = ifelse(AUC_table_display$P == 1, "1.0000",
                                           ifelse(AUC_table_display$P < 0.0001, "<0.0001",
                                                  paste0(formatC(round(AUC_table_display$P, digits = 4), format = 'f', flag='0', digits = 4))))
          AUC_table_display = AUC_table_display[,-6]
          names(AUC_table_display)[6] = "P"

          anno_tab = data.frame(t(AUC_table_display[,4:6]))
          anno_tab = cbind.data.frame(names = rownames(anno_tab),anno_tab)

          rownames(anno_tab) = NULL
          names(anno_tab) = c("",
                              "Worsening [solid area]",
                              "Improvement [striped area]")
        }


        if(bootstrap_ci == TRUE){

          if(!(0 < bootstrap_ci_alpha & bootstrap_ci_alpha < 1)){
            bootstrap_ci_alpha = 0.05
          }

          auc_tab_out2[] <- lapply(as.data.frame(t(anno_tab_dup)[2:3,], stringsAsFactors = FALSE), function(x) as.numeric(as.character(x)))

          worse_diff = auc_tab_out2[1,1] - auc_tab_out2[1,2]
          improv_diff = auc_tab_out2[2,1] - auc_tab_out2[2,2]

          worse_boots = stats::reshape(bladj_group_out_floor_inter_boot, timevar=arm_var, idvar=c("iteration"), direction="wide")
          worse_boots$diff = worse_boots[,2] - worse_boots[,3]
          worse_boots$resid = worse_boots$diff - worse_diff
          worse_delta = stats::quantile(worse_boots$resid, c(bootstrap_ci_alpha/2 , 1 - bootstrap_ci_alpha/2))

          improv_boots = stats::reshape(bladj_group_out_ceiling_inter_boot, timevar=arm_var, idvar=c("iteration"), direction="wide")
          improv_boots$diff = improv_boots[,2] - improv_boots[,3]
          improv_boots$resid = improv_boots$diff - improv_diff
          improv_delta = stats::quantile(improv_boots$resid, c(bootstrap_ci_alpha/2 , 1 - bootstrap_ci_alpha/2))

          diff_ci=rbind(paste0(round(worse_diff,round_dec)," (", round(worse_diff-worse_delta[2],round_dec),", ",round(worse_diff-worse_delta[1],round_dec),")"),
                        paste0(round(improv_diff,round_dec)," (", round(improv_diff-improv_delta[2],round_dec),", ",round(improv_diff-improv_delta[1],round_dec),")"))

          temp_i = data.frame(auc_tab_out2)
          rownames(temp_i) = NULL
          AUC_CI = cbind.data.frame(item = sub(pattern = sub('^.* ([[:alnum:]]+)$', ' \\1', ref_labs[ref_labs[,"name"]==item, "short_label"]),
                                               replacement = "", ref_labs[ref_labs[,"name"]==item, "short_label"]),
                                    type = sub('^.* ([[:alnum:]]+)$', '\\1', ref_labs[ref_labs[,"name"]==item, "short_label"]),
                                    AUC_type = sub(" .*$", "", rownames(auc_tab_out2)),
                                    temp_i,
                                    diff_ci)
          names(AUC_CI)[6] = paste0("difference (", 100*round(1-bootstrap_ci_alpha,2),"%CI)")
          anno_tab_boot = data.frame(t(AUC_CI[,4:6]))
          anno_tab_boot = cbind.data.frame(names = rownames(anno_tab_boot),anno_tab_boot)
          rownames(anno_tab_boot) = NULL
          names(anno_tab_boot) = c("",
                                   "Worsening [solid area]",
                                   "Improvement [striped area]")

          if(permute_tests == TRUE){
            anno_tab = rbind(anno_tab, anno_tab_boot[3,])
            rownames(anno_tab) = NULL
            AUC_table = cbind.data.frame(AUC_table, ci = AUC_CI[,6])
            names(AUC_table)[7] = paste0("difference (", 100*round(1-bootstrap_ci_alpha,2),"%CI)")

          } else {
            anno_tab = anno_tab_boot
            AUC_table = AUC_CI
          }
        }

        ribbon_dat = ribbon_plot[[i]]

        ribbon_pattern_dat0 = ribbon_dat[ribbon_dat[,item] <= ribbon_dat$bl_val, ]
        ribbon_pattern_dat0 = ribbon_pattern_dat0[order(ribbon_pattern_dat0[,arm_var], ribbon_pattern_dat0[,cycle_var]),]
        ribbon_pattern_dat1 = ribbon_pattern_dat0[ribbon_pattern_dat0[,item] <= ribbon_pattern_dat0$bl_val & ribbon_pattern_dat0[,arm_var] == name_arm1, ]
        ribbon_pattern_dat2 = ribbon_pattern_dat0[ribbon_pattern_dat0[,item] <= ribbon_pattern_dat0$bl_val & ribbon_pattern_dat0[,arm_var] == name_arm2, ]

        ribbon_pattern_dat = rbind.data.frame(ribbon_pattern_dat1[-1,], ribbon_pattern_dat2[-1,])

        item_title = ref_labs[ref_labs[,"name"]==item, "short_label"]

        df1 <- data.frame(x1 = baseline_val, x2 = cycle_limit, y1 = bl_arm1, y2 = bl_arm1, arm = name_arm1)
        df2 <- data.frame(x1 = baseline_val, x2 = cycle_limit, y1 = bl_arm2, y2 = bl_arm2, arm = name_arm2)

        df = rbind.data.frame(df1, df2)
        names(df)[5] = arm_var

        # -- Adjustments in vertical placement of the AUC table within the figure
        if(is.na(tab_ymin) | is.na(tab_ymax)){
          tab_ymin = 2.75 - (4 - y_limit)
          tab_ymax = 3.95 - (4 - y_limit)
        }

        # -- Force the fill colors to be the same for all plots
        plot_auc[,arm_var] = as.factor(plot_auc[,arm_var])
        ribbon_dat[,arm_var] = as.factor(ribbon_dat[,arm_var])

        # -- Solid cols
        ribbon_dat$arm_col = ifelse(ribbon_dat[,arm_var] == name_arm1, "black", "#E69F00")
        arm_fill_cols = unique(ribbon_dat$arm_col)

        # -- Pattern cols
        ribbon_pattern_dat$arm_col = ifelse(ribbon_pattern_dat[,arm_var] == name_arm1, "black", "#E69F00")
        arm_fill_pattern_cols = unique(ribbon_pattern_dat$arm_col)

        # -- Add note for permutation tests
        if(permute_tests == TRUE & bootstrap_ci == TRUE){
          caption_lab = paste0("Dashed horizontal line indicates arm-level baseline symptomatic adverse event level\n",
                               "Positive AUC represents increased symptomatic adverse events from baseline (Worsening)\n",
                               "Negative AUC represents decreased symptomatic adverse events from baseline (Improvement)\n",
                               paste0("P values reflect two-sided permutation tests for exchangeability with ", permute_n," replications\n"),
                               paste0("Condifence intervals constructed via bootstrap with ", bootstrap_n," replications"))
        } else if(permute_tests == TRUE){
          caption_lab = paste0("Dashed horizontal line indicates arm-level baseline symptomatic adverse event level\n",
                               "Positive AUC represents increased symptomatic adverse events from baseline (Worsening)\n",
                               "Negative AUC represents decreased symptomatic adverse events from baseline (Improvement)\n",
                               paste0("P values reflect two-sided permutation tests for exchangeability with ", permute_n," replications"))
        } else if(bootstrap_ci == TRUE){
          caption_lab = paste0("Dashed horizontal line indicates arm-level baseline symptomatic adverse event level\n",
                               "Positive AUC represents increased symptomatic adverse events from baseline (Worsening)\n",
                               "Negative AUC represents decreased symptomatic adverse events from baseline (Improvement)\n",
                               paste0("Condifence intervals constructed via bootstrap with ", bootstrap_n," replications"))
        } else {
          caption_lab = paste0("Dashed horizontal line indicates arm-level baseline symptomatic adverse event level\n",
                               "Positive AUC represents increased symptomatic adverse events from baseline (Worsening)\n",
                               "Negative AUC represents decreased symptomatic adverse events from baseline (Improvement)\n")
        }

        figure_i = ggplot2::ggplot(plot_auc, ggplot2::aes_string(x=cycle_var, y=item, group=arm_var, fill=arm_var)) +
          ggplot2::geom_line(ggplot2::aes_string(color=arm_var, linetype=arm_var)) +
          ggplot2::geom_point(ggplot2::aes_string(color=arm_var)) +

          # -- Worsening
          ggplot2::geom_ribbon(data = ribbon_dat[ribbon_dat[,item] >= ribbon_dat$bl_val & ribbon_dat[,arm_var] == name_arm1, ],
                      ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, fill=arm_var), alpha = 0.4, inherit.aes = FALSE) +
          ggplot2::geom_ribbon(data = ribbon_dat[ribbon_dat[,item] >= ribbon_dat$bl_val & ribbon_dat[,arm_var] == name_arm2, ],
                      ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, fill=arm_var), alpha = 0.3, inherit.aes = FALSE) +

          ggpattern::geom_ribbon_pattern(
            data = ribbon_pattern_dat,
            ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, pattern_fill = arm_var, pattern_color = arm_var),
            inherit.aes = FALSE,
            pattern_alpha = 0.4,
            fill            = NA,
            pattern         = 'stripe',
            pattern_spacing = 0.02,
            pattern_density = 0.5,
            show.legend = FALSE
          ) +

          # --- Add horizontal line indicating baseline value
          ggplot2::geom_segment(ggplot2::aes_string(x = "x1", y = "y1", xend = "x2", yend = "y2", colour = arm_var),
                                linetype="dashed",
                                data = df,
                                inherit.aes = FALSE) +

          ggplot2::scale_x_continuous(name=cycle_var, limits=c(baseline_val, cycle_limit), breaks = baseline_val:cycle_limit) +
          ggplot2::scale_y_continuous(name="Mean Score", limits=c(0, y_limit)) +
          ggplot2::scale_linetype_manual(values=c("solid", "solid")) +
          ggplot2::scale_color_manual(values=arm_fill_cols) +
          ggplot2::scale_fill_manual(values=arm_fill_cols) +
          ggpattern::scale_pattern_color_manual(values=arm_fill_pattern_cols) +
          ggpattern::scale_pattern_fill_manual(values=arm_fill_pattern_cols) +
          ggplot2::annotation_custom(gridExtra::tableGrob(t(anno_tab)), xmin=0.38*cycle_limit, xmax=0.8125*cycle_limit, ymin=tab_ymin, ymax=tab_ymin) +
          ggplot2::theme_bw() +
          ggplot2::ggtitle(ifelse(is.na(overwrite_title),
                         paste0(item_title, ": Baseline adjusted AUC"),
                         overwrite_title)) +
          ggplot2::labs(caption = caption_lab) +

          ggplot2::theme(legend.position = c(.15,.85),
                legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                legend.title = ggplot2::element_blank(),
                plot.caption = ggplot2::element_text(hjust = 0, size = 10))


      } else if(nrow(bl_dat)==3){
        bl_arm1 = bl_dat[1,item]
        name_arm1 = bl_dat[1,arm_var]

        bl_arm2 = bl_dat[2,item]
        name_arm2 = bl_dat[2,arm_var]

        bl_arm3 = bl_dat[3,item]
        name_arm3 = bl_dat[3,arm_var]

        plot_auc$bl_val = ifelse(plot_auc[,arm_var] == name_arm1, bl_arm1,
                                 ifelse(plot_auc[,arm_var] == name_arm2, bl_arm2, bl_arm3))

        # -- AUC table
        group_bl_floor_inter = bladj_group_out_floor_inter[,c(arm_var, item)]
        group_bl_ceiling_inter =  bladj_group_out_ceiling_inter[,c(arm_var, item)]

        group_bl_floor_inter = round(group_bl_floor_inter[,2], round_dec)
        group_bl_ceiling_inter = round(group_bl_ceiling_inter[,2], round_dec)

        anno_tab = cbind(data.frame(bl_dat[,1]),
                         group_bl_floor_inter,
                         group_bl_ceiling_inter)

        names(anno_tab) = c("",
                            "Worsening [solid area]",
                            "Improvement [striped area]")

        ribbon_dat = ribbon_plot[[i]]

        ribbon_pattern_dat0 = ribbon_dat[ribbon_dat[,item] <= ribbon_dat$bl_val, ]
        ribbon_pattern_dat0 = ribbon_pattern_dat0[order(ribbon_pattern_dat0[,arm_var], ribbon_pattern_dat0[,cycle_var]),]
        ribbon_pattern_dat1 = ribbon_pattern_dat0[ribbon_pattern_dat0[,item] <= ribbon_pattern_dat0$bl_val & ribbon_pattern_dat0[,arm_var] == name_arm1, ]
        ribbon_pattern_dat2 = ribbon_pattern_dat0[ribbon_pattern_dat0[,item] <= ribbon_pattern_dat0$bl_val & ribbon_pattern_dat0[,arm_var] == name_arm2, ]
        ribbon_pattern_dat3 = ribbon_pattern_dat0[ribbon_pattern_dat0[,item] <= ribbon_pattern_dat0$bl_val & ribbon_pattern_dat0[,arm_var] == name_arm3, ]

        ribbon_pattern_dat = rbind.data.frame(ribbon_pattern_dat1[-1,], ribbon_pattern_dat2[-1,], ribbon_pattern_dat3[-1,])

        item_title = ref_labs[ref_labs[,"name"]==item, "short_label"]

        df1 <- data.frame(x1 = baseline_val, x2 = cycle_limit, y1 = bl_arm1, y2 = bl_arm1, arm = name_arm1)
        df2 <- data.frame(x1 = baseline_val, x2 = cycle_limit, y1 = bl_arm2, y2 = bl_arm2, arm = name_arm2)
        df3 <- data.frame(x1 = baseline_val, x2 = cycle_limit, y1 = bl_arm3, y2 = bl_arm3, arm = name_arm3)

        df = rbind.data.frame(df1, df2, df3)
        names(df)[5] = arm_var

        # -- Adjustments in vertical placement of the AUC table within the figure
        if(is.na(tab_ymin) | is.na(tab_ymax)){
          tab_ymin = 2.75 - (4 - y_limit)
          tab_ymax = 3.95 - (4 - y_limit)
        }

        # -- Force the fill colors to be the same for all plots
        plot_auc[,arm_var] = as.factor(plot_auc[,arm_var])
        ribbon_dat[,arm_var] = as.factor(ribbon_dat[,arm_var])

        # -- Solid cols
        ribbon_dat$arm_col = ifelse(ribbon_dat[,arm_var] == name_arm1, "black",
                                    ifelse(ribbon_dat[,arm_var] == name_arm2, "#E69F00", "#0A7B83"))

        arm_fill_cols = unique(ribbon_dat$arm_col)

        # -- Pattern cols
        ribbon_pattern_dat$arm_col = ifelse(ribbon_pattern_dat[,arm_var] == name_arm1, "black",
                                            ifelse(ribbon_pattern_dat[,arm_var] == name_arm2, "#E69F00", "#0A7B83"))


        arm_fill_pattern_cols = unique(ribbon_pattern_dat$arm_col)

        figure_i = ggplot2::ggplot(plot_auc, ggplot2::aes_string(x=cycle_var, y=item, group=arm_var, fill=arm_var)) +
          ggplot2::geom_line(ggplot2::aes_string(color=arm_var, linetype=arm_var)) +
          ggplot2::geom_point(ggplot2::aes_string(color=arm_var)) +
          # -- Worsening
          ggplot2::geom_ribbon(data = ribbon_dat[ribbon_dat[,item] >= ribbon_dat$bl_val & ribbon_dat[,arm_var] == name_arm1, ],
                               ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, fill=arm_var), alpha = 0.4, inherit.aes = FALSE) +
          ggplot2::geom_ribbon(data = ribbon_dat[ribbon_dat[,item] >= ribbon_dat$bl_val & ribbon_dat[,arm_var] == name_arm2, ],
                               ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, fill=arm_var), alpha = 0.3, inherit.aes = FALSE) +
          ggplot2::geom_ribbon(data = ribbon_dat[ribbon_dat[,item] >= ribbon_dat$bl_val & ribbon_dat[,arm_var] == name_arm3, ],
                               ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, fill=arm_var), alpha = 0.3, inherit.aes = FALSE) +

          ggpattern::geom_ribbon_pattern(
            data = ribbon_pattern_dat,
            ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, pattern_fill = arm_var, pattern_color = arm_var),
            inherit.aes = FALSE,
            pattern_alpha = 0.4,
            fill            = NA,
            pattern         = 'stripe',
            pattern_spacing = 0.02,
            pattern_density = 0.5,
            show.legend = FALSE
          ) +

          # --- Add horizontal line indicating baseline value
          ggplot2::geom_segment(ggplot2::aes_string(x = "x1", y = "y1", xend = "x2", yend = "y2", colour = arm_var),
                                linetype="dashed",
                                data = df,
                                inherit.aes = FALSE) +

          ggplot2::scale_x_continuous(name=cycle_var, limits=c(baseline_val, cycle_limit), breaks = baseline_val:cycle_limit) +
          ggplot2::scale_y_continuous(name="Mean Score", limits=c(0, y_limit)) +
          ggplot2::scale_linetype_manual(values=c("solid", "solid", "solid")) +
          ggplot2::scale_color_manual(values=arm_fill_cols) +
          ggplot2::scale_fill_manual(values=arm_fill_cols) +
          ggpattern::scale_pattern_color_manual(values=arm_fill_pattern_cols) +
          ggpattern::scale_pattern_fill_manual(values=arm_fill_pattern_cols) +
          ggplot2::annotation_custom(gridExtra::tableGrob(t(anno_tab)), xmin=0.375*cycle_limit, xmax=0.8125*cycle_limit, ymin=tab_ymin, ymax=tab_ymin) +
          ggplot2::theme_bw() +
          ggplot2::ggtitle(ifelse(is.na(overwrite_title),
                         paste0(item_title, ": Baseline adjusted AUC"),
                         overwrite_title)) +
          ggplot2::labs(caption = paste0("Dashed horizontal line indicates arm-level baseline symptomatic adverse event level\n",
                                "Positive AUC represents increased symptomatic adverse events from baseline (Worsening)\n",
                                "Negative AUC represents decreased symptomatic adverse events from baseline (Improvement)\n")) +

          ggplot2::theme(legend.position = c(.15,.85),
                legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                legend.title = ggplot2::element_blank(),
                plot.caption = ggplot2::element_text(hjust = 0, size = 10))

        auc_tab_out = as.data.frame(t(anno_tab)[2:3,], stringsAsFactors = FALSE)
        names(auc_tab_out) = t(anno_tab)[1,]
        rownames(auc_tab_out) = rownames(t(anno_tab)[2:3,])
      } else if(nrow(bl_dat)==4){
        bl_arm1 = bl_dat[1,item]
        name_arm1 = bl_dat[1,arm_var]

        bl_arm2 = bl_dat[2,item]
        name_arm2 = bl_dat[2,arm_var]

        bl_arm3 = bl_dat[3,item]
        name_arm3 = bl_dat[3,arm_var]

        bl_arm4 = bl_dat[4,item]
        name_arm4 = bl_dat[4,arm_var]

        plot_auc$bl_val = ifelse(plot_auc[,arm_var] == name_arm1, bl_arm1,
                                 ifelse(plot_auc[,arm_var] == name_arm2, bl_arm2,
                                        ifelse(plot_auc[,arm_var] == name_arm3, bl_arm3, bl_arm4)))

        # -- AUC table
        group_bl_floor_inter = bladj_group_out_floor_inter[,c(arm_var, item)]
        group_bl_ceiling_inter =  bladj_group_out_ceiling_inter[,c(arm_var, item)]

        group_bl_floor_inter = round(group_bl_floor_inter[,2], round_dec)
        group_bl_ceiling_inter = round(group_bl_ceiling_inter[,2], round_dec)

        anno_tab = cbind(data.frame(bl_dat[,1]),
                         group_bl_floor_inter,
                         group_bl_ceiling_inter)

        names(anno_tab) = c("",
                            "Worsening [solid area]",
                            "Improvement [striped area]")

        ribbon_dat = ribbon_plot[[i]]

        ribbon_pattern_dat0 = ribbon_dat[ribbon_dat[,item] <= ribbon_dat$bl_val, ]
        ribbon_pattern_dat0 = ribbon_pattern_dat0[order(ribbon_pattern_dat0[,arm_var], ribbon_pattern_dat0[,cycle_var]),]
        ribbon_pattern_dat1 = ribbon_pattern_dat0[ribbon_pattern_dat0[,item] <= ribbon_pattern_dat0$bl_val & ribbon_pattern_dat0[,arm_var] == name_arm1, ]
        ribbon_pattern_dat2 = ribbon_pattern_dat0[ribbon_pattern_dat0[,item] <= ribbon_pattern_dat0$bl_val & ribbon_pattern_dat0[,arm_var] == name_arm2, ]
        ribbon_pattern_dat3 = ribbon_pattern_dat0[ribbon_pattern_dat0[,item] <= ribbon_pattern_dat0$bl_val & ribbon_pattern_dat0[,arm_var] == name_arm3, ]
        ribbon_pattern_dat4 = ribbon_pattern_dat0[ribbon_pattern_dat0[,item] <= ribbon_pattern_dat0$bl_val & ribbon_pattern_dat0[,arm_var] == name_arm4, ]

        ribbon_pattern_dat = rbind.data.frame(ribbon_pattern_dat1[-1,], ribbon_pattern_dat2[-1,], ribbon_pattern_dat3[-1,], ribbon_pattern_dat4[-1,])

        item_title = ref_labs[ref_labs[,"name"]==item, "short_label"]

        df1 <- data.frame(x1 = baseline_val, x2 = cycle_limit, y1 = bl_arm1, y2 = bl_arm1, arm = name_arm1)
        df2 <- data.frame(x1 = baseline_val, x2 = cycle_limit, y1 = bl_arm2, y2 = bl_arm2, arm = name_arm2)
        df3 <- data.frame(x1 = baseline_val, x2 = cycle_limit, y1 = bl_arm3, y2 = bl_arm3, arm = name_arm3)
        df4 <- data.frame(x1 = baseline_val, x2 = cycle_limit, y1 = bl_arm4, y2 = bl_arm4, arm = name_arm4)

        df = rbind.data.frame(df1, df2, df3, df4)
        names(df)[5] = arm_var

        # -- Adjustments in vertical placement of the AUC table within the figure
        if(is.na(tab_ymin) | is.na(tab_ymax)){
          tab_ymin = 2.75 - (4 - y_limit)
          tab_ymax = 3.95 - (4 - y_limit)
        }

        # -- Force the fill colors to be the same for all plots
        plot_auc[,arm_var] = as.factor(plot_auc[,arm_var])
        ribbon_dat[,arm_var] = as.factor(ribbon_dat[,arm_var])

        # -- Solid cols
        ribbon_dat$arm_col = ifelse(ribbon_dat[,arm_var] == name_arm1, "black",
                                    ifelse(ribbon_dat[,arm_var] == name_arm2, "#E69F00",
                                           ifelse(ribbon_dat[,arm_var] == name_arm3, "#0A7B83", "#CE4D45")))

        arm_fill_cols = unique(ribbon_dat$arm_col)

        # -- Pattern cols
        ribbon_pattern_dat$arm_col = ifelse(ribbon_pattern_dat[,arm_var] == name_arm1, "black",
                                            ifelse(ribbon_pattern_dat[,arm_var] == name_arm2, "#E69F00",
                                                   ifelse(ribbon_pattern_dat[,arm_var] == name_arm3, "#0A7B83", "#CE4D45")))


        arm_fill_pattern_cols = unique(ribbon_pattern_dat$arm_col)

        figure_i = ggplot2::ggplot(plot_auc, ggplot2::aes_string(x=cycle_var, y=item, group=arm_var, fill=arm_var)) +
          ggplot2::geom_line(ggplot2::aes_string(color=arm_var, linetype=arm_var)) +
          ggplot2::geom_point(ggplot2::aes_string(color=arm_var)) +
          # -- Worsening
          ggplot2::geom_ribbon(data = ribbon_dat[ribbon_dat[,item] >= ribbon_dat$bl_val & ribbon_dat[,arm_var] == name_arm1, ],
                      ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, fill=arm_var), alpha = 0.4, inherit.aes = FALSE) +
          ggplot2::geom_ribbon(data = ribbon_dat[ribbon_dat[,item] >= ribbon_dat$bl_val & ribbon_dat[,arm_var] == name_arm2, ],
                      ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, fill=arm_var), alpha = 0.3, inherit.aes = FALSE) +
          ggplot2::geom_ribbon(data = ribbon_dat[ribbon_dat[,item] >= ribbon_dat$bl_val & ribbon_dat[,arm_var] == name_arm3, ],
                      ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, fill=arm_var), alpha = 0.3, inherit.aes = FALSE) +
          ggplot2::geom_ribbon(data = ribbon_dat[ribbon_dat[,item] >= ribbon_dat$bl_val & ribbon_dat[,arm_var] == name_arm4, ],
                      ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, fill=arm_var), alpha = 0.3, inherit.aes = FALSE) +

          ggpattern::geom_ribbon_pattern(
            data = ribbon_pattern_dat,
            ggplot2::aes_string(x = cycle_var, ymin = "bl_val", ymax = item, pattern_fill = arm_var, pattern_color = arm_var),
            inherit.aes = FALSE,
            pattern_alpha = 0.4,
            fill            = NA,
            pattern         = 'stripe',
            pattern_spacing = 0.02,
            pattern_density = 0.5,
            show.legend = FALSE
          ) +

          # --- Add horizontal line indicating baseline value
          ggplot2::geom_segment(ggplot2::aes_string(x = "x1", y = "y1", xend = "x2", yend = "y2", colour = arm_var),
                       linetype="dashed",
                       data = df,
                       inherit.aes = FALSE) +

          ggplot2::scale_x_continuous(name=cycle_var, limits=c(baseline_val, cycle_limit), breaks = baseline_val:cycle_limit) +
          ggplot2::scale_y_continuous(name="Mean Score", limits=c(0, y_limit)) +
          ggplot2::scale_linetype_manual(values=c("solid", "solid", "solid", "solid")) +
          ggplot2::scale_color_manual(values=arm_fill_cols) +
          ggplot2::scale_fill_manual(values=arm_fill_cols) +
          ggpattern::scale_pattern_color_manual(values=arm_fill_pattern_cols) +
          ggpattern::scale_pattern_fill_manual(values=arm_fill_pattern_cols) +
          ggplot2::annotation_custom(gridExtra::tableGrob(t(anno_tab)), xmin=0.375*cycle_limit, xmax=0.8125*cycle_limit, ymin=tab_ymin, ymax=tab_ymin) +
          ggplot2::theme_bw() +
          ggplot2::ggtitle(ifelse(is.na(overwrite_title),
                         paste0(item_title, ": Baseline adjusted AUC"),
                         overwrite_title)) +
          ggplot2::labs(caption = paste0("Dashed horizontal line indicates arm-level baseline symptomatic adverse event level\n",
                                "Positive AUC represents increased symptomatic adverse events from baseline (Worsening)\n",
                                "Negative AUC represents decreased symptomatic adverse events from baseline (Improvement)\n")) +

          ggplot2::theme(legend.position = c(.15,.85),
                legend.background = ggplot2::element_rect(fill = "white", color = "black"),
                legend.title = ggplot2::element_blank(),
                plot.caption = ggplot2::element_text(hjust = 0, size = 10))

        auc_tab_out = as.data.frame(t(anno_tab)[2:3,], stringsAsFactors = FALSE)
        names(auc_tab_out) = t(anno_tab)[1,]
        rownames(auc_tab_out) = rownames(t(anno_tab)[2:3,])
      }
          # ------------------------------------------------------------------------------
          # --- Allow for permutation tests for the difference in AUC between 2 arms (iff 2 arms)
          # ------------------------------------------------------------------------------

      if(nrow(bl_dat) == 2 & permute_tests == TRUE){

        list_out[[i]] = list()
        list_out[[i]][[1]] = ref_labs[ref_labs[,"name"]==item, "short_label"]
        list_out[[i]][[2]] = figure_i
        list_out[[i]][[3]] = AUC_table
      } else{

        auc_tab_out[] <- lapply(auc_tab_out, function(x) as.numeric(as.character(x)))

        list_out[[i]] = list()
        list_out[[i]][[1]] = ref_labs[ref_labs[,"name"]==item, "short_label"]
        list_out[[i]][[2]] = figure_i
        list_out[[i]][[3]] = auc_tab_out
      }

    }

    ## -- Reference table
    ref_out = ref_labs[,c("index", "short_label")]
    rownames(ref_out) = NULL

    ## -- Reference table
    print(ref_out)

    ## -- Object return
    invisible(list_out)

}


