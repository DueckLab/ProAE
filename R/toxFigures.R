#' Create PRO-CTCAE severity frequency distribution figures for individual
#' survey items and composite scores
#'
#'
#' 	  Data format should be in 'long' format, where each PRO-CTCAE item is a
#' 	  variable/column.
#'
#' @param dsn A data.frame object with PRO-CTCAE data
#' @param id_var A character string.Name of ID variable differentiating each
#'   PRO-CTCAE survey/participant entered as a quoted string.
#' @param cycle_var A character string. Name of variable differentiating one
#'   longitudinal/repeated PRO-CTCAE survey from another, within an individual
#'   ID.
#' @param baseline_val A number indicating the expected baseline cycle/time
#'   point.
#' @param arm_var A character string. Name of arm variable differentiating
#'   treatment groups. Must be character or factor class. Overall AUC
#'   will be reported if no arm/grouping variable is provided. Defaults to
#'   \code{NA}.
#' @param plot_limit A number. Limit the number of cycles to be plotted up to
#'   and including a given cycle number. All available cycle time points are
#'   plotted if no cycle number is provided. Defaults to \code{NA}.
#' @param colors A number. Specify the coloring scheme of symptom scores within
#'   frequency bars. Options include: 1 = Blue and red color shading, 2 =
#'   qualitative color shades (color blind friendly), 3 = black and white.
#'   Defaults to 1.
#' @param bar_label A number. Label frequency bars with sample size (n) or percent
#'   shown on the y-axis. Label options include: \code{1} = sample size (n)
#'   within each cycle (symptom score 0 or higher), \code{2} = sample size (n)
#'   within each cycle with present symptoms (symptom score > 0), \code{3} =
#'   sample size (n) within each cycle with severe symptoms (symptom score >=
#'   3), \code{4} = percent of subjects within each cycle with present
#'   symptoms (symptom score > 0), \code{5} = percent of subjects within each
#'   cycle with severe symptoms (symptom score >= 3). No labels will be applied
#'   if not specified. Defaults to \code{NA}.
#' @param cycle_label Logical. Assign custom labels to cycles/time point. If
#'   \code{TRUE}, the \code{cycle_vals} and \code{cycle_labs} must also be specified.
#' @param cycle_vals Numeric column vector. Vector of values seen within the
#'   \code{cycle_var} variable. Must be same length of \code{cycle_labs}. Defaults
#'   to \code{NA}.
#' @param cycle_labs Character column vector. Vector of labels to be mapped to
#'   the associated \code{cycle_vals}. Must be same length of \code{cycle_vals}.
#'   Defaults to \code{NA}.
#' @param summary_only Logical. Only display the summary measures in figures /
#'   Suppress the individual time points from plotting. Defaults to
#'   \code{FALSE}.
#' @param cycles_only Logical. Only display the longitudinal time points in
#'   figures / Suppress the summary measures from plotting. Defaults to
#'   \code{TRUE}.
#' @param x_label A character string. Label for the x axis of the plot. Defaults
#'   to \code{"Randomized Treatment Assignment"} if \code{arm_var} is specified,
#'   defaults to \code{"Overall"} if not \code{arm_var} is specified.
#' @param y_label A character string. Label for the y axis of the plot. Defaults
#'   to \code{"Percent of Total Frequency"}.
#' @param x_lab_angle A integer between 0 and 360. Allows the user to rotate the
#'   x axis labels in order to fit long arm names (0 or 45 recommended).
#'   Defaults to \code{0}.
#' @param x_lab_vjust A number. A ggplot2 object option. Allows the user to
#'   vertically adjusts the x axis labels in order to fit arm names. Defaults to
#'   \code{1}.
#' @param x_lab_hjust A number. A ggplot2 object option. Allows the user to
#'   horizontally adjusts the x axis labels in order to fit arm names. Defaults
#'   to \code{0}.
#' @param footnote_break Logical. Add a line break to the footnote Defaults to
#'   \code{FALSE}.
#' @param suppress_legend Logical. Suppress the legend from appearing in figure.
#'    Defaults to \code{FALSE}.
#' @param add_item_title Logical. Add the symptom item name as a title to the
#'    figure. Defaults to \code{FALSE}.
#' @param summary_highlight Logical. Add black box around summary measure bar
#'    chart. Defaults to \code{FALSE}.
#' @return A list object. The returned object is a (k X 2) or (k x 3) nested
#'   list. Where k is the number of PRO-CTCAE item groups (e.g. pain, fatigue,
#'   nausea); list[[1 ... i ... k]]. For each list item there are 2 or 3
#'   elements. The 1st element of each list item is the name of the PRO-CTCAE
#'   item group returned as a string. The 2nd element is the PRO-CTCAE figure as
#'   a ggplot object. These objects can be modified as such.
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' fig_acute = toxFigures(dsn = ProAE::tox_acute,
#'  cycle_var = "Cycle",
#'  baseline_val = 1,
#'  arm_var = "arm",
#'  id_var = "id",
#'  x_lab_angle = -45,
#'  x_lab_vjust =  .3,
#'  x_lab_hjust = .2,
#'  colors = 2)
#' fig_acute[[1]]
#' }
#' @export

toxFigures = function(dsn,
                     id_var,
                     cycle_var,
                     baseline_val,
                     arm_var=NA,
                     plot_limit=NA,
                     colors=1,
                     bar_label=0,
                     cycle_label = FALSE,
                     cycle_vals = NA,
                     cycle_labs = NA,
                     summary_only=FALSE,
                     summary_highlight = FALSE,
                     cycles_only=TRUE,
                     x_lab_angle=0,
                     x_lab_vjust=1,
                     x_lab_hjust=0,
                     x_label = "Randomized Treatment Assignment",
                     y_label = "Percent of Total Frequency",
                     footnote_break = FALSE,
                     suppress_legend = FALSE,
                     add_item_title = FALSE){

  # ----------------------------------------------------------------
  # -- Checks 1/2
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

  if(exists("baseline_val")){
    if(!(is.numeric(baseline_val) | is.integer(baseline_val) | length(baseline_val)==1)){
      stop("param baseline_val must be provided as a single number, of class numeric or integer")
    }
  } else {stop("param baseline_val not provided")}

  ## -- Check for any duplicate individuals within cycles

  if(any(duplicated(dsn[,c(id_var, cycle_var)]))){
    stop(paste0("Duplicate observations were found within id_var and cycle_var combinations (", id_var, " and ", cycle_var, ")"))
  }

  ## -- Check for conflicts with figure options
  if(cycle_label == TRUE){
    if(any(is.na(cycle_vals) | is.na(cycle_labs))){
      stop("params cycle_vals and cycle_labs must be provided while using cycle_label")
    } else if(length(cycle_vals) != length(cycle_labs)){
      stop("params cycle_vals and cycle_labs must be the same length column vector")
    } else if(!is.numeric(cycle_vals) | !is.character(cycle_labs)){
      stop("params cycle_vals and cycle_labs must be numeric and character vectors, respectively")
    }
  }

  # ----------------------------------------------------------------
  # -- Defaults
  # ----------------------------------------------------------------

  if(is.na(arm_var)){
    arm_var = "overall_"
    if(x_label=="Randomized Treatment Assignment"){
      x_label = "Overall"
    }
  }

  if(summary_only==TRUE & cycles_only==TRUE){
    summary_only = TRUE
  }

  ## -------------------------------------------------------------
  ## --- Reference data sets
  ## -------------------------------------------------------------

  proctcae_vars[, ] = lapply(proctcae_vars[, ], as.character)

  # ================================================================
  # ==== ACCOMMODATION FOR 27A HAIR LOSS AMOUNT ====================
  # ================================================================
  ## -- Item 27 asks about hair loss and is referenced as an 'amount'
  ## -- However this item utilizes the Interference response bank.
  ## -- This is accommodated here by treating it as Interference then
  ## -- resulting plot is manually labeled 'Amount' BTL - 20AUG2020

  proctcae_vars[proctcae_vars$name=="PROCTCAE_27A_SCL",]$short_label = "Hair Loss Interference"

  # ================================================================
  # ================================================================
  # ================================================================


  # ----------------------------------------------------------------
  # -- Get existing PROCTCAE variables in dsn
  # ----------------------------------------------------------------

  ## - Individual items
  dsn_items0 = toupper(names(dsn)[toupper(names(dsn)) %in% proctcae_vars$name])
  dsn_items = dsn_items0[! dsn_items0 %in% as.character(proctcae_vars$name[proctcae_vars$fmt %in% c("yn_2_fmt", "yn_3_fmt", "yn_4_fmt")])]

  ## - Composites
  proctcae_vars_comp0 = proctcae_vars %>% dplyr::mutate_if(is.factor, as.character)

  proctcae_vars_comp0 = proctcae_vars_comp0[!proctcae_vars_comp0$name %in% as.character(proctcae_vars_comp0$name[proctcae_vars_comp0$fmt %in% c("yn_2_fmt",
                                                                                                                                                "yn_3_fmt",
                                                                                                                                                "yn_4_fmt",
                                                                                                                                                "gp5_fmt")]),]
  proctcae_vars_comp = c()
  proctcae_vars_comp$name = paste0(substr(proctcae_vars_comp0$name, 1, nchar(proctcae_vars_comp0$name)-5), "_COMP")
  proctcae_vars_comp$short_label = sub(proctcae_vars_comp0$short_label, pattern = " [[:alpha:]]*$", replacement = "")
  proctcae_vars_comp = unique(data.frame(proctcae_vars_comp, stringsAsFactors=FALSE))
  dsn_comps = toupper(names(dsn)[toupper(names(dsn)) %in% proctcae_vars_comp$name])

  # ----------------------------------------------------------------
  # -- Checks 2/2
  # ----------------------------------------------------------------

  ## -- Confirm there are available PRO-CTCAE variables within dsn with expected naming convention
  if(identical(dsn_items, character(0)) & identical(dsn_comps, character(0))){
    stop(paste0("No compatible ProAE variables found within dsn (", deparse(substitute(dsn)), ") meeting the expected naming convention"))
  }

  # ----------------------------------------------------------------
  # -- Cluster available items and composites together
  # ----------------------------------------------------------------

  ## -- Available PRO-CTCAE items and composites available in dsn
  if(any(dsn_items %in% "GP5")){
    refset = data.frame(name = c(dsn_items[!(dsn_items %in% "GP5")], dsn_comps, dsn_items[(dsn_items %in% "GP5")]))
    refset$number = ifelse(!refset$name %in% "GP5",
                           as.integer(gsub("[^0-9]", "", refset$name)),
                           999)
  } else {
    refset = data.frame(name = c(dsn_items, dsn_comps))
    refset$number = as.integer(gsub("[^0-9]", "", refset$name))
  }

  refset$rank = ifelse(grepl("COMP", refset$name), 4,
                       ifelse(grepl("A_", refset$name), 1,
                              ifelse(grepl("B_", refset$name), 2, 3)))

  refset = refset[order(refset$number,refset$rank),]
  refset = merge(x=refset,
                 y= data.frame(number = unique(refset$number),
                               group_rank = 1:length(unique(refset$number)),
                               stringsAsFactors=FALSE),
                 by="number", all.x=TRUE)

  ## -- All non-indicator PRO-CTCAE with possible composites

  reflabs0 = data.frame(name=proctcae_vars$name,
                        short_label = proctcae_vars$short_label,
                        lab = unlist(lapply(strsplit(proctcae_vars$short_label, split=" "), 1, FUN=utils::tail)),
                        stringsAsFactors=FALSE)

  reflabs0$group_lab = ifelse(reflabs0$lab=="Severity", trimws(gsub("Severity","", reflabs0$short_label)),
                              ifelse(reflabs0$lab=="Frequency", trimws(gsub("Frequency","", reflabs0$short_label)),
                                     trimws(gsub("Interference","", reflabs0$short_label))))
  reflabs0 = reflabs0[!grepl("_IND", reflabs0$name),]
  reflabs0$comp_name = paste0(substr(reflabs0$name,1,nchar(reflabs0$name)-5),"_COMP")
  reflabs1 = unique(reflabs0[c("group_lab", 'comp_name')])
  names(reflabs1) = c("group_lab","name")
  reflabs1$lab = "Composite"
  reflabs = rbind.data.frame(reflabs0[c("name", "group_lab", "lab")],
                             reflabs1[c("name", "group_lab", "lab")],
                             stringsAsFactors=FALSE)

  refset = merge(x=refset,
                 y=reflabs,
                 by= "name", all.x=TRUE)
  refset = refset[order(refset$number,refset$rank),]


  # ----------------------------------------------------------------
  # -- For each available ProAE variable grouping [i], process
  # ----------------------------------------------------------------

  list_out = vector(mode='list', max(refset$group_rank))

  # ----------------------------------------------------------------
  # --  PRO-CTCAE items
  # ----------------------------------------------------------------

  if(any(c(dsn_comps, dsn_items) %in% c(proctcae_vars_comp$name, proctcae_vars[proctcae_vars$name!="GP5",]$name))){
    for (i in 1:max(refset[refset$name!="GP5",]$group_rank)){

      # ----------------------------------------------------------------
      # -- For each PRO-CTCAE item within the grouping [j], process
      # ----------------------------------------------------------------
      group_i = refset[refset$group_rank==i,]
      plot_combined0 = data.frame()
      for (j in 1:NROW(group_i)){
        item = as.character(group_i$name[j])

        if(arm_var != "overall_"){
          dsn0 = stats::na.omit(dsn[,c(id_var, cycle_var, arm_var, item)])
        } else if (arm_var == "overall_"){
          dsn0 = stats::na.omit(dsn[,c(id_var, cycle_var, item)])
          dsn0$overall_ = "Overall"
        }


        # ----------------------------------------------------------------
        # -- Construct summary measures
        # ----------------------------------------------------------------
        ## -- Maximum grade post baseline
        max_post0 = stats::aggregate(dsn0[dsn0[,cycle_var]>baseline_val, item],
                                     by = list(dsn0[dsn0[,cycle_var]>baseline_val,id_var]),
                                     FUN=max)
        colnames(max_post0) = c(id_var, item)
        max_post0[,cycle_var] = "Maximum*"
        max_post1 = unique(merge(x=max_post0, y=dsn0[,c(id_var, arm_var)], by=id_var, all.x = TRUE))

        ## -- Baseline adjusted
        base_adj0 = dsn0[dsn0[,cycle_var]==baseline_val, c(id_var, item)]
        colnames(base_adj0)[2] = "base_val"

        base_adj1 = merge(x=dsn0, y=base_adj0, by=id_var, all.x=TRUE)
        base_adj2 = stats::aggregate(base_adj1[base_adj1[,cycle_var]>baseline_val, item],
                                     by = list(base_adj1[base_adj1[,cycle_var]>baseline_val,id_var]),
                                     FUN=max)
        colnames(base_adj2) = c(id_var, "max_post_bl")
        base_adj2 = merge(x=base_adj1, y=unique(base_adj2), by=id_var, all.x=TRUE)
        base_adj2$bl_adjusted = ifelse(!is.na(base_adj2$max_post_bl) & base_adj2$base_val >= base_adj2$max_post_bl,
                                       0, base_adj2$max_post_bl)
        base_adj3 = unique(base_adj2[!is.na(base_adj2$bl_adjusted),c(id_var, arm_var, "bl_adjusted")])
        base_adj3[,cycle_var] = "Adjusted**"
        colnames(base_adj3)[colnames(base_adj3) == "bl_adjusted"] = item


        dsn1 = rbind.data.frame(dsn0, max_post1, base_adj3)

        # ----------------------------------------------------------------
        # -- Create labels for the bars
        # ----------------------------------------------------------------

        if(bar_label==0){
          labs = dsn1[,c(cycle_var, arm_var)]
          labs = unique(labs[c(cycle_var, arm_var)])
          labs$bar_lab_opt = NA

        } else if(bar_label==1){
          labs = stats::aggregate(dsn1[dsn1[,item]>=0,item],
                                  by=list(dsn1[dsn1[,item]>=0,cycle_var],
                                          dsn1[dsn1[,item]>=0,arm_var]),
                                  FUN=length)
        } else if(bar_label %in% c(2,4)){
          labs = stats::aggregate(dsn1[dsn1[,item]>0,item],
                                  by=list(dsn1[dsn1[,item]>0,cycle_var],
                                          dsn1[dsn1[,item]>0,arm_var]),
                                  FUN=length)
        } else if(bar_label %in% c(3,5)){
          labs = stats::aggregate(dsn1[dsn1[,item]>=3,item],
                                  by=list(dsn1[dsn1[,item]>=3,cycle_var],
                                          dsn1[dsn1[,item]>=3,arm_var]),
                                  FUN=length)
        }
        names(labs) = c(cycle_var, arm_var, "bar_lab_opt")

        # -- Force zero counts
        combo_counts = data.frame(table(dsn1[,cycle_var], dsn1[,arm_var]), stringsAsFactors=FALSE)
        names(combo_counts) = c(cycle_var, arm_var, "total")
        combo_counts$flag = 1

        labs2 = merge(x=labs, y=combo_counts, by=c(cycle_var, arm_var), all.y = TRUE)

        if(any(is.na(labs2$bar_lab_opt))){
          labs2[is.na(labs2$bar_lab_opt),]$bar_lab_opt = 0
        }

        if(bar_label %in% c(4,5)){
          labs2[labs2$bar_lab_opt!=0,]$bar_lab_opt=round(labs2[labs2$bar_lab_opt!=0,]$bar_lab_opt/labs2[labs2$bar_lab_opt!=0,]$total,digits=2)*100
        }

        dsn2 = merge(dsn1, labs2, by = c(cycle_var, arm_var))

        # -- Only one label per cycle/arm combination
        high = by(dsn2, INDICES = list(dsn2[,cycle_var], dsn2[,arm_var]), FUN = utils::tail, n=1)
        high2 = do.call("rbind", as.list(high))
        row_keeps = as.integer(rownames(high2))
        dsn2[!(1:NROW(dsn2) %in% row_keeps),]$bar_lab_opt = NA

        # ----------------------------------------------------------------
        # -- Create labels for the cycle_var facets
        # ----------------------------------------------------------------
        dsn2$cycle_var_v_plot=factor(ifelse(dsn2[,cycle_var]=="Maximum*",
                                            "Maximum*",
                                            ifelse(dsn2[,cycle_var]=="Adjusted**", "Adjusted**",
                                                   paste0(cycle_var," ", dsn2[,cycle_var]))),
                                     levels=c(paste0(cycle_var," ", min(dsn0[,cycle_var]):max(dsn0[,cycle_var])),
                                              "Maximum*",
                                              "Adjusted**"))
        dsn2$lab = group_i$lab[j]
        dsn2$item = item
        colnames(dsn2)[4] = "grade"

        # -- Append
        plot_combined0 = rbind.data.frame(plot_combined0, dsn2, stringsAsFactors=FALSE)
      }

      plot_combined0$grade_item = ifelse(plot_combined0$lab!="Composite", plot_combined0$grade, NA)
      plot_combined0$grade_comp = ifelse(plot_combined0$lab=="Composite", plot_combined0$grade, NA)
      plot_combined0$item_lab = factor(plot_combined0$lab, levels=unique(plot_combined0$lab))

      # ----------------------------------------------------------------
      # -- Color palette options
      # ----------------------------------------------------------------

      if(colors==1){
        # -- Default
        item_col0 = "white"
        item_col1 = "#90B0D9"
        item_col2 = "#4D7EBF"
        item_col3 = "#13478C"
        item_col4 = "#142233"

        comp_col0 = "white"
        comp_col1 = "#E5BFC6"
        comp_col2 = "#D9576E"
        comp_col3 = "#99293D"

        colset1 = c(
          "0" = item_col0,
          "1" = item_col1,
          "2" = item_col2,
          "3" = item_col3,
          "4" = item_col4
        )
        colset2 = c(
          "0" = comp_col0,
          "1" = comp_col1,
          "2" = comp_col2,
          "3" = comp_col3
        )
      } else if(colors==2){
        # -- Qualitative (color-blind friendly)
        item_col0 = "white"
        item_col1 = "#0571B0"
        item_col2 = "#92C5DE"
        item_col3 = "#F4A582"
        item_col4 = "#CA0020"

        comp_col0 = "white"
        comp_col1 = "#A6DBA0"
        comp_col2 = "#C2A5CF"
        comp_col3 = "#7B3294"

        colset1 = c(
          "0" = item_col0,
          "1" = item_col1,
          "2" = item_col2,
          "3" = item_col3,
          "4" = item_col4
        )
        colset2 = c(
          "0" = comp_col0,
          "1" = comp_col1,
          "2" = comp_col2,
          "3" = comp_col3
        )
      } else if(colors==3){
        # -- Black and white
        item_col0 = "white"
        item_col1 = "#E6E6E6"
        item_col2 = "#BDBDBD"
        item_col3 = "#636363"
        item_col4 = "#000000"

        comp_col0 = "white"
        comp_col1 = "#BDBDBD"
        comp_col2 = "#636363"
        comp_col3 = "#000000"

        colset1 = c(
          "0" = item_col0,
          "1" = item_col1,
          "2" = item_col2,
          "3" = item_col3,
          "4" = item_col4
        )
        colset2 = c(
          "0" = comp_col0,
          "1" = comp_col1,
          "2" = comp_col2,
          "3" = comp_col3
        )
      }

      ## -- Line color adjustments
      colset1_line = colset1
      colset1_line[1] = "grey"

      colset2_line = colset2
      colset2_line[1] = "grey"

      ## -- Optional title
      title_i = refset[refset$group_rank==i,]$group_lab[1]

      plot_combined0$grade_item = ifelse(plot_combined0$lab!="Composite", plot_combined0$grade, NA)
      plot_combined0$grade_comp = ifelse(plot_combined0$lab=="Composite", plot_combined0$grade, NA)
      plot_combined0$item_lab = factor(plot_combined0$lab, levels=unique(plot_combined0$lab))

      # -- Legend components
      facets = factor(unique(plot_combined0$lab), levels = c("Frequency", "Severity", "Interference", "Composite"))

      foot_note = ""
      ptsize1="15pt"
      ptsize2="10pt"

      for (k in facets){
        if(k=="Frequency"){
          # freq0 = paste0("<span style='font-size:",ptsize1,"'>\U23CD </span>Never /   ")
          freq0 = paste0("<span style='font-size:",ptsize1,"'>\U25A1</span> Never /   ")
          freq1 = paste0("<span style='font-size:",ptsize2,";color:", item_col1, ";'>\U25A0</span> Rarely   /   ")
          freq2 = paste0("<span style='font-size:",ptsize2,";color:", item_col2, ";'>\U25A0</span> Occasionally   /   ")
          freq3 = paste0("<span style='font-size:",ptsize2,";color:", item_col3, ";'>\U25A0</span> Frequently   /   ")
          freq4 = paste0("<span style='font-size:",ptsize2,";color:", item_col4, ";'>\U25A0</span> Almost constantly")
          foot_note = paste0(foot_note, "**Frequency:** ", freq0, freq1, freq2, freq3, freq4, "<br>")
        } else if(k=="Severity"){
          # sev0 = paste0("<span style='font-size:",ptsize1,";'>\U23CD </span>None  /  ")
          sev0 = paste0("<span style='font-size:",ptsize1,";'>\U25A1</span> None  /  ")
          sev1 = paste0("<span style='font-size:",ptsize2,";color:", item_col1, "'>\U25A0</span>  Mild   /   ")
          sev2 = paste0("<span style='font-size:",ptsize2,";color:", item_col2, ";'>\U25A0</span> Moderate   /   ")
          sev3 = paste0("<span style='font-size:",ptsize2,";color:", item_col3, ";'>\U25A0</span> Severe   /   ")
          sev4 = paste0("<span style='font-size:",ptsize2,";color:", item_col4, ";'>\U25A0</span> Very severe")
          foot_note = paste0(foot_note, "**Severity:** ", sev0, sev1, sev2, sev3, sev4, "<br>")
        } else if(k=="Interference"){
          # int0 = paste0("<span style='font-size:",ptsize1,";'>\U23CD </span>Not at all  /  ")
          int0 = paste0("<span style='font-size:",ptsize1,";'>\U25A1</span> Not at all  /  ")
          int1 = paste0("<span style='font-size:",ptsize2,";color:", item_col1, ";'>\U25A0</span> A little bit   /   ")
          int2 = paste0("<span style='font-size:",ptsize2,";color:", item_col2, ";'>\U25A0</span> Somewhat   /   ")
          int3 = paste0("<span style='font-size:",ptsize2,";color:", item_col3, ";'>\U25A0</span> Quite a bit   /   ")
          int4 = paste0("<span style='font-size:",ptsize2,";color:", item_col4, ";'>\U25A0</span> Very much")

          # ================================================================
          # ==== ACCOMMODATION FOR 27A HAIR LOSs AMOUNT =====================
          # ================================================================
          if(refset[refset$group_rank==i,]$group_lab[1]=="Hair Loss"){
            foot_note = paste0(foot_note, "**Amount:** ", int0, int1, int2, int3, int4, "<br>")
          } else {
            foot_note = paste0(foot_note, "**Interference:** ", int0, int1, int2, int3, int4, "<br>")
          }
          # ================================================================
          # ================================================================
          # ================================================================

        } else if(k=="Composite"){
          # comp0 = paste0("<span style='font-size:",ptsize1,";'>\U23CD </span>0  /  ")
          comp0 = paste0("<span style='font-size:",ptsize1,";'>\U25A1</span> 0  /  ")
          comp1 = paste0("<span style='font-size:",ptsize2,";color:", comp_col1 ,";'>\U25A0</span> 1  /  ")
          comp2 = paste0("<span style='font-size:",ptsize2,";color:", comp_col2 ,";'>\U25A0</span> 2  /  ")
          comp3 = paste0("<span style='font-size:",ptsize2,";color:", comp_col3 ,";'>\U25A0</span> 3 ")
          foot_note = paste0(foot_note, "**Composite Score:** ", comp0, comp1, comp2, comp3, "<br>")
        }
      }

      foot1 = "\\*Maximum score reported post-baseline per patient.<br>"
      foot2a = "\\*\\*Maximum score reported post-baseline per patient when"
      if(summary_only == TRUE || footnote_break == TRUE){
        foot2b = " including only scores<br> that were worse than the patient's baseline score."
      } else {
        foot2b = " including only scores that were worse than the patient's baseline score."
      }

      # -- Label footnote addition
      label_foot = ""
      if(bar_label %in% c(2,4)){
        foot_grade = 1
      }
      if(bar_label %in% c(3,5)){
        foot_grade = 3
      }
      if(bar_label %in% c(2,3)){
        foot_symbol = "(n)"
        foot_type = "number"
        y_scale_lab = c("0","25","50","75","100", "n ")
      }
      if(bar_label %in% c(4,5)){
        foot_symbol = "(%)"
        foot_type = "percent"
        y_scale_lab = c("0","25","50","75","100", "% ")
      }

      if(bar_label==1){
        label_foot = "Column labels (n) show the number of subjects with an observed symptom score or grade.<br>"
        y_scale_lab = c("0","25","50","75","100", "n ")
      } else if (bar_label %in% c(2, 3, 4, 5)){
        label_foot = paste0("Column labels ", foot_symbol, " show the ", foot_type,
                            " of subjects with symptom score ",
                            foot_grade," or greater.<br>")
      }

      if(cycles_only==TRUE){
        foot_note = paste0(foot_note, label_foot)
      } else if(cycles_only==FALSE){
        foot_note = paste0(foot_note, label_foot, foot1, foot2a, foot2b)
      }

      if(suppress_legend == TRUE){
        foot_note = ""
      }

      ## ----------------------------------------------------------------
      ## -- Conditional ggplot object combinations
      ## ----------------------------------------------------------------

      # -- Restrict the cycles to be plotted
      if(!is.na(plot_limit)){
        plot_combined0$cycle_plot_lim = as.integer(gsub("[^0-9]", "", plot_combined0[,cycle_var]))
        plot_combined0=plot_combined0[(!is.na(plot_combined0$cycle_plot_lim)&plot_combined0$cycle_plot_lim<=plot_limit)|plot_combined0[,cycle_var] %in% c("Adjusted**","Maximum*"),]
      }

      # -- Add labels to cycles to be plotted
      if(cycle_label == TRUE){
        names(cycle_labs) = paste0(cycle_var," ", cycle_vals)
        plot_combined0$cycle_var_v_plot = dplyr::recode(plot_combined0[, "cycle_var_v_plot"], !!!cycle_labs)
      }

      # -- Only show summary measures
      if(summary_only==TRUE){
        plot_combined0 = plot_combined0[plot_combined0[,cycle_var] %in% c("Adjusted**", "Maximum*"),]
      }

      # -- Do not show summary measures
      if(cycles_only==TRUE){
        plot_combined0 = plot_combined0[!(plot_combined0[,cycle_var] %in% c("Adjusted**", "Maximum*")),]
      }

      names(plot_combined0)[names(plot_combined0) == arm_var] = "arm"
      if(arm_var=="overall_"){plot_combined0$arm=""}

      if(summary_highlight == TRUE){
        plot_combined0$fill_highlight = ifelse(plot_combined0$cycle_var_v_plot %in% c("Maximum*","Adjusted**"), "black",  "#ededed")
      } else {plot_combined0$fill_highlight = NA}

      # ================================================================
      # ==== ACCOMMODATION FOR 27A HAIR LOSs AMOUNT ====================
      # ================================================================
      if(refset[refset$group_rank==i,]$group_lab[1]=="Hair Loss"){
        levels(plot_combined0$item_lab) = c("Amount", "Composite", "Interference")
        plot_combined0[plot_combined0$lab=="Interference",]$item_lab = "Amount"
      }

      # ================================================================
      # ================================================================
      # ================================================================

      # ----------------------------------------------------------------
      # -- No labels requested
      # ----------------------------------------------------------------
      if(bar_label==0){

        # ----------------------------------------------------------------
        # --
        # ----------------------------------------------------------------
        if("Composite" %in% unique(plot_combined0$lab) &
           sum(c("Frequency", "Severity", "Interference", "Amount") %in% unique(plot_combined0$lab))>0){

          figure_i = ggplot2::ggplot(plot_combined0, ggplot2::aes(arm)) +
          ggplot2::geom_rect(colour=plot_combined0$fill_highlight, fill="white", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +

            ## -- Individual items
            ggplot2::scale_fill_manual(values = colset1,
                                       limits = names(colset1),
                                       guide = ggplot2::guide_legend(order = 2)) +
            ggplot2::scale_color_manual(values = colset1_line) +
            ggplot2::geom_bar(data = plot_combined0[!is.na(plot_combined0$grade_item),],
                              ggplot2::aes(fill = as.factor(grade_item), color=as.factor(grade_item)), position = "fill", width = .8) +

            ggnewscale::new_scale(new_aes = "fill") +
            ggnewscale::new_scale(new_aes = "color") +

            # ## -- Composite items
            ggplot2::scale_fill_manual(values = colset2,
                                       limits = names(colset2),
                                       guide = ggplot2::guide_legend(order = 3)) +
            ggplot2::scale_color_manual(values = colset2_line) +
            ggplot2::geom_bar(data = plot_combined0[!is.na(plot_combined0$grade_comp),],
                              ggplot2::aes(fill = as.factor(grade_comp), color=as.factor(grade_comp)), position = "fill", width = .8) +

            ggplot2:: scale_y_continuous(labels = c("0","25","50","75","100"),
                                         breaks = c(0, .25, .5, .75, 1),
                                         limits = c(0, 1.01)) +

            ggplot2::facet_grid(item_lab~cycle_var_v_plot) +
            ggplot2::xlab(x_label) +
            ggplot2::ylab(y_label) +

            ggplot2::labs(caption = foot_note) +

            ggplot2::theme(
              plot.caption = ggtext::element_markdown(hjust=0, lineheight = -20),
              plot.title.position = "plot",
              plot.caption.position =  "plot",
              axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
              axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
              axis.text.x = ggplot2::element_text(angle = x_lab_angle, vjust=x_lab_vjust, hjust=x_lab_hjust),
              legend.position="none",
              legend.box = "vertical",
              legend.direction = "horizontal",
              legend.box.just = "left",
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              panel.background = ggplot2::element_rect(color="grey", fill="white"),
              strip.background = ggplot2::element_rect(colour="grey", fill="#ededed"),
              strip.text.y = ggplot2::element_text(angle=90))
        }

        # ----------------------------------------------------------------
        # -- Individual items only
        # ----------------------------------------------------------------
        else if(!("Composite" %in% unique(plot_combined0$lab))){

          figure_i = ggplot2::ggplot(plot_combined0, ggplot2::aes(arm)) +

            ggplot2::geom_rect(colour=plot_combined0$fill_highlight, fill="white", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +

            ## -- Individual items
            ggplot2::scale_fill_manual(values = colset1,
                                       limits = names(colset1),
                                       guide = ggplot2::guide_legend(order = 2)) +
            ggplot2::scale_color_manual(values = colset1_line) +
            ggplot2::geom_bar(ggplot2::aes(fill = as.factor(grade_item), color=as.factor(grade_item)), position = "fill", width = .8) +

            ggplot2::scale_y_continuous(labels = c("0","25","50","75","100"),
                               breaks = c(0, .25, .5, .75, 1),
                               limits = c(0, 1.01)) +

            ggplot2::facet_grid(item_lab~cycle_var_v_plot) +
            ggplot2::xlab(x_label) +
            ggplot2::ylab(y_label) +

            ggplot2::labs(caption = foot_note) +

            ggplot2::theme(
              plot.caption = ggtext::element_markdown(hjust=0, lineheight = -20),
              plot.title.position = "plot",
              plot.caption.position =  "plot",
              axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
              axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
              axis.text.x = ggplot2::element_text(angle = x_lab_angle, vjust=x_lab_vjust, hjust=x_lab_hjust),
              legend.position="none",
              legend.box = "vertical",
              legend.direction = "horizontal",
              legend.box.just = "left",
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              panel.background = ggplot2::element_rect(color="grey", fill="white"),
              strip.background = ggplot2::element_rect(colour="grey", fill="#ededed"),
              strip.text.y = ggplot2::element_text(angle=90))

        }

        # ----------------------------------------------------------------
        # -- Composite items only
        # ----------------------------------------------------------------
        else if("Composite" %in% unique(plot_combined0$lab)&
                sum(c("Frequency", "Severity", "Interference", "Amount") %in% unique(plot_combined0$lab))==0){

          figure_i = ggplot2::ggplot(plot_combined0, ggplot2::aes(x=arm)) +
          ggplot2::geom_rect(colour=plot_combined0$fill_highlight, fill="white", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +

            ## -- Composite items
            ggplot2::scale_fill_manual(values = colset2,
                                       limits = names(colset2),
                                       guide = ggplot2::guide_legend(order = 3)) +
            ggplot2::scale_color_manual(values = colset2_line) +
            ggplot2::geom_bar(ggplot2::aes(fill = as.factor(grade_comp), color=as.factor(grade_comp)), position = "fill", width = .8) +

            ggplot2::scale_y_continuous(labels = c("0","25","50","75","100"),
                               breaks = c(0, .25, .5, .75, 1),
                               limits = c(0, 1.01)) +

            ggplot2::facet_grid(item_lab~cycle_var_v_plot) +
            ggplot2::xlab(x_label) +
            ggplot2::ylab(y_label) +

            ggplot2::labs(caption = foot_note) +

            ggplot2::theme(
              plot.caption = ggtext::element_markdown(hjust=0, lineheight = -20),
              plot.title.position = "plot",
              plot.caption.position =  "plot",
              axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
              axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
              axis.text.x = ggplot2::element_text(angle = x_lab_angle, vjust=x_lab_vjust, hjust=x_lab_hjust),
              legend.position="none",
              legend.box = "vertical",
              legend.direction = "horizontal",
              legend.box.just = "left",
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              panel.background = ggplot2::element_rect(color="grey", fill="white"),
              strip.background = ggplot2::element_rect(colour="grey", fill="#ededed"),
              strip.text.y = ggplot2::element_text(angle=90))

        }
      }

      # ----------------------------------------------------------------
      # -- Labels requested
      # ----------------------------------------------------------------
      else if(bar_label!=0){

        # ----------------------------------------------------------------
        # -- Individual and composite items
        # ----------------------------------------------------------------
        if("Composite" %in% unique(plot_combined0$lab) &
           sum(c("Frequency", "Severity", "Interference", "Amount") %in% unique(plot_combined0$lab))>0){

          figure_i = ggplot2::ggplot(plot_combined0, ggplot2::aes(x=arm)) +

          ggplot2::geom_rect(colour=plot_combined0$fill_highlight, fill="white", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +

            ## -- Individual items
            ggplot2::scale_fill_manual(values = colset1,
                                       limits = names(colset1),
                                       guide = ggplot2::guide_legend(order = 2)) +
            ggplot2::scale_color_manual(values = colset1_line) +
            ggplot2::geom_bar(data = plot_combined0[!is.na(plot_combined0$grade_item),],
                              ggplot2::aes(fill = as.factor(grade_item), color=as.factor(grade_item)), position = "fill", width = .8) +
            ggnewscale::new_scale(new_aes = "fill") +
            ggnewscale::new_scale(new_aes = "color") +


            ## -- Composite items
            ggplot2::scale_fill_manual(values = colset2,
                                       limits = names(colset2),
                                       guide = ggplot2::guide_legend(order = 3)) +
            ggplot2::scale_color_manual(values = colset2_line) +
            ggplot2::geom_bar(data = plot_combined0[!is.na(plot_combined0$grade_comp),],
                              ggplot2::aes(fill = as.factor(grade_comp), color=as.factor(grade_comp)), position = "fill", width = .8) +

            ## -- Customize
            ggplot2::geom_hline(yintercept=1.04, linetype="solid", color = "darkgrey") +


            ggplot2::geom_text(ggplot2::aes(y= 1.13, label=bar_lab_opt), size=3, na.rm=TRUE) +
            ggplot2::scale_y_continuous(labels = y_scale_lab,
                                        breaks = c(0, .25, .5, .75, 1, 1.13),
                                        limits = c(0, 1.14)) +

            ggplot2::facet_grid(item_lab~cycle_var_v_plot) +
            ggplot2::xlab(x_label) +
            ggplot2::ylab(y_label) +

            ggplot2::labs(caption = foot_note) +


            ggplot2::theme(
              plot.caption = ggtext::element_markdown(hjust=0, lineheight = -20),
              plot.title.position = "plot",
              plot.caption.position =  "plot",
              axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
              axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
              axis.text.x = ggplot2::element_text(angle = x_lab_angle, vjust=x_lab_vjust, hjust=x_lab_hjust),
              legend.position="none",
              legend.box = "vertical",
              legend.direction = "horizontal",
              legend.box.just = "left",
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              # panel.background = ggplot2::element_rect(color="grey", fill="white"),
              # strip.background = ggplot2::element_rect(colour="grey", fill="#ededed"),
              panel.background = ggplot2::element_rect(color="grey", fill="white"),
              strip.background = ggplot2::element_rect(colour="grey", fill="#ededed"),
              strip.text.y = ggplot2::element_text(angle=90))
        }

        # ----------------------------------------------------------------
        # -- Individual items only
        # ----------------------------------------------------------------
        else if(!("Composite" %in% unique(plot_combined0$lab))){

          figure_i = ggplot2::ggplot(plot_combined0, ggplot2::aes(x=arm)) +
          ggplot2::geom_rect(colour=plot_combined0$fill_highlight, fill="white", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +

            ## -- Individual items
            ggplot2::scale_fill_manual(values = colset1,
                                       limits = names(colset1),
                                       guide = ggplot2::guide_legend(order = 2)) +
            ggplot2::scale_color_manual(values = colset1_line) +
            ggplot2::geom_bar(ggplot2::aes(fill = as.factor(grade_item), color=as.factor(grade_item)), position = "fill", width = .8) +

            ## -- Customize
            ggplot2::geom_hline(yintercept=1.04, linetype="solid", color = "darkgrey") +

            ggplot2::geom_text(ggplot2::aes(y= 1.13, label=bar_lab_opt), size=3, na.rm=TRUE) +
            ggplot2::scale_y_continuous(labels = y_scale_lab,
                                        breaks = c(0, .25, .5, .75, 1, 1.13),
                                        limits = c(0, 1.14)) +

            ggplot2::facet_grid(item_lab~cycle_var_v_plot) +
            ggplot2::xlab(x_label) +
            ggplot2::ylab(y_label) +

            ggplot2::labs(caption = foot_note) +

            ggplot2::theme(
              plot.caption = ggtext::element_markdown(hjust=0, lineheight = -20),
              plot.title.position = "plot",
              plot.caption.position =  "plot",
              axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
              axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
              axis.text.x = ggplot2::element_text(angle = x_lab_angle, vjust=x_lab_vjust, hjust=x_lab_hjust),
              legend.position="none",
              legend.box = "vertical",
              legend.direction = "horizontal",
              legend.box.just = "left",
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              panel.background = ggplot2::element_rect(color="grey", fill="white"),
              strip.background = ggplot2::element_rect(colour="grey", fill="#ededed"),
              strip.text.y = ggplot2::element_text(angle=90))
        }

        # ----------------------------------------------------------------
        # -- Composite items only
        # ----------------------------------------------------------------
        else if("Composite" %in% unique(plot_combined0$lab)&
                sum(c("Frequency", "Severity", "Interference", "Amount") %in% unique(plot_combined0$lab))==0){

          figure_i = ggplot2::ggplot(plot_combined0, ggplot2::aes(x=arm)) +
          ggplot2::geom_rect(colour=plot_combined0$fill_highlight, fill="white", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +

            ## -- Composite items
            ggplot2::scale_fill_manual(values = colset2,
                                       limits = names(colset2),
                                       guide = ggplot2::guide_legend(order = 3)) +
            ggplot2::scale_color_manual(values = colset2_line) +
            ggplot2::geom_bar(ggplot2::aes(fill = as.factor(grade_comp), color=as.factor(grade_comp)), position = "fill", width = .8) +

            ## -- Customize
            ggplot2::geom_hline(yintercept=1.04, linetype="solid", color = "darkgrey") +

            ggplot2::geom_text(ggplot2::aes(y= 1.13, label=bar_lab_opt), size=3, na.rm=TRUE) +
            ggplot2::scale_y_continuous(labels = y_scale_lab,
                                        breaks = c(0, .25, .5, .75, 1, 1.13),
                                        limits = c(0, 1.14)) +

            ggplot2::facet_grid(item_lab~cycle_var_v_plot) +
            ggplot2::xlab(x_label) +
            ggplot2::ylab(y_label) +

            ggplot2::labs(caption = foot_note) +

            ggplot2::theme(
              plot.caption = ggtext::element_markdown(hjust=0, lineheight = -20),
              plot.title.position = "plot",
              plot.caption.position =  "plot",
              axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
              axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
              axis.text.x = ggplot2::element_text(angle = x_lab_angle, vjust=x_lab_vjust, hjust=x_lab_hjust),
              legend.position="none",
              legend.box = "vertical",
              legend.direction = "horizontal",
              legend.box.just = "left",
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              panel.background = ggplot2::element_rect(color="grey", fill="white"),
              strip.background = ggplot2::element_rect(colour="grey", fill="#ededed"),
              strip.text.y = ggplot2::element_text(angle=90))
        }
      }

      list_out[[i]] = list()
      list_out[[i]][[1]] = refset[refset$group_rank==i,]$group_lab[1]

      if (add_item_title == TRUE){
        list_out[[i]][[2]] = figure_i +
          ggplot2::ggtitle(refset[refset$group_rank==i,]$group_lab[1])
      } else {list_out[[i]][[2]] = figure_i}

    }
  }
  # ----------------------------------------------------------------
  # --  Other PROs
  # ----------------------------------------------------------------

  # ----------------------------------------------------------------
  # --  GP5
  # ----------------------------------------------------------------

  if(any(refset$name=="GP5")){
    i = refset[refset$name=="GP5",]$group_rank
    group_i = refset[refset$group_rank==i,]
    plot_combined0 = data.frame()
    for (j in 1:NROW(group_i)){
      item = as.character(group_i$name[j])

      if(arm_var != "overall_"){
        dsn0 = stats::na.omit(dsn[,c(id_var, cycle_var, arm_var, item)])
      } else if (arm_var == "overall_"){
        dsn0 = stats::na.omit(dsn[,c(id_var, cycle_var, item)])
        dsn0$overall_ = "Overall"
      }

      # ----------------------------------------------------------------
      # -- Construct summary measures
      # ----------------------------------------------------------------
      ## -- Maximum grade post baseline
      max_post0 = stats::aggregate(dsn0[dsn0[,cycle_var]>baseline_val, item],
                                   by = list(dsn0[dsn0[,cycle_var]>baseline_val,id_var]),
                                   FUN=max)
      colnames(max_post0) = c(id_var, item)
      max_post0[,cycle_var] = "Maximum*"
      max_post1 = unique(merge(x=max_post0, y=dsn0[,c(id_var, arm_var)], by=id_var, all.x = TRUE))

      ## -- Baseline adjusted
      base_adj0 = dsn0[dsn0[,cycle_var]==baseline_val, c(id_var, item)]
      colnames(base_adj0)[2] = "base_val"

      base_adj1 = merge(x=dsn0, y=base_adj0, by=id_var, all.x=TRUE)
      base_adj2 = stats::aggregate(base_adj1[base_adj1[,cycle_var]>baseline_val, item],
                                   by = list(base_adj1[base_adj1[,cycle_var]>baseline_val,id_var]),
                                   FUN=max)
      colnames(base_adj2) = c(id_var, "max_post_bl")
      base_adj2 = merge(x=base_adj1, y=unique(base_adj2), by=id_var, all.x=TRUE)
      base_adj2$bl_adjusted = ifelse(!is.na(base_adj2$max_post_bl) & base_adj2$base_val >= base_adj2$max_post_bl,
                                     0, base_adj2$max_post_bl)
      base_adj3 = unique(base_adj2[!is.na(base_adj2$bl_adjusted),c(id_var, arm_var, "bl_adjusted")])
      base_adj3[,cycle_var] = "Adjusted**"
      colnames(base_adj3)[colnames(base_adj3) == "bl_adjusted"] = item


      dsn1 = rbind.data.frame(dsn0, max_post1, base_adj3)
      # ----------------------------------------------------------------
      # -- Create labels for the bars
      # ----------------------------------------------------------------

      if(bar_label==0){
        labs = dsn1[,c(cycle_var, arm_var)]
        labs = unique(labs[c(cycle_var, arm_var)])
        labs$bar_lab_opt = NA

      } else if(bar_label==1){
        labs = stats::aggregate(dsn1[dsn1[,item]>=0,item],
                                by=list(dsn1[dsn1[,item]>=0,cycle_var],
                                        dsn1[dsn1[,item]>=0,arm_var]),
                                FUN=length)
      } else if(bar_label %in% c(2,4)){
        labs = stats::aggregate(dsn1[dsn1[,item]>0,item],
                                by=list(dsn1[dsn1[,item]>0,cycle_var],
                                        dsn1[dsn1[,item]>0,arm_var]),
                                FUN=length)
      } else if(bar_label %in% c(3,5)){
        labs = stats::aggregate(dsn1[dsn1[,item]>=3,item],
                                by=list(dsn1[dsn1[,item]>=3,cycle_var],
                                        dsn1[dsn1[,item]>=3,arm_var]),
                                FUN=length)
      }
      names(labs) = c(cycle_var, arm_var, "bar_lab_opt")

      # -- Force zero counts
      combo_counts = data.frame(table(dsn1[,cycle_var], dsn1[,arm_var]), stringsAsFactors=FALSE)
      names(combo_counts) = c(cycle_var, arm_var, "total")
      combo_counts$flag = 1

      labs2 = merge(x=labs, y=combo_counts, by=c(cycle_var, arm_var), all.y = TRUE)

      if(any(is.na(labs2$bar_lab_opt))){
        labs2[is.na(labs2$bar_lab_opt),]$bar_lab_opt = 0
      }

      if(bar_label %in% c(4,5)){
        labs2[labs2$bar_lab_opt!=0,]$bar_lab_opt=round(labs2[labs2$bar_lab_opt!=0,]$bar_lab_opt/labs2[labs2$bar_lab_opt!=0,]$total,digits=2)*100
      }

      dsn2 = merge(dsn1, labs2, by = c(cycle_var, arm_var))

      # -- Only one label per cycle/arm combination
      high = by(dsn2, INDICES = list(dsn2[,cycle_var], dsn2[,arm_var]), FUN = utils::tail, n=1)
      high2 = do.call("rbind", as.list(high))
      row_keeps = as.integer(rownames(high2))
      dsn2[!(1:NROW(dsn2) %in% row_keeps),]$bar_lab_opt = NA

      # ----------------------------------------------------------------
      # -- Create labels for the cycle_var facets
      # ----------------------------------------------------------------
      dsn2$cycle_var_v_plot=factor(ifelse(dsn2[,cycle_var]=="Maximum*",
                                          "Maximum*",
                                          ifelse(dsn2[,cycle_var]=="Adjusted**", "Adjusted**",
                                                 paste0(cycle_var," ", dsn2[,cycle_var]))),
                                   levels=c(paste0(cycle_var," ", min(dsn0[,cycle_var]):max(dsn0[,cycle_var])),
                                            "Maximum*",
                                            "Adjusted**"))
      dsn2$lab = group_i$lab[j]
      dsn2$item = item
      colnames(dsn2)[4] = "grade"

      # -- Append
      plot_combined0 = rbind.data.frame(plot_combined0, dsn2, stringsAsFactors=FALSE)
    }


    plot_combined0$grade_item = plot_combined0$grade
    plot_combined0$item_lab = factor(plot_combined0$lab, levels=unique(plot_combined0$lab))

    # ----------------------------------------------------------------
    # -- Color palette options
    # ----------------------------------------------------------------

    if(colors==1){
      # -- Default
      item_col0 = "white"
      item_col1 = "#90B0D9"
      item_col2 = "#4D7EBF"
      item_col3 = "#13478C"
      item_col4 = "#142233"

    } else if(colors==2){
      # -- Qualitative (color-blind friendly)
      item_col0 = "white"
      item_col1 = "#0571B0"
      item_col2 = "#92C5DE"
      item_col3 = "#F4A582"
      item_col4 = "#CA0020"

    } else if(colors==3){
      # -- Black and white
      item_col0 = "white"
      item_col1 = "#E6E6E6"
      item_col2 = "#BDBDBD"
      item_col3 = "#636363"
      item_col4 = "#000000"
    }

    colset1 = c(
      "0" = item_col0,
      "1" = item_col1,
      "2" = item_col2,
      "3" = item_col3,
      "4" = item_col4
    )

    ## -- Line color adjustments
    colset1_line = colset1
    colset1_line[1] = "grey"


    ## -- Optional title
    title_i = refset[refset$group_rank==i,]$group_lab[1]

    # -- Legend components
    facets = factor(unique(plot_combined0$lab))

    foot_note = ""
    ptsize1="15pt"
    ptsize2="10pt"

    freq0 = paste0("<span style='font-size:",ptsize1,"'>\U25A1</span> Not at all /   ")
    freq1 = paste0("<span style='font-size:",ptsize2,";color:", item_col1, ";'>\U25A0</span> A little bit   /   ")
    freq2 = paste0("<span style='font-size:",ptsize2,";color:", item_col2, ";'>\U25A0</span> Somewhat   /   ")
    freq3 = paste0("<span style='font-size:",ptsize2,";color:", item_col3, ";'>\U25A0</span> Quite a bit   /   ")
    freq4 = paste0("<span style='font-size:",ptsize2,";color:", item_col4, ";'>\U25A0</span> Very much")
    foot_note = paste0(foot_note, "**Symptom Bother:** ", freq0, freq1, freq2, freq3, freq4, "<br>")

    foot1 = "\\*Maximum bother reported post-baseline per patient.<br>"
    foot2a = "\\*\\*Maximum bother reported post-baseline per patient when"
    if(summary_only == TRUE || footnote_break == TRUE){
      foot2b = " including only scores<br> that were worse than the patient's baseline score."
    } else {
      foot2b = " including only scores that were worse than the patient's baseline score."
    }

    # -- Label footnote addition
    label_foot = ""
    if(bar_label %in% c(2,4)){
      foot_grade = 1
    }
    if(bar_label %in% c(3,5)){
      foot_grade = 3
    }
    if(bar_label %in% c(2,3)){
      foot_symbol = "(n)"
      foot_type = "number"
      y_scale_lab = c("0","25","50","75","100", "n ")
    }
    if(bar_label %in% c(4,5)){
      foot_symbol = "(%)"
      foot_type = "percent"
      y_scale_lab = c("0","25","50","75","100", "% ")
    }

    if(bar_label==1){
      label_foot = "Column labels (n) show the number of subjects with an observed bother score.<br>"
      y_scale_lab = c("0","25","50","75","100", "n ")
    } else if (bar_label %in% c(2, 3, 4, 5)){
      label_foot = paste0("Column labels ", foot_symbol, " show the ", foot_type,
                          " of subjects with bother score ",
                          foot_grade," or greater.<br>")
    }

    if(cycles_only==TRUE){
      foot_note = paste0(foot_note, label_foot)
    } else if(cycles_only==FALSE){
      foot_note = paste0(foot_note, label_foot, foot1, foot2a, foot2b)
    }

    if(suppress_legend == TRUE){
      foot_note = ""
    }

    ## ----------------------------------------------------------------
    ## -- Create ggplot object
    ## ----------------------------------------------------------------

    # -- Restrict the cycles to be plotted
    if(!is.na(plot_limit)){
      plot_combined0$cycle_plot_lim = as.integer(gsub("[^0-9]", "", plot_combined0[,cycle_var]))
      plot_combined0=plot_combined0[(!is.na(plot_combined0$cycle_plot_lim)&plot_combined0$cycle_plot_lim<=plot_limit)|plot_combined0[,cycle_var] %in% c("Adjusted**","Maximum*"),]
    }

    # -- Restrict the cycles to be plotted
    if(cycle_label == TRUE){
      names(cycle_labs) = paste0(cycle_var," ", cycle_vals)
      plot_combined0$cycle_var_v_plot = dplyr::recode(plot_combined0[, "cycle_var_v_plot"], !!!cycle_labs)
    }

    # -- Only show summary measures
    if(summary_only==TRUE){
      plot_combined0 = plot_combined0[plot_combined0[,cycle_var] %in% c("Adjusted**", "Maximum*"),]
    }

    # -- Do not show summary measures
    if(cycles_only==TRUE){
      plot_combined0 = plot_combined0[!(plot_combined0[,cycle_var] %in% c("Adjusted**", "Maximum*")),]
    }

    names(plot_combined0)[names(plot_combined0) == arm_var] = "arm"
    if(arm_var=="overall_"){plot_combined0$arm=""}

    if(summary_highlight == TRUE){
      plot_combined0$fill_highlight = ifelse(plot_combined0$cycle_var_v_plot %in% c("Maximum*","Adjusted**"), "black",  "#ededed")
    } else {plot_combined0$fill_highlight = NA}

    # ----------------------------------------------------------------
    # -- No labels requested
    # ----------------------------------------------------------------
    if(bar_label==0){

      figure_i = ggplot2::ggplot(plot_combined0, ggplot2::aes(arm)) +

        ggplot2::geom_rect(colour=plot_combined0$fill_highlight, fill="white", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +

        ## -- Individual items
        ggplot2::scale_fill_manual(values = colset1,
                                   limits = names(colset1),
                                   guide = ggplot2::guide_legend(order = 2)) +
        ggplot2::scale_color_manual(values = colset1_line) +
        ggplot2::geom_bar(ggplot2::aes(fill = as.factor(grade_item), color=as.factor(grade_item)), position = "fill", width = .8) +

        ggplot2::scale_y_continuous(labels = c("0","25","50","75","100"),
                                    breaks = c(0, .25, .5, .75, 1),
                                    limits = c(0, 1.01)) +

        ggplot2::facet_grid(.~cycle_var_v_plot) +
        ggplot2::xlab(x_label) +
        ggplot2::ylab(y_label) +

        ggplot2::labs(caption = foot_note) +

        ggplot2::theme(
          plot.caption = ggtext::element_markdown(hjust=0, lineheight = -20),
          plot.title.position = "plot",
          plot.caption.position =  "plot",
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
          axis.text.x = ggplot2::element_text(angle = x_lab_angle, vjust=x_lab_vjust, hjust=x_lab_hjust),
          legend.position="none",
          legend.box = "vertical",
          legend.direction = "horizontal",
          legend.box.just = "left",
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_rect(color="grey", fill="white"),
          strip.background = ggplot2::element_rect(colour="grey", fill="#ededed"),
          strip.text.y = ggplot2::element_text(angle=90))

    }

    # ----------------------------------------------------------------
    # -- Labels requested
    # ----------------------------------------------------------------
    else if(bar_label!=0){

      figure_i = ggplot2::ggplot(plot_combined0, ggplot2::aes(x=arm)) +
        ggplot2::geom_rect(colour=plot_combined0$fill_highlight, fill="white", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +

        ## -- Individual items
        ggplot2::scale_fill_manual(values = colset1,
                                   limits = names(colset1),
                                   guide = ggplot2::guide_legend(order = 2)) +
        ggplot2::scale_color_manual(values = colset1_line) +
        ggplot2::geom_bar(ggplot2::aes(fill = as.factor(grade_item), color=as.factor(grade_item)), position = "fill", width = .8) +

        ## -- Customize
        ggplot2::geom_hline(yintercept=1.04, linetype="solid", color = "darkgrey") +

        ggplot2::geom_text(ggplot2::aes(y= 1.13, label=bar_lab_opt), size=3, na.rm=TRUE) +
        ggplot2::scale_y_continuous(labels = y_scale_lab,
                                    breaks = c(0, .25, .5, .75, 1, 1.13),
                                    limits = c(0, 1.14)) +

        ggplot2::facet_grid(.~cycle_var_v_plot) +
        ggplot2::xlab(x_label) +
        ggplot2::ylab(y_label) +

        ggplot2::labs(caption = foot_note) +

        ggplot2::theme(
          plot.caption = ggtext::element_markdown(hjust=0, lineheight = -20),
          plot.title.position = "plot",
          plot.caption.position =  "plot",
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
          axis.text.x = ggplot2::element_text(angle = x_lab_angle, vjust=x_lab_vjust, hjust=x_lab_hjust),
          legend.position="none",
          legend.box = "vertical",
          legend.direction = "horizontal",
          legend.box.just = "left",
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_rect(color="grey", fill="white"),
          strip.background = ggplot2::element_rect(colour="grey", fill="#ededed"),
          strip.text.y = ggplot2::element_text(angle=90))
    }

    # --- add to output list object

    list_out[[i]] = list()
    list_out[[i]][[1]] = refset[refset$group_rank==i,]$group_lab[1]

    if (add_item_title == TRUE){
      list_out[[i]][[2]] = figure_i +
        ggplot2::ggtitle(refset[refset$group_rank==i,]$group_lab[1])
    } else {list_out[[i]][[2]] = figure_i}
  }


  ## -- Reference table for user to view the indexing of PRO-CTCAE item groups withing the list output
  for (i in 1:length(list_out)){
    if (i==1) {fig_tab = data.frame()}
    list_deat = data.frame(PRO_AE = list_out[[i]][[1]], item_index=i)
    fig_tab = rbind(fig_tab, list_deat)
  }

  ## -- Reference table
  print(fig_tab)

  ## -- Object return
  invisible(list_out)

}

