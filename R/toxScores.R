#' Re-code PRO-CTCAE text responses, apply a zero-imputation procedures, and
#' construct PRO-CTCAE composite grades.
#'
#'
#'    This function takes in a data frame set with PRO-CTCAE survey text
#'    fields/responses and returns a data frame with appropriate numerical
#'    re-coding. This function will accept 1 or up to all 124 PRO-CTCAE survey
#'    fields. All PRO-CTCAE variable names MUST conform to a pre-specified
#'    naming structure. PRO-CTCAE variable names are made up of FOUR components:
#'    1)'PROCTCAE', 2) number [1,2,3, ..., i, ..., 80], 3) 'A', 'B', or 'C'
#'    component of the i-th PRO-CTCAE field, 4) and 'SCL' (if severity,
#'    interference, or frequency) or 'IND' (if yes/no variable). Each component
#'    must be delimited by an underscore (_)
#'
#'    \itemize{
#'      \item[EX1] Question 1 of PRO-CTCAE should be: PROCTCAE_1A_SCL
#'      \item[EX2] Question 48 of PRO-CTCAE should be: PROCTCAE_48A_SCL,
#'      PROCTCAE_48B_SCL, PROCTCAE_48C_SCL
#'      \item[EX3]Question 73 of PRO-CTCAE should be: PROCTCAE_73A_IND
#'    }
#'
#'    This function also constructs PRO-CTCAE composite grades. Composite grade
#'    variables for respective PRO-CTCAE item groups are created and named as
#'    PROCTCAE_##_COMP.
#'
#'    \enumerate{
#'      \item https://healthcaredelivery.cancer.gov/pro-ctcae/pro-ctcae_english.pdf
#'      \item Ethan Basch, et al. Development of a Composite Scoring Algorithm
#'            for the National Cancer Institute's Patient-Reported Outcomes
#'            version of the Common Terminology Criteria for Adverse Events
#'            (PRO-CTCAE). ISOQOL 2019
#'      \item Basch E, et al. Composite Grading Algorithm for the National
#'            Cancer Institute’s Patient-Reported Outcomes version of the
#'            Common Terminology Criteria for Adverse Events (PRO-CTCAE).
#'            Clinical Trials 2020.
#'    }
#'
#' 	  Data format should be in 'long' format, where each PRO-CTCAE item is a
#' 	  variable/column.
#'
#' @param dsn A data.frame object with PRO-CTCAE data
#' @param reformat Reformat PRO-CTCAE text responses to numeric scores. Defaults
#'   to \code{FALSE}.
#' @param impute Apply zero-imputation where appropriate. Defaults to \code{FALSE}.
#' @param composites Construct composite grade using available PRO-CTCAE
#'   variables within \code{dsn}. Defaults to \code{FALSE}.
#' @param short_labels Add PRO-CTCAE short labels to available PRO-CTCAE
#'   variables within returned object
#' @return A data.frame object.
#' @examples
#' tox_acute_comp = toxScores(dsn = ProAE::tox_acute, composites = TRUE)
#' @export

toxScores = function(dsn,
                     reformat=FALSE,
                     impute=FALSE,
                     composites=FALSE,
                     short_labels=FALSE){

  # ----------------------------------------------------------------
  # -- Checks 1/1
  # ----------------------------------------------------------------

  ## -- Assign binding for data.frame variables used within subset function as global variables

  fmt_name = NULL

  ## -- Required parameters

  if(exists("dsn")){
    if(!is.data.frame(dsn)){
      stop("param dsn must be provided as a data.frame object")
    }
  } else {stop("param dsn not provided")}

  ## -------------------------------------------------------------
  ## --- Reference data sets
  ## -------------------------------------------------------------

  # - Do not evaluate yes/no PROCTCAE items in composite grading/imputation
  dsn_items = toupper(names(dsn)[toupper(names(dsn)) %in% proctcae_vars$name])
  dsn_items_comp = dsn_items[! dsn_items %in% as.character(proctcae_vars$name[proctcae_vars$fmt %in% c("yn_2_fmt", "yn_3_fmt", "yn_4_fmt")])]

  dsn_othrs = names(dsn)[!toupper(names(dsn)) %in% proctcae_vars$name]

  ## --------------------------------------------------------------------------
  ## --- Apply formats across those items found in dsn
  ## --------------------------------------------------------------------------

  if(reformat == TRUE){

    # -------------------------------------------------------------
    # --- Formatting
    # -------------------------------------------------------------

    fmts = c("int_5_fmt","sev_5_fmt", "sev_6_fmt", "sev_7_fmt", "frq_5_fmt", "frq_7_fmt", "yn_2_fmt", "yn_3_fmt", "yn_4_fmt")

    fmt_func = function(dat, fmt){
      if (fmt == "int_5_fmt"){
        as.numeric(dplyr::recode(toupper(dat),
                                 "NOT AT ALL" = "0",
                                 "A LITTLE BIT" = "1",
                                 "SOMEWHAT" = "2",
                                 "QUITE A BIT" = "3",
                                 "VERY MUCH" = "4"))
      }
      else if (fmt == "sev_5_fmt"){
        as.numeric(dplyr::recode(toupper(dat),
                                 "NONE" = "0",
                                 "MILD" = "1",
                                 "MODERATE" = "2",
                                 "SEVERE" = "3",
                                 "VERY SEVERE" = "4"))
      }
      else if (fmt == "sev_6_fmt"){
        as.numeric(dplyr::recode(toupper(dat),
                                 "NONE" = "0",
                                 "MILD" = "1",
                                 "MODERATE" = "2",
                                 "SEVERE" = "3",
                                 "VERY SEVERE" = "4",
                                 "NOT APPLICABLE" = NA))
      }
      else if (fmt == "sev_7_fmt"){
        as.numeric(dplyr::recode(toupper(dat),
                                 "NONE" = "0",
                                 "MILD" = "1",
                                 "MODERATE" = "2",
                                 "SEVERE" = "3",
                                 "VERY SEVERE" = "4",
                                 "NOT SEXUALLY ACTIVE" = NA,
                                 "PREFER NOT TO ANSWER" = NA))
      }
      else if (fmt == "frq_5_fmt"){
        as.numeric(dplyr::recode(toupper(dat),
                                 "NEVER" = "0",
                                 "RARELY" = "1",
                                 "OCCASIONALLY" = "2",
                                 "FREQUENTLY" = "3",
                                 "ALMOST CONSTANTLY" = "4"))
      }
      else if (fmt == "frq_7_fmt"){
        as.numeric(dplyr::recode(toupper(dat),
                                 "NEVER" = "0",
                                 "RARELY" = "1",
                                 "OCCASIONALLY" = "2",
                                 "FREQUENTLY" = "3",
                                 "ALMOST CONSTANTLY" = "4",
                                 "NOT SEXUALLY ACTIVE" = NA,
                                 "PREFER NOT TO ANSWER" = NA))
      }
      else if (fmt == "yn_2_fmt"){
        as.numeric(dplyr::recode(toupper(dat),
                                 "NO" = "0",
                                 "YES" = "1"))
      }
      else if (fmt == "yn_3_fmt"){
        as.numeric(dplyr::recode(toupper(dat),
                                 "NO" = "0",
                                 "YES" = "1",
                                 "NOT APPLICABLE" = NA))
      }
      else if (fmt == "yn_4_fmt"){
        as.numeric(dplyr::recode(toupper(dat),
                                 "NO" = "0",
                                 "YES" = "1",
                                 "NOT APPLICABLE" = NA,
                                 "PREFER NOT TO ANSWER" = NA))
      }
    }

    # Convert all PRO-CTCAE items to character to avoid factor issues
    dsn1 = dsn
    dsn1[dsn_items] = lapply(dsn[dsn_items], factor)

    dsn2 = as.data.frame(dsn1[dsn_othrs])
    for(i in fmts){
      fmt_vars = as.character(subset(proctcae_vars, fmt_name == i)$name[subset(proctcae_vars, fmt_name == i)$name %in% dsn_items])
      if (!identical(fmt_vars, character(0))){
        dsn_tmp = sapply(dsn1[fmt_vars], fmt = i, fmt_func)
        dsn2 = cbind(dsn2, dsn_tmp)
      }
    }
  } else {dsn2 = dsn}

  # --------------------------------------------------------------------------
  # --- Imputation
  # --------------------------------------------------------------------------

  # - Fills in a '0' score for PROCTCAE_##B/C fields if that PROCTCAE_##B/C is missing and the associated PROCTCAE_##A = 0.

  if (impute==TRUE){
    for(k in 1:7){
      for(i in as.numeric(strsplit(as.character(map_ref$qset[map_ref$rank==k]), " ")[[1]])){

        # ---------------------------------------------
        # -- Two-item groups (A/B)
        # ---------------------------------------------
        if(k %in% 4:6){
          # If all A/B items exist in the dataframe
          if(all(paste0("PROCTCAE_",i,c("A","B"),"_SCL") %in% dsn_items_comp)){

            tmp_logic = (!is.na(dsn2[paste0("PROCTCAE_",i,"A_SCL")]) & dsn2[paste0("PROCTCAE_",i,"A_SCL")]==0) & is.na(dsn2[paste0("PROCTCAE_",i,"B_SCL")])
            dsn2[paste0("PROCTCAE_",i,"B_SCL")][tmp_logic,] = 0
          }
        }

        # ---------------------------------------------
        # -- Three-item composite Grades
        # ---------------------------------------------
        else if(k == 7){
          # If all A/B/C variables are exist in the data frame
          if(all(paste0("PROCTCAE_",i,c("A","B"),"_SCL") %in% dsn_items_comp)){

            tmp_logic = (!is.na(dsn2[paste0("PROCTCAE_",i,"A_SCL")]) & dsn2[paste0("PROCTCAE_",i,"A_SCL")]==0) & is.na(dsn2[paste0("PROCTCAE_",i,"B_SCL")])
            dsn2[paste0("PROCTCAE_",i,"B_SCL")][tmp_logic,] = 0
          }

          if(all(paste0("PROCTCAE_",i,c("A","C"),"_SCL") %in% dsn_items_comp)){

            tmp_logic = (!is.na(dsn2[paste0("PROCTCAE_",i,"A_SCL")]) & dsn2[paste0("PROCTCAE_",i,"A_SCL")]==0) & is.na(dsn2[paste0("PROCTCAE_",i,"C_SCL")])
            dsn2[paste0("PROCTCAE_",i,"C_SCL")][tmp_logic,] = 0
          }
        }

      }
    }
  }

  # --------------------------------------------------------------------------
  # --- Composite Grading
  # --------------------------------------------------------------------------

  # - NOTE :
  # - If a freq item is included in the composite and is 0, then comp grade = 0
  # - If no freq item included, but there is a sev item in the composite and is 0, then comp grade = 0
  # - If no freq/sev item included (i.e. only an int item), and the int item is 0, then comp grade = 0

  if (composites==TRUE){
    for(k in 1:7){
      for(i in as.numeric(strsplit(as.character(map_ref$qset[map_ref$rank==k]), " ")[[1]])){

        # ---------------------------------------------
        # -- One-item composite grades
        # ---------------------------------------------
        if(k %in% 1:3){
          # If A item exists in the data.frame
          if(all(paste0("PROCTCAE_",i,c("A"),"_SCL") %in% dsn_items_comp)){

            comp_tab_tmp0 = comp_tab[comp_tab$rank==k,]

            # remove empty columns from composite ref table
            comp_tab_tmp = comp_tab_tmp0[,colSums(is.na(comp_tab_tmp0))<nrow(comp_tab_tmp0)]

            dsn2[paste0("PROCTCAE_",i,"_COMP")] = NA

            for(
              j in 1:nrow(dsn2)){

              # If row-wise A items are non-missing, construct composite
              if(all(!is.na(dsn2[c(paste0("PROCTCAE_",i,c("A"),"_SCL"))][j,]))){

                dsn2[paste0("PROCTCAE_",i,"_COMP")][j,] = comp_tab_tmp$comp[which(comp_tab_tmp[,1] == dsn2[paste0("PROCTCAE_",i,"A_SCL")][j,])]
              }
            }
          }
        }

        # ---------------------------------------------
        # -- Two-item composite grades
        # ---------------------------------------------
        else if(k %in% 4:6){

          # If all A/B items exisit in the dataframe
          if(all(paste0("PROCTCAE_",i,c("A","B"),"_SCL") %in% dsn_items_comp)){

            comp_tab_tmp0 = comp_tab[comp_tab$rank==k,]

            # remove empty columns from composite ref table
            comp_tab_tmp = comp_tab_tmp0[,colSums(is.na(comp_tab_tmp0))<nrow(comp_tab_tmp0)]

            dsn2[paste0("PROCTCAE_",i,"_COMP")] = NA

            for(j in 1:nrow(dsn2)){

              # If row-wise A/B items are non-missing, construct composite
              if(all(!is.na(dsn2[c(paste0("PROCTCAE_",i,c("A","B"),"_SCL"))][j,]))){

                # - If a freq item is included in the composite and is 0, then comp grade = 0
                # - If no freq item included, but there is a sev item in the composite and is 0, then comp grade = 0
                # - Note: For two-question item groups, question A is either a freq or sev item
                if(dsn2[paste0("PROCTCAE_",i,"A_SCL")][j,]==0){
                  dsn2[paste0("PROCTCAE_",i,"_COMP")][j,] = 0
                }
                else {
                  dsn2[paste0("PROCTCAE_",i,"_COMP")][j,] = comp_tab_tmp$comp[which(comp_tab_tmp[,1] == dsn2[paste0("PROCTCAE_",i,"A_SCL")][j,] &
                                                                                      comp_tab_tmp[,2] == dsn2[paste0("PROCTCAE_",i,"B_SCL")][j,])]
                }
              }
            }
          }
        }

        # ---------------------------------------------
        # -- Three-item composite grades
        # ---------------------------------------------
        else if(k == 7){

          # If all A/B/C variables are exisit in the dataframe
          if(all(paste0("PROCTCAE_",i,c("A","B","C"),"_SCL") %in% dsn_items_comp)){

            comp_tab_tmp = comp_tab[comp_tab$rank==k,]
            dsn2[paste0("PROCTCAE_",i,"_COMP")] = NA

            for(j in 1:nrow(dsn2)){

              # If all A/B/C variables are non-missing
              if(all(!is.na(dsn2[c(paste0("PROCTCAE_",i,c("A","B","C"),"_SCL"))][j,]))){

                # If a freq item is included in the composite and is 0, then comp grade = 0
                if(dsn2[paste0("PROCTCAE_",i,"A_SCL")][j,]==0){
                  dsn2[paste0("PROCTCAE_",i,"_COMP")][j,] = 0
                }
                else {
                  dsn2[paste0("PROCTCAE_",i,"_COMP")][j,] = comp_tab_tmp$comp[which(comp_tab_tmp$frq == dsn2[paste0("PROCTCAE_",i,"A_SCL")][j,] &
                                                                                      comp_tab_tmp$sev == dsn2[paste0("PROCTCAE_",i,"B_SCL")][j,] &
                                                                                      comp_tab_tmp$int == dsn2[paste0("PROCTCAE_",i,"C_SCL")][j,])]
                }
              }
            }
          }
        }

      }
    }
  }

  # --------------------------------------------------------------------------
  # --- Apply short labels for existing PRO-CTCAE items and composite grades
  # --------------------------------------------------------------------------

  if(short_labels==TRUE){

    proctcae_vars_tmp0 = as.data.frame(lapply(proctcae_vars[,-1], as.character))

    proctcae_vars_tmp = c()
    proctcae_vars_tmp$name = paste0(substr(proctcae_vars_tmp0$name, 1, nchar(proctcae_vars_tmp0$name)-5), "_COMP")
    proctcae_vars_tmp$short_label = sub(proctcae_vars_tmp0$short_label, pattern = " [[:alpha:]]*$", replacement = "")
    proctcae_vars_tmp = as.data.frame(proctcae_vars_tmp)

    proctcae_vars_tmp = unique(proctcae_vars_tmp)
    proctcae_labels0 = rbind.data.frame(proctcae_vars_tmp0, proctcae_vars_tmp)

    proctcae_labels1 = t(proctcae_labels0)
    colnames(proctcae_labels1) = proctcae_labels0$name
    rownames(proctcae_labels1) = NULL
    proctcae_labels = proctcae_labels1[-1,]

    Hmisc::label(dsn2) = as.list(proctcae_labels[match(names(dsn2), names(proctcae_labels))])
  }

  # -- Reorder dsn2 columns to match dsn (then composites in order if created here)
  if(length(colnames(dsn2)[!(names(dsn2) %in% names(dsn))]) > 0){
    col_ranks = data.frame(colnames_end = colnames(dsn2)[!(names(dsn2) %in% names(dsn))])
    col_ranks$colnames_num = as.numeric(gsub(".*PROCTCAE_\\s*|_COMP.*", "", col_ranks$colnames_end))
    col_end = levels(factor(col_ranks[order(col_ranks$colnames_num),]$colnames_end))
  }else{
    col_end = c()
  }

  dsn2 = dsn2[,c(colnames(dsn),col_end)]

  # --------------------------------------------------------------------------
  # --- Exit
  # --------------------------------------------------------------------------

  return(dsn2)

}
