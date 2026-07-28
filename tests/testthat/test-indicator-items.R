indicator_data <- data.frame(
  id = rep(1:4, each = 3),
  Cycle = rep(1:3, times = 4),
  arm = rep(c("A", "A", "B", "B"), each = 3),
  PROCTCAE_73A_IND = c(
    "NO", "YES", "YES",
    "NO", "NO", "YES",
    "YES", "YES", "YES",
    "NO", NA, "NO"
  ),
  stringsAsFactors = FALSE
)

test_that("toxScores reformats indicator items", {
  scored <- toxScores(indicator_data, reformat = TRUE)
  expect_equal(scored$PROCTCAE_73A_IND,
               c(0, 1, 1, 0, 0, 1, 1, 1, 1, 0, NA, 0))
})

test_that("toxScores does not create indicator composites", {
  scored <- toxScores(indicator_data, reformat = TRUE, composites = TRUE)
  expect_false("PROCTCAE_73_COMP" %in% names(scored))
})

test_that("toxTables reports 1+ and leaves 3+ blank for indicators", {
  scored <- toxScores(indicator_data, reformat = TRUE)
  tab <- toxTables(scored, id_var = "id", cycle_var = "Cycle",
                   baseline_val = 1, type = "max_post_bl")
  expect_equal(nrow(tab$individual), 1)
  expect_match(tab$individual$overall_pres, "3 \\(75%\\)")
  expect_identical(tab$individual$overall_sev, "")
})

test_that("toxTables supports indicator items by arm", {
  scored <- toxScores(indicator_data, reformat = TRUE)
  tab <- toxTables(scored, id_var = "id", cycle_var = "Cycle",
                   baseline_val = 1, type = "max_post_bl", arm_var = "arm")
  expect_true(all(c("A_pres", "B_pres", "A_sev", "B_sev") %in%
                    names(tab$individual)))
  expect_identical(tab$individual$A_sev, "")
  expect_identical(tab$individual$B_sev, "")
})

test_that("toxFigures returns a binary indicator figure", {
  scored <- toxScores(indicator_data, reformat = TRUE)
  figs <- toxFigures(scored, id_var = "id", cycle_var = "Cycle",
                     baseline_val = 1, cycles_only = FALSE)
  expect_length(figs, 1)
  expect_s3_class(figs[[1]][[2]], "ggplot")
  built <- ggplot2::ggplot_build(figs[[1]][[2]])
  expect_true(any(built$data[[1]]$fill != "white"))
})

test_that("toxFigures combines score and indicator figures", {
  scored <- toxScores(indicator_data, reformat = TRUE)
  scored$PROCTCAE_1A_SCL <- rep(c(0, 1, 2), 4)
  figs <- toxFigures(scored, id_var = "id", cycle_var = "Cycle",
                     baseline_val = 1)
  expect_gte(length(figs), 2)
  expect_true(all(vapply(figs, function(x) inherits(x[[2]], "ggplot"), logical(1))))
})
