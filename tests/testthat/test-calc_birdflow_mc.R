test_that("calc_birdflow_mc() produces consistent results", {

  bf <- BirdFlowModels::amewoo
  mc <- calc_birdflow_mc(bf, season = "prebreeding")
  expect_equal(mc, 0.394119668270)
  # Before bringing MC calculation into  calc_birdflow_mc():
  # 0.394119668270, 10.3 seconds, max total memory = 8255.3 Mb
  # After:
  # 0.394119678294, 463 milliseconds, max total memory   67.0 Mb

})
