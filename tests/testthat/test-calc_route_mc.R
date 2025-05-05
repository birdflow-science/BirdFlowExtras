test_that("mc on synthetic routes converges on BirdFlow mc", {

  skip("Slow test, always skipped") # 4 minutes is too long for now.

  bf <- BirdFlowModels::amewoo
  expect_no_error(rts <- route(bf, n = 10000, season = "postbreeding"))
  # takes about 4 minutes.  Smaller n didn't seem to speed it up
  t <- system.time(
    expect_no_error(mc <- calc_route_mc(rts, bf, season = "postbreeding"))
  )
  # mc is 0.2705

})


test_that("mc works on routes", {

  bf <- BirdFlowModels::amewoo
  set.seed(1)
  expect_no_error(rts <- route(bf, n = 500, season = "postbreeding"))

  expect_no_error(mc <- calc_route_mc(rts, bf, season = "postbreeding"))

  bf_mc <- calc_birdflow_mc(bf, season = "postbreeding")

  expect_equal(mc, bf_mc, tolerance = 0.1)

})
