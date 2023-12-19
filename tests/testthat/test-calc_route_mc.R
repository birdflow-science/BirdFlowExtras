test_that("mc on routes works", {

  skip_if(TRUE) # 4 minutes is too long for now.

  bf <- BirdFlowModels::amewoo
  expect_no_error(rts <- route(bf, n = 10000, season = "postbreeding"))
  # taskes about 4 minutes.  Smaller n didn't seem to speed it up
  t <- system.time(
    expect_no_error(mc <- calc_route_mc(rts, bf))
  )
  # mc = 0.2705

})
