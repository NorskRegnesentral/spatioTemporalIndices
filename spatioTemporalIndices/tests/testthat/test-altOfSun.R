test_that("altOfSun works correctly", {
  expect_equal(round(altOfSun(55, 6,28,4,70,20)$alt.of.sun,1), 24.5)
  expect_equal(round(altOfSun(55, 6,28,4,70,20)$sun.rise,1), 1.7)
})
