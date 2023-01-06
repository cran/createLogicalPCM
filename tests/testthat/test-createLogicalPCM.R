test_that("errors", {
  testthat::expect_error(
    createLogicalPCM(NA, c(1,2,3)),
    "The first parameter is mandatory"
  )

  testthat::expect_error(
    createLogicalPCM(3.25),
    "The first parameter has to be an integer"
  )

  testthat::expect_error(
    createLogicalPCM(3, c(1, "a", 2.5)),
    "The second parameter has to be a numeric vector"
  )

  testthat::expect_error(
    createLogicalPCM(3, c(1,2,.75, 3)),
    "The length of the second parameter has to be the same as the first parameter"
  )

})
