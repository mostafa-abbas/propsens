# Package skeleton: the package loads and exports its functions.

test_that("package skeleton loads and declares the MIT licence", {
  description <- utils::packageDescription("propsens")
  expect_equal(description$Package, "propsens")
  expect_match(description$License, "MIT")
})
