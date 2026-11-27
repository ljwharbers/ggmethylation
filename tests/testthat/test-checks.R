test_that(".check_number accepts a single number in range and returns it", {
  expect_identical(.check_number(3, "x"), 3)
  expect_identical(.check_number(0, "x", min = 0), 0)
  expect_identical(.check_number(0.5, "x", min = 0, max = 1), 0.5)
})

test_that(".check_number rejects non-numbers, vectors, NA and out-of-range values", {
  expect_error(.check_number("3", "x"), "`x` must be a single number\\.")
  expect_error(.check_number(c(1, 2), "x"), "`x` must be a single number")
  expect_error(.check_number(NA_real_, "x"), "`x` must be a single number")
  expect_error(.check_number(NULL, "x"), "`x` must be a single number")
  expect_error(.check_number(-1, "x", min = 0), "`x` must be a single number >= 0\\.")
  expect_error(.check_number(2, "x", max = 1), "`x` must be a single number <= 1\\.")
  expect_error(.check_number(2, "x", min = 0, max = 1),
               "`x` must be a single number in \\[0, 1\\]\\.")
})

test_that(".check_number handles infinite values and whole numbers", {
  expect_error(.check_number(Inf, "x"), "`x` must be a single number")
  expect_identical(.check_number(Inf, "x", allow_inf = TRUE), Inf)
  expect_identical(.check_number(4, "x", whole = TRUE), 4L)
  expect_identical(.check_number(Inf, "x", whole = TRUE, allow_inf = TRUE), Inf)
  expect_error(.check_number(2.5, "x", min = 1, whole = TRUE),
               "`x` must be a single whole number >= 1\\.")
})

test_that(".check_flag accepts only a single TRUE or FALSE", {
  expect_true(.check_flag(TRUE, "f"))
  expect_false(.check_flag(FALSE, "f"))
  expect_error(.check_flag(NA, "f"), "`f` must be TRUE or FALSE\\.")
  expect_error(.check_flag(c(TRUE, FALSE), "f"), "`f` must be TRUE or FALSE")
  expect_error(.check_flag("yes", "f"), "`f` must be TRUE or FALSE")
  expect_error(.check_flag(1, "f"), "`f` must be TRUE or FALSE")
})

test_that(".check_string accepts only a single non-empty string", {
  expect_identical(.check_string("a", "s"), "a")
  expect_error(.check_string("", "s"), "`s` must be a single non-empty character string")
  expect_error(.check_string(NA_character_, "s"), "`s` must be a single")
  expect_error(.check_string(c("a", "b"), "s"), "`s` must be a single")
  expect_error(.check_string(1, "s"), "`s` must be a single")
})

test_that(".check_colours accepts names and hex codes and keeps names", {
  expect_identical(.check_colours("red", "col", n = 1L), "red")
  cols = c("+" = "#4393C3", "-" = "#D6604D")
  expect_identical(.check_colours(cols, "col"), cols)
})

test_that(".check_colours rejects invalid colours, NAs and wrong lengths", {
  expect_error(.check_colours("notacolour", "col"),
               "`col` contains invalid colour\\(s\\): notacolour\\.")
  expect_error(.check_colours(c("red", NA), "col"), "`col` must be a character vector")
  expect_error(.check_colours(1, "col"), "`col` must be a character vector")
  expect_error(.check_colours(character(0), "col"), "`col` must be a character vector")
  expect_error(.check_colours(c("red", "blue"), "col", n = 1L),
               "`col` must be a single colour")
})
