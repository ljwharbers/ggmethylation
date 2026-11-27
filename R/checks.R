# Argument validators shared by every exported function. Each stops with a
# message that names the argument in backticks (call. = FALSE) and otherwise
# returns the checked value, so a call can be used as `x = .check_number(x, "x")`.

# A single non-missing number with `min <= x <= max`. Infinite values are
# rejected unless `allow_inf`. With `whole`, the value must be a whole number
# and is returned as an integer (an infinite value is returned unchanged).
.check_number = function(x, name, min = -Inf, max = Inf, whole = FALSE,
                         allow_inf = FALSE) {
  ok = is.numeric(x) && length(x) == 1L && !is.na(x) &&
    (allow_inf || is.finite(x)) && x >= min && x <= max
  if (ok && whole && is.finite(x)) ok = x == round(x)
  if (!ok) {
    bounds = if (is.finite(min) && is.finite(max)) {
      sprintf(" in [%s, %s]", format(min), format(max))
    } else if (is.finite(min)) {
      sprintf(" >= %s", format(min))
    } else if (is.finite(max)) {
      sprintf(" <= %s", format(max))
    } else {
      ""
    }
    stop(sprintf("`%s` must be a single %s%s.", name,
                 if (whole) "whole number" else "number", bounds), call. = FALSE)
  }
  if (whole && is.finite(x)) as.integer(x) else x
}

# A single TRUE or FALSE.
.check_flag = function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}

# A single non-missing, non-empty character string.
.check_string = function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("`%s` must be a single non-empty character string.", name),
         call. = FALSE)
  }
  x
}

# A character vector of colours R can draw (names or hex codes), without NAs.
# `n`, when given, is the required length.
.check_colours = function(x, name, n = NULL) {
  if (!is.character(x) || length(x) == 0L || anyNA(x) ||
        (!is.null(n) && length(x) != n)) {
    stop(sprintf("`%s` must be %s without NAs.", name,
                 if (is.null(n)) "a character vector of colours"
                 else if (n == 1L) "a single colour"
                 else sprintf("%d colours", n)), call. = FALSE)
  }
  valid = vapply(x, function(col) {
    !inherits(tryCatch(grDevices::col2rgb(col), error = function(e) e), "error")
  }, logical(1L))
  if (!all(valid)) {
    stop(sprintf("`%s` contains invalid colour(s): %s.", name,
                 paste(x[!valid], collapse = ", ")), call. = FALSE)
  }
  x
}
