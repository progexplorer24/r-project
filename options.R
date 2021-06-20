
#  Dwumianowe drzewo wyceny opcji  Cox-Ross-Rubinstein
crr <-
  function(type_flag = c("ce", "pe", "ca", "pa"), S, X, time, r, b, sigma, n) {
    type_flag <- type_flag[1]
    z <- NA
    if (type_flag == "ce" || type_flag == "ca") z <- +1
    if (type_flag == "pe" || type_flag == "pa") z <- -1
    if (is.na(z)) stop("Wprowadzono złą flagę. Wybierz jedną z opcji: ce|pe|ca|pa")

    # Parametry
    dt <- time / n
    u <- exp(sigma * sqrt(dt))
    d <- 1 / u
    p <- (exp(b * dt) - d) / (u - d)
    Df <- exp(-r * dt)

    # Iteration:
    OptionValue <- z * (S * u^(0:n) * d^(n:0) - X)
    OptionValue <- (abs(OptionValue) + OptionValue) / 2

    # European Option:
    if (type_flag == "ce" || type_flag == "pe") {
      for (j in seq(from = n - 1, to = 0, by = -1)) {
        for (i in 0:j) {
          OptionValue[i + 1] <-
            (p * OptionValue[i + 2] + (1 - p) * OptionValue[i + 1]) * Df
        }
      }
    }

    # American Option:
    if (type_flag == "ca" || type_flag == "pa") {
      for (j in seq(from = n - 1, to = 0, by = -1)) {
        for (i in 0:j) {
          OptionValue[i + 1] <- max(
            (z * (S * u^i * d^(abs(i - j)) - X)),
            (p * OptionValue[i + 2] + (1 - p) * OptionValue[i + 1]) * Df
          )
        }
      }
    }

    # Return Value:
    OptionValue[1]
  }

# -------------------------------------------------------------------------------------
#  Dwumianowe drzewo wyceny opcji  Jarrow Rudd
jr <-
  function(type_flag = c("ce", "pe", "ca", "pa"), S, X, time, r, b, sigma, n) {
    type_flag <- type_flag[1]
    z <- NA
    if (type_flag == "ce" || type_flag == "ca") z <- +1
    if (type_flag == "pe" || type_flag == "pa") z <- -1
    if (is.na(z)) stop("Wprowadzono złą flagę. Wybierz jedną z opcji: ce|pe|ca|pa")
    # Parameters:
    dt <- time / n
    u <- exp((b - sigma^2 / 2) * dt + sigma * sqrt(dt))
    d <- exp((b - sigma^2 / 2) * dt - sigma * sqrt(dt))
    p <- 1 / 2
    Df <- exp(-r * dt)

    # Iteration:
    OptionValue <- z * (S * u^(0:n) * d^(n:0) - X)
    OptionValue <- (abs(OptionValue) + OptionValue) / 2

    # European Option:
    if (type_flag == "ce" || type_flag == "pe") {
      for (j in seq(from = n - 1, to = 0, by = -1)) {
        for (i in 0:j) {
          OptionValue[i + 1] <-
            (p * OptionValue[i + 2] + (1 - p) * OptionValue[i + 1]) * Df
        }
      }
    }

    # American Option:
    if (type_flag == "ca" || type_flag == "pa") {
      for (j in seq(from = n - 1, to = 0, by = -1)) {
        for (i in 0:j) {
          OptionValue[i + 1] <- max(
            (z * (S * u^i * d^(abs(i - j)) - X)),
            (p * OptionValue[i + 2] + (1 - p) * OptionValue[i + 1]) * Df
          )
        }
      }
    }

    # Return Value:
    OptionValue[1]
  }

# -------------------------------------------------------------------------------------

tian <-
  function(type_flag = c("ce", "pe", "ca", "pa"), S, X, time, r, b, sigma, n) {
    # FUNCTION:

    # Check Flags:
    type_flag <- type_flag[1]
    if (type_flag == "ce" || type_flag == "ca") z <- +1
    if (type_flag == "pe" || type_flag == "pa") z <- -1

    # Parameters:
    dt <- time / n
    M <- exp(b * dt)
    V <- exp(sigma^2 * dt)
    u <- (M * V / 2) * (V + 1 + sqrt(V * V + 2 * V - 3))
    d <- (M * V / 2) * (V + 1 - sqrt(V * V + 2 * V - 3))
    p <- (M - d) / (u - d)
    Df <- exp(-r * dt)

    # Iteration:
    OptionValue <- z * (S * u^(0:n) * d^(n:0) - X)
    OptionValue <- (abs(OptionValue) + OptionValue) / 2

    # European Option:
    if (type_flag == "ce" || type_flag == "pe") {
      for (j in seq(from = n - 1, to = 0, by = -1)) {
        for (i in 0:j) {
          OptionValue[i + 1] <-
            (p * OptionValue[i + 2] + (1 - p) * OptionValue[i + 1]) * Df
        }
      }
    }

    # American Option:
    if (type_flag == "ca" || type_flag == "pa") {
      for (j in seq(from = n - 1, to = 0, by = -1)) {
        for (i in 0:j) {
          OptionValue[i + 1] <- max(
            (z * (S * u^i * d^(abs(i - j)) - X)),
            (p * OptionValue[i + 2] + (1 - p) * OptionValue[i + 1]) * Df
          )
        }
      }
    }


    OptionValue[1]
  }

# Example from Hull's book:
# Expected value: 4.488459
crr(
  type_flag = "pa", S = 50, X = 50,
  time = 5 / 12, r = 0.1, b = 0.1, sigma = 0.4, n = 5
)

# Another example
# Expected value: 4.92
crr(
  type_flag = "pa", S = 100, X = 95,
  time = 0.5, r = 0.08, b = 0.08, sigma = 0.3, n = 5
)

# Expected Value: 14.93
crr(
  type_flag = "ce", S = 100, X = 100,
  time = 1, r = 0.1, b = 0.1, sigma = 0.25, n = 50
)

# Expected Value: 14.94
jr(
  type_flag = "ce", S = 100, X = 100,
  time = 1, r = 0.1, b = 0.1, sigma = 0.25, n = 50
)

# Expected Value: 14.99
tian(
  type_flag = "ce", S = 100, X = 100,
  time = 1, r = 0.1, b = 0.1, sigma = 0.25, n = 50
)