
#  Dwumianowe drzewo wyceny opcji  Cox-Ross-Rubinstein
crr <-
  function(TypeFlag = c("ce", "pe", "ca", "pa"), S, X, Time, r, b, sigma, n) {

    # Check Flags:
    TypeFlag <- TypeFlag[1]
    z <- NA
    if (TypeFlag == "ce" || TypeFlag == "ca") z <- +1
    if (TypeFlag == "pe" || TypeFlag == "pa") z <- -1
    if (is.na(z)) stop("TypeFlag misspecified: ce|ca|pe|pa")

    # Parametry
    dt <- Time / n
    u <- exp(sigma * sqrt(dt))
    d <- 1 / u
    p <- (exp(b * dt) - d) / (u - d)
    Df <- exp(-r * dt)

    # Iteration:
    OptionValue <- z * (S * u^(0:n) * d^(n:0) - X)
    OptionValue <- (abs(OptionValue) + OptionValue) / 2

    # European Option:
    if (TypeFlag == "ce" || TypeFlag == "pe") {
      for (j in seq(from = n - 1, to = 0, by = -1)) {
        for (i in 0:j) {
          OptionValue[i + 1] <-
            (p * OptionValue[i + 2] + (1 - p) * OptionValue[i + 1]) * Df
        }
      }
    }

    # American Option:
    if (TypeFlag == "ca" || TypeFlag == "pa") {
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

jr <-
  function(TypeFlag = c("ce", "pe", "ca", "pa"), S, X, Time, r, b, sigma, n) {

    # Description:
    #   JR Modfication to the Binomial Tree Option

    # FUNCTION:

    # Check Flags:
    TypeFlag <- TypeFlag[1]
    if (TypeFlag == "ce" || TypeFlag == "ca") z <- +1
    if (TypeFlag == "pe" || TypeFlag == "pa") z <- -1

    # Parameters:
    dt <- Time / n
    # DW Bug Fix: r -> b
    u <- exp((b - sigma^2 / 2) * dt + sigma * sqrt(dt))
    d <- exp((b - sigma^2 / 2) * dt - sigma * sqrt(dt))
    # DW End of Bug Fix
    p <- 1 / 2
    Df <- exp(-r * dt)

    # Iteration:
    OptionValue <- z * (S * u^(0:n) * d^(n:0) - X)
    OptionValue <- (abs(OptionValue) + OptionValue) / 2

    # European Option:
    if (TypeFlag == "ce" || TypeFlag == "pe") {
      for (j in seq(from = n - 1, to = 0, by = -1)) {
        for (i in 0:j) {
          OptionValue[i + 1] <-
            (p * OptionValue[i + 2] + (1 - p) * OptionValue[i + 1]) * Df
        }
      }
    }

    # American Option:
    if (TypeFlag == "ca" || TypeFlag == "pa") {
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
  function(TypeFlag = c("ce", "pe", "ca", "pa"), S, X, Time, r, b, sigma, n) {
    # FUNCTION:

    # Check Flags:
    TypeFlag <- TypeFlag[1]
    if (TypeFlag == "ce" || TypeFlag == "ca") z <- +1
    if (TypeFlag == "pe" || TypeFlag == "pa") z <- -1

    # Parameters:
    dt <- Time / n
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
    if (TypeFlag == "ce" || TypeFlag == "pe") {
      for (j in seq(from = n - 1, to = 0, by = -1)) {
        for (i in 0:j) {
          OptionValue[i + 1] <-
            (p * OptionValue[i + 2] + (1 - p) * OptionValue[i + 1]) * Df
        }
      }
    }

    # American Option:
    if (TypeFlag == "ca" || TypeFlag == "pa") {
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
  TypeFlag = "pa", S = 50, X = 50,
  Time = 5 / 12, r = 0.1, b = 0.1, sigma = 0.4, n = 5
)

# Another example
# Expected value: 4.919211
crr(
  TypeFlag = "pa", S = 100, X = 95,
  Time = 0.5, r = 0.08, b = 0.08, sigma = 0.3, n = 5
)

# Expected Value: 14.93
crr(
  TypeFlag = "ce", S = 100, X = 100,
  Time = 1, r = 0.1, b = 0.1, sigma = 0.25, n = 50
)

jr(
  TypeFlag = "ce", S = 100, X = 100,
  Time = 1, r = 0.1, b = 0.1, sigma = 0.25, n = 50
)


tian(
  TypeFlag = "ce", S = 100, X = 100,
  Time = 1, r = 0.1, b = 0.1, sigma = 0.25, n = 50
)


TIANBinomialTreeOption(
  TypeFlag = "ce", S = 100, X = 100,
  Time = 1, r = 0.1, b = 0.1, sigma = 0.25, n = 50
)

JRBinomialTreeOption(
  TypeFlag = "ce", S = 100, X = 100,
  Time = 1, r = 0.1, b = 0.1, sigma = 0.25, n = 50
)