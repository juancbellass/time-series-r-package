#' Amplitude Scaling variation imposed
#' 
#' A time series distortion that modify the amplitude scale of a time series _S_.
#' @details $S1(t) = beta * S(t)$
#' @param S A time series with real values.
#' @param beta A constant that is multiplied to all point in the time series. 
#' @return The input time series with the amplitude scaling variation imposed. 
#' @export
AmpSca <- function(S, beta=1) {
  S1 = ts(beta * S)
  return(S1)
}

#' Amplitude Shift variation imposed
#' 
#' A time series distortion that shift the amplitude range of a time series _S_.
#' @details $S1(t) = beta + S(t)$
#' @param S A time series with real values.
#' @param beta A constant that is added to all points in the time series. 
#' @return The input time series with the amplitude shift variation imposed.
#' @export
AmpShift <- function(S, beta=1) {
  S1 = ts(beta + S)
  return(S1)
}

#' Time Scaling variation imposed
#' 
#' A time series distortion that modify the time length of a time series _S_.
#' @details $S1(t) = S(t*beta)$
#' @param S A time series with real values.
#' @param beta a natural constant. 
#' @return The input time series with the time scaling variation imposed.
#' @export
TimeSca <- function(S, beta=1) {
  S1 = ts(beta + S)
  return(S1)
}

#' Time Shift variation imposed
#' 
#' A time series distortion that shift the time in a time series _S_.
#' @details $S1(t) = S(t+beta)$
#' @param S A time series with real values.
#' @param beta A natural constant. 
#' @return The input time series with the time shift variation imposed.
#' @export
TimeShift <- function(S, beta=1) {
  S1 = ts(beta + S)
  return(S1)
}

#' Salt and pepper variation imposed
#' 
#' Similar to salt and pepper image noise but applying to time series. In this case, 
#' salt and pepper variations correspond to max(Y) and min(Y) respectively.
#' @param Y A time series with real values.
#' @param prop contamination percentage for the salt and pepper noise [0,1]. 
#' @return The input time series with the salt and pepper noise variation imposed.
#' @export
TS.SaltPepper <- function(Y, prop) {
  if (prop<0 | prop>1) {stop('epsi must be in the interval [0,1]')}
  n = length(Y) 
  salt = max(Y)
  pepper = min(Y)
  ts.bin = rbinom(n, 1, prop)
  Y2 = (abs(ts.bin-1)*Y) # unaltered data from the original time-series
  for (i in 1:n) {
    if (ts.bin[i] == 1) {ifelse(rbinom(1,1,.5) == 0, 
                                ts.bin[i] <- salt, ts.bin[i] <- pepper)}}
  return(Y2 + ts.bin)
}

#' Additive Gaussian Noise
#' 
#' Generates a matrix _C_ from an initial matrix _Y_ and a matrix _Nu_ of the same 
#' dimension. _Nu_ values follow an _N_(_mu_, _sigma_^2) distribution. The _prop_% of 
#' the elements of the array _Y_ are selected randomly. If Y_ij_ is selected 
#' then it is replaced by _Y_ij_ + _Nu_ij_. This procedure is called additive 
#' contamination of _prob_%.
#' @param Y The time-series (as numeric) to be contaminated. 
#' @param prob The probability of success of the binomial distribution. A number between 0 and 1.
#' @param mu The mean of the normal distribution.
#' @param sigma The standard deviation of the normal distribution. A real number.
#' @return The y matrix with additive noise.
#' @export
AddGaussianNoise <- function(y, prob, mu, sigma) {
  n = length(y)
  ts.bin = rbinom(n, 1, prob)
  ts.normal = rnorm(n, mean=mu, sd=sigma)
  return(y + (ts.bin * ts.normal)) 
}

#' Gaussian Noise
#' 
#' Bla bla bla
#' @param Y The time-series (as numeric) to be contaminated. 
#' @param prob The probability of success of the binomial distribution. A number between 0 and 1.
#' @param \mu The mean of the normal distribution.
#' @param \sigma The standard deviation of the normal distribution. A real number.
#' @return The Y matrix with additive noise.
#' @export
GaussianNoise <- function(y, prob, mu, sigma) {
  n = length(y)
  ts.bin = rbinom(n, 1, prob)
  ts.normal = rnorm(n, mean=mu, sd=sigma)
  return((abs(ts.bin-1)*y) + (ts.bin*ts.normal))
}

#' t-Student Noise
#' 
#' Bla bla bla
#' @param Y The time serie as *numeric* to be contaminated.
#' @param prob The probability of success of the binomial distribution. A number between 0 and 1.
#' @param dfr The degrees of freedom of the t-student distribution.
#' @return The Y matrix with t-Student noise.
#' @export
TStudentNoise <- function(y, prob, dfr) {
  n = length(y)
  ts.bin = rbinom(n, 1, prob)
  ts.tst = rt(n, df=dfr)
  return((abs(ts.bin-1)*y) + (ts.bin*ts.tst))
}

#' Contaminación por otro proceso autorregresivo.
#' 
#' Se genrea un vector llamado "ts.cont" a partir de un vector inicial "y", y otrvector 
#' "ts.ar"de la misma dimensión. Los valores de "ts.ar"son generados con un proceso AR con parámetros φ. 
#' El prob,100% de los valores de "y"son elegidos al azar. Si y_i es seleccionado, entonces 
#' es reemplazado por la entrada i del vector "ts.ar". Este procedimiento es llamado 
#' contaminación de reemplazo por otro proceso AR al prob,100%.
#' @param y The time serie as *numeric* to be contaminated.
#' @param prob The probability of success of the binomial distribution. A number between 0 and 1.
#' @param phi1
#' @return The Y matrix with autoregressive process noise.
#' @export
ReplARCont <- function(y, prob, phi1) {
  N = length(y)
  ts.bin = rbinom(N, 1, prob)
  ar1.sim = arima.sim(model = list(ar = phi1), n = N)
  return((abs(ts.bin-1)*y) + (ts.bin*ar1.sim))
}


