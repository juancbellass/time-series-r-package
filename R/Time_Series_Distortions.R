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
#' @param beta the factor that multiplies the minimum and maximum of the series to determine the value of salt and pepper respectively
#' @return The input time series with the salt and pepper noise variation imposed.
#' @export
TS.SaltPepper <- function(Y, prop, beta=1) {
  if (prop<0 | prop>1) stop('epsi must be in the interval [0,1]')
  n = length(Y) 
  salt = beta * max(Y)
  pepper = beta * min(Y)
  ts.bin = rbinom(n, 1, prop)
  salt_pepper = rbinom(n, 1, 0.5) * (salt - pepper) + pepper
  Y2 = (1 - ts.bin) * Y # unaltered data from the original time-series
  return(Y2 + ts.bin * salt_pepper)
}

#' Gaussian Noise
#' 
#' Generates a vector _ts.normal_ is created where every element of _ts.normal_ follow a 
#' N(mu, sigma2) distribution. The _prop_% of the elements of the array _Y_ are 
#' selected randomly. The _prop_% of the elements of the array _Y_ are selected 
#' randomly. If Y_i_ is selected  then it is replaced by _Y_i_ + _ts.normal_i_. 
#' This procedure is called additive contamination of _prob_%.
#' @param Y The time-series (as numeric) to be contaminated. 
#' @param prob The probability of success of the binomial distribution. A number between 0 and 1.
#' @param mu The mean of the normal distribution.
#' @param sigma The standard deviation of the normal distribution. A real number.
#' @return The y matrix with additive noise.
#' @export
GaussianNoise <- function(y, prob, mu, sigma) {
  if (prob < 0 || prob > 1) stop('prob must be in the interval [0,1]')
  ts.bin <- rbinom(length(y), 1, prob)
  ts.normal <- rnorm(sum(ts.bin), mean = mu, sd = sigma)
  y[ts.bin == 1] <- y[ts.bin == 1] + ts.normal
  return(y)
}

#' t-Student Noise
#' 
#' A vector ts.student is created where every element of _ts.student_ follow a 
#' t(df) distribution. The _prop_% of the elements of the array _Y_ are 
#' selected randomly. If Y_i_ is selected then it is replaced by _Y_i_ + _ts.student_i_.
#' @param Y The time serie as *numeric* to be contaminated.
#' @param prob The probability of success of the binomial distribution. A number between 0 and 1.
#' @param dfr The degrees of freedom of the t-student distribution.
#' @return The Y matrix with t-Student noise.
#' @export
TStudentNoise <- function(y, prob, dfr) {
  if (prob < 0 || prob > 1) stop('prob must be in the interval [0,1]')
  ts.bin <- rbinom(length(y), 1, prob)
  ts.tst <- rt(sum(ts.bin), df = dfr)
  y[ts.bin == 1] <- y[ts.bin == 1] + ts.tst
  return(y)
}

#' Contaminación por otro proceso autorregresivo.
#' 
#' Se genera un vector llamado "ts.cont" a partir de un vector inicial "y", y otrvector 
#' "ts.ar"de la misma dimensión. Los valores de "ts.ar"son generados con un proceso AR con parámetros φ. 
#' El prob,100% de los valores de "y"son elegidos al azar. Si y_i es seleccionado, entonces 
#' es reemplazado por la entrada i del vector "ts.ar". Este procedimiento es llamado 
#' contaminación de reemplazo por otro proceso AR al prob,100%.
#' @param y The time serie as *numeric* to be contaminated.
#' @param prob The probability of success of the binomial distribution. A number between 0 and 1.
#' @param phi1 The parameter model. A number in the interval (0,1).
#' @return The Y matrix with autoregressive process noise.
#' @export
ReplARCont <- function(y, prob=1, phi1) {
  if (prob < 0 || prob > 1) stop('prob must be in the interval [0,1]')
  if (phi1 <= 0 || phi1 >= 1) stop('phi1 must be in the interval (0,1)')
  ts.bin <- rbinom(length(y), 1, prob)
  ar1.sim <- arima.sim(model=list(ar=phi1), n=sum(ts.bin))
  y2 <- y * (1 - ts.bin)
  return(y2 + ar1.sim * ts.bin)
}