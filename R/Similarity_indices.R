#' FRECHET DISTANCE 
#' 
#' Computes the Frechet distance between two numerical trajectories.
#' @param s1 a numeric vector containing the first time series.
#' @param s2 a numeric vector containing the second time series.
#' @return a positive real number resulting from the the calculated distance between the pair of series.
#' @import SimilarityMeasures
#' @export
Frechets.Measure <- function(s1, s2){
  return(Frechet(t(as.matrix(s1)), t(as.matrix(s2))))
}

#' PEARSON MEASURE 
#' 
#' Computes the Pearson correlation between two numerical vectors and return the result of 1-abs(Person).
#' @param s1 a numeric vector containing the first time series.
#' @param s2 a numeric vector containing the second time series.
#' @return a positive real number between 0 and 1. Where 0 is the maximum similarity and 1 the maximum dissimilarity.
#' @import stats
#' @export
Pearsons.Measure <- function(s1, s2){return(1-abs(cor(s1, s2)))}


#' SSIMT INDEX
#' 
#' Computes SSIMT similarity index between two same length time series X and Y.
#' @param X a time series of length n with real values.
#' @param Y a time series of length n with real values.
#' @param FUN mean by default.
#' @param ... some other parameter required by the function FUN.
#' @return a real number between 1 and -1; where 1 means max similarity, -1 means opposite behavior and 0 no similarity. 
#' @export
SSIMT <- function(x, y, FUN=mean, ...){
  X = as.numeric(x)
  Y = as.numeric(y)
  #primero calculamos los elementos a utilizar para evitar carcularlos a cada rato
  #constante incluida para eliminar la inestabilidad que se da cuando las suma de las medias al cuadarado es muy proxima a cero
  C1 = 576e-6 # =(k1*L)^2  ; Por recomendacion del autor se toma k1=0.01 y L=2.4
  C2 = 5184e-6 # =(k2*L)^2 ; Por recomendacion del autor se toma k1=0.01 y L=2.4
  C3 = 2592e-6 # =c2/2
  muX = FUN(X, ...) 
  muY = FUN(Y, ...)
  sdX = sqrt(sum((X - muX)^2) / (length(X) - 1))
  sdY = sqrt(sum((Y - muY)^2) / (length(Y) - 1))
  vXY = sum((X - muX)*(Y - muY)) / (length(Y)-1)

  #relacion de luminancia
  l.xy = ((2 * muX * muY) + C1 ) / (muX^2 + muY^2 + C1)
  #relacion de contraste
  c.xy = ((2 * sdX * sdY) + C2 ) / (sdX^2 + sdY^2 + C2)
  #correlacion
  s.xy = (vXY + C3) / ((sdX * sdY) + C3)
  #calculamos la SSIM
  z = l.xy * c.xy * s.xy #asumimos los pesos a, b, c todos iguales a 1
  return(z)
}

#' TEMPORAL CORRELATION
#' 
#' Computes the Temporal Correlation of Chouakria and Nagabhushan between two same length time series S1 and S2.
#' @param s1 a numeric vector containing the first time series.
#' @param s2 a numeric vector containing the second time series.
#' @param k a equal to 3.1 by default.
#' @return CORT temporal correlation between S1 and S2. 
#' @references  for more details read  Chouakria Douzal 2003.
#' @export
CORT <- function(s1, s2){
  p = length(s1)
  suma = 0 #we initialize the sum for the euclidean distance
  suma1 = 0 #we initialize the s1 sum for the denominator for the CORT
  suma2 = 0 #we initialize the s2 sum for the denominator for the CORT
  numerador = 0 #we initialize the numerator sum for the CORT
  for (i in 1:p) {
    if (i < p) {
      numerador = numerador + ((s1[i+1] - s1[i]) * (s2[i+1] - s2[i])) #calculate the CORT numerator sum
      suma1 = suma1 + (s1[i+1] - s1[i])^2 #calculate both CORT denomiator sums
      suma2 = suma2 + (s2[i+1] - s2[i])^2
    }
  }
  return(numerador / (sqrt(suma1) * sqrt(suma2)))
}

#' TEMPORAL CORRELATION MEASURE
#' 
#' Computes the Temporal Correlation of Chouakria and Nagabhushan between two same length time series s1 and s2, and return the result of 1-abs(CORT).
#' @param s1 a numeric vector containing the first time series.
#' @param s2 a numeric vector containing the second time series.
#' @param k a equal to 3.1 by default.
#' @return a positive real number between 0 and 1. Where 0 is the maximum similarity and 1 the maximum dissimilarity.
#' @references  for more details read  Chouakria Douzal 2003.
#' @export
CORT.Measure <- function(s1, s2){1-abs(CORT(s1, s2))}


#' LONGEST COMMON SUBSEQUENCE DISTANCE MEASURE
#' 
#' Computes the 'Longest Common Sub sequence distance' between a pair of numeric time series, and return difference between the length of s1 and the LCSS distance.
#' @param s1 a numeric vector containing the first time series.
#' @param s2 a numeric vector containing the second time series.
#' @param e  a positive threshold value that defines the distance.
#' @return a positive real number between 0 and 1. Where 0 is the maximum similarity and 1 the maximum dissimilarity.
#' @references  for more details read the ''TSdist' package documentation
#' @import TSdist
#' @export
LCSS.Measure <- function(s1, s2, e=0.01){length(s1) - LCSSDistance(s1, s2, e)}


#' CHEBYSHEV MEASURE
#' 
#' Computes the Chebyshev distance using the 'chebyshev' function from the 'philentropy' package by setting the testNA parameter as False.
#' @param s1 a numeric vector containing the first time series.
#' @param s2 a numeric vector containing the second time series.
#' @return a positive real number. The computed distance between the pair of series.
#' @references for more details read the 'philentropy' package documentation.
#' @import philentropy
#' @export
Chebyshev.Measure <- function(s1, s2) {return(chebyshev(s1, s2, testNA = F))} 


#' EDITH DISTANCE WITH REAL PENALTY 
#' Computes the Edit Distance with Real Penalty between a pair of numeric time series using the 'ERPDistance' function from the 'TSdist' package.
#' We by setting the 'sigma' parameter as the 10% of the time series length, and the 'g' parameter equal to 0.1.
#' @param s1 a numeric vector containing the first time series.
#' @param s2 a numeric vector containing the second time series.
#' @param g the reference value used to penalize gaps.
#' @param sima a Sakoe-Chiba windowing constraint can be added by specifying a positive integer representing the window size.
#' @return the computed distance between the pair of series.
#' @references for more details read the 'TSdist' package documentation. 
#' @import TSdist
#' @export
ERP.Measure <- function(s1, s2) {
  win <- ceiling(length(s1)/.1)
  return(ERPDistance(s1, s2, g=.1, sigma=win))}


#' EDITH DISTANCE FOR REAL SEQUENCE
#' Computes the Edit Distance for Real Sequence between a pair of numeric time series using the 'EDRDistance' function from the 'TSdist' package.
#' We by setting the 'sigma' parameter as the 10% of the time series length, and the 'epsilon' parameter equal to 0.1.
#' @param s1 a numeric vector containing the first time series.
#' @param s2 a numeric vector containing the second time series.
#' @param epsilon a positive threshold value that defines the distance.
#' @param sima a Sakoe-Chiba windowing constraint can be added by specifying a positive integer representing the window size.
#' @return the computed distance between the pair of series.
#' @references for more details read the 'TSdist' package documentation. 
#' @import TSdist
#' @export
EDR.Measure <- function(s1, s2) {
  win <- ceiling(length(s1)/.1)
  return(EDRDistance(s1, s2, epsilon=.1, sigma=win))}


#' DYNAMIC TIME WARPING
#' Computes the Dynamic Time Warping distance between a pair of numeric time series using the 'DTWDistance' function from the 'TSdist' package.
#' We by setting the 'DTWDistance' function by using a  Sakoe-Chiba windowing constraint with a window equal to the 10% of the time series length.
#' @param s1 a numeric vector containing the first time series.
#' @param s2 a numeric vector containing the second time series.
#' @return the computed distance between the pair of series.
#' @references for more details read the 'TSdist' package documentation. 
#' @import TSdist
#' @export 
DTW.Measure <- function(s1, s2){
  win <- ceiling(length(s1)/.1)
  return(DTWDistance(s1, s2, window.type="sakoechiba", window.size=win))
}