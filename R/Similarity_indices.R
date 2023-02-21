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
      numerador = numerador + ((S1[i+1] - S1[i]) * (S2[i+1] - S2[i])) #calculate the CORT numerator sum
      suma1 = suma1 + (S1[i+1] - S1[i])^2 #calculate both CORT denomiator sums
      suma2 = suma2 + (S2[i+1] - S2[i])^2
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
#' Computes the Longest Common Sub sequence distance between a pair of numeric time series, and return difference between the length of s1 and the LCSS distance.
#' @param s1 a numeric vector containing the first time series.
#' @param s2 a numeric vector containing the second time series.
#' @param e  a positive threshold value that defines the distance.
#' @return a positive real number between 0 and 1. Where 0 is the maximum similarity and 1 the maximum dissimilarity.
#' @references  for more details read  TSdist documentation
#' @import TSdist
#' @export
LCSS.Measure <- function(s1, s2, e=0.01){length(s1) - LCSSDistance(s1, s2, e)}