#' EUCLIDEAN DISTANCE WITH NA.OMIT 
#' 
#' Compute point-to-point euclidean distance between two time series x and y.
#' @param x a time series with real values.
#' @param y a time series with real values and not necessarily the same length as x.
#' @return a positive real number result of the sum of point-to-point euclidean distance between x and y.
#' @import philentropy
#' @export
EUC.NA.OMIT <- function(x,y){
  return(euclidean(as.numeric(na.omit(x)),as.numeric(na.omit(y)),testNA = FALSE))
}

#' DYNAMIC TIME WARPING WITH NA.OMIT 
#' 
#' Compute dynamic time warping distance between two time series x and y using the "dtw" library.
#' @param x a time series with real values.
#' @param y a time series with real values and not necessarily the same length as x.
#' @return a positive real number result of dynamic time warping distance between x and y using dtw(x,y)$distance function.
#' @import dtw
#' @export
DTW.NA.OMIT <- function (x,y){
  return(dtw(as.numeric(na.omit(x)),as.numeric(na.omit(y)))$distance)
}

#' SSIMT INDEX
#' 
#' Compute SSIMT similarity index between two same length time series X and Y.
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
  s.xy = (vXY + C3) / ((vX * vY) + C3)
  #calculamos la SSIM
  z = l.xy * c.xy * s.xy #asumimos los pesos a, b, c todos iguales a 1
  return(z)
}

#' SSIMT INDEX WITH NA.OMIT
#' 
#' Compute SSIMT similarity index between two same length time series X and Y.
#' @param X a time series of length n with real values.
#' @param Y a time series of length n with real values.
#' @return a real number between 1 and -1; where 1 means max similarity, -1 means opposite behavior and 0 no similarity. 
#' @export
SSIMT.NA.OMIT <- function(x,y){
  return(SSIMT(as.numeric(na.omit(x)),as.numeric(na.omit(y))))
}

#' CHOUAKRIA INDEX
#' 
#' Compute D index of Chouakria and Nagabhushan similarity index between two same length time series S1 and S2.
#' @param S1 a time series of length n with real values.
#' @param S2 a time series of length n with real values.
#' @param k a equal to 3.1 by default.
#' @return D a positive real number plus 0; where 0 means max similarity.
#' @return CORT temporal correlation between S1 and S2. 
#' @return dE point-to-point euclidean distance between S1 and S2.
#' @references  for more details read Chouakria and Nagabhushan 2007.
#' @export
CHOU <- function(x, y, FUN, k=3.1){
  S1 = as.numeric(x)
  S2 = as.numeric(y)
  p = length(S1)
  suma = 0 #inicializamos la suma para la distancia euclidea
  suma1 = 0 #inicializamos la sumatoria S1 del denominador para el CORT
  suma2 = 0 #inicializamos la sumatoria S2 del denominador para el CORT
  numerador = 0 #inicializamos la sumatoria del numerador para el CORT
  for (i in 1:p) {
    if (i < p) {
      numerador = numerador + ((S1[i+1] - S1[i]) * (S2[i+1] - S2[i])) #calculamos la sumatoria del numerador del CORT
      suma1 = suma1 + (S1[i+1] - S1[i])^2 #calculamos las dos sumatorias del denominador del CORT
      suma2 = suma2 + (S2[i+1] - S2[i])^2
    }
  }
  if (missing(FUN)) {
    dE = euclidean(S1,S2,testNA = FALSE)
  } else {
    dE = FUN(S1,S2)
  }
  CORT = numerador / (sqrt(suma1) * sqrt(suma2))
  #finalmente calculamos la dissimilaridad
  D = (2 / (1 + exp(k * CORT))) * dE
  y = as.numeric(D)
  return(y)
}

#' CHOUAKRIA INDEX WITH NA.OMIT
#' 
#' Compute D index of Chouakria and Nagabhushan similarity index between two same length time series S1 and S2.
#' @param S1 a time series of length n with real values.
#' @param S2 a time series of length n with real values.
#' @param k a equal to 3.1 by default.
#' @return D a positive real number plus 0; where 0 means max similarity.
#' @return CORT temporal correlation between S1 and S2. 
#' @return dE point-to-point euclidean distance between S1 and S2.
#' @references  for more details read Chouakria and Nagabhushan 2007.
#' @export
CHOU.NA.OMIT <- function(x,y, K=3.1){
  return(CHOU(as.numeric(na.omit(x)),as.numeric(na.omit(y)),k=K)[1])
}

#' TEMPORAL CORRELATION
#' 
#' Compute the Temporal Correlation of Chouakria and Nagabhushan between two same length time series S1 and S2.
#' @param S1 a time series of length n with real values.
#' @param S2 a time series of length n with real values.
#' @param k a equal to 3.1 by default.
#' @return CORT temporal correlation between S1 and S2. 
#' @references  for more details read  Chouakria Douzal 2003.
#' @export
CORT <- function(x, y){
  S1 = as.numeric(x)
  S2 = as.numeric(y)
  p = length(S1)
  suma = 0 #inicializamos la suma para la distancia euclidea
  suma1 = 0 #inicializamos la sumatoria S1 del denominador para el CORT
  suma2 = 0 #inicializamos la sumatoria S2 del denominador para el CORT
  numerador = 0 #inicializamos la sumatoria del numerador para el CORT
  for (i in 1:p) {
    if (i < p) {
      numerador = numerador + ((S1[i+1] - S1[i]) * (S2[i+1] - S2[i])) #calculamos la sumatoria del numerador del CORT
      suma1 = suma1 + (S1[i+1] - S1[i])^2 #calculamos las dos sumatorias del denominador del CORT
      suma2 = suma2 + (S2[i+1] - S2[i])^2
    }
  }
  return(numerador / (sqrt(suma1) * sqrt(suma2)))
}