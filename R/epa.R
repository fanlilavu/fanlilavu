#' Epanechnicov Kern
#'
#' Die Funktion epa() berechnet die Epanechnikov-Kernfunktion.
#'
#' @author Ingo Lange (\emph{ilange@mail.uni-mannheim.de})
#' @param x Reelle Zahl, fuer welche die Epanechnikov-Funktion berechnet wird.
#' @export


epa <- function(x){
  (1-x^2)*(abs(x)<1)
}
