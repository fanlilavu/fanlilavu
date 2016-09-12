#' Uniform Kern
#'
#' Die Funktion uni() berechnet die Uniform-Kernfunktion.
#'
#' @author Ingo Lange (\emph{ilange@mail.uni-mannheim.de})
#' @param x Reelle Zahl, fuer welche die Uniform-Funktion berechnet wird.
#' @export



uni <- function(x){
  0.5*(abs(x)<=1)
}
