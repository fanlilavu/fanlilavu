% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lavu.R
\name{lavu}
\alias{lavu}
\title{Teststatistik Lavergne und Vuong}
\usage{
lavu(Y, X1, X2, type = "epanechnikov", const = 1,
  bandwidth = "rule_of_thumb", simplify.var.est = T)
}
\arguments{
\item{Y}{nx1 Datastream, Liste oder Matrix, welche die abhaengigen Variablen enthaelt.}

\item{X1}{nxq1 Datastream oder Matrix, welche die q1 Variablen enthaelt, welche in der Nullhypothese als signifikant angesehen werden.}

\item{X2}{nxq2 Datastream oder Matrix, welche die q2 Variablen enthaelt, welche in der Nullhypothese als insignifikant angesehen werden.}

\item{type}{Der Kern, welcher fuer die Kernfunktionen verwendendet werden soll. Es stehen der Gauss-Kern (\emph{="gaussian"}), der Epanechnikov-Kern (\emph{="epanechnikov"}) sowie der Uniform-Kern (\emph{="uniform"}) zur Verfuegung.}

\item{const}{Eine Konstanten mit der die Bandweiten multipliziert werden. Falls nicht weiter spezifiziert, werden die Bandbreiten nicht veraendert. Der Sinn dieser Konstanten ist die Untersuchung der Auswirkung der gew?hlten Bandbreite.}

\item{bandwidth}{Gibt an, ob die Bandbreite automatisch mithilfe des npregbw-Befehls berechnet werden soll (\emph{automatic}). Die Standardeinstellung ist die Verwendung einer Daumenregel \emph{"rule_of_thumb"}.}

\item{simplify.var.est}{Wird die boolsche Variable auf TRUE gesetzt, wird die von Lavergne und Vuong (2000: 581) angegebene Vereinfachung f?r die Varianz genutzt. Dies ist aufgrund der extremen Laufzeitverl?ngerung f?r FALSE sehr zu empfehlen.}
}
\value{
Gibt den Wert der Teststatistik nach Lavergne und Vuong zurueck.
}
\description{
Die Funktion lavu() berechnet die Teststatistik, welche Lavergne und Vuong (2000) in ihrer Arbeit
"Nonparametric Significance Testing" vorgestellt haben. Lavergne, P. & Vunong Q.: Nonparametric
Significance Testing. (\emph{Econometric Theory}), 16(4), 576-601
}
\author{
Ingo Lange (\emph{ilange@mail.uni-mannheim.de})
}

