#' Teststatistik Lavergne und Vuong
#'
#' Die Funktion lavu() berechnet die Teststatistik, welche Lavergne und Vuong (2000) in ihrer Arbeit
#' "Nonparametric Significance Testing" vorgestellt haben. Lavergne, P. & Vunong Q.: Nonparametric
#' Significance Testing. (\emph{Econometric Theory}), 16(4), 576-601
#'
#' @author Ingo Lange (\emph{ilange@mail.uni-mannheim.de})
#' @param Y nx1 Datastream, Liste oder Matrix, welche die abhaengigen Variablen enthaelt.
#' @param X1 nxq1 Datastream oder Matrix, welche die q1 Variablen enthaelt, welche in der Nullhypothese als signifikant angesehen werden.
#' @param X2 nxq2 Datastream oder Matrix, welche die q2 Variablen enthaelt, welche in der Nullhypothese als insignifikant angesehen werden.
#' @param type Der Kern, welcher fuer die Kernfunktionen verwendendet werden soll. Es stehen der Gauss-Kern (\emph{="gaussian"}), der Epanechnikov-Kern (\emph{="epanechnikov"}) sowie der Uniform-Kern (\emph{="uniform"}) zur Verfuegung.
#' @param const Eine Konstanten mit der die Bandweiten multipliziert werden. Falls nicht weiter spezifiziert, werden die Bandbreiten nicht veraendert. Der Sinn dieser Konstanten ist die Untersuchung der Auswirkung der gew?hlten Bandbreite.
#' @param bandwidth Gibt an, ob die Bandbreite automatisch mithilfe des npregbw-Befehls berechnet werden soll (\emph{automatic}). Die Standardeinstellung ist die Verwendung einer Daumenregel \emph{"rule_of_thumb"}.
#' @param simplify.var.est Wird die boolsche Variable auf TRUE gesetzt, wird die von Lavergne und Vuong (2000: 581) angegebene Vereinfachung f?r die Varianz genutzt. Dies ist aufgrund der extremen Laufzeitverl?ngerung f?r FALSE sehr zu empfehlen.
#' @return Gibt den Wert der Teststatistik nach Lavergne und Vuong zurueck.
#' @export


lavu <- function(Y,X1,X2,type="epanechnikov", const=1, bandwidth="rule_of_thumb",simplify.var.est=T){

  X1m <- as.matrix(X1)  # Umwandlung der ggf. als Dataframe eingelesenen Daten
  X2m <- as.matrix(X2)  # in Matrizen
  Ym  <- as.matrix(Y)

  q1 <- ncol(X1m)  # Festlegen der Dimension q1 von W
  q2 <- ncol(X2m)  # Festlegen der Dimension q2 von Z
  n <- nrow(X1m)   # Anzahl der Daten n

  if(nrow(X2m)!=n||(nrow(Ym)!=n)){
    stop("number of rows are not equal")}   # Warnmeldung und Abbruch, falls die Zeilen
  # von W,Z und Y nicht übereinstimmen

  X  <- cbind(X1,X2)        # X=(W,Z)
  Xm  <- as.matrix(X)       # Matrixumwandlung

  # Je nachdem, welche Einstellung gew?hlt wurde, wird eine Bandbreite nach Daumenregel oder mit npregbw
  bw <- NULL
  if(bandwidth=="automatic"){       # Automatische Bandbreitenwahl f?r den type mit npregbw
    bw <- npregbw(ydat=Y,xdat=X, ckertype = type)
    bw <- const*bw$bw               # Exraktion des Bandbreitenvektors
  }else{
    if(bandwidth=="rule_of_thumb"){   # Daumenregel
      if(type=="gaussian"){             # (i) f?r Gau?-Kern
        for(i in 1:(q1+q2)){
          bw[i] <- const*(4/3)^(1/5)*sd(Xm[,i])*n^(-1/5)
        }
      }else{
        if(type=="epanechnikov"){       # (i) f?r den Epanechnikov-Kern
          for(i in 1:(q1+q2)){
            bw[i] <- const*sd(Xm[,i])*n^(-1/5)
          }
        }
      }
    }
  }

  if(type=="epanechnikov") { # Berechnung f?r den Epanechnikov-Kern
    V_n <- 0   # Der Wert Vn wird zu Beginn auf Null gesetzt und die Summe wird dann in
    # der Schleife ?ber (i,j,k,l) in (I,J,K,L) bebildet

    I <- (1:n)          # I ist 1:n
    for (i in I){
      J <- I[!(I%in%i)]   # J ist 1:n au?er i
      for(j in J){
        K <- J[!(J%in%j)]    # K ist 1:n au?er i,j
        for (k in K){
          L <- K[!(K%in%k)]    # L ist 1:n au?er i,j,l
          for(l in L){
            V_n <- V_n + (Ym[i,]-Ym[k,])*(Ym[j,]-Ym[l,])*prod(epa((X1m[i,]-X1m[k,])%*%diag(1/bw[1:q1],q1,q1)))*prod(epa((X1m[j,]-X1m[l,])%*%diag(1/bw[1:q1],q1,q1)))*prod(epa((Xm[i,]-Xm[j,])%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))
          }
        }
      }
    }
    V_n <- V_n/(factorial(n)/factorial(n-4)*prod(bw[1:q1])^2*prod(bw[1:(q1+q2)]))


    if(simplify.var.est==T){   # Berechnung der vereinfachten Varianz
      fhat_w  <- NULL # Definition der leeren Vektoren f?r: - Dichtesch?tzer
      Yhat    <- NULL #                                     - Sch?tzer f?r y
      u_hat   <- NULL #                                     - Sch?tzer f?r u
      # Berechnung der Sch?tzer
      for (i in 1:n){
        K_ij_w      <- rowProds(epa((matrix(rep(X1m[i,],n),nrow=n,byrow = T)-X1m)[(1:n)[!((1:n)%in%i)],]%*%diag(1/bw[1:q1],q1,q1)))
        fhat_w[i]   <- sum(K_ij_w)/((n-1)*prod(bw[1:q1]))
        Yhat[i]     <- ifelse(fhat_w[i]!=0,sum(K_ij_w*Ym[(1:n)[!((1:n)%in%i)],])/((n-1)*prod(bw[1:q1])*fhat_w[i]),0)
        u_hat[i]    <- Y[i]-Yhat[i]
      }

      sig2 <- 0 # Die gesch?tzte Varianz wird auf Null gesetzt und dann in den Schleifen berechnet
      I <- (1:n) # I enth?lt 1:n
      for (i in I){
        J <- I[!(I%in%i)] # J enth?lt 1:n au?er i
        for(j in J){
          sig2 <- sig2 + u_hat[i]^2*u_hat[j]^2*fhat_w[i]^2*fhat_w[j]^2*prod(epa((Xm[i,]-Xm[j,])%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))^2
        }
      }
      sig2 <- 2*sig2/(n*(n-1)*prod(bw[1:(q1+q2)]))        # Nach der Bildung der Summe werden die Operationen vorgenommen, die aus der Summe gezogen werden d?rfen
    }else{
      if(simplify.var.est==F){   # Berechnung der komplexeren Varianz
        sig2 <- 0   # Der Wert sig2 wird zu Beginn auf Null gesetzt und die Summe wird dann in
        # der Schleife ?ber (i,j,k,l) in (I,J,K,L) bebildet

        I <- (1:n)          # I ist 1:n
        for (i in I){
          J <- I[!(I%in%i)]   # J ist 1:n au?er i
          for(j in J){
            K <- J[!(J%in%j)]    # K ist 1:n au?er i,j
            for (k in K){
              L <- K[!(K%in%k)]    # L ist 1:n au?er i,j,l
              for(l in L){
                O <- L[!(L%in%l)]    # L ist 1:n au?er i,j,l
                for(o in O){
                  P <- O[!(O%in%o)]    # L ist 1:n au?er i,j,l
                  for(p in P){
                    sig2 <- sig2 + (Ym[i,]-Ym[k,])*(Ym[i,]-Ym[o,])*(Ym[j,]-Ym[l,])*(Ym[j,]-Ym[p,])*
                      prod(epa((X1m[i,]-X1m[k,])%*%diag(1/bw[1:q1],q1,q1)))*
                      prod(epa((X1m[i,]-X1m[o,])%*%diag(1/bw[1:q1],q1,q1)))*
                      prod(epa((X1m[j,]-X1m[l,])%*%diag(1/bw[1:q1],q1,q1)))*
                      prod(epa((X1m[j,]-X1m[p,])%*%diag(1/bw[1:q1],q1,q1)))*
                      prod(epa((Xm[i,]-Xm[j,])%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))^2
                  }
                }
              }
            }
          }
        }
        sig2 <- 2*sig2/(factorial(n)/factorial(n-6)*prod(bw[1:q1])*prod(bw[1:(q1+q2)]))
      }else{
        stop("simplify.var.est is not TRUE or FALSE")       #Fehlermeldung bei falscher Eingabe
      }
    }
  }
  if(type=="gaussian") { # Berechnung f?r den Gau?-Kern
    V_n <- 0   # Der Wert Vn wird zu Beginn auf Null gesetzt und die Summe wird dann in
    # der Schleife ?ber (i,j,k,l) in (I,J,K,L) bebildet

    I <- (1:n)          # I ist 1:n
    for (i in I){
      J <- I[!(I%in%i)]   # J ist 1:n au?er i
      for(j in J){
        K <- J[!(J%in%j)]    # K ist 1:n au?er i,j
        for (k in K){
          L <- K[!(K%in%k)]    # L ist 1:n au?er i,j,l
          for(l in L){
            V_n <- V_n + (Ym[i,]-Ym[k,])*(Ym[j,]-Ym[l,])*prod(dnorm((X1m[i,]-X1m[k,])%*%diag(1/bw[1:q1],q1,q1)))*prod(dnorm((X1m[j,]-X1m[l,])%*%diag(1/bw[1:q1],q1,q1)))*prod(dnorm((Xm[i,]-Xm[j,])%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))
          }
        }
      }
    }
    V_n <- V_n/(factorial(n)/factorial(n-4)*prod(bw[1:q1])^2*prod(bw[1:(q1+q2)]))


    if(simplify.var.est==T){   # Berechnung der vereinfachten Varianz
      fhat_w  <- NULL # Definition der leeren Vektoren f?r: - Dichtesch?tzer
      Yhat    <- NULL #                                     - Sch?tzer f?r y
      u_hat   <- NULL #                                     - Sch?tzer f?r u
      # Berechnung der Sch?tzer
      for (i in 1:n){
        K_ij_w      <- rowProds(dnorm((matrix(rep(X1m[i,],n),nrow=n,byrow = T)-X1m)[(1:n)[!((1:n)%in%i)],]%*%diag(1/bw[1:q1],q1,q1)))
        fhat_w[i]   <- sum(K_ij_w)/((n-1)*prod(bw[1:q1]))
        Yhat[i]     <- ifelse(fhat_w[i]!=0,sum(K_ij_w*Ym[(1:n)[!((1:n)%in%i)],])/((n-1)*prod(bw[1:q1])*fhat_w[i]),0)
        u_hat[i]    <- Y[i]-Yhat[i]
      }

      sig2 <- 0 # Die gesch?tzte Varianz wird auf Null gesetzt und dann in den Schleifen berechnet
      I <- (1:n) # I enth?lt 1:n
      for (i in I){
        J <- I[!(I%in%i)] # J enth?lt 1:n au?er i
        for(j in J){
          sig2 <- sig2 + u_hat[i]^2*u_hat[j]^2*fhat_w[i]^2*fhat_w[j]^2*prod(dnorm((Xm[i,]-Xm[j,])%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))^2
        }
      }
      sig2 <- 2*sig2/(n*(n-1)*prod(bw[1:(q1+q2)]))        # Nach der Bildung der Summe werden die Operationen vorgenommen, die aus der Summe gezogen werden d?rfen
    }else{
      if(simplify.var.est==F){   # Berechnung der komplexeren Varianz
        sig2 <- 0   # Der Wert sig2 wird zu Beginn auf Null gesetzt und die Summe wird dann in
        # der Schleife ?ber (i,j,k,l) in (I,J,K,L) bebildet

        I <- (1:n)          # I ist 1:n
        for (i in I){
          J <- I[!(I%in%i)]   # J ist 1:n au?er i
          for(j in J){
            K <- J[!(J%in%j)]    # K ist 1:n au?er i,j
            for (k in K){
              L <- K[!(K%in%k)]    # L ist 1:n au?er i,j,l
              for(l in L){
                O <- L[!(L%in%l)]    # L ist 1:n au?er i,j,l
                for(o in O){
                  P <- O[!(O%in%o)]    # L ist 1:n au?er i,j,l
                  for(p in P){
                    sig2 <- sig2 + (Ym[i,]-Ym[k,])*(Ym[i,]-Ym[o,])*(Ym[j,]-Ym[l,])*(Ym[j,]-Ym[p,])*
                      prod(dnorm((X1m[i,]-X1m[k,])%*%diag(1/bw[1:q1],q1,q1)))*
                      prod(dnorm((X1m[i,]-X1m[o,])%*%diag(1/bw[1:q1],q1,q1)))*
                      prod(dnorm((X1m[j,]-X1m[l,])%*%diag(1/bw[1:q1],q1,q1)))*
                      prod(dnorm((X1m[j,]-X1m[p,])%*%diag(1/bw[1:q1],q1,q1)))*
                      prod(dnorm((Xm[i,]-Xm[j,])%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))^2
                  }
                }
              }
            }
          }
        }
        sig2 <- 2*sig2/(factorial(n)/factorial(n-6)*prod(bw[1:q1])*prod(bw[1:(q1+q2)]))
      }else{
        stop("simplify.var.est is not TRUE or FALSE")       #Fehlermeldung bei falscher Eingabe
      }
    }
  }
  if(type=="uniform") {# Berechnung f?r den Uniform-Kern
    V_n <- 0   # Der Wert Vn wird zu Beginn auf Null gesetzt und die Summe wird dann in
    # der Schleife ?ber (i,j,k,l) in (I,J,K,L) bebildet

    I <- (1:n)          # I ist 1:n
    for (i in I){
      J <- I[!(I%in%i)]   # J ist 1:n au?er i
      for(j in J){
        K <- J[!(J%in%j)]    # K ist 1:n au?er i,j
        for (k in K){
          L <- K[!(K%in%k)]    # L ist 1:n au?er i,j,l
          for(l in L){
            V_n <- V_n + (Ym[i,]-Ym[k,])*(Ym[j,]-Ym[l,])*prod(uni((X1m[i,]-X1m[k,])%*%diag(1/bw[1:q1],q1,q1)))*prod(uni((X1m[j,]-X1m[l,])%*%diag(1/bw[1:q1],q1,q1)))*prod(uni((Xm[i,]-Xm[j,])%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))
          }
        }
      }
    }
    V_n <- V_n/(factorial(n)/factorial(n-4)*prod(bw[1:q1])^2*prod(bw[1:(q1+q2)]))


    if(simplify.var.est==T){   # Berechnung der vereinfachten Varianz
      fhat_w  <- NULL # Definition der leeren Vektoren f?r: - Dichtesch?tzer
      Yhat    <- NULL #                                     - Sch?tzer f?r y
      u_hat   <- NULL #                                     - Sch?tzer f?r u
      # Berechnung der Sch?tzer
      for (i in 1:n){
        K_ij_w      <- rowProds(uni((matrix(rep(X1m[i,],n),nrow=n,byrow = T)-X1m)[(1:n)[!((1:n)%in%i)],]%*%diag(1/bw[1:q1],q1,q1)))
        fhat_w[i]   <- sum(K_ij_w)/((n-1)*prod(bw[1:q1]))
        Yhat[i]     <- ifelse(fhat_w[i]!=0,sum(K_ij_w*Ym[(1:n)[!((1:n)%in%i)],])/((n-1)*prod(bw[1:q1])*fhat_w[i]),0)
        u_hat[i]    <- Y[i]-Yhat[i]
      }

      sig2 <- 0 # Die gesch?tzte Varianz wird auf Null gesetzt und dann in den Schleifen berechnet
      I <- (1:n) # I enth?lt 1:n
      for (i in I){
        J <- I[!(I%in%i)] # J enth?lt 1:n au?er i
        for(j in J){
          sig2 <- sig2 + u_hat[i]^2*u_hat[j]^2*fhat_w[i]^2*fhat_w[j]^2*prod(uni((Xm[i,]-Xm[j,])%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))^2
        }
      }
      sig2 <- 2*sig2/(n*(n-1)*prod(bw[1:(q1+q2)]))        # Nach der Bildung der Summe werden die Operationen vorgenommen, die aus der Summe gezogen werden d?rfen
    }else{
      if(simplify.var.est==F){   # Berechnung der komplexeren Varianz
        sig2 <- 0   # Der Wert sig2 wird zu Beginn auf Null gesetzt und die Summe wird dann in
        # der Schleife ?ber (i,j,k,l) in (I,J,K,L) bebildet

        I <- (1:n)          # I ist 1:n
        for (i in I){
          J <- I[!(I%in%i)]   # J ist 1:n au?er i
          for(j in J){
            K <- J[!(J%in%j)]    # K ist 1:n au?er i,j
            for (k in K){
              L <- K[!(K%in%k)]    # L ist 1:n au?er i,j,l
              for(l in L){
                O <- L[!(L%in%l)]    # L ist 1:n au?er i,j,l
                for(o in O){
                  P <- O[!(O%in%o)]    # L ist 1:n au?er i,j,l
                  for(p in P){
                    sig2 <- sig2 + (Ym[i,]-Ym[k,])*(Ym[i,]-Ym[o,])*(Ym[j,]-Ym[l,])*(Ym[j,]-Ym[p,])*
                      prod(uni((X1m[i,]-X1m[k,])%*%diag(1/bw[1:q1],q1,q1)))*
                      prod(uni((X1m[i,]-X1m[o,])%*%diag(1/bw[1:q1],q1,q1)))*
                      prod(uni((X1m[j,]-X1m[l,])%*%diag(1/bw[1:q1],q1,q1)))*
                      prod(uni((X1m[j,]-X1m[p,])%*%diag(1/bw[1:q1],q1,q1)))*
                      prod(uni((Xm[i,]-Xm[j,])%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))^2
                  }
                }
              }
            }
          }
        }
        sig2 <- 2*sig2/(factorial(n)/factorial(n-6)*prod(bw[1:q1])*prod(bw[1:(q1+q2)]))
      }else{
        stop("simplify.var.est is not TRUE or FALSE")       #Fehlermeldung bei falscher Eingabe
      }
    }
  }

  T_V  <- n*sqrt(prod(bw[1:(q1+q2)]))*V_n/sqrt(sig2)  # Berechnung der Teststatistik
  return(T_V)                                         # Ausgabe der Teststatistik
}
