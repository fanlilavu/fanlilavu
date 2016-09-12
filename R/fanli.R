#' Teststatistik Fan und Li
#'
#' Die Funktion fanli() berechnet die Teststatistik, welche Fan und Li (1996) in ihrer Arbeit
#' "Consistent Model Specification Tests: Omitted Variables and Semiparamtric Functional Forms"
#' vorgestellt haben: Fan Y. & Li, Q (1996): Consistent Model Specification Tests: Omitted Variables and Semiparamtric
#'  Functional Forms. \emph{Econometrica}, 64 (4), 868-890
#'
#' @author Ingo Lange (\emph{ilange@mail.uni-mannheim.de})
#' @param Y nx1 Datastream, Liste oder Matrix, welche die abhaengigen Variablen enthaelt.
#' @param X1 nxq1 Datastream oder Matrix, welche die q1 Variablen enthaelt, welche in der Nullhypothese als signifikant angesehen werden.
#' @param X2 nxq2 Datastream oder Matrix, welche die q2 Variablen enthaelt, welche in der Nullhypothese als insignifikant angesehen werden.
#' @param type Der fuer die Kernfunktionen zu verwendende Kern. Es stehen neben dem Gauss-Kern (\emph{="gaussian"}), der Epanechnikov-Kern (\emph{="epanechnikov"}) und der Uniform-Kern (\emph{="uniform"}) zur Verfuegung.
#' @param const Konstanten mit der die Bandweiten multipliziert werden sollen. Falls nicht weiter spezifiziert, werden die Bandbreiten nicht veraendert. Der Sinn dieser Konstanten ist die Untersuchung der Sensitivitaet der Bandbreitenwahl.
#' @param bandwidth Gibt an, ob die Bandbreite automatisch berechnet werden soll (\emph{"automatic"}). Fuer individuelle Bandbreitenwahl waehle \emph{"manually"}.
#' @param return.bw Bei der Standardeinstellung (\emph{FALSE}) wird lediglich der Wert der Teststatistik zurück gegeben.
#' @param a Gibt die Bandbreite a an, welche fuer die auf der Variablenmenge X1 basierenden Kernfunktion verwendet werden soll. Wird nur fuer \emph{bandwidth="manually"} verwendet.
#' @param h Gibt die Bandbreite h an, welche fuer die auf der Gesamtvariablenmenge basierende Kernfunktion gewaehlt werden soll. Wird nur fuer \emph{bandwidth="manually"} verwendet.
#' @return Gibt den Wert der Teststatistik zurueck.
#' @export


fanli <- function(Y,X1,X2,type="gaussian", const=1, bandwidth="automatic",return.bw=F,a=NA,h=NA){
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

  if(bandwidth=="automatic"|bandwidth=="rule_of_thumb"){

  bw <- NULL
  if(bandwidth=="automatic"){
    bw <- npregbw(ydat=Y,xdat=X, ckertype = type)
    bw <- const*bw$bw               # Exraktion des Bandbreitenvektors
  }else{
    if(bandwidth=="rule_of_thumb"){
      if(type=="gaussian"){
        for(i in 1:(q1+q2)){
          bw[i] <- const*(4/3)^(1/5)*sd(Xm[,i])*n^(-1/5)
        }
      }else{
        if(type=="epanechnikov"){
          for(i in 1:(q1+q2)){
            bw[i] <- const*sd(Xm[,i])*n^(-1/5)
          }
        }
      }
    }
  }
  fhat_w  <- NULL # Definition der leeren Vektoren für: - Dichteschätzer
  Yhat    <- NULL #                                     - Schätzer für y
  u_hat   <- NULL #                                     - Schätzer für u
  if(type=="gaussian"){
    # f?r i=1,...,n werden fhat_w_i, Yhat_i, u_hat_i berechnet und in den leeren Vektoren gespeichert
    for (i in 1:n){
      K_ij_w      <- rowProds(dnorm((matrix(rep(X1m[i,],n),nrow=n,byrow = T)-X1m)[(1:n)[!((1:n)%in%i)],]%*%diag(1/bw[1:q1],q1,q1)))
      fhat_w[i]   <- sum(K_ij_w)/((n-1)*prod(bw[1:q1]))
      Yhat[i]     <- sum(K_ij_w*Ym[(1:n)[!((1:n)%in%i)],])/((n-1)*prod(bw[1:q1])*fhat_w[i])
      u_hat[i]    <- Y[i]-Yhat[i]
    }
    # Berechnung von I_n^a:
    Ina <- 0  # Der Wert von I_n^a wird auf Null gesetzt und die Summe ?ber alle n in einer Schleife gerechnet
    for (i in 1:n){
      K_ij      <- rowProds(dnorm((matrix(rep(Xm[i,],n),nrow=n,byrow = T)-Xm)[(1:n)[!((1:n)%in%i)],]%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))
      Ina <- Ina+sum((u_hat[i]*fhat_w[i])*(u_hat[(1:n)[!((1:n)%in%i)]]*fhat_w[(1:n)[!((1:n)%in%i)]])*K_ij)
    }
    Ina <- Ina/(n*(n-1)*prod(bw[1:(q1+q2)]))  # Die Summe wird durch (n*(n-1)*PROD[bw]) geteilt, wodurch der entg?ltige Wert berechnet ist

    # Berechnung des Sch?tzers f?r die Standartabweichung:
    sig2a <- 0 # Der Wert von I_n^a wird auf Null gesetzt und die Summe ?ber alle n in einer Schleife gerechnet
    for (i in 1:n){
      K_ij      <- rowProds(dnorm((matrix(rep(Xm[i,],n),nrow=n,byrow = T)-Xm)[(1:n)[!((1:n)%in%i)],]%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))
      sig2a <- sig2a+sum((u_hat[i]*fhat_w[i])^2*(u_hat[(1:n)[!((1:n)%in%i)]]*fhat_w[(1:n)[!((1:n)%in%i)]])^2*K_ij)
    }
    sig2a <- sig2a/(n*(n-1)*prod(bw[1:(q1+q2)])*(2*sqrt(pi))^(q1+q2))  # Die Summe wird durch (n*(n-1)*PROD[bw]) geteilt, wodurch der entg?ltige Wert berechnet ist


    Ta <- n*sqrt(prod(bw[1:(q1+q2)]))*Ina/(sqrt(2*sig2a))
  }
  else{
    if(type=="epanechnikov"){
      # fuer i=1,...,n werden fhat_w_i, Yhat_i, u_hat_i berechnet und in den leeren Vektoren gespeichert
      for (i in 1:n){
        K_ij_w      <- rowProds(epa((matrix(rep(X1m[i,],n),nrow=n,byrow = T)-X1m)[(1:n)[!((1:n)%in%i)],]%*%diag(1/bw[1:q1],q1,q1)))
        fhat_w[i]   <- sum(K_ij_w)/((n-1)*prod(bw[1:q1]))
        Yhat[i]     <- ifelse(fhat_w[i]!=0,sum(K_ij_w*Ym[(1:n)[!((1:n)%in%i)],])/((n-1)*prod(bw[1:q1])*fhat_w[i]),0)
        u_hat[i]    <- Y[i]-Yhat[i]
      }
      # Berechnung von I_n^a:
      Ina <- 0  # Der Wert von I_n^a wird auf Null gesetzt und die Summe ?ber alle n in einer Schleife gerechnet
      for (i in 1:n){
        K_ij      <- rowProds(epa((matrix(rep(Xm[i,],n),nrow=n,byrow = T)-Xm)[(1:n)[!((1:n)%in%i)],]%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))
        Ina <- Ina+sum((u_hat[i]*fhat_w[i])*(u_hat[(1:n)[!((1:n)%in%i)]]*fhat_w[(1:n)[!((1:n)%in%i)]])*K_ij)
      }
      Ina <- Ina/(n*(n-1)*prod(bw[1:(q1+q2)]))  # Die Summe wird durch (n*(n-1)*PROD[bw]) geteilt, wodurch der entg?ltige Wert berechnet ist

      # Berechnung des Schaetzers fuer die Standartabweichung:
      sig2a <- 0 # Der Wert von I_n^a wird auf Null gesetzt und die Summe ?ber alle n in einer Schleife gerechnet
      for (i in 1:n){
        K_ij      <- rowProds(epa((matrix(rep(Xm[i,],n),nrow=n,byrow = T)-Xm)[(1:n)[!((1:n)%in%i)],]%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))
        sig2a <- sig2a+sum((u_hat[i]*fhat_w[i])^2*(u_hat[(1:n)[!((1:n)%in%i)]]*fhat_w[(1:n)[!((1:n)%in%i)]])^2*K_ij)
      }
      sig2a <- (sig2a*(9/15)^(q1+q2))/(n*(n-1)*prod(bw[1:(q1+q2)]))  # Die Summe wird durch (n*(n-1)*PROD[bw]) geteilt, wodurch der entg?ltige Wert berechnet ist


      Ta <- n*sqrt(prod(bw[1:(q1+q2)]))*Ina/(sqrt(2*sig2a))
    }
    else{
      if(type=="uniform"){
        # f?r i=1,...,n werden fhat_w_i, Yhat_i, u_hat_i berechnet und in den leeren Vektoren gespeichert
        for (i in 1:n){
          K_ij_w      <- rowProds(uni((matrix(rep(X1m[i,],n),nrow=n,byrow = T)-X1m)[(1:n)[!((1:n)%in%i)],]%*%diag(1/bw[1:q1],q1,q1)))
          fhat_w[i]   <- sum(K_ij_w)/((n-1)*prod(bw[1:q1]))
          Yhat[i]     <- ifelse(fhat_w[i]!=0,sum(K_ij_w*Ym[(1:n)[!((1:n)%in%i)],])/((n-1)*prod(bw[1:q1])*fhat_w[i]),0)
          u_hat[i]    <- Y[i]-Yhat[i]
        }
        # Berechnung von I_n^a:
        Ina <- 0  # Der Wert von I_n^a wird auf Null gesetzt und die Summe ?ber alle n in einer Schleife gerechnet
        for (i in 1:n){
          K_ij      <- rowProds(uni((matrix(rep(Xm[i,],n),nrow=n,byrow = T)-Xm)[(1:n)[!((1:n)%in%i)],]%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))
          Ina <- Ina+sum((u_hat[i]*fhat_w[i])*(u_hat[(1:n)[!((1:n)%in%i)]]*fhat_w[(1:n)[!((1:n)%in%i)]])*K_ij)
        }
        Ina <- Ina/(n*(n-1)*prod(bw[1:(q1+q2)]))  # Die Summe wird durch (n*(n-1)*PROD[bw]) geteilt, wodurch der entg?ltige Wert berechnet ist

        # Berechnung des Sch?tzers f?r die Standartabweichung:
        sig2a <- 0 # Der Wert von I_n^a wird auf Null gesetzt und die Summe ?ber alle n in einer Schleife gerechnet
        for (i in 1:n){
          K_ij      <- rowProds(uni((matrix(rep(Xm[i,],n),nrow=n,byrow = T)-Xm)[(1:n)[!((1:n)%in%i)],]%*%diag(1/bw[1:(q1+q2)],(q1+q2),(q1+q2))))
          sig2a <- sig2a+sum((u_hat[i]*fhat_w[i])^2*(u_hat[(1:n)[!((1:n)%in%i)]]*fhat_w[(1:n)[!((1:n)%in%i)]])^2*K_ij)
        }
        sig2a <- (sig2a*0.5^(q1+q2))/(n*(n-1)*prod(bw[1:(q1+q2)]))  # Die Summe wird durch (n*(n-1)*PROD[bw]) geteilt, wodurch der entg?ltige Wert berechnet ist


        Ta <- n*sqrt(prod(bw[1:(q1+q2)]))*Ina/(sqrt(2*sig2a))
      }
      else{
        stop("kernel is not gaussian,epanechnikov or uniform")}   # Warnmeldung und Abbruch, falls die Zeilen
      # von W,Z und Y nicht übereinstimmen
    }
  }
  if(return.bw==F){
    return(Ta)}else{
      if(return.bw==T){
        return(c(Ta,bw))}else{
          stop("return.bw is not boolean")
        }
    }

  }


  if(bandwidth=="manually"){
    a <- rep(a,q1)
    h <- rep(h,(q1+q2))
    fhat_w  <- NULL # Definition der leeren Vektoren für: - Dichteschätzer
    Yhat    <- NULL #                                     - Schätzer für y
    u_hat   <- NULL #                                     - Schätzer für u
    if(type=="gaussian"){
      # f?r i=1,...,n werden fhat_w_i, Yhat_i, u_hat_i berechnet und in den leeren Vektoren gespeichert
      for (i in 1:n){
        K_ij_w      <- rowProds(dnorm((matrix(rep(X1m[i,],n),nrow=n,byrow = T)-X1m)[(1:n)[!((1:n)%in%i)],]%*%diag(1/a,q1,q1)))
        fhat_w[i]   <- sum(K_ij_w)/((n-1)*prod(a))
        Yhat[i]     <- sum(K_ij_w*Ym[(1:n)[!((1:n)%in%i)],])/((n-1)*prod(a)*fhat_w[i])
        u_hat[i]    <- Y[i]-Yhat[i]
      }
      # Berechnung von I_n^a:
      Ina <- 0  # Der Wert von I_n^a wird auf Null gesetzt und die Summe ?ber alle n in einer Schleife gerechnet
      for (i in 1:n){
        K_ij      <- rowProds(dnorm((matrix(rep(Xm[i,],n),nrow=n,byrow = T)-Xm)[(1:n)[!((1:n)%in%i)],]%*%diag(1/h,(q1+q2),(q1+q2))))
        Ina <- Ina+sum((u_hat[i]*fhat_w[i])*(u_hat[(1:n)[!((1:n)%in%i)]]*fhat_w[(1:n)[!((1:n)%in%i)]])*K_ij)
      }
      Ina <- Ina/(n*(n-1)*prod(h))  # Die Summe wird durch (n*(n-1)*PROD[bw]) geteilt, wodurch der entg?ltige Wert berechnet ist

      # Berechnung des Sch?tzers f?r die Standartabweichung:
      sig2a <- 0 # Der Wert von I_n^a wird auf Null gesetzt und die Summe ?ber alle n in einer Schleife gerechnet
      for (i in 1:n){
        K_ij      <- rowProds(dnorm((matrix(rep(Xm[i,],n),nrow=n,byrow = T)-Xm)[(1:n)[!((1:n)%in%i)],]%*%diag(1/h,(q1+q2),(q1+q2))))
        sig2a <- sig2a+sum((u_hat[i]*fhat_w[i])^2*(u_hat[(1:n)[!((1:n)%in%i)]]*fhat_w[(1:n)[!((1:n)%in%i)]])^2*K_ij)
      }
      sig2a <- sig2a/(n*(n-1)*prod(h)*(2*sqrt(pi))^(q1+q2))  # Die Summe wird durch (n*(n-1)*PROD[bw]) geteilt, wodurch der entg?ltige Wert berechnet ist


      Ta <- n*sqrt(prod(h))*Ina/(sqrt(2*sig2a))
    }
    else{
      if(type=="epanechnikov"){
        # fuer i=1,...,n werden fhat_w_i, Yhat_i, u_hat_i berechnet und in den leeren Vektoren gespeichert
        for (i in 1:n){
          K_ij_w      <- rowProds(epa((matrix(rep(X1m[i,],n),nrow=n,byrow = T)-X1m)[(1:n)[!((1:n)%in%i)],]%*%diag(1/a,q1,q1)))
          fhat_w[i]   <- sum(K_ij_w)/((n-1)*prod(a))
          Yhat[i]     <- ifelse(fhat_w[i]!=0,sum(K_ij_w*Ym[(1:n)[!((1:n)%in%i)],])/((n-1)*prod(a)*fhat_w[i]),0)
          u_hat[i]    <- Y[i]-Yhat[i]
        }
        # Berechnung von I_n^a:
        Ina <- 0  # Der Wert von I_n^a wird auf Null gesetzt und die Summe ?ber alle n in einer Schleife gerechnet
        for (i in 1:n){
          K_ij      <- rowProds(epa((matrix(rep(Xm[i,],n),nrow=n,byrow = T)-Xm)[(1:n)[!((1:n)%in%i)],]%*%diag(1/h,(q1+q2),(q1+q2))))
          Ina <- Ina+sum((u_hat[i]*fhat_w[i])*(u_hat[(1:n)[!((1:n)%in%i)]]*fhat_w[(1:n)[!((1:n)%in%i)]])*K_ij)
        }
        Ina <- Ina/(n*(n-1)*prod(h))  # Die Summe wird durch (n*(n-1)*PROD[bw]) geteilt, wodurch der entg?ltige Wert berechnet ist

        # Berechnung des Schaetzers fuer die Standartabweichung:
        sig2a <- 0 # Der Wert von I_n^a wird auf Null gesetzt und die Summe ?ber alle n in einer Schleife gerechnet
        for (i in 1:n){
          K_ij      <- rowProds(epa((matrix(rep(Xm[i,],n),nrow=n,byrow = T)-Xm)[(1:n)[!((1:n)%in%i)],]%*%diag(1/h,(q1+q2),(q1+q2))))
          sig2a <- sig2a+sum((u_hat[i]*fhat_w[i])^2*(u_hat[(1:n)[!((1:n)%in%i)]]*fhat_w[(1:n)[!((1:n)%in%i)]])^2*K_ij)
        }
        sig2a <- (sig2a*(9/15)^(q1+q2))/(n*(n-1)*prod(h))  # Die Summe wird durch (n*(n-1)*PROD[bw]) geteilt, wodurch der entg?ltige Wert berechnet ist


        Ta <- n*sqrt(prod(h))*Ina/(sqrt(2*sig2a))
      }
      else{
        if(type=="uniform"){
          # f?r i=1,...,n werden fhat_w_i, Yhat_i, u_hat_i berechnet und in den leeren Vektoren gespeichert
          for (i in 1:n){
            K_ij_w      <- rowProds(uni((matrix(rep(X1m[i,],n),nrow=n,byrow = T)-X1m)[(1:n)[!((1:n)%in%i)],]%*%diag(1/a,q1,q1)))
            fhat_w[i]   <- sum(K_ij_w)/((n-1)*prod(a))
            Yhat[i]     <- ifelse(fhat_w[i]!=0,sum(K_ij_w*Ym[(1:n)[!((1:n)%in%i)],])/((n-1)*prod(a)*fhat_w[i]),0)
            u_hat[i]    <- Y[i]-Yhat[i]
          }
          # Berechnung von I_n^a:
          Ina <- 0  # Der Wert von I_n^a wird auf Null gesetzt und die Summe ?ber alle n in einer Schleife gerechnet
          for (i in 1:n){
            K_ij      <- rowProds(uni((matrix(rep(Xm[i,],n),nrow=n,byrow = T)-Xm)[(1:n)[!((1:n)%in%i)],]%*%diag(1/h,(q1+q2),(q1+q2))))
            Ina <- Ina+sum((u_hat[i]*fhat_w[i])*(u_hat[(1:n)[!((1:n)%in%i)]]*fhat_w[(1:n)[!((1:n)%in%i)]])*K_ij)
          }
          Ina <- Ina/(n*(n-1)*prod(h))  # Die Summe wird durch (n*(n-1)*PROD[bw]) geteilt, wodurch der entg?ltige Wert berechnet ist

          # Berechnung des Sch?tzers f?r die Standartabweichung:
          sig2a <- 0 # Der Wert von I_n^a wird auf Null gesetzt und die Summe ?ber alle n in einer Schleife gerechnet
          for (i in 1:n){
            K_ij      <- rowProds(uni((matrix(rep(Xm[i,],n),nrow=n,byrow = T)-Xm)[(1:n)[!((1:n)%in%i)],]%*%diag(1/h,(q1+q2),(q1+q2))))
            sig2a <- sig2a+sum((u_hat[i]*fhat_w[i])^2*(u_hat[(1:n)[!((1:n)%in%i)]]*fhat_w[(1:n)[!((1:n)%in%i)]])^2*K_ij)
          }
          sig2a <- (sig2a*0.5^(q1+q2))/(n*(n-1)*prod(h))  # Die Summe wird durch (n*(n-1)*PROD[bw]) geteilt, wodurch der entg?ltige Wert berechnet ist


          Ta <- n*sqrt(prod(h))*Ina/(sqrt(2*sig2a))
        }
        else{
          stop("kernel is not gaussian,epanechnikov or uniform")}   # Warnmeldung und Abbruch, falls die Zeilen
        # von W,Z und Y nicht übereinstimmen
      }
    }
    if(return.bw==F){
      return(Ta)}else{
        if(return.bw==T){
          return(c(Ta,bw))}else{
            stop("return.bw is not boolean")
          }
      }

  }
}
