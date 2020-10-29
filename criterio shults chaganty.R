### MATRICES DE CORRELACIÓN DE TRABAJO

#Intercambiable
exc_cor <- function(n,rho){
  corr = matrix(rho,nrow=n,ncol=n)
  for (i in 1:n){
    for (j in 1:n){
      if (i==j){
        corr[i,j]=1
      }
    }
  }
  return(corr)
}
#AR(1)
ar1_cor <- function(n, rho){
  exponent <- abs(matrix(1:n-1, nrow=n, ncol=n, byrow=TRUE)-(1:n-1))
  return(rho^exponent)
}


### CRITERIO SHULTS CHAGANTY

shults_chaganty <- function(estructura,modelo,dataset){
  datos=as.matrix.data.frame(dataset)
  p=ncol(datos)
  #Creo 3 columnas que serán: mu, S=y-mu, varianza=mu(1-mu)
  #mu
  mu=modelo$fitted.values 
  datos=cbind(datos,mu)
  #S=y-mu (y tiene que estar en la columna 2 del dataset)
  S=datos[,2]-datos[,p+1]
  datos=cbind(datos,S)
  #varianza=mu(1-mu)
  varianza=mu*(1-mu)
  datos=cbind(datos,varianza)
  
  colnames(datos)<-NULL
  rownames(datos)<-NULL
  
  id=unique(datos[,1]) #individuos (id tiene que estar en la columna 1 del dataset)
  
  SC_individuo=rep(0,times=max(id))
  m=p+3 #número de columnas de la nueva tabla
  
  alpha=modelo$geese$alpha
  
  for (i in id){
    matriz_individuo = datos[datos[,1]==i,] #en cada iteración, coge los datos del individuo i
    t=nrow(matriz_individuo)
    #Determinar S
    S_trans=matrix(c(matriz_individuo[,m-1]),nrow=1) #matriz 1xt
    S=matrix(c(matriz_individuo[,m-1]),ncol=1) #matriz tx1
    #Determinar A
    var=matriz_individuo[,m]
    A=matrix(0,nrow=t,ncol=t)
    for (j in 1:t){
      A[j,j]=var[j]
    }
    A_sqrt=A^(1/2)
    #Determinar R
    if (estructura=="AR1"){R=ar1_cor(t,alpha)}
    if (estructura=="EXC"){R=exc_cor(t,alpha)}
    if (estructura=="IND"){R=diag(t)}
    #Determinar V=sqrt(A)*R*sqrt(A)
    V=(A_sqrt%*%R%*%A_sqrt)
    V_inv=solve(V)
    
    #Determinar SC de cada invividuo
    SC_valor=as.numeric(S_trans%*%V_inv%*%S)
    SC_individuo[i]=SC_valor
  }
  SC=sum(SC_individuo)
  return(SC)
}

