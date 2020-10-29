set.seed(1)

### FUNCIONES AUXILIARES

alpha_iniciales <- function(p,corr){
  matriz=matrix(0,nrow=3,ncol=3)
  for (i in 1:3){
    for (j in 1:3){
      q = 1-p
      matriz[i,j]=log(1+corr[i,j]*(q[i]*(1/p[i])*q[j]*(1/p[j]))^(1/2))
    }
  }
  return(matriz) 
}

posicion <- function(valor,matriz){  
  matriz[lower.tri(matriz)]=0
  posicion = c()
  for (i in 1:3){
    for (j in 1:3){
      if (matriz[i,j]==valor){
        posicion=c(i,j)
      }
    }
  }
  return (posicion)
}

minimo <- function(matriz){ #Devuelve el valor mínimo entre los elementos >0 de una matriz
  valores_positivos=matriz[matriz>0] #Tomo los valores positivos
  minimo=sort(valores_positivos)[1] #Los ordeno de menor a mayor y cojo el primero 
  return(minimo) 
}

library(gtools)
permutaciones <- function(vector){ #Devuelve todas las permutaciones (con reemplazamiento) del conjunto de valores {1,2,3}
  permutaciones=matrix(0,ncol=2)
  long=length(vector)
  if (long==1){ 
    permutaciones=matrix(c(vector,vector),ncol=2)
  }
  if (long==2){ 
    posiciones=permutations(long,2,repeats.allowed=TRUE) #Las posiciones del vector a permutar
    for (i in 1:nrow(posiciones)){ 
      permutaciones=rbind(permutaciones,c(vector[posiciones[i,1]],vector[posiciones[i,2]]))
    }
    permutaciones=permutaciones[-1,] #Elimino la primera fila de ceros
  }
  if (long==3){
    permutaciones=permutations(long,2,repeats.allowed=TRUE)
  }
  return(permutaciones)
}


### FUNCIÓN PRINCIPAL  

algoritmo_park_parte1 <- function(p,corr){
  #En betas, S y matrices_alphas se almacenan los resultados en cada iteración l
  betas=c()
  S=list()
  matrices_alphas=list()
  
  matriz_alpha=alpha_iniciales(p,corr) #Esta matriz_alpha se actualiza en cada l=l+1
  matrices_alphas[[1]]=alpha_iniciales(p,corr) #en l=0

  l=1 #contador
  while (sum(c(matriz_alpha)>0)){ #Mientras que todos los elementos de alpha sean distintos de 0
    
    #PASO 1: Determinar beta_l y S_l
    alpha_minimo=minimo(matriz_alpha) #Es el mínimo alpha >0
    betas=c(betas,alpha_minimo) #alpha_minimo se almacena en betas
    posicion_beta=posicion(alpha_minimo,matriz_alpha) #Saco su posición en matriz_alpha
    r=posicion_beta[1]
    s=posicion_beta[2]
    if (sum(matriz_alpha[r,r]==0|matriz_alpha[s,s]==0)>0){
      return(print("Se para el algoritmo"))
    }else{
      S_inicial=unique(c(r,s)) 
      S_l=S_inicial #La posición de beta_l va al conjunto S_l
      for (k in 1:3){ #k=1,2,3 porque hemos hay 3 periodos de tiempo
        cont=0 #Este contador +1 si algún elemento es 0
        long=length(S_l) #El tamaño de S_l es de 1 o 2
        for (i in 1:long){
          if (matriz_alpha[i,k]==0){cont=cont+1} #Para cada k, se prueba con cada elemento de S_l -> Ver si en esa posicion es =0
        }
        if (cont==0){S_l=union(S_l,k)} #Si ninguno de los valores observados es =0, entonces k se almacenará en S_l (junto con S_inicial)
      }
      S[[l]]=S_l #El conjunto S_l al final de la etapa se guarda en S
      
      #PASO 2: Actualizar las alphas
      posiciones_actualizar=permutaciones(S_l) #Las posiciones de los alpha que se van a actualizar
      for (i in 1:nrow(posiciones_actualizar)){
        if (matriz_alpha[posiciones_actualizar[i,1],posiciones_actualizar[i,2]]>0){
          matriz_alpha[posiciones_actualizar[i,1],posiciones_actualizar[i,2]]=matriz_alpha[posiciones_actualizar[i,1],posiciones_actualizar[i,2]]-alpha_minimo
        }
      }
    }
    matrices_alphas[[l+1]]=matriz_alpha #Guardo la matriz_alpha del paso l en matrices_alphas
    l=l+1
  }
  return(list(betas,S,matrices_alphas))
}

algoritmo_park_parte2 <- function(p,corr,m){ #m es el número de variables Poisson de la simulación
  #1a PARTE: Determinamos la distribución conjunta (Y1,Y2,Y3) en función de las X1,...,X6
  Y1_indicesX=c()
  Y2_indicesX=c()
  Y3_indicesX=c()
  S=algoritmo_park_parte1(p,corr)[[2]]
  for (i in 1:length(S)){
    numero_elementos=length(S[[i]])
    for (j in 1:numero_elementos){
      if (S[[i]][j]==1){Y1_indicesX=append(Y1_indicesX,i)}
      if (S[[i]][j]==2){Y2_indicesX=append(Y2_indicesX,i)}
      if (S[[i]][j]==3){Y3_indicesX=append(Y3_indicesX,i)}
    }
  }

  #2a PARTE: Generamos m muestras de Z simulando las X
  Z=matrix(0,nrow=m,ncol=3)
  #Saco las betas
  betas=algoritmo_park_parte1(p,corr)[[1]]
  for (i in 1:m){
    poisson=c()
    for (j in 1:length(betas)){
      poisson[j]=rpois(1,betas[j])
    }
    #Hago la suma de los valores Poisson 
    Y1_valor=sum(poisson[Y1_indicesX])
    Y2_valor=sum(poisson[Y2_indicesX])
    Y3_valor=sum(poisson[Y3_indicesX])
    #Zi es la función indicatriz de que Yi=0   
    Z_valor=c()
    if (Y1_valor==0){Z_valor[1]=1}else{Z_valor[1]=0}
    if (Y2_valor==0){Z_valor[2]=1}else{Z_valor[2]=0}
    if (Y3_valor==0){Z_valor[3]=1}else{Z_valor[3]=0}
    Z[i,]=Z_valor
  }
  
  #3a PARTE: Calculamos las probabilidades de cada posible resultado de (Z1,Z2,Z3)
  posiciones=permutations(2,3,repeats.allowed=TRUE) #Las posiciones del vector a permutar
  vector=c(0,1)
  resultados_Z=matrix(nrow=nrow(posiciones),ncol=3) #Será una matriz con todas las posibles combinaciones de Z
  for (i in 1:nrow(posiciones)){ 
    resultados_Z[i,]=c(vector[posiciones[i,1]],vector[posiciones[i,2]],vector[posiciones[i,3]])
  }
  resultados_Z=cbind(resultados_Z,rep(0,times=nrow(resultados_Z))) #Añado una columna donde se irán contando las frecuencias
  
  for (i in 1:nrow(Z)){
    for (j in 1:nrow(resultados_Z)){
      if (sum(Z[i,]==resultados_Z[j,(1:3)])==3){ #Cojo la muestra i de Z y busco su fila en resultados_Z -> Cuando coincida, +1
        resultados_Z[j,4]=resultados_Z[j,4]+1
      }
    }
  }
  resultados_Z[,4]=resultados_Z[,4]/m #prob=frecuencia/tamaño muestra total
  
  prob=resultados_Z[,4]
  prob_acumuladas=c(prob[1])
  for (i in 2:nrow(resultados_Z)){
    prob_acumuladas[i]=sum(prob[1:i])
  }
  prob_Z=cbind(resultados_Z,prob_acumuladas)
  colnames(prob_Z)<-c("visita_1","visita_2","visita_3","prob","prob_acum")
  return(prob_Z)
}

algoritmo_park_final <- function(p,corr,m,n){
  prob_Z=algoritmo_park_parte2(p,corr,m)
  prob_acum=prob_Z[,5]
  datos_correlacionados=matrix(nrow=n,ncol=3)
  for (i in 1:n){
    u=runif(1,0,1)
    if (u<=prob_acum[1]){
      datos_correlacionados[i,]=prob_Z[1,(1:3)]
    }else{
      for (j in 2:nrow(prob_Z)){
        if (prob_acum[j-1]<u & u<=prob_acum[j]){
          datos_correlacionados[i,]=prob_Z[j,(1:3)]
        }
      }
    }
  }
  return(datos_correlacionados)
}


### Datos del paper de Park para comprobar

# Matriz de correlación desestructurada para t=3
uns_corr <- function(a12,a13,a23){
  corr = matrix(c(1,a12,a13,a12,1,a23,a13,a23,1),nrow=3,ncol=3)
  return(corr)
}

p = c(0.9,0.8,0.7)
corr=uns_corr(0.1,0.5,0.5)

betas=algoritmo_park_parte1(p,corr)[[1]]
S=algoritmo_park_parte1(p,corr)[[2]]
matrices_alphas=algoritmo_park_parte1(p,corr)[[3]]

prob=algoritmo_park_parte2(p,corr,m=10000)
