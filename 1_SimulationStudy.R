tiempos=c(1,2,3) #j=1,2,3
betas_verdaderas=c(0.1,0.2,0.2,0.1) #(beta0,beta1,beta2,beta3)

#Medias de GRUPO 1 (x_i=0)
log_ods_1=c()
for (i in 1:length(tiempos)){
  log_ods_1[i]=betas_verdaderas[1]+betas_verdaderas[3]*tiempos[i]
}
mu1=exp(log_ods_1)/(1+exp(log_ods_1))

#Medias de GRUPO 2 (x_i=1)
log_ods_2=c()
for (i in 1:length(tiempos)){
  log_ods_2[i]=betas_verdaderas[1]+betas_verdaderas[2]+(betas_verdaderas[3]+betas_verdaderas[4])*tiempos[i]
}
mu2=exp(log_ods_2)/(1+exp(log_ods_2))


### FUNCIÓN AUXILIAR: convierte las muestras obtenidas en algoritmo_park_final en formato de datos de panel
datos_correlacionados <- function(mu1,mu2,corr,n1,n2){
  t=length(tiempos)
  filas=(n1+n2)*t #número de muestras
  datos_correlacionados=matrix(0,nrow=filas,ncol=4) #4 columnas: id, y, visita, grupo
  #GRUPO 1
  muestra_grupo1=algoritmo_park_final(mu1,corr,10000,n1) #m=10000 iteraciones para simular las X
  for (i in 1:n1){
    datos_correlacionados[(3*i)-2,]=c(i,muestra_grupo1[i,1],1,0) #visita=1, grupo=0
    datos_correlacionados[(3*i)-1,]=c(i,muestra_grupo1[i,2],2,0) #visita=2, grupo=0
    datos_correlacionados[(3*i),]=c(i,muestra_grupo1[i,3],3,0) #visita=3, grupo=0
  }
  #GRUPO 2
  muestra_grupo2=algoritmo_park_final(mu2,corr,10000,n2)
  ultima_fila_grupo1=n1*3
  for (j in 1:n2){
    datos_correlacionados[ultima_fila_grupo1+((3*j)-2),]=c(n1+j,muestra_grupo2[j,1],1,1) #visita=1, grupo=1
    datos_correlacionados[ultima_fila_grupo1+((3*j)-1),]=c(n1+j,muestra_grupo2[j,2],2,1) #visita=2, grupo=1
    datos_correlacionados[ultima_fila_grupo1+((3*j)),]=c(n1+j,muestra_grupo2[j,3],3,1) #visita=3, grupo=1
  }
  datos_correlacionados_df=as.data.frame(datos_correlacionados)
  colnames(datos_correlacionados_df)=c("id","y","visita","grupo")
  return(datos_correlacionados_df)
}

library(geepack)
### FUNCIÓN PRINCIPAL
simulacion <- function(mu1,mu2,n1,n2,alpha,corstr,iteraciones){
  #Matriz ECM de betas
  matriz_ECM_ind=matrix(nrow=iteraciones,ncol=4) #4 columnas porque hay 4 betas
  matriz_ECM_gee=matrix(nrow=iteraciones,ncol=4)
  
  #Matriz error estándar de betas
  matriz_error_ind=matrix(nrow=iteraciones,ncol=4)
  matriz_error_gee=matrix(nrow=iteraciones,ncol=4)
  
  if (corstr=="ar1"){corr=ar1_cor(3,alpha)}
  if (corstr=="exchangeable"){corr=exc_cor(3,alpha)}
  
  for (i in 1:iteraciones){
    print(i)
    datos=datos_correlacionados(mu1,mu2,corr,n1,n2)
    #Independiente
    modelo_ind = geeglm(y ~ grupo + visita + grupo*visita, family=binomial, data=datos, id=id, corstr="independence")
    matriz_ECM_ind[i,]=(c(modelo_ind$coefficients)-betas_verdaderas)^2
    matriz_error_ind[i,]=(summary(modelo_ind)$coef[,2])
    
    #GEE con corstr
    modelo_gee = geeglm(y ~ grupo + visita + grupo*visita, family=binomial, data=datos, id=id, corstr=corstr)
    matriz_ECM_gee[i,]=(c(modelo_gee$coefficients)-betas_verdaderas)^2
    matriz_error_gee[i,]=(summary(modelo_gee)$coef[,2])
  }
  
  #Función que detecta NaN
  id_sin_NaN <- function(matriz){
    id_sin_NaN=c()
    for (i in 1:nrow(matriz)){
      if (sum(is.na(matriz[i,]))==0){id_sin_NaN=append(id_sin_NaN,i)}   
    }
    return(id_sin_NaN)
  }
  
  #Identifico las filas sin NaN en el error estándar
  id_ind = id_sin_NaN(matriz_error_ind)
  id_gee = id_sin_NaN(matriz_error_gee)
  conteo_NaN_id = iteraciones-length(id_ind)
  conteo_NaN_gee = iteraciones-length(id_gee)
  
  #Me quedo con las iteraciones que no tienen NaN
  matriz_ECM_ind=matriz_ECM_ind[id_ind,]
  matriz_ECM_gee=matriz_ECM_gee[id_gee,]
  matriz_error_ind=matriz_error_ind[id_ind,]
  matriz_error_gee=matriz_error_gee[id_gee,]
  
  ##ECM de betas
  ECM_ind=sqrt(apply(matriz_ECM_ind, 2, mean))
  ECM_gee=sqrt(apply(matriz_ECM_gee, 2, mean))
  
  ##Promedio de error estándar
  Promedio_error_ind=apply(matriz_error_ind, 2, mean)
  Promedio_error_gee=apply(matriz_error_gee, 2, mean)
  
  return(list(ECM_ind,ECM_gee,Promedio_error_ind,Promedio_error_gee,c(conteo_NaN_id,conteo_NaN_gee)))
}

