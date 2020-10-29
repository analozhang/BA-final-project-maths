library(readr)
madras <- read_table2("http://faculty.washington.edu/heagerty/Books/AnalysisLongitudinal/madras.data", col_names = FALSE, col_types = cols(X1 = col_integer(), X2 = col_integer(), X3 = col_integer(), X4 = col_integer(), X5 = col_integer(), X6 = col_integer(), X7 = col_integer()))
names(madras) = c("id", "symptom", "month", "age", "gender", "monthxage", "monthxgender")
madras = as.data.frame(madras)

id=madras[madras[,3]==11,][,1] #id de los individuos con observaciones completas
filas_id=c()
for (i in id){
  for (j in 1:nrow(madras)){
    if(madras[j,1]==i){filas_id=append(filas_id,j)}
  }
}
madras=madras[filas_id,]
madras_df=as.data.frame(madras) #dataset balanceado (12 observaciones por individuo)

library(geepack)

#Ajustamos por GEE
#1)Independiente
m.ind <- geeglm(symptom ~ month + age + gender + monthxage + monthxgender, data=madras_df, id=id, family=binomial, corstr="independence")
summary(m.ind)
#2)Intercambiable
m.ex <- geeglm(symptom ~ month + age + gender + monthxage + monthxgender, data=madras_df, id=id, family=binomial, corstr="exchangeable")
summary(m.ex)
#3)AR(1)
m.ar1 <- geeglm(symptom ~ month + age + gender + monthxage + monthxgender, data=madras_df, id=id, family=binomial, corstr="ar1")
summary(m.ar1)

#Criterios de selección de la matriz de correlación de trabajo
#QIC y CIC
QIC(m.ind)[c(1,4)]
QIC(m.ex)[c(1,4)]
QIC(m.ar1)[c(1,4)]
#Shults-Chaganty
shults_chaganty("IND",m.ind,madras_df)
shults_chaganty("EXC",m.ex,madras_df)
shults_chaganty("AR1",m.ar1,madras_df)

