library("bnlearn")

#se crea el modelo teorico
sachs.modelstring <-
  paste("[PKC][PKA|PKC][Raf|PKC:PKA][Mek|PKC:PKA:Raf][Erk|Mek:PKA][Akt|Erk:PKA][P38|PKC:PKA][Jnk|PKC:PKA][Plcg][PIP3|Plcg][PIP2|Plcg:PIP3]")
dag.sachs <- model2network(sachs.modelstring)

#se crea modelo por parte del archivo entregado
data = read.table("SACHS10k.csv", header = TRUE,sep=",", colClasses = "factor")
res = hc(data)
#se comparan ambos modelos
unlist(compare(dag.sachs, res))
plot(res)

#se comparan los 3 metodos
set.seed(6)
res = hc(data)
unlist(compare(dag.sachs, res))

set.seed(6)
res = mmhc(data)
unlist(compare(dag.sachs, res))

set.seed(6)
res = mmpc(data)
unlist(compare(dag.sachs, res))

#se realiza el cambio de columnas en el modelo
data = read.table("SACHS10o.csv", header = TRUE,sep=",", colClasses = "factor")

#se comparan los 3 metodos
set.seed(6)
res = hc(data)
unlist(compare(dag.sachs, res))

set.seed(6)
res = mmhc(data)
unlist(compare(dag.sachs, res))

set.seed(6)
res = mmpc(data)
unlist(compare(dag.sachs, res))

#al saber que hc entrega los mejores resultados se intenta variar la semilla y el restart
maximoTp = 0
for(seed in 1:500){
  for(i in 1:100){
    set.seed(seed)
    res = hc(data,restart = i)
    tph = unlist(compare(dag.sachs, res))["tp"]
    if(maximoTp < tph){
      maximoTp = tph
      print(paste("tp: ",maximoTp," i: ",i," seed: ",seed, "HC"))
    }
  }
}

# [1] "tp:  11  i:  1  seed:  1 HC"
# [1] "tp:  12  i:  8  seed:  1 HC"
# [1] "tp:  13  i:  16  seed:  2 HC"
# [1] "tp:  14  i:  36  seed:  2 HC"
# [1] "tp:  15  i:  100  seed:  2 HC"
# [1] "tp:  16  i:  35  seed:  6 HC"

#se intenta buscar un rango en el archivo que sea funcional y nos entrege el mismo resultado de prediccion
rango = 100
maximoTp = 0
for(i in 1:100){
  set.seed(6)
  res = hc(data[1:rango,],restart = 35)
  tph = unlist(compare(dag.sachs, res))["tp"]
  if(tph == 16){
    maximoTp = tph
    print(paste("tp: ",maximoTp," rango: ",rango, "HC"))
  }
  rango = rango + 100
}
# [1] "tp:  3  rango:  100 HC"
# [1] "tp:  4  rango:  200 HC"
# [1] "tp:  5  rango:  300 HC"
# [1] "tp:  7  rango:  500 HC"
# [1] "tp:  8  rango:  900 HC"
# [1] "tp:  9  rango:  1000 HC"
# [1] "tp:  10  rango:  1600 HC"
# [1] "tp:  11  rango:  1800 HC"
# [1] "tp:  13  rango:  3400 HC"
# [1] "tp:  16  rango:  3600 HC"

set.seed(6)
res = hc(data,restart = 35)
unlist(compare(dag.sachs, res))
score(res,data)
plot(res)


#comparar metodos
data = read.table("SACHS10o.csv", header = TRUE,sep=",", colClasses = "factor")
set.seed(6)
res = hc(data)
unlist(compare(dag.sachs, res))
score(res,data)

set.seed(6)
res = mmhc(data)
unlist(compare(dag.sachs, res))
score(res,data)

set.seed(6)
res = mmpc(data)
unlist(compare(dag.sachs, res))
score(res,data)

#ontener tablas de progagacion de la evidencia
set.seed(6)
res = hc(data ,restart = 35)
unlist(compare(dag.sachs, res))

#tenemos 2 modelos y vamos a obtener algunas tablas de progagacion
#de la evidencia

fittedbn = bn.fit(res, data = data)
print(fittedbn$Mek)

fittedbn = bn.fit(dag.sachs, data = data)
print(fittedbn$Mek)

fittedbn = bn.fit(res, data = data)
print(fittedbn$Plcg)

fittedbn = bn.fit(dag.sachs, data = data)
print(fittedbn$Plcg)

#realizar consultas
for (i in names(data))
  levels(data[, i]) = c("LOW", "AVG", "HIGH")

fittedbn = bn.fit(res, data = data)
fittedbn2 = bn.fit(dag.sachs, data = data)


set.seed(6)
cpquery(fittedbn, event = (Jnk == "LOW"), evidence = (PKA=="HIGH" & PKC=="HIGH"))

set.seed(6)
cpquery(fittedbn, event = (Jnk =="AVG"), evidence = (PKA=="HIGH" & PKC=="HIGH"))

set.seed(6)
cpquery(fittedbn, event = (PKA=="HIGH" & PKC=="HIGH"), evidence = (Jnk == "LOW"))

set.seed(6)
cpquery(fittedbn, event = (PKA=="HIGH" & PKC=="HIGH"), evidence = (Jnk =="AVG"))


set.seed(6)
cpquery(fittedbn, event = (Erk=="LOW"), evidence = (PKC =="HIGH"))

set.seed(6)
cpquery(fittedbn, event = (Erk=="HIGH"), evidence = (PKC =="HIGH"))

set.seed(6)
cpquery(fittedbn, event = (Erk=="AVG"), evidence = (PKC =="HIGH"))

set.seed(6)
cpquery(fittedbn, event = (PKC =="LOW"), evidence = ((Erk =="AVG") | (Erk =="LOW") | (Erk =="HIGH")))
