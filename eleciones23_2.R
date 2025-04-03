library(combinat)
#PP       "A"137
#UPN       "B"1
#ERC       "C"7
#Cca       "D"1
#Vox       "E"33
#EAJ/PNV   "F"5
#Junts     "G"7
#BNG       "H"1
#Sumar     "I"31
#PSOE      "J"121
#EH Bildu  "K"6
y <- c("A","B","C","D","E","F","G","H","I","J","K")
escaños <- c(137,1,7,1,33,5,7,1,31,121,6)
complejo <- c("A","B","C","D","E","F","G","H","I","J","K",
              "HC","DI","GB","DB","KJ","HG","HK","GA",
              "DC","AB","CI","DF","HJ","KI","JI",
              "KC","JB","DJ","DK","HI","FB","CJ","KG",
              "DA","FC","FA","FJ","HD","DKC","HDJ","HDI",
              "GAB","DFJ","HDK","DAB","DKJ","HCI","HKC",
              "HCJ","DFC","DJI","DKI","DCI","HKI","FJB",
              "KCJ","HJI","FAB","KCI","HKJ","DCJ","DJB",
              "HDC","FCJ","KJI","CJI","HKG","DFA","DFB",
              "HDCJ","HDKJ","KCJI","HDJI","HDKI","HCJI",
              "DKJI","DFAB","HDKC","DFJB","DCJI","DKCI",
              "DKCJ","HKJI","HDCI","DFCJ","HKCJ","HKCI",
              "HKCJI","HDKJI","HDCJI","HDKCI","HDKCJ","DKCJI","HDKCJI")

L <- c()
for (i in complejo){  
  i <- strsplit(i, "")[[1]]
  i <- sort(i)
  i <- paste(i,collapse = "")
  L <-c (L,i)
}
complejo <- L
votos <- 0


###############################################################################
#                                                                             #
#                  Calculo de los vecino de una coalicion                     #
#                                                                             #
###############################################################################

vecinos <- function(x){
  aristas <- c()
  for(simplice in complejo){
    simplice <- strsplit(simplice, "")[[1]]
    if(length(simplice)==2){
      arista <- paste(simplice,collapse = "")
      aristas <- c(aristas,arista)
    }
  }
  vecinos <- c()
  vertices <- strsplit(x, "")[[1]]
  for(arista in aristas){
    arista <- strsplit(arista, "")[[1]]
    for(vertice in vertices){
      if(vertice %in% arista){
        vecino <- setdiff(arista,vertice)
        if(vecino %in% vertices | vecino %in% vecinos){next}
        vecinos <- c(vecinos,vecino)
      }
    }
  }
  return(vecinos)
}

###############################################################################
#                                                                             #
#                  Calculo de las componentes conexas de K                    #
#                                                                             #
###############################################################################

#componentes_conexas <- function(x,complejo){
componentes_conexas <- complejo
nuevas_componentes <- c()
for(componente in complejo){
  T <- setdiff(complejo,componente)
  for(conectado in T){
    componente <- paste(componente,collapse = "")
    componente <- strsplit(componente, "")[[1]]
    conectado <- strsplit(conectado, "")[[1]]
    if(any(componente %in% conectado)){
      componente_conexa <- union(componente,conectado)
      componente_conexa <- sort(componente_conexa)
      componente_conexa <- paste(componente_conexa,collapse = "")
      if(componente_conexa %in% componentes_conexas){next}else{
        componentes_conexas <- c(componentes_conexas,componente_conexa)
        nuevas_componentes <- c(nuevas_componentes,componente_conexa)}
    }
  }
}
cont <- 0
while(cont < 4){
  for(componente in componentes_conexas){
    for(conectado in nuevas_componentes){
      componente <- paste(componente,collapse = "")
      componente <- strsplit(componente, "")[[1]]
      conectado <- strsplit(conectado, "")[[1]]
      if(any(componente %in% conectado)){
        componente_conexa <- union(componente,conectado)
        componente_conexa <- sort(componente_conexa)
        componente_conexa <- paste(componente_conexa,collapse = "")
        ifelse(componente_conexa %in% componentes_conexas,next,componentes_conexas <- c(componentes_conexas,componente_conexa))
      }
    }
  }
  #componente_conexa <- strsplit(componente_conexa, "")[[1]]
  #cont <- length(componente_conexa)+1
  cont <- cont + 1
}
#   return(componentes_conexas)
# }

###############################################################################
#                                                                             #
#         Calculo de las particiones maximales de una coalicion               #
#                                                                             #
###############################################################################

f_particion_maximal <- function(x){
  maximales <- c()
  resto <- complejo
  for(simplice in complejo){
    resto <- setdiff(resto,simplice)
    contenidos <- 0
    for(simplice2 in resto){
      simplice1 <- strsplit(simplice, "")[[1]]
      simplice2 <- strsplit(simplice2, "")[[1]]
      if(all(simplice1 %in% simplice2)){contenidos <- contenidos + 1}
    }
    if(contenidos == 0){maximales <- c(maximales,simplice)}
  }
  
  posibles_intersecciones <- c()
  x <- strsplit(x, "")[[1]]
  for(max in maximales){
    max <- strsplit(max, "")[[1]]
    intersec <- intersect(max,x)
    if(length(intersec)==0){next}
    interseccion <- paste(intersec,collapse = "")
    posibles_intersecciones <- c(posibles_intersecciones, interseccion)
  }
  
  intersecciones <- c()
  for(int in posibles_intersecciones){
    resto <- posibles_intersecciones
    resto <- setdiff(resto,int)
    contenidos <- 0
    for(int2 in resto){
      int1 <- strsplit(int, "")[[1]]
      int2 <- strsplit(int2, "")[[1]]
      if(all(int1 %in% int2)){contenidos <- contenidos + 1}
    }
    if(contenidos == 0){intersecciones <- c(intersecciones, int)}
  }
  
  cuasi_particiones_maximales <- c()
  cuasi_particiones_max <- c()
  fijos <- c()
  i <- 0
  indices <- c()
  for(fijo in intersecciones){
    i = i + 1
    fijos <- c(fijos,fijo)
    cuasi_particion_maximal <- c()
    resto <- intersecciones
    resto <- setdiff(resto,fijo)
    for(fijo2 in resto){
      fijo1 <- strsplit(fijo, "")[[1]]
      fijo2 <- strsplit(fijo2, "")[[1]]
      interseccion <- intersect(fijo1,fijo2)
      cuasi_particion <- setdiff(fijo2,interseccion)
      cuasi_particion <- paste(cuasi_particion,collapse = "")
      cuasi_particion_maximal <- c(cuasi_particion_maximal,cuasi_particion)
    }
    cuasi_particion_maximal <- unique(cuasi_particion_maximal)
    cuasi_particion_maximal_opt <- c()
    for(cpm in cuasi_particion_maximal){
      resto <- cuasi_particion_maximal
      resto <- setdiff(resto,cpm)
      resto <- unique(resto)
      if(length(resto)==0){
        cuasi_particion_maximal_opt <- c(cuasi_particion_maximal_opt, cpm)
        next
      }
      for(cpm2 in resto){
        contenidos <- 0
        cpm1 <- strsplit(cpm, "")[[1]]
        cpm2 <- strsplit(cpm2, "")[[1]]
        if(all(cpm1 %in% cpm2)){contenidos <- contenidos + 1}
      }
      if(contenidos == 0){
        cuasi_particion_maximal_opt <- c(cuasi_particion_maximal_opt, cpm)
      }
    }
    cuasi_particion_maximal_opt <- sort(cuasi_particion_maximal_opt)
    indices <- c(indices,rep(i,length(cuasi_particion_maximal_opt)))
    cuasi_particiones_max <- c(cuasi_particiones_max,cuasi_particion_maximal_opt)
  }
  cuasi_particiones_maximales <- split(cuasi_particiones_max,indices)
  
  particiones_maximales <- c()
  indices <- c()
  l <- 0
  for(i in 1:length(cuasi_particiones_maximales)){
    if(length(cuasi_particiones_maximales[[i]])==1){
      l = l + 1
      particion_maximal <- c()
      particion_maximal <- c(particion_maximal,fijos[i],cuasi_particiones_maximales[[i]])
      indices <- c(indices,rep(l,length(particion_maximal)))
      particiones_maximales <- c(particiones_maximales,particion_maximal)
      next
    }else{permutaciones <- permn(cuasi_particiones_maximales[[i]])}
    for(h in 1:(length(permutaciones)/2)){
      l = l + 1
      particion_maximal <- c()
      particion_maximal <- c(particion_maximal,fijos[i])
      acumula <- c()
      for(j in 1:length(permutaciones[[h]])){
        particion <- c()
        if(j == 1){
          acumula1 <- strsplit(permutaciones[[h*2]][j], "")[[1]]
          acumula <- c(acumula,acumula1)
          acumula1 <- paste(acumula1,collapse = "")
          particion_maximal <- c(particion_maximal,acumula1)
        }
        if(j==2){
          acumula <- c(acumula,permutaciones[[h*2]][j])
          acumulado <- unlist(strsplit(as.character(acumula),""))
          acumula <- acumulado[acumulado != ""]
          acumula <- unique(acumula)
          cpmm <- permutaciones[[h*2]][j-1]
          cpmm <- strsplit(cpmm, "")[[1]]
          cpmm2 <- permutaciones[[h*2]][j]
          cpmm2 <- strsplit(cpmm2, "")[[1]]
          interseccion <- intersect(cpmm,cpmm2)
          particion <- setdiff(cpmm2,interseccion)
          if(length(particion)==0){next}
          particion <- paste(particion,collapse = "")
          particion_maximal <- c(particion_maximal,particion)
        }
        if(j>2){
          acumula <- c(acumula,permutaciones[[h*2]][j])
          acumulado <- unlist(strsplit(as.character(acumula),""))
          acumula <- acumulado[acumulado != ""]
          acumula <- unique(acumula)
          cpmm <- permutaciones[[h*2]][j-1]
          cpmm <- strsplit(cpmm, "")[[1]]
          interseccion <- intersect(acumula,cpmm)
          resta <- union(interseccion,cpmm)
          particion <- setdiff(cpmm,resta)
          if(length(particion)==0){next}
          particion <- paste(particion,collapse = "")
          particion_maximal <- c(particion_maximal,particion)  
        }
      }
      indices <- c(indices,rep(l,length(particion_maximal)))
      particiones_maximales <- c(particiones_maximales,particion_maximal)
    }
  }
  particiones_maximales <- split(particiones_maximales,indices)
  
  #Quitar Duplicados
  ordenar_convertir <- function(x) toString(sort(x))
  duplicados_inversos <- duplicated(sapply(particiones_maximales,ordenar_convertir))
  particiones_maximales <- particiones_maximales[!duplicados_inversos]
  
  #Funcion el cuadrado de |S|
  #   suma_maxima <- c()
  #   for(i in particiones_maximales){
  #     suma <- 0
  #     for(j in i){
  #       j <- strsplit(j, "")[[1]]
  #       suma <- suma + length(j)^2
  #     }
  #     suma_maxima <- c(suma_maxima,suma)
  #   }
  #   return(max(suma_maxima))
  
  #Funcion el |S| * maximo
  #   funcion_particion <- c()
  #   for(i in particiones_maximales){
  #     suma <- 0
  #     for(j in i){
  #       j <- strsplit(j, "")[[1]]
  #       vector <- c()
  
  #       for(k in j){
  #         contador <- 0
  #         for(g in complejo){
  #           contador <- contador + 1
  #           if(k==g){
  #             break
  #           }
  #         }
  #         vector <- c(vector,y[contador])
  #       }
  #       funcion <- max(vector) * length(vector)
  #       suma <- suma + funcion
  #     }
  #     funcion_particion <- c(funcion_particion,suma)
  #   }
  #Funcion |S|*sum(S)*10
  #   funcion_particion <- c()
  #   for(i in particiones_maximales){
  #     suma <- 0
  #     for(j in i){
  #       j <- strsplit(j, "")[[1]]
  #       vector <- c()
  #       for(k in j){
  #         contador <- 0
  #         for(g in complejo){
  #           contador <- contador + 1
  #           if(k==g){
  #             break
  #           }
  #         }
  #         vector <- c(vector,y[contador])
  #       }
  #       funcion <- sum(vector)*10
  #       suma <- suma + funcion
  #     }
  #     funcion_particion <- c(funcion_particion,suma)
  #   }
  #Funcion caracteristica proporcional a las areas
  funcion_particion <- c()
  for(i in particiones_maximales){
    suma <- 0
    for(j in i){
      j <- strsplit(j, "")[[1]]
      vector <- c()
      for(k in j){
        contador <- 0
        for(g in complejo){
          contador <- contador + 1
          if(k==g){
            break
          }
        }
        vector <- c(vector,y[contador])
        votos <- c(votos,escaños[contador])
      }
      funcion <- (sum(votos)*length(vector))^(1/3)
      suma <- suma + funcion
    }
    funcion_particion <- c(funcion_particion,suma)
  }
  return(max(funcion_particion))
}

###############################################################################
#                                                                             #
#                            Calculo del K-valor                              #
#                                                                             #
###############################################################################

f_conexa <- function(x){
  if(x %in% complejo){
    vector <- c()
    x <- strsplit(x, "")[[1]]
    for(k in x){
      contador <- 0
      for(g in complejo){
        contador <- contador + 1
        if(k==g){
          break
        }
      }
      vector <- c(vector,y[contador])
      votos <- c(votos,escaños[contador])
    }
    salida <- (sum(votos)*length(vector))^(1/3)
    return(salida)
  }else{f_particion_maximal(x)}
}

# f_conexa <- function(x){
#   if(x %in% complejo){
#     vector <- c()
#     x <- strsplit(x, "")[[1]]
#     salida <- length(x)^2
#     return(salida)
#   }else{f_particion_maximal(x)}
# }

k_valor <- function(x){
  jugadores <- complejo[1:length(x)]
  GKV <- c()
  for(i in jugadores){
    print(paste("Se esta calculando el jugador = ",i))
    KV <- 0
    vecindad <- c(i)
    vecindad <- c(vecindad,vecinos(i))
    componentes <- c()
    for(cc in componentes_conexas){
      cc <- strsplit(cc, "")[[1]]
      if(any(vecindad %in% cc)){
        cc <- paste(cc,collapse = "")
        componentes <- c(componentes,cc)
      }
    }
    for(j in componentes){
      j <- strsplit(j, "")[[1]]
      if(i %in% j){
        c <- paste(j,collapse = "")
        formula <- f_conexa(c)*(factorial(length(j)-1)*factorial(length(vecinos(c))))/factorial(length(vecinos(c))+length(j))
        KV <- KV + formula
      }else{
        c <- paste(j,collapse = "")
        formula <- f_conexa(c)*(factorial(length(j))*factorial(length(vecinos(c))-1))/factorial(length(vecinos(c))+length(j))
        KV <- KV - formula
      }
    }
    K_por_jugador <- round(KV,2)
    GKV <- c(GKV,K_por_jugador)
    print(paste("Su valor es = ",K_por_jugador))
  }
  return(GKV)
}

k_valor <- k_valor(y)

normalized_k_valor <- function(x){
  V_n <- sum(escaños)
  NKV <- (V_n/sum(k_valor))*k_valor
  return(NKV)
}


round(normalized_k_valor(y),4)
round(normalized_k_valor(y)/sum(escaños),4)
