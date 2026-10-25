##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)
# Carga de datos
datos <- read.table("archivo-entrada.txt", header = TRUE, sep = "\t")

# Inspección visual de la tabla
head(datos)
tail(datos)
dim(datos)
### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)
# Cálculo de la media y desviación estándar para Control y Condición
datos$Media_Control <- apply(datos[, 2:4], 1, mean)
datos$SD_Control    <- apply(datos[, 2:4], 1, sd)

datos$Media_Condicion <- apply(datos[, 5:7], 1, mean)
datos$SD_Condicion    <- apply(datos[, 5:7], 1, sd)
### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).
# Histogramas de las medias en una sola imagen (2 gráficos lado a lado)
par(mfrow = c(1, 2))

hist(datos$Media_Control, 
     main = "Media del Control", 
     xlab = "Expresión (TPM)", 
     col = "lightblue", 
     border = "darkblue")

hist(datos$Media_Condicion, 
     main = "Media de la Condición", 
     xlab = "Expresión (TPM)", 
     col = "orange", 
     border = "brown")

par(mfrow = c(1, 1))
### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)
# Boxplot comparativo con leyenda
boxplot(datos$Media_Control, datos$Media_Condicion,
        names = c("Control", "Condición"),
        col = c("lightblue", "orange"),
        main = "Comparación de Expresión entre Control y Condición",
        ylab = "Valores de TPM")

legend("topleft", 
       legend = c("Control", "Condición"), 
       fill = c("lightblue", "orange"))
### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.
# Test de Shapiro-Wilk para evaluar la normalidad
shapiro.test(datos$Media_Control)
shapiro.test(datos$Media_Condicion)

# Test de comparativa de medias no paramétrico (Wilcoxon) para cada gen
pvalues <- apply(datos, 1, function(fila) {
  ctrl <- as.numeric(fila[2:4])
  cond <- as.numeric(fila[5:7])
  test <- wilcox.test(ctrl, cond, exact = FALSE)
  return(test$p.value)
})

# Guardar el p-value en una nueva columna de la tabla original
datos$pvalue <- pvalues
### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)
# Aplicamos el corte de p-value < 0.1 indicado por la profesora
# Justificación: Con p-value < 0.05 apenas se obtienen genes debido al bajo tamaño muestral (3 réplicas).
# Siguiendo el aviso de la profesora, ajustamos el corte a p-value < 0.1, obteniendo 462 genes significativos (226 sobreexpresados y 236 reprimidos).
corte_pvalue <- 0.1
deg <- datos[!is.na(datos$pvalue) & datos$pvalue < corte_pvalue, ]

cat("Genes significativos con p-value <", corte_pvalue, ":", nrow(deg), "\n")

# Conteo de genes sobreexpresados y reprimidos
sobreexpresados <- sum(deg$Media_Condicion > deg$Media_Control)
reprimidos      <- sum(deg$Media_Condicion < deg$Media_Control)

cat("Genes sobreexpresados en Condición:", sobreexpresados, "\n")
cat("Genes reprimidos en Condición:", reprimidos, "\n")

# Listado de los nombres de los genes diferencialmente expresados
deg$Gene.Name
### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)
# Construcción de la matriz con las réplicas de los genes significativos
matriz_heatmap <- as.matrix(deg[, 2:7])
rownames(matriz_heatmap) <- deg$Gene.Name

# Generación del Heatmap
heatmap(matriz_heatmap, 
        scale = "row", 
        col = heat.colors(256),
        main = "Heatmap de Genes Diferencialmente Expresados (p < 0.1)")
