# Guion de sustentación — Proyecto 1

Duración objetivo: **18–20 minutos** (≈ 1 min por diapo, los divisores son de 20 segundos). Lo que está en *cursiva* es lo que se dice; lo que está en **Por qué** es la intención de la diapo, para que quien la presente sepa qué defender si preguntan.

Sugerencia de reparto: César → diapos 1–7 (problema, datos, método), Yesenia → 8–13 (modelo y correlación espacial), Juan Pablo → 14–20 (validación, predicción, conclusiones). Cada quien cambia en un divisor, así el cambio de voz coincide con el cambio de parte.

---

## 1 · Portada

*Buenas. Somos César, Yesenia y Juan Pablo. Nuestro proyecto es sobre la lluvia semanal en el Valle del Cauca, pero la pregunta de fondo no es "cuánto llueve", sino algo un poco más interesante: cuando uno ya tiene una regresión decente para explicar la lluvia, ¿cuánto más se gana si además le pone atención a dónde están los puntos, o sea, a la correlación espacial? Eso fue lo que medimos.*

*Todo está hecho en R base; el único paquete que usamos es terra, para leer los mapas.*

**Por qué:** en 20 segundos se planta la pregunta del enunciado con palabras simples. El mapa de la derecha ya muestra el contraste (>200 mm en la costa vs 30–50 en el valle) sin que haya que explicarlo todavía. Mencionar R base y terra desde el inicio evita la pregunta "¿usaron gstat?" más adelante.

---

## 2 · Problema

*Miren el mapa. Eso es una sola semana de lluvia. En la punta de la izquierda, hacia el Pacífico, caen más de 200 milímetros. Ciento cincuenta kilómetros hacia el oriente, en el valle, caen entre 30 y 50. Es cinco a siete veces menos en una distancia que uno recorre en carro en dos horas.*

*¿Por qué pasa eso? La masa de aire húmedo viene del Pacífico, choca con la cordillera Occidental y descarga ahí. Al valle llega ya seca. Ese contraste tan brutal es lo que hace del Valle un buen caso para nuestra pregunta.*

*La estrategia que seguimos se llama kriging universal, y la idea es muy sencilla: la lluvia en un punto es dos cosas sumadas. Una parte que se explica con "dónde estoy" —longitud, latitud, altura— y que la modelamos con una regresión. Y una parte que sobra, el residuo, que ya no depende de una fórmula sino de cómo se parece un punto a sus vecinos. Esa segunda parte es la que se modela con kriging.*

**Por qué:** la fórmula Z(s) = m(s) + ε(s) es el esqueleto de toda la presentación. Si el público entiende "tendencia + residuo" aquí, entiende todo lo demás. No hay que leer la fórmula: hay que decir "dos cosas sumadas".

---

## 3 · Divisor I — Datos y método

*Vamos a arrancar por lo menos glamuroso pero más importante: los datos. Porque nos encontramos tres cosas que, si uno no las arregla antes, contaminan todo lo que viene después.*

**Por qué:** un divisor sirve para respirar y para anunciar. Los tres números (688, −76 867, S44) son un adelanto de lo que viene.

---

## 4 · Auditoría de datos

*Primera cosa: los archivos de CHIRPS no marcan los huecos como "sin dato", sino con un número negativo enorme: menos setenta y seis mil. Si uno calcula un promedio sin filtrarlos, sale basura. Lo filtramos con una línea: todo lo negativo se vuelve NA. Después de eso el mínimo pasa a 2,77 milímetros, que ya es un valor físico.*

*Segunda cosa: los rásteres son un rectángulo, y el Valle no es un rectángulo. De los 1 443 píxeles, solo 688 caen dentro del departamento. Ese 688 es el número que usamos como base para todo. Esto parece obvio, pero si uno calcula la cobertura de la altitud contra el rectángulo completo, sale que solo cubre el 47 % y uno la descarta, cuando en realidad cubre el 100 %.*

*Tercera cosa, y esta sí nos sorprendió: el enunciado dice que temperatura y radiación son comparables píxel a píxel con la lluvia. Geométricamente sí, pero en contenido no. NASA POWER tiene una resolución de medio grado, diez veces más gruesa que CHIRPS. Al remuestrear, solo quedaron los centros de las celdas originales: la temperatura tiene dato en el 18 % del área, y la radiación en el 4 %. Para la radiación son literalmente tres números para todo el Valle. Eso no es una variable continua, es un escalón costa-valle-cordillera, que ya lo captura la longitud. Así que las descartamos y trabajamos con altitud y coordenadas, que tienen cobertura completa y no nos obligan a inventar ningún dato.*

**Por qué:** esta diapo es la que muestra criterio. El profesor sabe que las capas de POWER tienen huecos; que el grupo lo haya detectado y justificado la exclusión vale más que haber metido temperatura a la fuerza. Si preguntan "¿por qué no interpolaron los huecos?": porque habría sido arriesgado y habría requerido funciones que no se vieron en clase.

---

## 5 · Elección de la semana

*Teníamos 16 años por 52 semanas. Elegimos trabajar sobre la climatología, o sea el promedio de los 16 años semana a semana, y de esas 52 semanas tomamos la más lluviosa: la semana 44, primeros días de noviembre, con 92 milímetros en promedio. La más seca, la semana 3, tiene 29. Se ve el ciclo bimodal típico de la zona: dos picos, abril-mayo y octubre-noviembre, que es la Zona de Convergencia Intertropical pasando dos veces por año.*

*¿Por qué el promedio y no un año concreto? Porque lo que nos preguntan es la estructura espacial de la lluvia, no cómo llovió en 2015. Si tomamos un año suelto, el residuo mezcla estructura espacial con la anomalía de ese año —El Niño, La Niña— y el rango que después estima el semivariograma deja de significar algo claro.*

**Por qué:** anticipa la pregunta "¿y por qué no usaron los 16 años?". La respuesta corta: promediar quita el ruido interanual y deja limpia la parte espacial, que es la que se modela.

---

## 6 · Método en seis pasos

*Esta figura es el proyecto completo en una sola imagen, y la vamos a recorrer rápido porque cada panel es una sección que viene después.*

*Panel 1: el dato, la lluvia de la semana 44 con 70 estaciones que sorteamos. Panel 2: lo que la regresión logra explicar; se ve que reproduce el gradiente de occidente a oriente pero de forma muy suave. Panel 3: lo que sobra, el residuo. Y fíjense en algo: no es ruido salpicado al azar, son manchas. Manchas rojas donde el modelo se queda corto, manchas azules donde se pasa. Eso es correlación espacial. Panel 4: la medimos con el índice de Moran. Panel 5: con el semivariograma vemos cuánta correlación hay y hasta qué distancia llega, unos 56 kilómetros. Panel 6: con eso, el kriging calcula para cada punto qué peso darle a cada estación; el círculo grande es el vecino más cercano, con peso 0,5.*

*La predicción final es la suma: tendencia más kriging del residuo.*

**Por qué:** es la diapo "mapa de ruta". Si alguien se pierde después, puede volver mentalmente a estos seis paneles. No hay que explicar las ecuaciones de la derecha; están para que se vea que cada paso tiene su fórmula y que está en el informe.

---

## 7 · Análisis exploratorio

*Sorteamos 70 estaciones dentro del departamento, con semilla fija para que sea reproducible, y proyectamos las coordenadas a kilómetros. Eso no es un detalle estético: el semivariograma mide distancia, y en grados un grado de longitud y uno de latitud no miden lo mismo.*

*La lluvia en esas 70 estaciones es muy asimétrica. La media es 91 pero la mediana es 64, un 42 % de diferencia. La cola larga hacia la derecha es la vertiente pacífica. Shapiro-Wilk rechaza normalidad sin ninguna duda.*

*Y aquí viene algo clave. Uno podría pensar "no es normal, hay que sacarle logaritmo". Pero la pregunta correcta no es si la lluvia es normal, sino si lo son los residuales después de quitar la tendencia. Eso lo respondemos en dos diapos.*

*Una cosa más: ¿por qué 70 estaciones y no las 688 celdas? Dos razones. Costo: el kriging invierte una matriz para cada punto que predice y en cada vuelta de la validación. Y sentido: si uso las 688 celdas ya tengo el mapa completo, no queda nada que interpolar. Las 70 estaciones simulan una red de observación real, y así podemos evaluar el kriging contra un valor que sí conocemos.*

**Por qué:** las correlaciones de la tabla siembran la sorpresa de la siguiente diapo: altitud con −0,65 "parece" fuerte. Y la frase "no le sacamos logaritmo" es una de las decisiones que más conviene dejar clara, porque va contra el reflejo habitual.

---

## 8 · Divisor II — Modelo y correlación espacial

*Ahora sí, el modelo. Y les adelanto los tres resultados de esta parte: la regresión explica el 84 % de la varianza, la altitud cambia de signo cuando uno controla la posición, y lo que sobra tiene un Moran de 0,79. Vamos uno por uno.*

**Por qué:** anunciar los tres números antes de mostrarlos hace que el público los esté esperando.

---

## 9 · La tendencia y la confusión

*Probamos dos modelos anidados: solo coordenadas, y coordenadas más altitud. El test F dice que la altitud aporta, con p de 0,03, así que se queda. El modelo final explica el 84,5 % de la varianza.*

*Pero el resultado interesante no es ese. Miren el coeficiente de la altitud: es positivo. Cuando calculamos la correlación simple en la diapo anterior era negativa, −0,65, y uno diría "a más altura, menos lluvia". Resulta que eso era mentira. Lo que esa correlación decía en realidad era "más al oriente, menos lluvia", porque en el Valle la costa es baja y está al occidente, y la cordillera es alta y está al oriente. Altitud y longitud están correlacionadas en 0,75. Cuando uno controla la posición, la altitud pasa a ser positiva: a más altura, más lluvia. Es el efecto orográfico, que es lo que la física dice que debería pasar.*

*El gráfico de la izquierda lo muestra dividiendo por franjas de longitud: en el Pacífico la pendiente es negativa, porque ahí la nube descarga al nivel del mar; en el centro y el oriente es positiva. El gráfico de la derecha, el de variable añadida, muestra la pendiente ya limpia de posición, y es exactamente el coeficiente de la tabla.*

**Por qué:** esta es la diapo con más "aha" de toda la presentación. Es un caso de manual de confusión, y conviene decir la palabra: *confounding*. Si el profesor pregunta por qué no metieron una interacción altitud × longitud: porque se salía de lo visto en clase, y lo declaramos como limitación al final.

---

## 10 · Diagnóstico de residuales

*Aquí cerramos la duda que dejamos abierta. Los residuales del modelo sí son normales: Shapiro da p de 0,17, no se rechaza. La asimetría que veíamos en la lluvia era la tendencia, no la distribución. Quitando la posición, desaparece sola. Por eso no hacía falta ningún logaritmo, y el kriging queda justificado como predictor lineal óptimo.*

*Pero hay un problema. Miren el tercer panel, residual contra ajustado: hace una U. Donde el modelo predice poco, se queda corto; donde predice medio, se pasa; donde predice mucho, se vuelve a quedar corto. Eso es curvatura sistemática: una recta no puede doblarse. Y en el mapa de la derecha se ve lo mismo en el espacio: azul en el centro del valle, donde sobreestima; rojo en los dos flancos, donde subestima.*

*Ese mapa nos dice dos cosas. Que la tendencia necesita un término cuadrático, y que hay correlación espacial en el residuo. Las dos las vamos a atacar.*

**Por qué:** la U y las manchas son la evidencia de que la tendencia lineal en longitud está mal especificada. Esta diapo justifica la cuadrática que aparece más adelante; sin ella, la cuadrática parecería un capricho.

---

## 11 · Índice de Moran

*Esta es la pregunta central del enunciado y se responde con lo que dice el material del curso: aplicar Moran a los residuales de la regresión. Porque mínimos cuadrados supone errores independientes, y si los residuales están correlacionados espacialmente, ese supuesto se viola.*

*Programamos Moran y Geary a mano, con la fórmula de clase, con una matriz de pesos binaria por distancia. Y como el umbral de distancia cambia el resultado, en vez de reportar un número reportamos tres. Con vecinos a 20 kilómetros, el Moran de los residuales es 0,79; bajo aleatoriedad se esperaría −0,014. La significancia la evaluamos por permutación: barajamos los residuales 999 veces y ninguna de las 999 llegó al valor observado. Geary apunta en la misma dirección, así que no es un error de cálculo.*

*En castellano: después de que la regresión se llevó el 84 %, lo que sobra sigue teniendo muchísima estructura espacial. Un modelo de regresión solo deja información sobre la mesa. Ese es el hueco que va a llenar el kriging.*

**Por qué:** es la respuesta directa a "evaluar el papel de la correlación espacial". El diagrama de Moran se lee así: cada punto es una estación, horizontal es su residual, vertical es el promedio de sus vecinos, y la pendiente de la recta es I. Rojos = alto rodeado de alto; azules = bajo rodeado de bajo.

---

## 12 · Correlograma

*Este es el resultado más informativo de esta parte. Un correlograma es el Moran calculado a distintas distancias. Uno esperaría una curva que baja hasta cero y se queda ahí. Aquí no: es una onda. Positiva hasta unos 48 kilómetros, negativa alrededor de los 89, y vuelve a subir hacia los 149.*

*Cada tramo se lee sobre el mapa de manchas. Los primeros 48 kilómetros son dependencia genuina: los vecinos se parecen. Los 89 kilómetros son media anchura del departamento: el centro azul contra un flanco rojo, signos opuestos. Y los 149 son los dos flancos reencontrándose, ambos rojos.*

*Esto tiene una consecuencia técnica: el lóbulo negativo dice que a gran distancia los residuales no son estacionarios, todavía queda tendencia sin modelar. Y los modelos de semivariograma asumen que la curva crece hasta una meseta y se queda ahí. Eso solo se cumple hasta unos 50 kilómetros, así que ahí cortamos el semivariograma. El libro del curso respalda esto: recomienda limitar la distancia máxima cuando hay comportamiento errático a distancias mayores.*

**Por qué:** justifica el corte a 50 km, que de otro modo parecería arbitrario. Y conecta el número con el mapa: la onda no es un artefacto matemático, es la geografía del Valle.

---

## 13 · Semivariograma

*El diagnóstico anterior nos obligó a revisar la tendencia. Llevamos dos en paralelo, la lineal y una con término cuadrático en longitud, y dejamos que decidiera la evidencia.*

*A cada semivariograma empírico le ajustamos los tres modelos de clase: esférico, exponencial y gaussiano, por mínimos cuadrados con optim. Pero el criterio para aceptar un ajuste no es qué tan bonito pega la curva, sino un supuesto teórico: bajo estacionariedad, la meseta debe aproximar la varianza del proceso. Si la meseta supera claramente la varianza, la curva sigue subiendo y eso significa que quedó tendencia sin modelar.*

*Miren la tabla. Con la tendencia lineal, las tres mesetas superan la varianza entre 1,3 y 3 veces: no estacionario. Con la cuadrática, el esférico da una razón de 1,02 y el gaussiano 0,96: coherentes. En el gráfico se ve clarísimo: a la izquierda la curva sube sin parar y cruza la línea punteada de la varianza; a la derecha se aplana justo sobre ella hacia los 40 kilómetros.*

*Así que la cuadrática no es una preferencia nuestra, es lo que los datos piden.*

**Por qué:** el profesor valora que la decisión lineal/cuadrática salga de un criterio del libro (meseta ≈ varianza) y no de "se veía mejor". Los parámetros del modelo final (pepita 0, meseta 161,7, rango 56 km) aparecen en la diapo 16.

---

## 14 · Divisor III — Validación y predicción

*Última parte. Ya tenemos seis combinaciones posibles, dos tendencias por tres modelos. Ahora hay que elegir una, y les adelanto que no elegimos la de menor error. Y después, con la elegida, reconstruimos el mapa completo desde 70 puntos.*

**Por qué:** "no elegimos la de menor error" genera curiosidad y prepara el argumento de los tres filtros.

---

## 15 · Validación cruzada

*Hicimos validación dejando una estación afuera por vez: con las otras 69 reestimamos la regresión y el sistema de kriging, predecimos la que quedó afuera, y repetimos 70 veces. Para las seis combinaciones.*

*Primero, la columna que responde la pregunta del enunciado en términos de predicción: "RMSE sin kriging" contra "con kriging". Con la tendencia cuadrática, pasar de 14 a 8,9 milímetros es una mejora del 37 %, y esa mejora sale únicamente de aprovechar la correlación espacial del residuo. No agregamos ninguna variable.*

*Segundo, dos trampas. El modelo gaussiano ajusta muy bien el semivariograma y sin embargo destroza la predicción: R² negativo, de −8,9. La razón es que su matriz de covarianzas queda casi singular, mil veces peor condicionada que la del esférico, y los pesos explotan. Moraleja: ajustar mejor el variograma no es predecir mejor.*

*La otra trampa: la fila con menor RMSE, 7,09, es lineal con exponencial. Pero su meseta triplica la varianza —no estacionario— y su sd(z) es 0,69, que significa que la varianza de kriging está mal calibrada: el mapa de incertidumbre exageraría la precisión en un 45 %. Ahí el kriging está compensando una tendencia mal especificada; es un mérito prestado.*

*Por eso elegimos con tres filtros a la vez: estacionariedad, calibración y estabilidad. Exactamente una combinación pasa los tres.*

**Por qué:** es la diapo más densa. Si el tiempo aprieta, lo que no se puede omitir es: (1) 37 % de mejora solo por correlación espacial, (2) el menor RMSE no gana, (3) los tres filtros. Sd(z) explicado en una frase: "si la varianza que el kriging dice tener es creíble, el error dividido por esa desviación debería tener desviación 1".

---

## 16 · Modelo seleccionado

*El modelo final es tendencia cuadrática con semivariograma esférico: pepita cero, meseta 161,7 y rango 56 kilómetros, que es la distancia a partir de la cual dos puntos dejan de parecerse. Rinde RMSE de 8,86, MAE de 5,5, R² de validación 0,97 y sd(z) de 1,06.*

*Los diagnósticos de la izquierda: observado contra predicho pegado a la diagonal, residuos de validación sin patrón, y el Q-Q del error estandarizado sobre la recta. Eso último es lo que nos permite confiar en el mapa de incertidumbre que viene después.*

*Y el mapa de la derecha es la confirmación visual: comparen con el mapa de manchas de hace unas diapos. Aquí los colores ya no forman manchas. El kriging absorbió la estructura espacial que la regresión había dejado suelta. Cumplió su función.*

**Por qué:** cierra el arco "manchas → kriging → sin manchas". Es la demostración gráfica de que la correlación espacial se aprovechó.

---

## 17 · Predicción

*Con el modelo elegido predijimos las 688 celdas del departamento y comparamos contra el valor real de CHIRPS.*

*El número que vale es el de las 620 celdas que no tenían estación: RMSE de 10,8 milímetros y R² de 0,95, reconstruyendo el mapa completo a partir de 70 puntos. El rango predicho va de 34,6 a 223 y el observado de 34,9 a 225; la media predicha es 88,9 contra 90,1 observada, un sesgo del 1,4 %.*

*Fíjense que reportamos también las 68 celdas con estación, con RMSE de 2, pero no lo celebramos. Con pepita cero, el kriging está obligado a pasar exacto por los puntos observados. Ese número mide una propiedad del método, no su capacidad de predecir.*

*En los mapas: izquierda el observado, derecha el reconstruido. Se recupera el gradiente, la posición del máximo en el Pacífico y las magnitudes. El estimado es un poco más suave, que es exactamente lo que se espera de un interpolador óptimo en media cuadrática.*

**Por qué:** la honestidad de separar "celdas con estación" de "sin estación" es un punto a favor. Y la frase "más suave, como se espera" anticipa la crítica "el mapa predicho se ve borroso".

---

## 18 · Incertidumbre y descomposición

*Dos mapas para cerrar el análisis. El de la izquierda es la desviación estándar de kriging, y se comporta como manda la teoría: mínima sobre las estaciones, 1,1 milímetros, y máxima en los bordes y en los huecos sin muestrear, 12,6 en el filo oriental y la punta norte. Este mapa solo tiene sentido porque sd(z) dio 1,06; con el modelo de menor RMSE habría mentido.*

*Y el de la derecha es, literalmente, la respuesta gráfica al enunciado. El panel A es lo que aporta la tendencia determinística, el B es lo que agrega el kriging del residuo. Miren el panel B: no es ruido. Tiene estructura espacial organizada, manchas coherentes. Esa estructura es exactamente "el papel de la correlación espacial" que había que evaluar. Es la mitad del mapa que la regresión no sabía dibujar.*

**Por qué:** el panel B es el cierre conceptual: uno puede *ver* la correlación espacial como un mapa. Si solo se pudiera mostrar una figura para responder el enunciado, sería esta.

---

## 19 · Conclusiones

*Para cerrar, en una frase: la correlación espacial no es un adorno del modelo, es la mitad del modelo. La regresión explica el 84 %, lo que sobra tiene un Moran de 0,79, aprovecharlo con kriging baja el error un 37 %, y sobre el mapa completo se reconstruyen 620 celdas nunca vistas con R² de 0,95 desde 70 puntos.*

*Y hay tres cosas que nos parecen más valiosas que las cifras. Uno: la altitud parece secar el territorio cuando se mira sola, y en realidad lo humedece cuando se controla la posición; es un caso de confusión de manual. Dos: ajustar mejor el variograma no es predecir mejor; el gaussiano nos lo demostró. Tres: el menor RMSE no siempre es la mejor decisión; la combinación más precisa era inadmisible por no estacionaria y mal calibrada, y elegirla habría producido un mapa de incertidumbre engañoso.*

**Por qué:** se repiten los cuatro números clave a propósito; es lo que el jurado se lleva. Los tres hallazgos son los que diferencian este trabajo de uno que "solo corrió kriging".

---

## 20 · Limitaciones

*Y cuatro limitaciones que preferimos declarar antes que esconder. Las 688 celdas vienen de un producto grillado, no de estaciones reales, así que nuestro R² mide capacidad de reconstruir el campo CHIRPS, no destreza frente a observaciones independientes. El efecto de la altitud no es homogéneo, negativo en el Pacífico y positivo en el interior; lo correcto sería una interacción altitud por longitud, que omitimos por quedarnos dentro de lo visto en clase. Los residuales no son estacionarios más allá de 50 kilómetros, lo que nos obligó a cortar el variograma y deja fuera la estructura de gran escala. Y el análisis es de una sola semana; en temporada seca los parámetros del variograma bien podrían cambiar.*

*El código es un solo archivo de R, reproducible, en el repositorio. Gracias. Quedamos atentos a preguntas.*

**Por qué:** declarar limitaciones uno mismo desarma preguntas hostiles. Cada una lleva su "por qué" implícito (fuera de lo visto en clase, decisión respaldada por el libro, etc.).

---

## Preguntas probables y respuesta corta

| Pregunta | Respuesta |
|---|---|
| ¿Por qué no usaron gstat / geoR? | El enunciado exige R; todo se programó en base R para que cada paso fuera transparente. terra solo para leer rásteres. |
| ¿Por qué descartaron temperatura y radiación? | Cobertura real 18 % y 4 % dentro del departamento por la resolución de 0,5° de POWER. Radiación = 3 valores. Meterlas obligaba a interpolar huecos. |
| ¿Por qué 70 estaciones y no todas las celdas? | Costo (inversión de matriz por punto y por iteración) y sentido: con las 688 no queda nada que interpolar. Simula una red de observación. |
| ¿Por qué no transformaron la lluvia? | Los residuales del modelo con posición ya son normales (p = 0,17). La asimetría era la tendencia. |
| ¿Por qué cuadrática y no lineal? | Tres evidencias: la U en residual vs ajustado, las manchas en el mapa, y la meseta que no aparece con la lineal (razón meseta/varianza > 1,3 en los tres modelos). |
| ¿Por qué cortar a 50 km? | El correlograma cruza cero ahí y después es negativo: no estacionario a gran distancia. El libro recomienda limitar la distancia en ese caso. |
| ¿Por qué no eligieron el de menor RMSE? | Meseta triplica la varianza y sd(z) = 0,69: mapa de incertidumbre engañoso. El kriging estaría compensando una tendencia mal especificada. |
| ¿Qué es sd(z)? | Desviación del error estandarizado (obs − pred) / σ_kriging. Si la varianza de kriging es creíble, ≈ 1. |
| ¿Isotropía? | Se asumió: el semivariograma se calculó contra distancia, no contra dirección. Está declarado como supuesto. |
| ¿Qué pasaría en otra semana? | Los parámetros del variograma (rango, meseta) probablemente cambian en temporada seca. Está en limitaciones. |
