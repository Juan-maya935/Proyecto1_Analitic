# Proyecto 1 — Analítica de Datos

**Modelamiento de la precipitación semanal en el Valle del Cauca y evaluación del papel de la correlación espacial.**

Universidad Autónoma de Occidente · Ingeniería de Datos e Inteligencia Artificial
Asignatura: **Analítica de Datos** · Semestre **2026-2S**
Profesor: **Johann A. Ospina** — jaospina@uao.edu.co

| | |
|---|---|
| **Entrega** | domingo **13 de septiembre de 2026, 23:59** |
| **Plataforma** | Moodle — [Project 1](https://campus.uaovirtual.edu.co/mod/assign/view.php?id=529073) |
| **Peso** | 16,67 % de la nota del curso |
| **Modalidad** | Grupal, con sustentación |
| **Lenguaje** | R (obligatorio, 100 %) |

---

## 1. Integrantes

> ⚠️ **Obligatorio:** el enunciado exige que *todos* los integrantes estén relacionados en el documento entregado. Completar esta tabla antes de la entrega.

| Nombre completo | Código | Correo institucional |
|---|---|---|
| César Armando Reyes Oliveros | 2236379 | |
| Yesenia Díaz | 2231783 | |
| Juan Pablo Maya | 2236377 | |

---

## 2. El problema

El departamento del Valle del Cauca combina zonas con regímenes climáticos marcadamente distintos: la llanura del valle geográfico, la costa Pacífica y las cordilleras Occidental y Central. Esa heterogeneidad lo convierte en un caso de estudio adecuado para **evaluar el papel de la correlación espacial** en el modelamiento de variables ambientales.

**Objetivo:** modelar la **precipitación semanal** en el Valle del Cauca en función de las covariables ambientales disponibles.

**Área de estudio:** el límite administrativo oficial del departamento (fuente **GADM**). Todos los rasters ya vienen recortados y enmascarados a ese límite y comparten la misma grilla espacial.

### Alcance exigido

El grupo debe **definir y justificar con desarrollo teórico paso a paso**, usando las metodologías vistas en clase, y cubrir:

1. **Análisis Exploratorio de Datos**
2. **Modelación**
3. **Validación**
4. **Predicción**

---

## 3. Datos

### 3.1 Tabla oficial del enunciado

| Variable | Fuente | Resolución temporal | Periodo | Archivo |
|---|---|---|---|---|
| Precipitación | CHIRPS | Semanal (52 semanas/año) | 2010–2025 | `chirps_semanal_valle_<año>.tif` |
| Temperatura a 2 metros | NASA POWER | Semanal (52 semanas/año) | 2010–2025 | `power_semanal_valle_<año>.tif` |
| Radiación solar incidente | NASA POWER | Semanal (52 semanas/año) | 2010–2025 | `power_semanal_valle_<año>.tif` |
| Altitud | SRTM | Estática | — | `altitud_valle.tif` |
| Climatología semanal de precipitación | CHIRPS | Semanal (promedio histórico) | 2010–2025 | `chirps_climatologia_semanal_valle.tif` |
| Climatología semanal de temperatura | NASA POWER | Semanal (promedio histórico) | 2010–2025 | `power_temp_climatologia_semanal_valle.tif` |
| Climatología semanal de radiación | NASA POWER | Semanal (promedio histórico) | 2010–2025 | `power_radiacion_climatologia_semanal_valle.tif` |

> Cada archivo `.tif` contiene **52 bandas, una por semana ISO** del año correspondiente. Temperatura y radiación están remuestreadas a la resolución espacial de CHIRPS, por lo que **todas las capas son directamente comparables píxel a píxel**.

### 3.2 Geometría de la grilla (verificada sobre los archivos)

| Propiedad | Valor |
|---|---|
| CRS | WGS 84 (EPSG:4326) |
| Tamaño de píxel | 0,05° × 0,05° (≈ 5,5 km) |
| Dimensiones | 37 columnas × 39 filas = **1 443 píxeles/banda** |
| Origen (esquina sup. izq.) | lon −77,55 · lat 5,00 |
| Extensión | lon [−77,55, −75,70] · lat [3,05, 5,00] |
| Bandas por archivo anual | 52 (semanas ISO) |

### 3.3 Contenido de `data/datos_proyecto_1.zip`

17,8 MB comprimidos → **229 archivos, 22,9 MB** descomprimidos.

```
datos_proyecto_1/
├── imagenes_semanales/          # 36 archivos — los que pide el enunciado
│   ├── chirps_semanal_valle_2010.tif … _2025.tif        (16 archivos, 52 bandas c/u)
│   ├── power_semanal_valle_2010.tif … _2025.tif         (16 archivos, 52 bandas c/u)
│   ├── altitud_valle.tif                                (estático, SRTM)
│   ├── chirps_climatologia_semanal_valle.tif
│   ├── power_temp_climatologia_semanal_valle.tif
│   └── power_radiacion_climatologia_semanal_valle.tif
└── cache_raster/                # 193 archivos — insumos crudos
    ├── chirps_diario_2010_01.tif … chirps_diario_2025_12.tif   (192 mensuales de datos diarios)
    └── gadm41_COL_1_pk.rds      # límite administrativo GADM nivel 1 (Colombia)
```

`cache_raster/` es material de origen: **para el trabajo se usa `imagenes_semanales/`**. El `.rds` de GADM es el que sirve para dibujar y enmascarar el límite del departamento en R.

### 3.4 Cómo descomprimir

```bash
unzip data/datos_proyecto_1.zip -d data/
# → data/datos_proyecto_1/imagenes_semanales/...
```

Los `.tif` descomprimidos están en `.gitignore`: **no se versionan**, se regeneran del zip.

---

## 4. Condiciones de entrega (del enunciado)

**Formato**

- Documento en **PDF**, con asunto o nombre de archivo: **`Proyecto 1 AnalíticaDeDatos`**.
- **Máximo 10 páginas**, incluyendo las referencias bibliográficas.
- **Una sola columna.** No se permite formato a dos columnas.
- ❌ No se aceptan documentos escritos a mano, trabajos desarrollados **exclusivamente en Google Colab**, ni documentos elaborados en **R Markdown**.

**Código**

- Adjuntar **un único archivo de R, en archivo independiente**, con todo el código, **completamente reproducible**.
- El desarrollo debe realizarse **completamente en R**.
- Utilizar **únicamente los scripts, procedimientos y metodologías trabajados en clase**.
- Librerías adicionales permitidas **solo** para: elaboración de gráficos y carga de mapas o archivos de Excel. Para lo demás, herramientas y funciones vistas en clase.

**Resultados**

- Presentar resultados mediante **tablas, gráficos e indicadores** que respondan de forma clara y fundamentada a los planteamientos.
- **Todos los resultados deben interpretarse** según el contexto del problema. *No es suficiente con presentar tablas, gráficos o indicadores sin su respectiva interpretación.*

**Sustentación**

- Todos los integrantes deben estar **relacionados en el documento**.
- Se selecciona **aleatoriamente** a un integrante para sustentar. Si no está presente sin excusa debidamente justificada, **su calificación es 0.0** y se sortea a otro integrante.

**IA generativa**

- Se espera desarrollo autónomo, sin depender de IA generativa.
- Si se usa, hay que adjuntar **un documento de Word independiente** con **los prompts utilizados** e indicando **claramente la herramienta** (ChatGPT, Claude, Gemini, DeepSeek, etc.).

**Plazo y canal**

- ❌ No se aceptan trabajos después de la fecha y hora establecidas.
- ❌ No se aceptan trabajos enviados por medios distintos a la plataforma institucional.

---

## 5. Estructura del repositorio

```
.
├── README.md                                    # este archivo
├── Cuaderno_Proyecto1.ipynb                     # bitácora: qué se hizo y por qué
├── Proyecto 1 AnalíticaDeDatos.pdf              # ← entregable 1: el documento
├── .gitignore
├── docs/
│   └── Proyecto1_AnaliticaDeDatos_2026-2S.pdf   # enunciado original del profesor
├── data/
│   └── datos_proyecto_1.zip                     # dataset original, sin modificar
├── R/
│   └── proyecto1.R                              # ← entregable 2: el .R reproducible
├── informe/
│   ├── informe.html                             # fuente del PDF
│   └── generar_pdf.sh                           # lo imprime a PDF
└── resultados/                                  # figuras y tablas generadas por el script
```

**Regenerar el PDF** después de tocar el texto o las figuras:

```bash
Rscript R/proyecto1.R      # si cambiaron las figuras
informe/generar_pdf.sh     # reimprime "Proyecto 1 AnalíticaDeDatos.pdf"
```

El documento **no** está hecho en R Markdown, que el enunciado prohíbe: es HTML
con reglas `@page` de CSS, impreso por el motor del navegador. Una sola columna,
A4, 9 páginas de las 10 permitidas.

**Regla:** `docs/` y `data/` no se tocan — son los originales tal como los entregó el profesor. Todo lo que produzcamos vive en `R/` y `resultados/`.

---

## 6. Entregables

- [x] `Proyecto 1 AnalíticaDeDatos.pdf` — **9 páginas**, una columna, listo en la raíz
- [x] `proyecto1.R` — script único y reproducible (`R/proyecto1.R`), verificado desde un clon limpio
- [x] Todos los integrantes listados en el documento
- [ ] *(solo si se usó IA)* `prompts_IA.docx` — prompts y herramienta empleada
- [ ] Correos institucionales en la tabla de integrantes (sección 1)
- [ ] Resolver con el grupo qué versión de `R/proyecto1.R` queda en `main`
- [ ] Subido a Moodle antes del **13-sep-2026 23:59**

---

## 7. Material de clase de referencia

Insumos del curso aplicables a este proyecto (geoestadística y datos espaciales):

- `Practic_Geostatistics.R`
- `Practice Geostatistics/Ejemplo1–4_Geoestadística.R`
- `Script_Spatial_Analytics.R`, `Script_Class5.R`
- `Lecture_SpatialStatistics.pdf`
- `Class_5_Handle_Spatial_Data.pdf`
- `Patrones_Puntuales.pdf`, `Presentacion_Modelling_ppp.pdf`
- `Geoestadística_Datos_Espaciales_Espacio-Temporales_y_Funcionales.pdf`
- `The_Climate_Hazards_Infrared_Precipitation…pdf` (artículo fuente de CHIRPS)

Además hay **dos asesorías grabadas** del curso, del **3 y 4 de septiembre de 2026**, dedicadas a este proyecto.
