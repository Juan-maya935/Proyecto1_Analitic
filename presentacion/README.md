# Presentación — sustentación del Proyecto 1

20 diapositivas (16:9) construidas a partir del informe (`../Proyecto 1 AnalíticaDeDatos.pdf`) y de las figuras generadas por `R/proyecto1.R` (`../resultados/`).

| Archivo | Uso |
|---|---|
| `Proyecto1_diapositivas.pdf` | Deck completo, 20 páginas 1280×720. Para Moodle o presentar sin internet. |
| `pdf_por_hoja/diapo_NN.pdf` | Cada diapositiva como PDF independiente. |
| `index.html` | Versión interactiva: abrir en Chrome, F11, navegar con ← → / espacio. |
| `diapositivas.html` | Misma página sin envoltorio `<html>` (fuente para publicación web). |
| `img/`, `fonts/` | Figuras copiadas de `resultados/` y tipografías embebidas (Archivo, Source Serif 4, JetBrains Mono). |

Regenerar el PDF desde el HTML:

```bash
google-chrome --headless=new --disable-gpu --no-pdf-header-footer \
  --virtual-time-budget=12000 --allow-file-access-from-files \
  --print-to-pdf=Proyecto1_diapositivas.pdf "file://$PWD/index.html"
pdfseparate Proyecto1_diapositivas.pdf pdf_por_hoja/diapo_%02d.pdf
```
