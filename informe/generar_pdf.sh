#!/usr/bin/env bash
# Genera el PDF del informe a partir de informe.html.
#
# No se usa R Markdown porque el enunciado lo prohibe explicitamente. El HTML
# se imprime con el motor de Chrome, que respeta las reglas @page de CSS: A4,
# una sola columna y margenes fijos.
#
#   uso:  ./generar_pdf.sh
set -euo pipefail

AQUI="$(cd "$(dirname "$0")" && pwd)"
RAIZ="$(dirname "$AQUI")"
SALIDA="$RAIZ/Proyecto 1 AnalíticaDeDatos.pdf"

NAVEGADOR=""
for c in google-chrome chromium chromium-browser brave-browser; do
  if command -v "$c" >/dev/null 2>&1; then NAVEGADOR="$c"; break; fi
done
[ -n "$NAVEGADOR" ] || { echo "No se encontro Chrome/Chromium." >&2; exit 1; }

"$NAVEGADOR" --headless=new --disable-gpu --no-pdf-header-footer \
  --virtual-time-budget=25000 \
  --print-to-pdf="$SALIDA" "file://$AQUI/informe.html"

PAGS=$(pdfinfo "$SALIDA" 2>/dev/null | awk '/^Pages:/{print $2}')
echo "PDF generado: $SALIDA  (${PAGS:-?} paginas)"
[ "${PAGS:-0}" -le 10 ] || echo "AVISO: el enunciado permite maximo 10 paginas."
