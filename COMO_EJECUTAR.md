# Taller 4 — cómo ejecutarlo con sus propios datos

## 1. Instalar (una vez)
```
pip install -r requirements.txt
```

## 2. Datos que pide la actividad
**Ejercicio 2 — imagen del JWST (≥ 512 × 512 px, FITS o monocromática).** El enunciado exige JWST; el cuaderno trae el Hubble XDF sólo como respaldo para comprobar que el código corre.
1. Descarguen una imagen pública de campo profundo del JWST (por ejemplo, desde el portal MAST del STScI, https://mast.stsci.edu, buscando una observación del JWST/NIRCam y bajando el producto de imagen combinada `*_i2d.fits`; o una imagen en PNG/TIF desde la galería de Webb). Los datos públicos no necesitan cuenta. Si el profesor indicó un enlace concreto para «NASA / JWST High-Res Deep Field Images», usen ese.
2. Guárdenla en la carpeta `data/`.
3. Ejecuten (esto regenera el cuaderno **y los números de las discusiones** con su imagen):
   ```
   python construir.py --imagen data/su_imagen_i2d.fits
   ```
   Si la imagen mide más de 2048 px por lado se usa un recorte de 2048 × 2048 (el centrado, por defecto). **Para elegir la región o reducir la imagen se usan opciones de `construir.py`, no se edita la celda del cuaderno** (se regenera cada vez y se perdería el cambio):
   ```
   python construir.py --imagen data/su_imagen.tif --reducir 3          # promedia en bloques 3×3 antes de recortar
   python construir.py --imagen data/su_imagen.tif --origen 5000 4000   # esquina (fila, columna) del recorte, sobre la imagen ya reducida
   ```
   `--reducir N` conviene cuando las estrellas salen anchas: el cuaderno avisa si el FWHM de las fuentes supera ≈ 4.5 px (una estrella a la escala de píxel natural mide 1–3 px) y sugiere el factor. El filtro del enunciado (corte en 1/8 ciclos/px, L = 2) está pensado para fuentes de ese tamaño.

**Ejercicio 1 — ECG oficial de PhysioNet (opcional).** Por defecto se usa una copia del registro 208 de la MIT-BIH incluida en SciPy. Para la descarga oficial:
```
python construir.py --imagen data/su_imagen_i2d.fits --ecg physionet
```
(requiere internet; los datos de PhysioNet MIT-BIH son abiertos). Si falla la descarga, vuelvan a `--ecg local`.

## 3. Qué hacer después (importante)
`construir.py` actualiza **todas las cifras** del texto, pero **no puede reescribir los juicios**. Al terminar imprime una lista de las afirmaciones cualitativas que ya no se cumplen con sus datos (también aparecen como un recuadro ⚠️ dentro del cuaderno). Cada afirmación avisada hay que reescribirla a la vista de las tablas y figuras. Si no hay avisos, el texto sigue siendo coherente con sus resultados; revísenlo de todos modos antes de entregar.

Además: pongan los nombres de los integrantes en la portada y añadan `--html` si quieren también la versión HTML.

## 4. Archivos
* `Taller4_FFT_Wavelets.ipynb` — cuaderno ejecutado.
* `construir.py` — regenera, ejecuta y verifica el cuaderno.
* `build_notebook.py`, `nb_ej1.py`, `nb_ej2.py`, `disc_ej1.py`, `disc_ej2.py`, `disc_util.py`, `verificar.py`, `ej3_mapa.py` — generadores (código y textos).
* `data/`, `figuras/`, `resultados.json` — datos, figuras y cifras de la última ejecución.

## 5. Problemas frecuentes
* **`ValueError: <COMPRESSION.LZW: 5> requires the 'imagecodecs' package`** (TIFF comprimido): ejecute `pip install imagecodecs` y vuelva a lanzar `construir.py`. La versión actual de `nb_ej2.py` además intenta leerlo con Pillow si falta ese paquete.
* **`RuntimeError: No se encontraron fuentes compactas y aisladas…`** (imágenes de color estiradas, muy sobremuestreadas o con nebulosa intensa, como IC 348): **corregido** en la versión actual de `nb_ej2.py`. Si hay menos de 10 fuentes que cumplan el criterio estricto (área ≤ 12 px, aisladas), el cuaderno detecta sobre la imagen sin fondo suave y relaja poco a poco el área (hasta 80 px) y el aislamiento (hasta 65 %). Imprime qué criterio usó, lo guarda en `resultados.json`, lo menciona en la discusión de la Pregunta A y añade un aviso ⚠️. **Miren siempre la figura «fuentes de referencia»** (celda posterior a la detección): las 150 marcas deben caer sobre estrellas, no sobre filamentos de la nebulosa. Con una imagen de nebulosa casi todo el «cielo vacío» está en la zona oscura; para un análisis más parecido a un campo profundo, recorten esa zona con `--origen FILA COL`.
* **Aviso «el X % de los píxeles es negro puro»**: la imagen trae bordes vacíos (típico de mosaicos). Falsean la mediana (fondo) y los umbrales; usen `--origen FILA COL` para quedarse con una zona con señal.
* **`UnicodeEncodeError: 'charmap' codec can't encode character '\u2192'`** (Windows, al guardar `resultados.json`): Windows escribe por defecto en cp1252. **Corregido**: `build_notebook.py` y `construir.py` leen y escriben ahora siempre en UTF-8. Si les aparece algo parecido en otro punto, ejecuten antes `$env:PYTHONUTF8=1` (PowerShell) y repitan.
* **`RuntimeWarning: Proactor event loop does not implement add_reader...`** (Windows): es sólo un aviso de Jupyter; no impide la ejecución.
* **Memoria o lentitud con imágenes grandes**: baje `MAX_LADO` (p. ej. 1024) en la celda de carga del Ejercicio 2.
* Si `construir.py` falla en una celda, el mensaje muestra el código de esa celda y el error; la celda de carga del Ejercicio 2 es la primera que lee su imagen.
