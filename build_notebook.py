# -*- coding: utf-8 -*-
"""Ensambla y ejecuta el notebook del taller. Uso: python build_notebook.py [--final]"""
import json, re, sys
from pathlib import Path
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

from nb_ej1 import EJ1
from nb_ej2 import EJ2

import argparse
_ap = argparse.ArgumentParser()
_ap.add_argument("--final", action="store_true")
_ap.add_argument("--imagen", default=None, help="ruta a un FITS (JWST) o PNG/TIF para el Ejercicio 2")
_ap.add_argument("--ecg", default="local", choices=["local", "physionet", "scipy"])
_ap.add_argument("--reducir", type=int, default=1, help="promedia la imagen en bloques N×N antes de recortar (p. ej. 3 o 4 si las estrellas salen muy anchas)")
_ap.add_argument("--origen", type=int, nargs=2, metavar=("FILA", "COL"), default=None, help="esquina superior izquierda del recorte (sobre la imagen ya reducida)")
ARGS = _ap.parse_args()
FINAL = ARGS.final
R = json.loads(Path("resultados.json").read_text(encoding="utf-8")) if (FINAL and Path("resultados.json").exists()) else {}

cells = []
md = lambda s: cells.append(("md", s))
code = lambda s: cells.append(("code", s))

# ───────────────────────────── Portada ─────────────────────────────
if ARGS.imagen:
    fila_img = f"Imagen proporcionada: `{Path(ARGS.imagen).name}` (ver la celda de carga)"
    par_img = "(ii) en el Ejercicio 2 se usa la imagen indicada en `RUTA_IMAGEN` (la celda de carga describe la lectura y el preprocesado)."
    if ARGS.reducir > 1:
        par_img = par_img[:-1] + f" Antes de recortar se promedió la imagen en bloques de {ARGS.reducir}×{ARGS.reducir} píxeles (`REDUCIR`), para que las estrellas tengan un ancho comparable a la escala del filtro del enunciado."
else:
    fila_img = "Campo profundo del Hubble (**sustituto**; el enunciado pide JWST)"
    par_img = ("(ii) en el Ejercicio 2 se usa, como **sustituto**, el *Hubble eXtreme Deep Field* de *scikit-image*; **no es una imagen de JWST**, que es lo que pide el enunciado. "
               "Para usar JWST se asigna su ruta (FITS `*_i2d.fits` o PNG/TIF) a `RUTA_IMAGEN` en la celda del Ejercicio 2.")
if ARGS.ecg == "physionet":
    fila_ecg = "ECG real, MIT-BIH registro 208 (descarga oficial de PhysioNet)"
    par_ecg = "(i) el ECG es el registro 208 de la MIT-BIH, descargado del repositorio oficial de PhysioNet con `wfdb`;"
else:
    fila_ecg = "ECG real, MIT-BIH registro 208 (fragmento de 5 min incluido en SciPy)"
    par_ecg = "(i) el ECG es el registro 208 de la MIT-BIH tal como lo empaqueta SciPy (360 Hz, 5 min); para usar la descarga oficial de PhysioNet se pone `ECG_FUENTE = \"physionet\"`;"

md(f'''# Taller: Reducción de dimensionalidad, análisis espectral y transformadas integrales
**Asignatura:** Métodos de Machine Learning en Ciencias e Ingenierías (106023C) · **Semestre:** 2026A
**Temas:** SVD, PCA, Fourier (FFT) y Wavelets (DWT)
**Integrantes:** Juan Esteban Martinez Herrera · Juan David Cuartas
**Instituciones:** Universidad del Quindío · Universidad del Valle

| Ejercicio | Contenido | Datos |
|---|---|---|
| 1 | Compresión y reconstrucción 1D: FFT vs. DWT (en parejas) | {fila_ecg} |
| 2 | Aislamiento de fuentes puntuales 2D: 2D-DWT vs. 2D-FFT | {fila_img} |
| 3 | Mapa mental sistémico (5 nodos) | — |

**Cómo leer este cuaderno.** Todo el código se ejecuta de principio a fin sin intervención (no hay aleatoriedad). Las cifras citadas en las discusiones salen de las celdas de código anteriores; los resultados completos se guardan en `resultados.json`.

**Datos.** {par_ecg} {par_img}''')

code(r'''import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle
import pywt
from scipy.fft import fft, ifft, fft2, ifft2, fftshift, ifftshift, fftfreq
from scipy.signal import find_peaks
from scipy import ndimage as ndi
from skimage.color import rgb2gray
from skimage.feature import peak_local_max

plt.rcParams.update({"figure.dpi": 100, "savefig.dpi": 150, "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False})
DATA, FIG = Path("data"), Path("figuras"); FIG.mkdir(exist_ok=True)
RES = {}                                    # cifras clave del taller → resultados.json
print("NumPy", np.__version__, "| PyWavelets", pywt.__version__)''')

# ───────────────────────────── Ejercicio 1 ─────────────────────────────
cells.extend(EJ1)
if FINAL:
    from disc_ej1 import disc1
    from verificar import afirmaciones_fallidas
    _fallos = afirmaciones_fallidas(R)
    def _aviso(lista):
        if lista:
            md("> ⚠️ **Revisar antes de entregar.** Con los datos usados, estas afirmaciones del texto de discusión **no se cumplen**; hay que reescribirlas a la vista de las tablas y figuras:\n>\n" + "\n".join(f"> * {t}" for t in lista))
    _aviso([t for t in _fallos if t.startswith("Ej. 1")])
    md(disc1(R))

# ───────────────────────────── Ejercicio 2 ─────────────────────────────
cells.extend(EJ2)
if FINAL:
    from disc_ej2 import disc2A, disc2B, disc2C
    _aviso([t for t in _fallos if t.startswith(("Ej. 2", "Pregunta"))])
    md(disc2A(R)); md(disc2B(R)); md(disc2C(R))

# ───────────────────────────── Ejercicio 3 ─────────────────────────────
md(r'''---
# Ejercicio 3 · Mapa mental sistémico
Mapa conceptual jerárquico (raíz → 5 nodos temáticos → subconceptos) con **8 enlaces explícitos y etiquetados entre los cinco nodos** (flechas discontinuas). Se dibuja con *matplotlib* con diseño manual; el motor de disposición automática de Mermaid escalonaba los cinco nodos al ligarlos entre sí y perdía la jerarquía.''')

mapa_src = Path("ej3_mapa.py").read_text(encoding="utf-8")
mapa_src = mapa_src.replace('matplotlib.use("Agg")\n', "")
mapa_src = mapa_src.split('if __name__ == "__main__":')[0].rstrip() + "\n\n\nfig_mapa = dibujar_mapa(FIG / \"mapa_mental\", mostrar=True)\n"
code(mapa_src)

md(r'''**Enlaces entre los nodos (lectura del mapa).**

| Enlace | Relación | Justificación |
|---|---|---|
| ① → ② | la proyección elige un subespacio = elegir una base | Proyectar ortogonalmente sobre un subespacio de dimensión $k$ equivale a fijar una base ortonormal y quedarse con $k$ coordenadas; la base puede aprenderse de los datos (SVD/PCA) o ser fija (Fourier, DCT, wavelets). |
| ② → ③ | soporte global vs. compacto | Entre las bases fijas, Fourier tiene soporte global (frecuencia exacta, sin información temporal) y las wavelets tienen soporte compacto (compromiso tiempo–frecuencia, multirresolución). |
| ③ → ④ | localización ⇒ sparsity | Si la base está localizada como la señal (transitorios, bordes, fuentes puntuales), pocos coeficientes concentran casi toda la energía y el truncamiento es eficaz. |
| ② → ④ | base adaptada ⇒ decaimiento rápido | La base adaptada (KLT/SVD) concentra la energía de forma óptima: $\sigma_i$ y $\lvert c_k\rvert$ decaen lo más rápido posible (Eckart–Young). |
| ④ → ① | $k$ por varianza explicada | El número de componentes $k$ se fija con la varianza explicada acumulada (p. ej. ≥ 95 %), cerrando el ciclo de la reducción de dimensionalidad. |
| ④ → ⑤ | el umbral $\lambda$ | El umbral (duro/suave) decide cuánto ruido se elimina (denoising) y cuánto se comprime; con $L_1$ se promueve sparsity (el umbral suave es el operador proximal de $\lambda\lVert\cdot\rVert_1$). |
| ① → ⑤ | SVD/PCA ⇒ eigenfaces, POD | La SVD/PCA da las características (eigenfaces) y los modos POD; DMD añade dinámica sobre la SVD. |
| ③ → ⑤ | señales no estacionarias | La localización tiempo–frecuencia es la que permite analizar y comprimir señales no estacionarias (ECG, ondas gravitacionales, espectros) e imágenes. |''')

# ───────────────────────────── Cierre ─────────────────────────────
if FINAL:
    from disc_ej1 import conclusiones
    md(conclusiones(R))

md(r'''---
### Referencias
* Moody, G. B. y Mark, R. G. (2001). *The impact of the MIT-BIH Arrhythmia Database*. IEEE Eng. in Medicine and Biology 20(3), 45–50. · Goldberger, A. L. et al. (2000). *PhysioBank, PhysioToolkit, and PhysioNet*. Circulation 101(23), e215–e220.
* Mallat, S. (2009). *A Wavelet Tour of Signal Processing* (3.ª ed.). Academic Press.
* Brunton, S. L. y Kutz, J. N. (2019). *Data-Driven Science and Engineering*. Cambridge University Press.
* Donoho, D. L. y Johnstone, I. M. (1994). *Ideal spatial adaptation by wavelet shrinkage*. Biometrika 81(3), 425–455.
* Material del curso: *Overview on the Fourier Analysis* (A. Argüelles y C. Galíndez, Quantum Optics Group) y *2D transforms and image processing*.''')

code(r'''Path("resultados.json").write_text(json.dumps(RES, indent=1, ensure_ascii=False, default=float), encoding="utf-8")   # UTF-8 explícito (Windows usa cp1252 por defecto)
print("Resultados guardados en resultados.json:", ", ".join(sorted(k.split("_")[0] for k in RES)) and f"{len(RES)} entradas")''')

# ───────────────────────────── Parámetros de datos ─────────────────────────────
def _inyectar(tipo, txt):
    if tipo != "code": return txt
    if ARGS.imagen:
        txt = txt.replace("RUTA_IMAGEN = None ", f'RUTA_IMAGEN = {json.dumps(ARGS.imagen)} ', 1)
    if ARGS.reducir > 1:
        txt = txt.replace("REDUCIR = 1 ", f"REDUCIR = {ARGS.reducir} ", 1)
    if ARGS.origen:
        txt = txt.replace("ORIGEN_RECORTE = None ", f"ORIGEN_RECORTE = ({ARGS.origen[0]}, {ARGS.origen[1]}) ", 1)
    txt = txt.replace('ECG_FUENTE = "local"', f'ECG_FUENTE = "{ARGS.ecg}"', 1)
    return txt
cells = [(t, _inyectar(t, s)) for t, s in cells]

# ───────────────────────────── Construcción y ejecución ─────────────────────────────
nb = new_notebook()
for tipo, txt in cells:
    nb.cells.append(new_markdown_cell(txt) if tipo == "md" else new_code_cell(txt))
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb.metadata["language_info"] = {"name": "python"}
nbformat.write(nb, "Taller4_FFT_Wavelets.ipynb")
print(f"Notebook escrito: {len(nb.cells)} celdas ({'FINAL' if FINAL else 'pasada 1, sin discusiones'})" + (f" | avisos: {len(_fallos)}" if FINAL else ""))
