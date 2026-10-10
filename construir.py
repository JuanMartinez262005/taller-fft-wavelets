#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regenera, ejecuta y verifica el cuaderno del taller con sus propios datos.

Uso (desde esta carpeta):
    python construir.py                                   # Hubble de respaldo + ECG local
    python construir.py --imagen data/mi_jwst_i2d.fits    # imagen JWST (FITS o PNG/TIF)
    python construir.py --imagen data/mi_jwst.png --ecg physionet --html

Hace: ensamblar → ejecutar → (cifras a resultados.json) → ensamblar con las discusiones → ejecutar → verificar.
Requiere: pip install -r requirements.txt  (jupyter/nbconvert incluidos)
"""
import argparse, json, shutil, subprocess, sys
from pathlib import Path

import nbformat

for _f in (sys.stdout, sys.stderr):
    try: _f.reconfigure(errors="replace")
    except Exception: pass

NB = "Taller4_FFT_Wavelets.ipynb"
ap = argparse.ArgumentParser()
ap.add_argument("--imagen", default=None)
ap.add_argument("--ecg", default="local", choices=["local", "physionet", "scipy"])
ap.add_argument("--reducir", type=int, default=1, help="promedia la imagen en bloques N×N antes de recortar (p. ej. 3 ó 4 si las estrellas salen muy anchas)")
ap.add_argument("--origen", type=int, nargs=2, metavar=("FILA", "COL"), default=None, help="esquina superior izquierda del recorte (sobre la imagen ya reducida)")
ap.add_argument("--html", action="store_true", help="exportar también a HTML")
a = ap.parse_args()
extra = (["--imagen", a.imagen] if a.imagen else []) + ["--ecg", a.ecg]
if a.reducir > 1: extra += ["--reducir", str(a.reducir)]
if a.origen: extra += ["--origen", str(a.origen[0]), str(a.origen[1])]


def run(cmd):
    print("→", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run(cmd, check=True)


def ejecutar():
    run([sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", "--inplace",
         "--ExecutePreprocessor.timeout=3600", NB])


run([sys.executable, "build_notebook.py", *extra]); ejecutar()                  # pasada 1: calcula las cifras
shutil.copy("resultados.json", "resultados_pasada1.json")
run([sys.executable, "build_notebook.py", "--final", *extra]); ejecutar()       # pasada 2: añade las discusiones con esas cifras

nb = nbformat.read(NB, as_version=4)
errores = [i for i, c in enumerate(nb.cells) if c.cell_type == "code" and any(o.output_type == "error" for o in c.get("outputs", []))]
leer = lambda f: json.loads(Path(f).read_text(encoding="utf-8"))
igual = leer("resultados.json") == leer("resultados_pasada1.json")
avisos = [c.source for c in nb.cells if c.cell_type == "markdown" and c.source.startswith("> ⚠️")]
print("\n================ RESULTADO ================")
print(f"Celdas: {len(nb.cells)} | con error: {errores or 'ninguna'} | cifras idénticas entre pasadas: {igual}")
if avisos:
    print("\n⚠️  Hay afirmaciones del texto que NO se cumplen con sus datos (aparecen en el cuaderno como aviso):")
    for t in avisos: print(t)
else:
    print("Las afirmaciones cualitativas de las discusiones se cumplen con estos datos.")
if a.html:
    run([sys.executable, "-m", "jupyter", "nbconvert", "--to", "html", NB])
print("Listo:", NB)
sys.exit(1 if errores or not igual else 0)
