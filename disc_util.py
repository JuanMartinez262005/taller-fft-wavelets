# -*- coding: utf-8 -*-
"""Rellena plantillas Markdown con cifras de resultados.json: @@ruta>sub>clave|formato@@ (formato opcional)."""
import re


def fill(plantilla, R):
    def rep(m):
        ruta, _, fmt = m.group(1).partition("|")
        v = R
        for k in ruta.split(">"):
            v = v[k]
        return format(v, fmt) if fmt else str(v)
    return re.sub(r"@@(.+?)@@", rep, plantilla)


def lim_img(R):
    """Frase sobre la naturaleza de la imagen usada, según cómo se cargó."""
    if "SUSTITUTO" in str(R.get("ej2_fuente", "")):
        return "imagen de Hubble de 8 bits ya estirada (sustituto: no es JWST ni flujo lineal)"
    return "imagen normalizada por percentiles (0.5–99.95 %): los MSE son relativos, no están en unidades de flujo"


def desc_fuentes(R):
    """Cómo llamar a las fuentes de referencia según el criterio de selección que se aplicó de verdad."""
    c = R.get("ej2_criterio_fuentes")
    if not c or c["modo"] == "estricto":
        return "fuentes compactas y aisladas"
    if c["modo"] == "adaptativo" and not c.get("relajado", True):
        return "fuentes compactas y aisladas (detectadas sobre la imagen residual, por ser el fondo estructurado)"
    if c["modo"] == "adaptativo":
        return (f"fuentes de referencia con criterio **relajado** (área a mitad de máximo ≤ {c['area_max']} px, mediana {c['area_mediana']:.0f} px, "
                f"frente a ≤ 12 px del criterio estricto, que sólo cumplían {c['n_estricto']}): son más extensas que un punto")
    return ("máximos brillantes **sin filtro de forma** (el criterio estricto no encontró fuentes suficientes): pueden no ser puntuales ni estar aisladas")


def fwhm_fuentes(R):
    """FWHM equivalente (px) de las fuentes de referencia, a partir del área mediana a mitad de máximo (círculo)."""
    c = R.get("ej2_criterio_fuentes")
    return 2 * (c["area_mediana"] / 3.141592653589793) ** 0.5 if c else None


def vs_fft(q, tol=0.05):
    """Frase sobre un MSE comparado con el de la FFT; q = MSE_método / MSE_FFT."""
    if q > 1 + tol:
        return (f"{q:.0f}" if q >= 10 else f"{q:.2f}") + " veces el MSE de la FFT (**peor**)"
    if q < 1 - tol:
        return f"{1 / q:.2f} veces menos MSE que la FFT (**mejor**)"
    return f"prácticamente igual a la FFT (cociente de MSE {q:.3f})"
