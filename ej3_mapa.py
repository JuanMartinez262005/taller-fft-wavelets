"""Mapa mental sistémico (Ejercicio 3): jerarquía raíz -> 5 nodos -> 3 subconceptos + enlaces cruzados etiquetados."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.path import Path as MplPath

NODOS = [
    dict(id=1, color="#2563eb", fondo="#eff6ff", titulo="① Reducción de dimensionalidad",
         hijos=["Lineales: PCA · SVD · MDS clásico\nproyección ortogonal $\\hat{x}=U_kU_k^{T}x$",
                "No lineales: Kernel-PCA · Isomap\nt-SNE · UMAP · autoencoders",
                "Eckart–Young: $X_k$ es la mejor\naproximación de rango $k$"]),
    dict(id=2, color="#16a34a", fondo="#f0fdf4", titulo="② Bases adaptativas vs. rígidas",
         hijos=["Adaptativas (dependen de los datos):\nSVD · PCA · Karhunen–Loève",
                "Rígidas (fijas, universales):\nFourier · DCT · Wavelets",
                "Compromiso: óptima en MSE vs.\ncoste, entrenamiento, universalidad"]),
    dict(id=3, color="#d97706", fondo="#fffbeb", titulo="③ Localización espacio/tiempo–frecuencia",
         hijos=["Fourier: frecuencia pura, soporte global\n$\\Delta\\omega\\to 0$, $\\Delta t\\to\\infty$",
                "STFT / Gabor: ventana fija,\nresolución t–f uniforme",
                "Wavelets: multirresolución, soporte compacto\n$\\Delta t\\,\\Delta\\omega\\geq 1/2$ (Heisenberg–Gabor)"]),
    dict(id=4, color="#db2777", fondo="#fdf2f8", titulo="④ Truncamiento y sparsity",
         hijos=["Varianza explicada:\n$\\sum_{i\\leq k}\\sigma_i^2\\,/\\,\\sum_i\\sigma_i^2\\geq 95\\%$",
                "Umbralizado de energía: duro (top-$M$)\no suave $S_\\lambda(w)$",
                "Norma $L_2$ vs. $L_1$: $L_2$ ↔ Parseval, MSE;\n$L_1$ promueve sparsity"]),
    dict(id=5, color="#7c3aed", fondo="#f5f3ff", titulo="⑤ Aplicaciones en ML y Física",
         hijos=["Denoising y compresión:\numbralizado wavelet, JPEG, JPEG2000",
                "Extracción de características:\neigenfaces, rasgos wavelet, espectrogramas",
                "Modos propios: POD (SVD de snapshots)\ny DMD (frecuencia y tasa de crecimiento)"]),
]

# (origen, destino, altura del arco, etiqueta)
ENLACES = [
    (1, 2, 0.65, "la proyección elige un subespacio\n= elegir una base"),
    (2, 3, 0.65, "Fourier: soporte global\nWavelets: soporte compacto"),
    (3, 4, 0.65, "localización ⇒ pocos coeficientes\ngrandes (sparsity)"),
    (4, 5, 0.65, "el umbral $\\lambda$ controla\ndenoising y compresión"),
    (2, 4, 1.45, "base adaptada ⇒ $\\sigma_i$ y $|c_k|$\ndecaen rápido"),
    (3, 5, 1.45, "señales no estacionarias\n(ECG, ondas gravitacionales)"),
    (4, 1, 2.25, "$k$ = nº de modos con\nvarianza explicada ≥ 95 %"),
    (1, 5, 3.05, "SVD / PCA ⇒ eigenfaces, POD"),
]


def dibujar_mapa(ruta_base=None, mostrar=False):
    W, Hh = 22.2, 11.7
    fig, ax = plt.subplots(figsize=(W, Hh))
    ax.set_xlim(0, W); ax.set_ylim(0, Hh); ax.set_aspect("equal"); ax.axis("off")

    ancho_main, alto_main = 3.9, 0.85
    paso = 4.2
    cx = {n["id"]: 2.4 + paso*(n["id"]-1) for n in NODOS}
    y_main = 5.55                      # centro de la fila de nodos principales
    y_top = y_main + alto_main/2       # borde superior
    y_bus = y_top + 3.75               # barra horizontal desde la raíz
    y_root = y_bus + 0.95

    def caja(x, y, w, h, fc, ec, lw=1.6, z=3):
        p = FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle="round,pad=0.02,rounding_size=0.13",
                           fc=fc, ec=ec, lw=lw, zorder=z)
        ax.add_patch(p)

    # raíz
    caja(W/2, y_root, 9.6, 1.05, "#1f2937", "#111827", lw=2)
    ax.text(W/2, y_root+0.17, "Representación compacta de datos y señales", color="white",
            ha="center", va="center", fontsize=15, fontweight="bold", zorder=4)
    ax.text(W/2, y_root-0.2, "$x\\approx\\sum_k c_k\\,\\varphi_k$ : elegir la base $\\varphi$ y truncar a $M$ términos",
            color="white", ha="center", va="center", fontsize=12.5, zorder=4)

    # barra raíz -> nodos principales
    ax.plot([W/2, W/2], [y_root-0.525, y_bus], color="#9ca3af", lw=1.3, zorder=1)
    ax.plot([cx[1], cx[5]], [y_bus, y_bus], color="#9ca3af", lw=1.3, zorder=1)
    for n in NODOS:
        ax.add_patch(FancyArrowPatch((cx[n["id"]], y_bus), (cx[n["id"]], y_top), arrowstyle="-|>",
                                     mutation_scale=14, color="#9ca3af", lw=1.3, zorder=1, shrinkA=0, shrinkB=0))

    # nodos principales + hijos
    alto_hijo, sep_hijo = 1.05, 0.3
    ancho_hijo = 3.4
    for n in NODOS:
        x = cx[n["id"]]
        caja(x, y_main, ancho_main, alto_main, n["fondo"], n["color"], lw=2.6)
        ax.text(x, y_main, n["titulo"], ha="center", va="center", fontsize=12.2, fontweight="bold",
                color=n["color"], zorder=4)
        x_tronco = x - ancho_main/2 + 0.2
        y_prev = y_main - alto_main/2
        y_fin = None
        for k, h in enumerate(n["hijos"]):
            yc = y_main - alto_main/2 - 0.55 - alto_hijo/2 - k*(alto_hijo+sep_hijo)
            xl = x + ancho_main/2 - ancho_hijo/2
            caja(xl, yc, ancho_hijo, alto_hijo, n["fondo"], n["color"], lw=1.4)
            ax.text(xl, yc, h, ha="center", va="center", fontsize=10.6, color="#111827", zorder=4, linespacing=1.35)
            ax.add_patch(FancyArrowPatch((x_tronco, yc), (xl-ancho_hijo/2, yc), arrowstyle="-|>",
                                         mutation_scale=11, color=n["color"], lw=1.5, zorder=2, shrinkA=0, shrinkB=0))
            y_fin = yc
        ax.plot([x_tronco, x_tronco], [y_main-alto_main/2, y_fin], color=n["color"], lw=1.5, zorder=2)

    # enlaces cruzados: puntos de anclaje en el borde superior con anidamiento sin cruces
    color = {n["id"]: n["color"] for n in NODOS}
    anclas = {}
    for i in cx:
        izq, der = [], []
        for e, (a, b, h, _) in enumerate(ENLACES):
            otro = b if a == i else a if b == i else None
            if otro is None: continue
            (izq if otro < i else der).append((abs(otro-i), e))
        izq.sort()                      # más cercano primero  -> más a la izquierda
        der.sort(reverse=True)          # más lejano primero   -> más cerca del centro
        slots_izq = [-1.45, -0.95, -0.5][: len(izq)] if len(izq) > 1 else [-0.75]
        slots_der = [0.5, 0.95, 1.45][: len(der)] if len(der) > 1 else [0.75]
        # orden izquierda->derecha: cercanos izq (más a la izq) ... lejanos izq (cerca del centro)
        izq_pos = sorted(slots_izq)
        der_pos = sorted(slots_der)
        for (_, e), s in zip(izq, izq_pos): anclas[(i, e)] = cx[i] + s
        for (_, e), s in zip(der, der_pos): anclas[(i, e)] = cx[i] + s

    for e, (a, b, h, txt) in enumerate(ENLACES):
        x0, x1 = anclas[(a, e)], anclas[(b, e)]
        y0 = y_top + 0.04
        verts = [(x0, y0), ((x0+x1)/2, y0+2*h), (x1, y0)]
        path = MplPath(verts, [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3])
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=15, lw=2.0,
                                     linestyle=(0, (5, 3)), color=color[a], zorder=5, shrinkA=0, shrinkB=0))
        xm, ym = (x0+x1)/2, y0+h
        ax.text(xm, ym, txt, ha="center", va="center", fontsize=9.8, color=color[a], zorder=6, linespacing=1.3,
                bbox=dict(boxstyle="round,pad=0.28", fc="white", ec=color[a], lw=1.1))

    ax.text(0.3, 0.3, "Flechas continuas: jerarquía (contiene a)   ·   Flechas discontinuas: relaciones entre nodos (color = nodo de origen)",
            fontsize=10.5, color="#4b5563", ha="left", va="center")
    fig.tight_layout(pad=0.3)
    if ruta_base:
        for ext in ("png", "svg", "pdf"):
            fig.savefig(f"{ruta_base}.{ext}", dpi=150 if ext == "png" else None, facecolor="white")
    if mostrar:
        plt.show()
    return fig


if __name__ == "__main__":
    dibujar_mapa("figuras/mapa_mental")
    print("ok")
