# -*- coding: utf-8 -*-
"""Celdas del Ejercicio 2 (aislamiento de fuentes puntuales 2D: DWT vs FFT)."""

EJ2 = []
md = lambda s: EJ2.append(("md", s))
code = lambda s: EJ2.append(("code", s))

md(r'''---
# Ejercicio 2 · Aislamiento de fuentes puntuales en 2D (campo profundo)

**Objetivo.** Aislar las fuentes compactas (estrellas y galaxias distantes) de una imagen de campo profundo con dos herramientas: **2D-DWT** (umbralizado suave de los detalles, aproximación anulada) y **2D-FFT** (filtro pasa-altas gaussiano), y compararlas.

> **Sobre los datos.** La actividad pide una imagen del **JWST** (≥ 512 × 512 px, FITS o monocromática). La celda siguiente acepta un FITS de JWST (archivos `*_i2d.fits`, que guardan la imagen en la extensión `SCI`) o una imagen PNG/TIF/JPG, asignando su ruta a `RUTA_IMAGEN`; todos los umbrales se recalculan a partir de los datos. Si `RUTA_IMAGEN = None` se usa como **sustituto** el *Hubble eXtreme Deep Field* de *scikit-image* (872 × 1000 px, luminancia de 8 bits ya estirada, **no es JWST**); sirve para comprobar que el código funciona, pero no cumple lo que pide el enunciado.''')

code(r'''RUTA_IMAGEN = None        # ← p. ej. "data/jwst_nircam_f200w_i2d.fits" o "data/jwst_deep_field.png"; None = Hubble XDF (sólo de respaldo)
MAX_LADO = 2048           # si la imagen es mayor, se usa un recorte de MAX_LADO × MAX_LADO (memoria y tiempo)
ORIGEN_RECORTE = None     # (fila, columna) de la esquina superior izquierda del recorte (ya sobre la imagen reducida); None = recorte centrado
REDUCIR = 1               # factor entero: promedia la imagen en bloques REDUCIR × REDUCIR antes de recortar (1 = no hacer nada). Útil si las estrellas salen muy anchas

def _gris_reducido(d, f=1):
    """Luminancia 2D (float64) de una imagen en gris, RGB/RGBA o con canales primero; con f > 1 la promedia en bloques f × f.
    Se procesa por BANDAS de filas: una imagen de ~1 Gpíxel no cabe entera en float64 (≈ 9 GB sólo en gris y ≈ 26 GB si se convierte el RGB completo)."""
    import warnings
    while d.ndim > 3: d = d[0]                                           # pila de páginas/canales extra
    if d.ndim == 3:
        if d.shape[-1] in (3, 4):   pass                                 # RGB / RGBA
        elif d.shape[0] in (3, 4):  d = np.moveaxis(d[:3], 0, -1)        # canales primero
        else:                       d = d[..., 0] if d.shape[-1] <= 2 else d[0]   # varias páginas: la primera
    f = max(1, int(f))
    h, w = (d.shape[0] // f) * f, (d.shape[1] // f) * f
    if (h // f) * (w // f) > 3e8:
        raise MemoryError(f"La imagen reducida mediría {h // f}×{w // f} px ({h * w / f**2 / 1e9:.2f} Gpíxeles): no cabe en memoria como decimales de 64 bits. "
                          f"Use un factor mayor: python construir.py --imagen ... --reducir {max(2, int(np.ceil((h * w / 1.5e8) ** 0.5)))} (o más).")
    sal = np.empty((h // f, w // f))
    paso = f * max(1, 4_000_000 // (w * f))                              # ≈ 4 Mpíxeles por banda
    for i in range(0, h, paso):
        b = d[i:min(i + paso, h), :w]
        g = rgb2gray(b[..., :3]) if b.ndim == 3 else np.asarray(b, dtype=float)
        if f > 1:
            g = g.reshape(g.shape[0] // f, f, w // f, f)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")                          # bloques enteros de NaN → NaN (se rellenan después)
                g = np.nanmean(g, axis=(1, 3)) if np.isnan(g).any() else g.mean(axis=(1, 3))
        sal[i // f:i // f + g.shape[0]] = g
    return sal

def _leer_raster(ruta):
    """PNG/TIF/JPG → array. Un TIFF comprimido (LZW, etc.) necesita 'imagecodecs'; si falta, se prueba con Pillow."""
    from skimage import io
    try:
        return io.imread(ruta)
    except (ValueError, ImportError, RuntimeError) as e1:
        from PIL import Image
        try:
            with Image.open(ruta) as im:
                if im.mode in ("P", "PA"): im = im.convert("RGB")
                return np.asarray(im)
        except Exception as e2:
            raise RuntimeError(f"No se pudo leer {Path(ruta).name}.\n  scikit-image: {e1}\n  Pillow: {e2}\n"
                               "Si es un TIFF comprimido, instale el decodificador:  pip install imagecodecs") from e2

def leer_imagen(ruta, f=1):
    """Lee un FITS (JWST: extensión 'SCI') o una imagen PNG/TIF/JPG y devuelve un array 2D float (promediado en bloques f × f si f > 1)."""
    ruta = Path(ruta)
    if ruta.suffix.lower() in (".fits", ".fit", ".fz"):
        from astropy.io import fits
        with fits.open(ruta) as hdul:
            cand = [h for h in hdul if h.data is not None and np.ndim(h.data) >= 2]     # la HDU primaria de JWST está vacía
            if not cand:
                raise ValueError(f"{ruta.name}: ninguna extensión contiene una imagen 2D")
            h = next((h for h in cand if h.name.upper() == "SCI"), cand[0])
            d = h.data; nombre = f"{ruta.name} [extensión {h.name or 'primaria'}]"
            while d.ndim > 2: d = d[0]
            d = _gris_reducido(d, f)
    else:
        d = _gris_reducido(_leer_raster(ruta), f)
        nombre = f"{ruta.name} (imagen, luminancia)"
    return d, nombre

def cargar_imagen(ruta=None):
    lo = hi = None
    if ruta:
        d, nombre = leer_imagen(ruta, REDUCIR)
        if REDUCIR > 1: nombre += f" · promediada en bloques de {REDUCIR}×{REDUCIR}"
        malos = ~np.isfinite(d)
        s = int(np.ceil(np.sqrt(d.size / 1e7))) if d.size > 3e7 else 1    # imágenes enormes: las estadísticas globales se estiman con una submuestra
        if malos.any():                                                   # JWST: NaN fuera del campo cubierto
            d[malos] = np.nanmedian(d[::s, ::s])
        lo, hi = np.percentile(d[::s, ::s], [0.5, 99.95]); hi = max(hi, lo + 1e-12)
        fuente = f"{nombre} · normalizada con percentiles 0.5–99.95 → [0, 1] · {100*malos.mean():.1f} % de píxeles no finitos rellenados"
        del malos
    else:
        ruta_h = DATA / "hubble_deep_field_rgb.npy"
        if ruta_h.exists():
            rgb = np.load(ruta_h)
        else:                                                            # misma imagen, descargada por scikit-image (necesita red y pooch)
            from skimage import data as skdata
            rgb = skdata.hubble_deep_field()
        d = rgb2gray(rgb)
        fuente = "SUSTITUTO: Hubble eXtreme Deep Field (scikit-image), luminancia en [0, 1] — NO es JWST"
        if REDUCIR > 1:
            d = _gris_reducido(d, REDUCIR); fuente += f" · promediada en bloques de {REDUCIR}×{REDUCIR}"
    if min(d.shape) > MAX_LADO or max(d.shape) > MAX_LADO:                # recorte central si es muy grande
        y0, x0 = ORIGEN_RECORTE if ORIGEN_RECORTE else (max(0, (d.shape[0] - MAX_LADO) // 2), max(0, (d.shape[1] - MAX_LADO) // 2))
        d = d[y0:y0 + MAX_LADO, x0:x0 + MAX_LADO]; fuente += f" · recorte {d.shape[0]}×{d.shape[1]} desde (fila {y0}, col {x0})"
    if lo is not None:                                                    # se normaliza sólo el recorte (mismos valores que normalizar todo y recortar)
        d = np.clip((d - lo) / (hi - lo), 0, 1)
    H, W = (d.shape[0] // 4) * 4, (d.shape[1] // 4) * 4                  # múltiplos de 4 → L = 2 exacto
    if min(H, W) < 512: print(f"AVISO: la imagen mide {H}×{W}; el enunciado pide al menos 512×512")
    negro = np.mean(d[:H, :W] <= 0)
    if ruta and negro > 0.03:                                             # bordes vacíos (mosaicos): falsean el fondo y la normalización
        print(f"AVISO: el {100*negro:.1f} % de los píxeles es negro puro (¿bordes vacíos del mosaico?). Use --origen FILA COL "
              "para recortar una zona con señal; si no, el fondo (mediana) y los umbrales quedan falseados.")
    return d[:H, :W], fuente

img, fuente = cargar_imagen(RUTA_IMAGEN)
Hd, Wd = img.shape
MED = np.median(img)                                              # nivel del fondo
SIG = np.median(np.abs(img - MED)) / 0.6745                       # σ del fondo (MAD)
print(f"{fuente}\nTamaño: {Hd} × {Wd} px | rango [{img.min():.3f}, {img.max():.3f}] | fondo = {MED:.4f} ± {SIG:.4f} (MAD)")
RES.update(ej2_shape=[Hd, Wd], ej2_fondo=float(MED), ej2_sig_fondo=float(SIG), ej2_fuente=fuente, ej2_reducir=int(REDUCIR))''')

md(r'''### Parte 1 · Adquisición e inspección (escala logarítmica, `magma`)''')

code(r'''fig, ax = plt.subplots(1, 2, figsize=(14, 6.3))
ax[0].imshow(img, cmap="gray", vmin=0, vmax=1); ax[0].set_title("Lineal")
im = ax[1].imshow(np.log10(img + 1e-3), cmap="magma"); ax[1].set_title(r"Escala logarítmica: $\log_{10}(I + 10^{-3})$")
plt.colorbar(im, ax=ax[1], fraction=.046, pad=.02, label=r"$\log_{10}$ intensidad")
for a in ax: a.axis("off")
plt.tight_layout(); plt.show()''')

md(r'''### Parte 2 · 2D-DWT a nivel $L=2$ (`db2`, ortogonal y de soporte muy corto)
Se usa `mode="periodization"` (transformada ortogonal exacta). Subbandas: aproximación $cA_2$ y detalles $(cH_i, cV_i, cD_i)$, $i=1,2$ — horizontal, vertical y diagonal.''')

code(r'''WAV2, L2 = "db2", 2
coefs2 = pywt.wavedec2(img, WAV2, level=L2, mode="periodization")
cA2, (cH2, cV2, cD2), (cH1, cV1, cD1) = coefs2
print("cA2:", cA2.shape, "| nivel 2:", cH2.shape, "| nivel 1:", cH1.shape)

def mostrar_detalle(ax, w, titulo):
    a = np.percentile(np.abs(w), 99.5)
    ax.imshow(w, cmap="RdBu_r", vmin=-a, vmax=a, interpolation="nearest"); ax.set_title(titulo, fontsize=10); ax.axis("off")

fig = plt.figure(figsize=(15, 8.2)); gs = gridspec.GridSpec(2, 4, wspace=.06, hspace=.14)
a = fig.add_subplot(gs[0, 0]); a.imshow(np.log10(np.clip(cA2 / 2**L2, 0, None) + 1e-3), cmap="magma"); a.set_title(r"$cA_2$ (aproximación, log)", fontsize=10); a.axis("off")
for k, (w, n) in enumerate([(cH2, "cH2"), (cV2, "cV2"), (cD2, "cD2")], start=1): mostrar_detalle(fig.add_subplot(gs[0, k]), w, f"{n}  (detalle nivel 2)")
for k, (w, n) in enumerate([(cH1, "cH1"), (cV1, "cV1"), (cD1, "cD1")]): mostrar_detalle(fig.add_subplot(gs[1, k]), w, f"{n}  (detalle nivel 1)")
a = fig.add_subplot(gs[1, 3]); a.axis("off")
a.text(0, .5, "Convención:\ncH: bordes horizontales\ncV: bordes verticales\ncD: diagonales\n\nDetalles en escala simétrica\n(percentil 99.5 de |w|)", va="center", fontsize=10)
plt.show()''')

md(r'''### Parte 3 · 2D-FFT centrada: espectro de magnitud logarítmica y espectro de fase''')

code(r'''F = fftshift(fft2(img))                                          # espectro centrado
fxc, fyc = fftshift(fftfreq(Wd)), fftshift(fftfreq(Hd))          # ciclos/píxel
ext = [fxc[0], fxc[-1], fyc[-1], fyc[0]]
fig, ax = plt.subplots(1, 2, figsize=(14, 6))
i0 = ax[0].imshow(np.log1p(np.abs(F)), cmap="magma", extent=ext); ax[0].set_title(r"$\log(1+|F(u,v)|)$")
plt.colorbar(i0, ax=ax[0], fraction=.046, pad=.02)
i1 = ax[1].imshow(np.angle(F), cmap="twilight", extent=ext, vmin=-np.pi, vmax=np.pi); ax[1].set_title(r"$\arg F(u,v)$ [rad]")
plt.colorbar(i1, ax=ax[1], fraction=.046, pad=.02)
for a in ax: a.set(xlabel="u [ciclos/px]", ylabel="v [ciclos/px]")
plt.tight_layout(); plt.show()''')

md(r'''La energía se concentra en torno a $(u,v)=(0,0)$ (estructura de gran escala y fondo) y decae radialmente; la "cruz" sobre los ejes procede de la discontinuidad entre bordes opuestos de la imagen (la DFT supone periodicidad). **La fase**, aparentemente ruidosa, es la que codifica *dónde* está cada fuente: un objeto puntual en $(x_0,y_0)$ aporta un factor $e^{-i2\pi(ux_0+vy_0)}$ a todas las frecuencias.''')

md(r'''### Parte 4 · Aislamiento con wavelets: $cA_2=0$ + umbralizado suave $S_\lambda(w)=\operatorname{sign}(w)\max(|w|-\lambda,0)$

**Elección de $\lambda$ (justificada, no arbitraria).** El ruido del fondo se estima con el estimador robusto MAD sobre el detalle diagonal más fino, $\hat\sigma=\operatorname{mediana}(|cD_1|)/0.6745$, y se aplica el umbral universal de Donoho–Johnstone, $\lambda=\hat\sigma\sqrt{2\ln n}$ ($n$ = número de píxeles). Anular $cA_2$ elimina las componentes de baja frecuencia ($< 1/8$ ciclos/px aprox.): el fondo suave y el interior de objetos extensos.''')

code(r'''def umbral_suave(w, lam): return np.sign(w) * np.maximum(np.abs(w) - lam, 0.0)
def umbral_duro(w, lam):  return w * (np.abs(w) > lam)

def aislar_wavelet(imagen, lam, wavelet=WAV2, nivel=L2, tipo="suave", coefs=None):
    """cA = 0 + umbralizado de los detalles + reconstrucción 2D-IDWT."""
    if coefs is None:
        coefs = pywt.wavedec2(imagen, wavelet, level=nivel, mode="periodization")
    u = umbral_suave if tipo == "suave" else umbral_duro
    nuevo = [np.zeros_like(coefs[0])] + [tuple(u(d, lam) for d in det) for det in coefs[1:]]
    return pywt.waverec2(nuevo, wavelet, mode="periodization")

SIGMA_W = np.median(np.abs(cD1)) / 0.6745                         # ruido del fondo en la escala de los detalles
LAM = SIGMA_W * np.sqrt(2 * np.log(img.size))                     # umbral universal (VisuShrink)
f_wav = aislar_wavelet(img, LAM, coefs=coefs2)
det_todos = np.concatenate([d.ravel() for det in coefs2[1:] for d in det])
sobrev = np.mean(np.abs(det_todos) > LAM)
print(f"σ̂ (MAD de cD1) = {SIGMA_W:.4f}   |   λ = σ̂·√(2 ln n) = {LAM:.4f}   ({LAM/SIGMA_W:.2f} σ̂)")
print(f"Coeficientes de detalle que sobreviven al umbral: {100*sobrev:.2f} %  ({int(sobrev*det_todos.size)} de {det_todos.size})")
RES.update(ej2_sigma_w=float(SIGMA_W), ej2_lambda=float(LAM), ej2_sobrev_pct=float(100*sobrev))''')

md(r'''### Parte 5 · Aislamiento con 2D-FFT: filtro pasa-altas gaussiano
$H(u,v)=1-\exp\!\left(-\dfrac{D^2(u,v)}{2D_0^2}\right)$, con $D$ la distancia al centro del espectro desplazado (en ciclos/píxel).

**Elección de $D_0$ (comparación en igualdad de condiciones).** Anular $cA_2$ en la DWT equivale a cortar por debajo de $f_c=1/8$ ciclos/px (borde de la banda de aproximación del nivel 2). Se fija $D_0$ para que el filtro gaussiano valga **0.5 justo en $f_c$**: $H(f_c)=0.5\Rightarrow D_0=f_c/\sqrt{2\ln 2}\approx 0.106$ ciclos/px.''')

code(r'''FC = 1/8                                                    # ciclos/px: borde de la aproximación del nivel 2
D0 = FC / np.sqrt(2 * np.log(2))                            # H(FC) = 0.5

def filtro_fft(imagen, tipo="gauss", D0=D0):
    """2D-FFT centrada → filtro pasa-altas H(u,v) → transformada inversa (parte real)."""
    Fc = fftshift(fft2(imagen))
    u, v = fftshift(fftfreq(imagen.shape[1])), fftshift(fftfreq(imagen.shape[0]))
    Dd = np.hypot(*np.meshgrid(u, v))                       # D(u,v) [ciclos/px]
    Hh = 1 - np.exp(-Dd**2 / (2 * D0**2)) if tipo == "gauss" else (Dd > D0).astype(float)   # "ideal": corte abrupto
    return np.real(ifft2(ifftshift(Fc * Hh)))

f_fft   = filtro_fft(img, "gauss")                          # el filtro pedido en el enunciado
f_ideal = filtro_fft(img, "ideal")                          # variante de CORTE ABRUPTO (se usa en la pregunta B)
print(f"D0 = {D0:.4f} ciclos/px  →  H(f_c = 1/8) = {1-np.exp(-FC**2/(2*D0**2)):.2f}")
RES["ej2_D0"] = float(D0)''')

md(r'''### Parte 6 · Evaluación comparativa — panel 2 × 2
Original (escala log) · $f_{wavelet}$ · $f_{fft}$ · $|f_{wavelet}-f_{fft}|$. Las dos imágenes filtradas se muestran con **la misma escala lineal** (los valores negativos —halos— se recortan a negro en la visualización; se cuantifican en la tabla de abajo).''')

code(r'''vmax = np.percentile(np.maximum(np.maximum(f_wav, f_fft), 0), 99.5)
fig, ax = plt.subplots(2, 2, figsize=(14, 12.2))
ax[0, 0].imshow(np.log10(img + 1e-3), cmap="magma"); ax[0, 0].set_title("Original (escala logarítmica)")
ax[0, 1].imshow(f_wav, cmap="magma", vmin=0, vmax=vmax); ax[0, 1].set_title(r"$f_{wavelet}$: DWT db2 L=2, $cA_2=0$, umbral suave")
ax[1, 0].imshow(f_fft, cmap="magma", vmin=0, vmax=vmax); ax[1, 0].set_title(r"$f_{fft}$: pasa-altas gaussiano ($D_0$ = %.3f ciclos/px)" % D0)
ax[1, 1].imshow(np.abs(f_wav - f_fft), cmap="inferno", vmin=0, vmax=vmax); ax[1, 1].set_title(r"$|f_{wavelet}-f_{fft}|$")
for a in ax.ravel(): a.axis("off")
plt.tight_layout(); plt.savefig(FIG / "ej2_panel_2x2.png", bbox_inches="tight", dpi=110); plt.show()''')

md(r'''#### Evidencia cuantitativa para la discusión
Se seleccionan **fuentes de referencia** (máximos locales que destacan del fondo) y se mide, para cada método: (i) **posición** (¿el máximo local queda a ≤ 1 px del original?), (ii) **retención del pico**, (iii) **correlación de forma** (recorte 11×11 vs. original), (iv) **RMS en cielo vacío** (≥ 14 px de cualquier objeto) y (v) **halo**: media y mínimo del anillo 4–10 px alrededor de la fuente, normalizados por su amplitud (negativo ⇒ sub-oscilación).

**Criterio de selección.** *Estricto* (el habitual): fuentes **compactas** (área a mitad de máximo ≤ 12 px) y **aisladas** (en el borde de una ventana de 15 × 15 px la señal cae por debajo del 35 % de la amplitud); el umbral de detección baja 15σ → 10σ → 6σ → 4σ hasta reunir ≥ 30. **Si la imagen no tiene al menos 10 fuentes así** (típico de imágenes muy sobremuestreadas, con nebulosidad intensa o muy pobladas) el cuaderno cambia a un criterio *adaptativo*: detecta sobre la imagen sin el fondo suave y **relaja poco a poco** el área (hasta 80 px) y el aislamiento (hasta el 65 %). Imprime cuál acabó usando, lo guarda en `resultados.json` y la discusión lo menciona: con el criterio relajado las fuentes ya no son «puntuales» en sentido estricto. **Comprueben siempre la figura de selección** (siguiente celda): las fuentes elegidas deben ser estrellas o galaxias compactas, no filamentos de nebulosa.''')

code(r'''def _medir(imagen, y, x0):
    """Amplitud sobre el fondo local, área a mitad de máximo y máximo del borde de la ventana 15×15 centrada en (y, x0)."""
    p = imagen[y-7:y+8, x0-7:x0+8]; fondo = np.median(p)
    amp = imagen[y-1:y+2, x0-1:x0+2].max() - fondo
    lab, _ = ndi.label(p > fondo + 0.5 * amp)
    area = np.sum(lab == lab[7, 7]) if lab[7, 7] > 0 else 999
    borde = np.r_[p[0], p[-1], p[:, 0], p[:, -1]]
    return amp, area, borde.max(), fondo

def residuo(imagen):
    """Imagen sin la nebulosidad ni los gradientes suaves (σ = 20 px) y su σ robusto: permite detectar en campos con fondo estructurado."""
    bgs = ndi.gaussian_filter(imagen, 1); res = bgs - ndi.gaussian_filter(bgs, 20)
    return res, np.median(np.abs(res - np.median(res))) / 0.6745

def fuentes_puntuales(imagen, n_max=150, area_max=12, minimas=30, suficientes=10):
    """Fuentes de referencia → (posiciones, amplitudes, criterio aplicado).
    1) ESTRICTO: compactas (área a mitad de máximo ≤ area_max) y aisladas (borde de la ventana < 35 % de la amplitud); el umbral
       de detección baja 15σ → 10σ → 6σ → 4σ hasta reunir 'minimas'. Es lo que se usa si salen al menos 'suficientes'.
    2) ADAPTATIVO (imágenes sobremuestreadas, con nebulosidad o muy pobladas): detección sobre la imagen residual y relajación
       gradual del área y del aislamiento.
    3) SIN FILTRO DE FORMA (último recurso, si lo anterior da < 5): los máximos más brillantes."""
    # --- 1) criterio estricto
    bgs = ndi.gaussian_filter(imagen, 1)
    for k in (15, 10, 6, 4):
        cand = peak_local_max(bgs, min_distance=8, threshold_abs=MED + k*SIG, exclude_border=14)
        info = []
        for y, x0 in cand:
            amp, area, bmax, fondo = _medir(imagen, y, x0)
            if area <= area_max and bmax < fondo + 0.35 * amp:
                info.append((y, x0, amp, area))
        if len(info) >= minimas: break
    n_estricto = len(info)
    crit = dict(modo="estricto", area_max=int(area_max), aislamiento=0.35, k_sigma=int(k))
    if n_estricto < suficientes:
        # --- 2) criterio adaptativo
        res, sig_res = residuo(imagen)
        cand = peak_local_max(res, min_distance=8, threshold_abs=4*sig_res, exclude_border=14, num_peaks=4000)
        med = [(y, x0, *_medir(imagen, y, x0), res[y, x0] / sig_res) for y, x0 in cand]     # y, x0, amp, área, borde, fondo, nivel (σ)
        def elegir(am, ais, kk):
            return [(y, x0, amp, area) for y, x0, amp, area, bmax, fondo, nivel in med
                    if amp > 0 and area <= am and bmax < fondo + ais * amp and nivel >= kk]
        elegido = None
        for am, ais in ((area_max, .35), (25, .45), (50, .55), (80, .65)):               # primero se relaja la forma ...
            for kk in (15, 10, 6, 4):                                                    # ... y, dentro de cada forma, el umbral
                sel = elegir(am, ais, kk)
                if len(sel) >= minimas: elegido = (sel, am, ais, kk); break
            if elegido: break
        if elegido is None: elegido = (elegir(80, .65, 4), 80, .65, 4)                   # ninguno llega a 'minimas': el más laxo
        n_adapt = len(elegido[0])
        if max(n_estricto, n_adapt) >= 5:
            if n_adapt >= n_estricto:
                info = elegido[0]
                crit = dict(modo="adaptativo", area_max=int(elegido[1]), aislamiento=float(elegido[2]), k_sigma=int(elegido[3]))
        else:                                                                            # --- 3) sin filtro de forma
            info = [(y, x0, amp, area) for y, x0, amp, area, bmax, fondo, nivel in med if amp > 0 and nivel >= 4]
            crit = dict(modo="sin filtro", area_max=None, aislamiento=None, k_sigma=4)
        crit["n_candidatas"] = int(len(med))
    if len(info) == 0:
        raise RuntimeError("No se encontró ninguna fuente brillante; revise la imagen (¿está vacía o es constante?).")
    info = np.array(info, dtype=float); info = info[np.argsort(info[:, 2])[::-1]][:n_max]
    crit.update(n=int(len(info)), area_mediana=float(np.median(info[:, 3])), n_estricto=int(n_estricto))
    crit["relajado"] = bool(crit["modo"] == "sin filtro" or (crit["modo"] == "adaptativo" and (crit["area_max"] > area_max or crit["aislamiento"] > 0.35)))   # ¿se aflojó el criterio de FORMA?
    return info[:, :2].astype(int), info[:, 2], crit

fuentes, amps, crit = fuentes_puntuales(img)
if crit["modo"] == "estricto":
    obj4, obj15 = img > MED + 4*SIG, img > MED + 15*SIG                    # caso habitual: el fondo es plano
else:                                                                     # fondo estructurado: los objetos se definen sobre la imagen residual
    res_, sig_res_ = residuo(img); obj4, obj15 = res_ > 4*sig_res_, res_ > 15*sig_res_
cielo = ~ndi.binary_dilation(obj4, iterations=14)                         # cielo vacío
for it in (8, 4, 2):                                                      # campo muy poblado: reducir la distancia de seguridad
    if cielo.mean() >= 0.02: break
    cielo = ~ndi.binary_dilation(obj4, iterations=it)
    print(f"AVISO: poco cielo vacío; distancia de seguridad reducida a {it} px")
if cielo.sum() < 1000:
    cielo = ~obj4; print("AVISO: casi no hay cielo vacío; se usa todo píxel fuera de los objetos detectados")
cerca = ndi.binary_dilation(obj15, iterations=10) & ~cielo                # alrededor de fuentes brillantes
yy_, xx_ = np.ogrid[-10:11, -10:11]; rr_ = np.hypot(yy_, xx_); anillo = (rr_ >= 4) & (rr_ <= 10)
etq = {"estricto": "compactas y aisladas", "sin filtro": "SIN filtro de forma",
       "adaptativo": "con criterio RELAJADO" if crit["relajado"] else "compactas y aisladas (detectadas sobre la imagen residual)"}[crit["modo"]]
print(f"{len(fuentes)} fuentes {etq} | criterio: {crit['modo']}"
      + (f" (área ≤ {crit['area_max']} px, aislamiento < {crit['aislamiento']:.0%}, ≥ {crit['k_sigma']}σ)" if crit["area_max"] else "")
      + f" | área mediana = {crit['area_mediana']:.0f} px | cielo vacío = {100*cielo.mean():.1f} % de los píxeles")
if crit["modo"] == "adaptativo" and not crit["relajado"]:
    print("NOTA: el fondo es estructurado (nebulosa o gradientes), así que las fuentes se detectaron sobre la imagen residual; la forma exigida es la misma (área ≤ 12 px, aislamiento < 35 %).")
if crit["relajado"]:
    print(f"AVISO: con el criterio estricto (área ≤ 12 px, aislamiento < 35 %) sólo había {crit['n_estricto']} fuentes; las de referencia son más "
          "extensas que un punto. Compruebe en la figura siguiente que son fuentes reales y no filamentos de nebulosa.")
fwhm_ = 2 * np.sqrt(crit["area_mediana"] / np.pi)                          # FWHM equivalente de una fuente circular con esa área a mitad de máximo
if fwhm_ > 4.5:
    k_ = int(max(2, round(fwhm_ / 2.2)))
    print(f"AVISO: las fuentes miden ≈ {fwhm_:.1f} px de FWHM (una estrella a la escala de píxel natural mide 1–3 px). El filtro del enunciado (corte en 1/8 ciclos/px, L = 2) "
          f"elimina buena parte de su energía, así que las comparaciones del Ejercicio 2 valdrán sólo para esa escala. Para fuentes más finas: python construir.py --imagen ... --reducir {REDUCIR * k_}")
RES["ej2_criterio_fuentes"] = crit
RES["ej2_n_fuentes"] = int(len(fuentes))''')

code(r'''# Comprobación visual de la selección: posiciones sobre la imagen y 12 recortes (de la más brillante a la más tenue)
sel_g = np.linspace(0, len(fuentes) - 1, 12).astype(int)
fig = plt.figure(figsize=(15, 7)); gs = gridspec.GridSpec(4, 6, figure=fig, wspace=.05, hspace=.25)
a0 = fig.add_subplot(gs[:, :3]); a0.imshow(np.log10(img + 1e-3), cmap="magma")
a0.scatter(fuentes[:, 1], fuentes[:, 0], s=70, facecolors="none", edgecolors="cyan", linewidths=.7)
a0.set_title(f"{len(fuentes)} fuentes de referencia (criterio {crit['modo']})"); a0.axis("off")
for j, i in enumerate(sel_g):
    a = fig.add_subplot(gs[j // 3, 3 + j % 3]); y, x0 = fuentes[i]
    a.imshow(img[y-10:y+11, x0-10:x0+11], cmap="magma", interpolation="nearest"); a.set_title(f"#{i + 1} · amp. {amps[i]:.2f}", fontsize=8); a.axis("off")
plt.savefig(FIG / "ej2_fuentes_seleccionadas.png", bbox_inches="tight", dpi=110); plt.show()''')

code(r'''def evaluar_fuentes(f):
    desp, ret, cor, anil, lob = [], [], [], [], []
    for (y, x0), a in zip(fuentes, amps):
        sub = f[y-2:y+3, x0-2:x0+3]; dy, dx = np.unravel_index(np.argmax(sub), sub.shape)
        desp.append(np.hypot(dy-2, dx-2)); ret.append(sub.max() / a)
        cor.append(np.corrcoef(img[y-5:y+6, x0-5:x0+6].ravel(), f[y-5:y+6, x0-5:x0+6].ravel())[0, 1])
        pr = f[y-10:y+11, x0-10:x0+11]; anil.append(pr[anillo].mean() / a); lob.append(pr[anillo].min() / a)
    desp = np.array(desp)
    return {"posición ≤1 px": np.mean(desp <= 1), "retención del pico": np.mean(ret), "correlación de forma": np.mean(cor),
            "RMS en cielo vacío": np.sqrt(np.mean(f[cielo]**2)), "halo medio / amp.": np.mean(anil), "halo mínimo / amp.": np.mean(lob)}

metodos = {"Wavelet db2 L=2 · umbral suave (enunciado)": f_wav,
           "Wavelet db2 L=2 · umbral duro": aislar_wavelet(img, LAM, coefs=coefs2, tipo="duro"),
           "FFT pasa-altas gaussiano (enunciado)": f_fft,
           "FFT pasa-altas ideal (corte abrupto)": f_ideal}
tabla_A = pd.DataFrame({k: evaluar_fuentes(v) for k, v in metodos.items()}).T
display(tabla_A.style.format("{:.3f}"))
print(f"(Referencia: RMS del ruido de fondo de la imagen original = {SIG:.4f})")
RES["ej2_tabla_A"] = {k: {c: float(v) for c, v in fila.items()} for k, fila in tabla_A.iterrows()}
RES["ej2_n_fuentes"] = int(len(fuentes))''')

md(r'''**¿Depende la conclusión de los parámetros elegidos?** Se repite la evaluación variando el umbral de la wavelet ($\lambda=2,3,4\hat\sigma$ y el universal $5.2\hat\sigma$) y la frecuencia de corte del filtro gaussiano ($f_c=1/16,\,1/8,\,1/4$ ciclos/px, con $D_0=f_c/\sqrt{2\ln 2}$).''')

code(r'''filas_s = {}
for k in (2, 3, 4, LAM / SIGMA_W):
    filas_s[f"Wavelet suave · λ = {k:.1f} σ̂"] = evaluar_fuentes(aislar_wavelet(img, k * SIGMA_W, coefs=coefs2))
for fc in (1/16, 1/8, 1/4):
    filas_s[f"FFT gaussiano · f_c = 1/{round(1/fc)}"] = evaluar_fuentes(filtro_fft(img, "gauss", D0=fc / np.sqrt(2 * np.log(2))))
tabla_S = pd.DataFrame(filas_s).T[["posición ≤1 px", "retención del pico", "correlación de forma", "RMS en cielo vacío", "halo mínimo / amp."]]
display(tabla_S.style.format("{:.3f}"))
RES["ej2_tabla_sens"] = {k: {c: float(v) for c, v in fila.items()} for k, fila in tabla_S.iterrows()}''')

code(r'''# Recortes y perfiles 1D a través de tres fuentes (brillante, intermedia, tenue)
sel = [0, len(fuentes)//4, len(fuentes)//2]
mt = {"original": img, "wavelet (suave)": f_wav, "FFT gaussiano": f_fft, "FFT ideal": f_ideal}
fig, ax = plt.subplots(3, 5, figsize=(17, 9.6), gridspec_kw={"width_ratios": [1, 1, 1, 1, 1.7]})
for r, i in enumerate(sel):
    y, x0 = fuentes[i]; ya, xa = np.clip(y-20, 0, Hd-41), np.clip(x0-20, 0, Wd-41)
    for c, (nom, f) in enumerate(mt.items()):
        ax[r, c].imshow(f[ya:ya+41, xa:xa+41], cmap="magma", vmin=0 if c else 0, vmax=img[y-2:y+3, x0-2:x0+3].max()); ax[r, c].axis("off")
        if r == 0: ax[r, c].set_title(nom, fontsize=10)
    ax[r, 4].axhline(0, color="0.7", lw=.8)
    for nom, f, col in (("original", img, "k"), ("wavelet (suave)", f_wav, "C3"), ("FFT gaussiano", f_fft, "C0"), ("FFT ideal", f_ideal, "C2")):
        ax[r, 4].plot(np.arange(-20, 21), f[y, x0-20:x0+21] if x0-20 >= 0 and x0+21 <= Wd else np.zeros(41), color=col, lw=1.4, label=nom)
    ax[r, 4].set(xlabel="píxeles respecto al pico", ylabel="intensidad"); ax[r, 4].grid(alpha=.25)
    if r == 0: ax[r, 4].legend(fontsize=8); ax[r, 4].set_title("Perfil horizontal por el pico")
plt.tight_layout(); plt.savefig(FIG / "ej2_perfiles_fuentes.png", bbox_inches="tight"); plt.show()''')

md(r'''### Pregunta B · ¿Por qué la FFT produce oscilaciones (Gibbs) y las wavelets no? — experimentos

**Experimento 1 — respuesta al impulso.** Se inyecta una fuente puntual (gaussiana de $\sigma=1$ px, amplitud 0.9) en una imagen vacía y se mide qué hace cada filtro con ella: forma del perfil, número de cambios de signo y **radio que contiene el 99 % de la energía de la respuesta** (extensión espacial del artefacto).''')

code(r'''cy, cx = Hd // 2, Wd // 2
yy, xx = np.mgrid[0:Hd, 0:Wd]
src = 0.9 * np.exp(-((yy - cy)**2 + (xx - cx)**2) / (2 * 1.0**2))            # fuente puntual sintética sobre fondo nulo (los cambios de signo se cuentan donde |f| > 1 % del pico)
resp = {"FFT gaussiano": filtro_fft(src, "gauss"), "FFT ideal (corte abrupto)": filtro_fft(src, "ideal"),
        "Wavelet db2 L=2 (umbral suave)": aislar_wavelet(src, LAM)}
rr = np.hypot(yy - cy, xx - cx)
rint = np.round(rr).astype(int).ravel()                                      # radio entero de cada píxel

filas, curvas = [], {}
for nom, r in resp.items():
    perf = r[cy, cx:cx+60]; pico = perf[0]
    sg = perf[np.abs(perf) > 0.01 * abs(pico)]
    cruces = int(np.sum(np.diff(np.sign(sg)) != 0))
    e = np.cumsum(np.bincount(rint, weights=(r.ravel()**2))); e /= e[-1]; curvas[nom] = e      # energía acumulada vs. radio
    filas.append({"filtro": nom, "pico / entrada": pico / 0.9, "lóbulo negativo mín. / pico": r.min() / pico,
                  "cambios de signo": cruces, "radio 99 % energía [px]": float(np.searchsorted(e, 0.99))})
tabla_B = pd.DataFrame(filas).set_index("filtro"); display(tabla_B.style.format("{:.3f}", subset=["pico / entrada", "lóbulo negativo mín. / pico", "radio 99 % energía [px]"]))
RES["ej2_tabla_B"] = {k: {c: float(v) for c, v in fila.items()} for k, fila in tabla_B.iterrows()}

fig, ax = plt.subplots(1, 2, figsize=(13, 4.2))
for (nom, r), col in zip(resp.items(), ["C0", "C2", "C3"]):
    ax[0].plot(np.arange(-40, 41), r[cy, cx-40:cx+41], color=col, lw=1.6, label=nom)
    ax[1].semilogx(np.arange(len(curvas[nom]))[1:], curvas[nom][1:], color=col, lw=1.8, label=nom)
ax[0].plot(np.arange(-40, 41), src[cy, cx-40:cx+41], "k--", lw=1, label="fuente de entrada"); ax[0].axhline(0, color="0.7", lw=.8)
ax[0].set(xlabel="píxeles", ylabel="respuesta", title="Respuesta al impulso (perfil por el centro)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.25)
ax[1].axhline(.99, color="0.7", ls=":"); ax[1].set(xlabel="radio [px]", ylabel="energía acumulada (fracción)", title="¿Hasta dónde se extiende el artefacto?", xlim=(1, 300), ylim=(0.7, 1.005)); ax[1].legend(fontsize=8, loc="lower right"); ax[1].grid(alpha=.25, which="both")
plt.tight_layout(); plt.savefig(FIG / "ej2_respuesta_impulso.png", bbox_inches="tight"); plt.show()''')

md(r'''**Experimento 1b — la DWT no es invariante a traslaciones.** La rejilla diádica de la DWT con $L=2$ se repite cada 4 px, y el resultado anterior corresponde a una sola posición de la fuente (el centro, múltiplo de 4). Se repite el experimento para las 16 posiciones relativas posibles $(\delta_y,\delta_x)\in\{0,1,2,3\}^2$. La FFT, al ser un filtro lineal invariante, da siempre el mismo resultado.''')

code(r'''def respuesta_desplazada(dy, dx, metodo):
    s_ = 0.9 * np.exp(-((yy - (cy + dy))**2 + (xx - (cx + dx))**2) / (2 * 1.0**2))
    r_ = aislar_wavelet(s_, LAM) if metodo == "wavelet" else filtro_fft(s_, "gauss")
    return r_.max() / 0.9, r_.min() / r_.max()

filas_d = {}
for metodo in ("wavelet", "fft"):
    v = np.array([respuesta_desplazada(dy, dx, metodo) for dy in range(4) for dx in range(4)])
    filas_d["Wavelet db2 L=2" if metodo == "wavelet" else "FFT gaussiano"] = {
        "pico/entrada · mín.": v[:, 0].min(), "pico/entrada · máx.": v[:, 0].max(),
        "lóbulo neg./pico · mejor": v[:, 1].max(), "lóbulo neg./pico · peor": v[:, 1].min()}
tabla_D = pd.DataFrame(filas_d).T; display(tabla_D.style.format("{:.3f}"))
RES["ej2_traslacion"] = {k: {c: float(v) for c, v in fila.items()} for k, fila in tabla_D.iterrows()}''')

md(r'''**Experimento 2 — el fenómeno de Gibbs en una dimensión.** Suma parcial de la serie de Fourier de un escalón unidad (sólo armónicos impares hasta $N$): $S_N(x)=\tfrac12+\tfrac2\pi\sum_{k\ \mathrm{impar}\le N}\tfrac{\sin kx}{k}$. El sobreimpulso **no desaparece** al aumentar $N$.''')

code(r'''from scipy.special import sici

def sobreimpulso_gibbs(Nterm, npts=6001):
    xg = np.linspace(1e-4, 3*np.pi/Nterm, npts)
    k = np.arange(1, Nterm + 1, 2)[:, None]
    S = 0.5 + (2/np.pi) * np.sum(np.sin(k * xg[None, :]) / k, axis=0)
    return S.max() - 1                                       # fracción del salto (=1)

teo = sici(np.pi)[0] / np.pi - 0.5                           # (1/π)·Si(π) − 1/2
tab_g = pd.DataFrame({"N (armónicos)": [11, 31, 101, 301, 1001],
                      "sobreimpulso [% del salto]": [100*sobreimpulso_gibbs(n) for n in (11, 31, 101, 301, 1001)]})
display(tab_g.style.format({"sobreimpulso [% del salto]": "{:.3f}"}).hide(axis="index"))
print(f"Límite teórico (1/π)·Si(π) − 1/2 = {100*teo:.3f} % del salto  →  no depende de N")
RES["ej2_gibbs"] = dict(tabla={int(n): float(v) for n, v in zip(tab_g.iloc[:, 0], tab_g.iloc[:, 1])}, limite_pct=float(100*teo))''')

md(r'''### Pregunta C · Compresión reteniendo sólo el 5 % de los coeficientes
Se conservan los coeficientes de mayor módulo (FFT: conjugados por pares; DWT: umbral global) y se mide el **MSE en el dominio de la imagen**. Para no sacar conclusiones de un único caso se barren la fracción retenida y la **profundidad de la descomposición** ($L=2$ es la del enunciado).''')

code(r'''def rec_fft_top(imagen, frac):
    Fq = fft2(imagen); E = np.abs(Fq)**2; M = int(frac * imagen.size)
    thr = np.partition(E.ravel(), -M)[-M]
    keep = E >= thr; ky, kx = np.nonzero(keep); keep[(-ky) % Hd, (-kx) % Wd] = True       # simetría hermítica
    return np.real(ifft2(np.where(keep, Fq, 0))), int(keep.sum())

def rec_dwt_top(imagen, frac, wavelet, nivel):
    arr, sl = pywt.coeffs_to_array(pywt.wavedec2(imagen, wavelet, level=nivel, mode="periodization"))
    M = int(frac * imagen.size); thr = np.partition(np.abs(arr).ravel(), -M)[-M]
    arr_c = np.where(np.abs(arr) >= thr, arr, 0)
    return pywt.waverec2(pywt.array_to_coeffs(arr_c, sl, output_format="wavedec2"), wavelet, mode="periodization"), int(np.count_nonzero(arr_c))

configs = [("FFT", None, None), ("db2 L=2 (enunciado)", "db2", 2), ("db4 L=2", "db4", 2), ("bior2.2 L=2", "bior2.2", 2),
           ("db2 L=3", "db2", 3), ("db2 L=5", "db2", 5), ("db4 L=5", "db4", 5), ("bior2.2 L=5", "bior2.2", 5)]
fracs2 = [0.005, 0.01, 0.02, 0.05, 0.10, 0.20]
mse_tab = pd.DataFrame(index=[f"{100*f:g} %" for f in fracs2], columns=[c[0] for c in configs], dtype=float)
rec5 = {}
for fr in fracs2:
    for nom, w, L in configs:
        rec, Mq = rec_fft_top(img, fr) if w is None else rec_dwt_top(img, fr, w, L)
        mse_tab.loc[f"{100*fr:g} %", nom] = np.mean((img - rec)**2)
        if abs(fr - 0.05) < 1e-9: rec5[nom] = rec
display(mse_tab.style.format("{:.3e}"))
RES["ej2_mse_tab"] = {idx: {c: float(v) for c, v in fila.items()} for idx, fila in mse_tab.iterrows()}

fig, ax = plt.subplots(figsize=(8.5, 4.6))
for (nom, w, L), col in zip(configs, ["C0", "C3", "C1", "C6", "C2", "C4", "C5", "C8"]):
    ax.semilogy([100*f for f in fracs2], mse_tab[nom], "o-", color=col, label=nom, lw=2.2 if nom.startswith(("FFT", "db2 L=2")) else 1.3)
ax.axvline(5, color="0.6", ls=":"); ax.set(xlabel="coeficientes retenidos [% de N]", ylabel="MSE", title="MSE vs. fracción retenida"); ax.legend(fontsize=8); ax.grid(alpha=.3, which="both")
plt.tight_layout(); plt.savefig(FIG / "ej2_mse_vs_fraccion.png", bbox_inches="tight"); plt.show()
print("Tamaño de la banda de aproximación como fracción de N:", {f"L={L}": f"{100/4**L:.2f} %" for L in (2, 3, 5)})''')

code(r'''# ¿Dónde cae el error al 5 %? (cielo vacío / cerca de fuentes brillantes / total) y píxeles con valor negativo (físicamente imposible)
filas = {}
for nom in ["FFT", "db2 L=2 (enunciado)", "db2 L=5"]:
    e = img - rec5[nom]
    filas[nom] = {"RMS error en cielo": np.sqrt(np.mean(e[cielo]**2)), "RMS error cerca de fuentes": np.sqrt(np.mean(e[cerca]**2)),
                  "RMS total": np.sqrt(np.mean(e**2)), "mín. reconstruido": rec5[nom].min(), "píxeles < 0 [%]": 100*np.mean(rec5[nom] < 0)}
tabla_C = pd.DataFrame(filas).T; display(tabla_C.style.format("{:.4f}"))

def composicion_retenidos(imagen, frac, wavelet, nivel):
    """¿Qué parte del presupuesto de coeficientes se gasta en la banda de aproximación cA_L?"""
    arr, sl = pywt.coeffs_to_array(pywt.wavedec2(imagen, wavelet, level=nivel, mode="periodization"))
    M = int(frac * imagen.size); thr = np.partition(np.abs(arr).ravel(), -M)[-M]
    keep = np.abs(arr) >= thr; ka = keep[sl[0]]
    return dict(banda_aprox_pct_de_N=100 * ka.size / arr.size, retenidos_en_aprox_pct=100 * ka.sum() / keep.sum(), aprox_retenida_pct=100 * ka.mean())
comp = {f"{w} L={L}": composicion_retenidos(img, 0.05, w, L) for w, L in [("db2", 2), ("bior2.2", 2), ("db2", 3), ("db2", 5)]}
display(pd.DataFrame(comp).T.rename(columns={"banda_aprox_pct_de_N": "cA_L ocupa [% de N]", "retenidos_en_aprox_pct": "retenidos que son de cA_L [%]", "aprox_retenida_pct": "cA_L conservada [%]"}).style.format("{:.1f}"))
RES["ej2_composicion_5pct"] = comp
RES["ej2_tabla_C"] = {k: {c: float(v) for c, v in fila.items()} for k, fila in tabla_C.iterrows()}

# Recorte alrededor de la fuente más brillante: reconstrucciones y error absoluto al 5 %
y, x0 = fuentes[0]; ya, xa = np.clip(y-48, 0, Hd-96), np.clip(x0-48, 0, Wd-96)
fig, ax = plt.subplots(2, 4, figsize=(15, 7.8))
for c, (nom, im_) in enumerate([("original", img), ("FFT (5 %)", rec5["FFT"]), ("db2 L=2 (5 %)", rec5["db2 L=2 (enunciado)"]), ("db2 L=5 (5 %)", rec5["db2 L=5"])]):
    ax[0, c].imshow(np.arcsinh(im_[ya:ya+96, xa:xa+96] / 0.05), cmap="magma", vmin=-1.5, vmax=3.7); ax[0, c].set_title(nom); ax[0, c].axis("off")
    if c: ax[1, c].imshow(np.abs(img - im_)[ya:ya+96, xa:xa+96], cmap="inferno", vmin=0, vmax=0.2); ax[1, c].set_title(f"|error|  (MSE = {np.mean((img-im_)**2):.2e})", fontsize=9)
    ax[1, c].axis("off")
plt.tight_layout(); plt.savefig(FIG / "ej2_compresion_5pct.png", bbox_inches="tight"); plt.show()''')
