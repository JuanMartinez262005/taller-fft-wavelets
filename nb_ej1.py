# -*- coding: utf-8 -*-
"""Celdas del Ejercicio 1 (compresión 1D: FFT vs DWT). Cada celda: ("md"|"code", texto)."""

EJ1 = []
md = lambda s: EJ1.append(("md", s))
code = lambda s: EJ1.append(("code", s))

md(r'''---
# Ejercicio 1 · Compresión y reconstrucción 1D: FFT vs. DWT
*(en parejas: Integrante 1 → Partes 1 y 2 · Integrante 2 → Parte 3)*

**Señal elegida.** Electrocardiograma real de la *MIT-BIH Arrhythmia Database* (registro 208, 360 Hz). Es una señal **no estacionaria** por construcción: complejos QRS muy estrechos (transitorios de ~80 ms) sobre ondas P/T lentas y deriva de línea base, con latidos ectópicos intercalados.

**Marco de trabajo.**

* **DFT/FFT:** $\hat x_F[k]=\sum_{n=0}^{N-1}x[n]\,e^{-i2\pi kn/N}$, con Parseval $\;\sum_n|x[n]|^2=\tfrac1N\sum_k|\hat x_F[k]|^2$.
* **DWT** multirresolución con `db4` y $L=5$: $x[n]=\sum_k A_5[k]\phi_{5,k}[n]+\sum_{j=1}^{5}\sum_k D_j[k]\psi_{j,k}[n]$. Se usa `mode="periodization"`, que la hace una transformada **ortogonal exacta** de $N$ coeficientes (con Parseval, igual que la FFT normalizada).
* **Criterio de compresión (idéntico para ambos métodos):** conservar los $M$ coeficientes de mayor energía $|c|^2$ hasta acumular el **95 %** de la energía total y fijar el resto a cero.
* **Convención de conteo en la FFT.** Como $x$ es real, $\hat x_F[N-k]=\overline{\hat x_F[k]}$; los coeficientes se conservan **por pares** $(k,N-k)$ para que la reconstrucción sea real. Cada par equivale a un número complejo (2 reales), de modo que $M_F$ cuenta, en la práctica, tantos números reales como $M_W$ en la DWT y la comparación es justa.
* **Métricas:** $CR=N/M$ y $MSE=\frac1N\sum_n (x[n]-x_{rec}[n])^2$.

> **Observación previa (Parseval).** Como ambas bases son ortonormales y ambos métodos descartan ≈5 % de la energía, **el MSE será prácticamente el mismo por construcción**. La diferencia entre métodos aparece en **cuántos coeficientes (M) hacen falta** y en **dónde cae el error**.''')

code(r'''# ── Datos: ECG real, MIT-BIH Arrhythmia Database (registro 208, canal MLII, 360 Hz) ──
ECG_FUENTE = "local"      # "local": data/ecg_mitbih208.npy (copia de SciPy) · "physionet": descarga oficial con wfdb (pip install wfdb) · "scipy"
FS = 360.0                # Hz

def cargar_ecg(fuente):
    if fuente == "physionet":
        import wfdb                                               # datos públicos: no requiere cuenta ni credenciales
        reg = wfdb.rdrecord("208", pn_dir="mitdb")               # registro oficial completo (30 min, 650 000 muestras)
        return reg.p_signal[:, 0].astype(float), float(reg.fs), "PhysioNet · MIT-BIH 208 (descarga oficial con wfdb)"
    ruta_ecg = DATA / "ecg_mitbih208.npy"
    if fuente == "local" and ruta_ecg.exists():
        return np.load(ruta_ecg), 360.0, "MIT-BIH 208 · fragmento de 5 min (copia de SciPy incluida en data/)"
    from scipy import datasets                                    # misma señal que empaqueta SciPy (necesita red y pooch)
    return datasets.electrocardiogram(), 360.0, "MIT-BIH 208 · fragmento de 5 min (scipy.datasets)"

ecg, FS, ORIGEN_ECG = cargar_ecg(ECG_FUENTE)
RES["ej1_fuente"] = ORIGEN_ECG

N = 8192                                     # N ≥ 4096; potencia de 2 (≈ 22.8 s)
x = ecg[:N].astype(float)
t = np.arange(N) / FS
print(f"{ORIGEN_ECG}\nSeñal disponible: {len(ecg)} muestras ({len(ecg)/FS/60:.1f} min a {FS:.0f} Hz)")
print(f"Segmento de trabajo: N = {N} muestras ({N/FS:.2f} s) | media = {x.mean():+.3f} mV | σ = {x.std():.3f} mV")

fig, ax = plt.subplots(figsize=(13, 2.8))
ax.plot(t, x, "k", lw=.7)
ax.set(xlabel="t [s]", ylabel="ECG [mV]", title=f"x[n]: ECG MIT-BIH, N = {N}")
ax.grid(alpha=.25); plt.tight_layout(); plt.show()''')

code(r'''# ── Utilidades de compresión (mismo criterio de energía para FFT y DWT) ──
def m_para_energia(E_ordenada, frac):
    """Menor M tal que los M coeficientes de mayor energía acumulan >= frac de la energía total."""
    acum = np.cumsum(E_ordenada) / E_ordenada.sum()
    return int(np.searchsorted(acum, frac) + 1)

def comprimir_fft(x, frac=None, M=None):
    """FFT -> conservar los M coeficientes de mayor |X_k|^2 (o los necesarios para 'frac' de la energía)."""
    N = len(x)
    X = fft(x)                               # X[k] = Σ x[n] e^{-i2πkn/N}
    E = np.abs(X)**2
    orden = np.argsort(E)[::-1]
    if M is None:
        M = m_para_energia(E[orden], frac)
    keep = np.zeros(N, bool); keep[orden[:M]] = True
    idx = np.flatnonzero(keep); keep[(-idx) % N] = True     # señal real: se conservan k y N-k juntos
    return np.where(keep, X, 0), keep

def comprimir_dwt(x, frac=None, M=None, wavelet="db4", nivel=5):
    """DWT multirresolución -> umbral global por energía. Devuelve (coefs_comprimidos, arr, arr_comprimido, slices)."""
    coefs = pywt.wavedec(x, wavelet, level=nivel, mode="periodization")      # [A5, D5, D4, D3, D2, D1]
    arr, sl = pywt.coeffs_to_array(coefs)
    E = arr**2
    orden = np.argsort(E)[::-1]
    if M is None:
        M = m_para_energia(E[orden], frac)
    arr_c = np.zeros_like(arr); arr_c[orden[:M]] = arr[orden[:M]]
    return pywt.array_to_coeffs(arr_c, sl, output_format="wavedec"), arr, arr_c, sl''')

md(r'''### Parte 1 · Integrante 1 — FFT: número mínimo de coeficientes para el 95 % de la energía''')

code(r'''X = fft(x)                                   # x̂_F[k]
E_total = np.sum(x**2)
print(f"Parseval:  Σ|x[n]|² = {E_total:.4f}   |   (1/N) Σ|x̂_F[k]|² = {np.sum(np.abs(X)**2)/N:.4f}")

Xc, keep_F = comprimir_fft(x, frac=0.95)     # vector comprimido  x̃_F  (N-M_F coeficientes puestos a cero)
M_F = int(keep_F.sum())
E_ret_F = np.sum(np.abs(Xc)**2) / N / E_total
f_Hz = fftfreq(N, 1/FS)
f_max_F = np.abs(f_Hz[keep_F]).max()
dc_pct = 100 * N * x.mean()**2 / E_total
print(f"M_F = {M_F} de {N} coeficientes ({100*M_F/N:.2f} %)  →  energía retenida = {100*E_ret_F:.2f} %")
print(f"Frecuencia máxima conservada: {f_max_F:.1f} Hz   |   la componente DC aporta el {dc_pct:.1f} % de la energía total")
RES.update(ej1_M_F=M_F, ej1_Eret_F=100*E_ret_F, ej1_fmax_F=float(f_max_F), ej1_dc_pct=float(dc_pct))''')

md(r'''### Parte 2 · Integrante 1 — DWT multirresolución (db4, $L=5$)''')

code(r'''WAV, NIVEL = "db4", 5
coefs_c, arr, arr_c, slices = comprimir_dwt(x, frac=0.95, wavelet=WAV, nivel=NIVEL)
M_W = int(np.count_nonzero(arr_c))
E_ret_W = np.sum(arr_c**2) / np.sum(arr**2)
print(f"Parseval DWT:  Σ|x[n]|² = {E_total:.4f}   |   Σ coef² = {np.sum(arr**2):.4f}")
print(f"M_W = {M_W} de {N} coeficientes ({100*M_W/N:.2f} %)  →  energía retenida = {100*E_ret_W:.2f} %")

# Anatomía de la descomposición: banda de frecuencia, energía y coeficientes retenidos por subbanda
nombres = ["A5"] + [f"D{j}" for j in range(NIVEL, 0, -1)]
coefs_full = pywt.wavedec(x, WAV, level=NIVEL, mode="periodization")
filas_sb = []
for nom, c, cq in zip(nombres, coefs_full, coefs_c):
    j = NIVEL if nom == "A5" else int(nom[1])
    lo, hi = (0.0, FS/2**(NIVEL+1)) if nom == "A5" else (FS/2**(j+1), FS/2**j)
    filas_sb.append({"subbanda": nom, "banda [Hz]": f"{lo:.1f}–{hi:.1f}", "n.º coef.": len(c),
                  "energía [%]": 100*np.sum(c**2)/E_total, "retenidos": int(np.count_nonzero(cq))})
tabla_dwt = pd.DataFrame(filas_sb)
display(tabla_dwt)
RES.update(ej1_M_W=M_W, ej1_Eret_W=100*E_ret_W,
           ej1_niveles={r["subbanda"]: dict(energia=float(r["energía [%]"]), retenidos=int(r["retenidos"])) for r in filas_sb})

# "Entrega" de los vectores comprimidos al Integrante 2
np.savez(DATA / "ej1_comprimidos.npz", Xc=Xc, arr_c=arr_c)''')

md(r'''### Parte 3 · Integrante 2 — Reconstrucción, residuos, CR y MSE''')

code(r'''# ── Integrante 2: recibe los vectores comprimidos y aplica las transformadas inversas ──
datos = np.load(DATA / "ej1_comprimidos.npz")
x_F = np.real(ifft(datos["Xc"]))                                                          # iFFT
x_W = pywt.waverec(pywt.array_to_coeffs(datos["arr_c"], slices, output_format="wavedec"),
                   WAV, mode="periodization")                                              # iDWT
e_F, e_W = x - x_F, x - x_W                                                                # residuos e[n] = x[n] − x_rec[n]

# Señal original, reconstrucciones y residuos en una misma figura (columna izq.: completa · derecha: zoom de 3 s)
fig = plt.figure(figsize=(15, 9.5))
gs = gridspec.GridSpec(3, 2, height_ratios=[1.6, 1, 1], width_ratios=[2.3, 1.5], hspace=.12, wspace=.1)
lim_e = 1.08 * max(np.abs(e_F).max(), np.abs(e_W).max())
for col, (xl, tit) in enumerate([((0, N/FS), f"Señal completa ({N/FS:.1f} s)"), ((0.3, 3.3), "Zoom de 3 s")]):
    sel = (t >= xl[0]) & (t <= xl[1])
    a0 = fig.add_subplot(gs[0, col]); a1 = fig.add_subplot(gs[1, col], sharex=a0); a2 = fig.add_subplot(gs[2, col], sharex=a0)
    a0.plot(t, x, "k", lw=1.0, label="x[n] original")
    a0.plot(t, x_F, color="C0", lw=.9, label=f"x_rec,F  (FFT, M = {M_F})")
    a0.plot(t, x_W, color="C3", lw=.9, label=f"x_rec,W  (DWT {WAV}, L = {NIVEL}, M = {M_W})")
    a1.plot(t, e_F, color="C0", lw=.8); a2.plot(t, e_W, color="C3", lw=.8)
    ymin = min(x[sel].min(), x_F[sel].min(), x_W[sel].min()); ymax = max(x[sel].max(), x_F[sel].max(), x_W[sel].max())
    a0.set_xlim(*xl); a0.set_ylim(ymin - .06*(ymax-ymin), ymax + .22*(ymax-ymin))
    for a in (a1, a2): a.set_ylim(-lim_e, lim_e)
    a0.set_title(tit); a2.set_xlabel("t [s]")
    plt.setp(a0.get_xticklabels(), visible=False); plt.setp(a1.get_xticklabels(), visible=False)
    if col == 0:
        a0.legend(loc="upper right", ncol=3, fontsize=9); a0.set_ylabel("ECG [mV]")
        a1.set_ylabel(r"$e_F[n]=x-x_{rec,F}$ [mV]"); a2.set_ylabel(r"$e_W[n]=x-x_{rec,W}$ [mV]")
    for a in (a0, a1, a2): a.grid(alpha=.25)
plt.savefig(FIG / "ej1_reconstruccion.png", bbox_inches="tight"); plt.show()''')

code(r'''def metricas_1d(x, xr, M, nombre):
    e = x - xr; mse = np.mean(e**2)
    return {"método": nombre, "M": M, "CR = N/M": len(x)/M, "MSE [mV²]": mse,
            "RMSE [mV]": np.sqrt(mse), "SNR [dB]": 10*np.log10(np.sum(x**2)/np.sum(e**2))}

tabla_1d = pd.DataFrame([metricas_1d(x, x_F, M_F, "FFT (95 % energía)"),
                         metricas_1d(x, x_W, M_W, f"DWT {WAV} L={NIVEL} (95 % energía)")]).set_index("método")
display(tabla_1d)

# Comprobación por Parseval: MSE = (energía descartada)/N  (ambas transformadas son ortonormales)
mse_teo_F = np.sum(np.abs(X[~keep_F])**2) / N**2
mse_teo_W = np.sum(arr[arr_c == 0]**2) / N
print(f"MSE medido vs. predicho por Parseval →  FFT: {tabla_1d.iloc[0, 2]:.6f} vs {mse_teo_F:.6f}   |   DWT: {tabla_1d.iloc[1, 2]:.6f} vs {mse_teo_W:.6f}")
RES.update(ej1_CR_F=N/M_F, ej1_CR_W=N/M_W, ej1_MSE_F=float(tabla_1d.iloc[0, 2]), ej1_MSE_W=float(tabla_1d.iloc[1, 2]),
           ej1_SNR_F=float(tabla_1d.iloc[0, 4]), ej1_SNR_W=float(tabla_1d.iloc[1, 4]))''')

md(r'''#### ¿Dónde cae el error? Fidelidad en los complejos QRS (transitorios)
Para valorar la preservación de **eventos transitorios** se detectan los picos R y se mide, en ventanas de ±60 ms alrededor de cada uno: (i) cuánto de la **amplitud del pico R** conserva cada método, y (ii) el **RMS del error** dentro y fuera de esas ventanas.''')

code(r'''def picos_R(x, umbral=0.8, sep=0.25):
    pk, _ = find_peaks(np.abs(x - np.median(x)), height=umbral, distance=int(sep*FS))
    return pk

def fidelidad_transitorios(x, xF, xW, semi=0.06):
    w = int(semi * FS); pk = picos_R(x); pk = pk[(pk > w) & (pk < len(x) - w)]
    med = np.median(x); mask = np.zeros(len(x), bool)
    for p in pk: mask[p-w:p+w+1] = True
    out = {"n_picos": len(pk), "frac_QRS": mask.mean()}
    for nom, xr in (("F", xF), ("W", xW)):
        e = x - xr
        out[f"retR_{nom}"]    = np.mean(np.abs(xr[pk] - med) / np.abs(x[pk] - med))   # retención de la amplitud del pico R
        out[f"rmsQRS_{nom}"]  = np.sqrt(np.mean(e[mask]**2))                          # RMS del error dentro de ventanas QRS
        out[f"rmsResto_{nom}"] = np.sqrt(np.mean(e[~mask]**2))                        # RMS del error fuera de ellas
    return out

d0 = fidelidad_transitorios(x, x_F, x_W)
display(pd.DataFrame({"FFT": [d0["retR_F"], d0["rmsQRS_F"], d0["rmsResto_F"]],
                      "DWT": [d0["retR_W"], d0["rmsQRS_W"], d0["rmsResto_W"]]},
                     index=["retención amplitud pico R", "RMS error en ventanas QRS [mV]", "RMS error fuera de QRS [mV]"]))
RES["ej1_seg0"] = {k: float(v) for k, v in d0.items()}
print(f"({d0['n_picos']} picos R detectados; las ventanas QRS ocupan el {100*d0['frac_QRS']:.0f} % de la señal)")''')

md(r'''#### Robustez: los mismos cálculos en los 13 segmentos no solapados del registro (N = 8192)
Para que las conclusiones no dependan del segmento elegido, se repite todo en cada segmento. Además se compara la FFT **con el mismo número de coeficientes que usa la DWT** (a igual $M$).''')

code(r'''filas = []
for s in range(0, len(ecg) - N + 1, N):
    xs = ecg[s:s+N].astype(float)
    Xs, kF = comprimir_fft(xs, frac=0.95);  xF_s = np.real(ifft(Xs))
    cW, a_s, ac_s, sl_s = comprimir_dwt(xs, frac=0.95, wavelet=WAV, nivel=NIVEL)
    xW_s = pywt.waverec(cW, WAV, mode="periodization")
    MF_s, MW_s = int(kF.sum()), int(np.count_nonzero(ac_s))
    Xq, kq = comprimir_fft(xs, M=MW_s);     xFq = np.real(ifft(Xq))          # FFT con (≈) el mismo M que la DWT
    filas.append(dict(seg=s//N, M_F=MF_s, M_W=MW_s, M_Feq=int(kq.sum()),
                      MSE_F=np.mean((xs-xF_s)**2), MSE_W=np.mean((xs-xW_s)**2), MSE_F_igualM=np.mean((xs-xFq)**2),
                      **fidelidad_transitorios(xs, xF_s, xW_s)))
rob = pd.DataFrame(filas)
rob["CR_F"], rob["CR_W"] = N/rob.M_F, N/rob.M_W
display(rob[["seg", "M_F", "M_W", "CR_F", "CR_W", "MSE_F", "MSE_W", "retR_F", "retR_W", "rmsQRS_F", "rmsQRS_W", "rmsResto_F", "rmsResto_W"]]
        .style.format(precision=4).hide(axis="index"))

n = len(rob)
resumen = dict(
    n_seg=n, M_F_mean=rob.M_F.mean(), M_F_std=rob.M_F.std(), M_W_mean=rob.M_W.mean(), M_W_std=rob.M_W.std(),
    CR_F_mean=rob.CR_F.mean(), CR_W_mean=rob.CR_W.mean(), wins_M=int((rob.M_W < rob.M_F).sum()),
    ratio_M=float((rob.M_F / rob.M_W).mean()),
    mseratio_igualM=float((rob.MSE_F_igualM / rob.MSE_W).mean()), wins_igualM=int((rob.MSE_W < rob.MSE_F_igualM).sum()),
    Meq_over_MW=float((rob.M_Feq / rob.M_W).mean()),
    mse_gap_max=float((np.abs(rob.MSE_F - rob.MSE_W) / rob.MSE_F).max()),
    retR_F=rob.retR_F.mean(), retR_W=rob.retR_W.mean(), wins_retR=int((rob.retR_W > rob.retR_F).sum()),
    rmsQRS_F=rob.rmsQRS_F.mean(), rmsQRS_W=rob.rmsQRS_W.mean(), wins_QRS=int((rob.rmsQRS_W < rob.rmsQRS_F).sum()),
    rmsResto_F=rob.rmsResto_F.mean(), rmsResto_W=rob.rmsResto_W.mean(), wins_resto=int((rob.rmsResto_W < rob.rmsResto_F).sum()))
RES["ej1_robustez"] = {k: float(v) for k, v in resumen.items()}
print(f"M (95 % energía):  FFT = {resumen['M_F_mean']:.0f} ± {resumen['M_F_std']:.0f}   |   DWT = {resumen['M_W_mean']:.0f} ± {resumen['M_W_std']:.0f}   →  la DWT necesita menos coeficientes en {resumen['wins_M']}/{n} segmentos")
print(f"CR medio: FFT = {resumen['CR_F_mean']:.1f}   |   DWT = {resumen['CR_W_mean']:.1f}")
print(f"Diferencia relativa máxima de MSE entre métodos a igual energía: {100*resumen['mse_gap_max']:.2f} %  (Parseval)")
print(f"A igual M (FFT usa ≈{resumen['Meq_over_MW']:.2f}·M_W):  MSE_FFT / MSE_DWT = {resumen['mseratio_igualM']:.2f}  (DWT menor en {resumen['wins_igualM']}/{n})")
print(f"Retención de la amplitud del pico R:  FFT = {resumen['retR_F']:.3f}   |   DWT = {resumen['retR_W']:.3f}   (DWT mayor en {resumen['wins_retR']}/{n})")
print(f"RMS del error en ventanas QRS [mV]:    FFT = {resumen['rmsQRS_F']:.4f}  |   DWT = {resumen['rmsQRS_W']:.4f}   (DWT menor en {resumen['wins_QRS']}/{n})")
print(f"RMS del error fuera de QRS [mV]:       FFT = {resumen['rmsResto_F']:.4f}  |   DWT = {resumen['rmsResto_W']:.4f}   (DWT menor en {resumen['wins_resto']}/{n})")''')

md(r'''#### Compromiso compresión–error: variando el porcentaje de energía retenida
Se barre el umbral de energía (80 %–99.9 %) en el segmento de trabajo y se compara el número de coeficientes y el MSE.''')

code(r'''fracs = np.array([0.80, 0.90, 0.95, 0.98, 0.99, 0.995, 0.999])
curva = []
for fr in fracs:
    Xq, kq = comprimir_fft(x, frac=fr); cq, aa, aq, ss = comprimir_dwt(x, frac=fr, wavelet=WAV, nivel=NIVEL)
    curva.append(dict(energía=fr, M_F=int(kq.sum()), M_W=int(np.count_nonzero(aq)),
                      MSE_F=np.mean((x - np.real(ifft(Xq)))**2),
                      MSE_W=np.mean((x - pywt.waverec(cq, WAV, mode="periodization"))**2)))
curva = pd.DataFrame(curva)
curva["CR_F"], curva["CR_W"] = N/curva.M_F, N/curva.M_W
display(curva.style.format({"energía": "{:.3f}", "MSE_F": "{:.2e}", "MSE_W": "{:.2e}", "CR_F": "{:.1f}", "CR_W": "{:.1f}"}).hide(axis="index"))

from matplotlib.ticker import FuncFormatter
fmt = FuncFormatter(lambda v, _: f"{v:g}")
fig, ax = plt.subplots(1, 2, figsize=(12, 4))
ax[0].semilogy(100*curva.energía, curva.M_F, "o-", color="C0", label="FFT"); ax[0].semilogy(100*curva.energía, curva.M_W, "s-", color="C3", label=f"DWT {WAV}")
ax[0].set(xlabel="energía retenida [%]", ylabel="coeficientes M", title="Coeficientes necesarios"); ax[0].legend(); ax[0].grid(alpha=.3, which="both")
ax[1].loglog(curva.CR_F, curva.MSE_F, "o-", color="C0", label="FFT"); ax[1].loglog(curva.CR_W, curva.MSE_W, "s-", color="C3", label=f"DWT {WAV}")
ax[1].set(xlabel="razón de compresión CR = N/M", ylabel="MSE [mV²]", title="Error vs. compresión"); ax[1].legend(); ax[1].grid(alpha=.3, which="both")
for a_ in ax: a_.yaxis.set_major_formatter(fmt)
ax[1].xaxis.set_major_formatter(fmt); ax[1].set_xticks([4, 6, 10, 20, 40, 60])
plt.tight_layout(); plt.savefig(FIG / "ej1_curva_CR_MSE.png", bbox_inches="tight"); plt.show()

# Mismo MSE ≈ mismo tramo de energía, pero con menos coeficientes: CR de la DWT respecto a la FFT por nivel de energía
RES["ej1_curva"] = {f"{100*r.energía:.1f}": dict(M_F=int(r.M_F), M_W=int(r.M_W)) for r in curva.itertuples()}''')

md(r'''#### ¿Qué coeficientes se conservan? Frecuencia (FFT) frente a tiempo-escala (DWT)
Esta figura conecta el resultado con el principio de incertidumbre: la FFT solo informa de *qué frecuencias* hay; la DWT dice además *cuándo*.''')

code(r'''fig, axs = plt.subplots(1, 2, figsize=(14, 4.3), gridspec_kw={"width_ratios": [1, 1.5]})
# (a) FFT: espectro y coeficientes conservados
mitad = slice(0, N//2)
P = 10 * np.log10(np.abs(X)**2 / N + 1e-12)
axs[0].plot(f_Hz[mitad], P[mitad], color="0.75", lw=.6, label="todos los coeficientes")
kk = keep_F[mitad]
axs[0].plot(f_Hz[mitad][kk], P[mitad][kk], "o", color="C0", ms=3, label=f"conservados (M_F = {M_F}, contando k y N−k)")
axs[0].set(xlim=(0, 60), xlabel="frecuencia [Hz]", ylabel=r"$|\hat x_F[k]|^2/N$ [dB]", title="(a) FFT: conserva líneas espectrales en TODA la señal")
axs[0].legend(loc="upper right", fontsize=8); axs[0].grid(alpha=.25)

# (b) DWT: dónde (en el tiempo) viven los coeficientes conservados, por subbanda
for tp in picos_R(x): axs[1].axvline(tp/FS, color="0.85", lw=.6, zorder=0)
for i, (nom, cq) in enumerate(zip(nombres, coefs_c)):
    j = NIVEL if nom == "A5" else int(nom[1])
    kidx = np.flatnonzero(cq)
    axs[1].plot((kidx + .5) * 2**j / FS, np.full(len(kidx), i), "|", color="C3", ms=11, mew=1.3)
axs[1].set_yticks(range(len(nombres))); axs[1].set_yticklabels([f"{r['subbanda']}  ({r['banda [Hz]']} Hz)" for r in filas_sb], fontsize=8)
axs[1].set(xlim=(0, N/FS), xlabel="t [s]", title=f"(b) DWT (M_W = {M_W}): D3–D4 sólo aparecen en los complejos QRS; D2 y D1 se descartan (líneas grises = picos R)")
axs[1].grid(alpha=.25, axis="x")
plt.tight_layout(); plt.savefig(FIG / "ej1_coeficientes_conservados.png", bbox_inches="tight"); plt.show()

# ¿Cuántos coeficientes conservados de cada subbanda caen dentro de ventanas QRS (±60 ms)?
w = int(0.06 * FS); pkq = picos_R(x); pkq = pkq[(pkq > w) & (pkq < N - w)]
mask_qrs = np.zeros(N, bool)
for p in pkq: mask_qrs[p-w:p+w+1] = True
alin = {}
for nom, cq in zip(nombres, coefs_c):
    j = NIVEL if nom == "A5" else int(nom[1]); kidx = np.flatnonzero(cq)
    if len(kidx):
        pos = np.clip(np.round((kidx + .5) * 2**j).astype(int), 0, N - 1)
        alin[nom] = dict(retenidos=int(len(kidx)), en_QRS_pct=float(100 * mask_qrs[pos].mean()))
print(f"Las ventanas QRS ocupan el {100*mask_qrs.mean():.0f} % del tiempo. Coeficientes conservados que caen dentro de ellas:")
for nom, d in alin.items(): print(f"  {nom}: {d['en_QRS_pct']:5.1f} %  ({d['retenidos']} conservados)")
RES["ej1_alineacion"] = dict(frac_tiempo_QRS_pct=float(100 * mask_qrs.mean()), **alin)''')

code(r'''# Esquema de teselaciones tiempo–frecuencia (todas con celdas de área Δt·Δω constante ≥ 1/2)
def tesela(ax, celdas, titulo, colores=None):
    for i, (x0, y0, w, h) in enumerate(celdas):
        ax.add_patch(Rectangle((x0, y0), w, h, fc=(colores[i] if colores else "none"), ec="k", lw=.8, alpha=.85))
    ax.set(xlim=(0, 1), ylim=(0, 1), title=titulo, xlabel="tiempo"); ax.set_xticks([]); ax.set_yticks([])

fig, axs = plt.subplots(1, 4, figsize=(14, 3.2))
tesela(axs[0], [(i/16, 0, 1/16, 1) for i in range(16)], r"Serie temporal" "\n" r"($\Delta t$ mínimo, $\Delta\omega$ máximo)"); axs[0].set_ylabel("frecuencia")
tesela(axs[1], [(0, i/16, 1, 1/16) for i in range(16)], r"Fourier" "\n" r"($\Delta\omega\to 0$, $\Delta t\to T$)")
tesela(axs[2], [(i/4, k/4, 1/4, 1/4) for i in range(4) for k in range(4)], "STFT / Gabor\n(ventana fija)")
cel, col = [], []
paleta = ["#fde68a", "#fdba74", "#fca5a5", "#c4b5fd", "#93c5fd"]
for j, (lo, hi) in enumerate([(1/2, 1), (1/4, 1/2), (1/8, 1/4), (1/16, 1/8)], start=1):
    for i in range(2**(4-j)):
        cel.append((i/2**(4-j), lo, 1/2**(4-j), hi-lo)); col.append(paleta[j-1])
cel.append((0, 0, 1, 1/16)); col.append(paleta[4])
tesela(axs[3], cel, r"DWT (multirresolución)" "\n" r"($\Delta t\propto 2^{j}$, $\Delta\omega\propto 2^{-j}$)", col)
plt.tight_layout(); plt.savefig(FIG / "ej1_teselaciones.png", bbox_inches="tight"); plt.show()''')
