# -*- coding: utf-8 -*-
from disc_util import fill, lim_img, desc_fuentes, fwhm_fuentes, vs_fft

W = "Wavelet db2 L=2 · umbral suave (enunciado)"
G = "FFT pasa-altas gaussiano (enunciado)"
I = "FFT pasa-altas ideal (corte abrupto)"


def disc2A(R):
    R = dict(R)
    A = R["ej2_tabla_A"]
    R["d"] = dict(
        lim_img=lim_img(R), desc_fuentes=desc_fuentes(R),
        cielo=A[G]["RMS en cielo vacío"] / A[W]["RMS en cielo vacío"],
        ruido_fft=100 * A[G]["RMS en cielo vacío"] / R["ej2_sig_fondo"],
        ruido_wav=100 * A[W]["RMS en cielo vacío"] / R["ej2_sig_fondo"],
        lam_k=R["ej2_lambda"] / R["ej2_sigma_w"],
    )
    # --- frases que dependen de los números (no se asume el resultado del caso Hubble) ---
    d = R["d"]; r = d["cielo"]; Wa, Ga = A[W], A[G]
    d["fondo_txt"] = ("un **fondo mucho más limpio**" if r >= 2 else "un **fondo algo más limpio**" if r >= 1.15
                      else "un fondo **de limpieza similar**" if r >= 1 / 1.15 else "un fondo **menos limpio**")
    d["fondo_gana"] = f"aquí gana la wavelet, por un factor ≈ {r:.1f}" if r > 1 else f"aquí **no** gana la wavelet (factor {r:.2f})"
    fw = fwhm_fuentes(R)
    dif_pos = Ga["posición ≤1 px"] - Wa["posición ≤1 px"]
    if min(Wa["posición ≤1 px"], Ga["posición ≤1 px"]) >= 0.85:
        d["pos_txt"] = "ambos métodos localizan bien las fuentes; el filtro gaussiano " + ("algo mejor" if dif_pos < 0.10 else "mejor")
    else:
        d["pos_txt"] = ("**ninguno de los dos localiza el máximo con precisión en todas las fuentes**"
                        + (f" (con fuentes de ≈ {fw:.1f} px de FWHM, un desplazamiento de 1–2 px del máximo es fácil)" if fw and fw > 4.5 else "")
                        + "; el filtro gaussiano " + ("algo mejor" if 0 <= dif_pos < 0.10 else "mejor" if dif_pos >= 0.10 else "peor"))
    d["pico_txt"] = f"conserva en media el {Wa['retención del pico']:.0%} del pico (la FFT gaussiana, el {Ga['retención del pico']:.0%})"

    S = R["ej2_tabla_sens"]; wk = [k for k in S if k.startswith("Wavelet")]; fk = [k for k in S if k.startswith("FFT")]
    w_corr, w_sky = max(S[k]["correlación de forma"] for k in wk), max(S[k]["RMS en cielo vacío"] for k in wk)
    w_sky_min = min(S[k]["RMS en cielo vacío"] for k in wk)
    f_corr_min, f_sky_min = min(S[k]["correlación de forma"] for k in fk), min(S[k]["RMS en cielo vacío"] for k in fk)
    k4, k16 = "FFT gaussiano · f_c = 1/4", "FFT gaussiano · f_c = 1/16"
    lam_rango = rf"$\lambda\in[2,\,{d['lam_k']:.1f}]\,\hat\sigma$"
    if w_corr < f_corr_min and w_sky < f_sky_min:
        d["param_txt"] = (rf"**¿Es un artefacto de los parámetros?** No: la tabla de sensibilidad muestra que, para {lam_rango}, la wavelet queda siempre con correlación ≤ {w_corr:.2f} y RMS de cielo ≤ {w_sky:.4f}, "
                          rf"mientras que el filtro gaussiano con $f_c\in[1/16,\,1/4]$ obtiene correlaciones de {S[k4]['correlación de forma']:.2f} a {S[k16]['correlación de forma']:.2f} pero RMS de cielo ≥ {f_sky_min:.4f}. "
                          rf"Hay un compromiso en la FFT: bajando $f_c$ mejora la forma ({S[k16]['correlación de forma']:.2f} con $f_c=1/16$) pero deja más ruido en el cielo.")
    elif w_corr < f_corr_min:
        cmp_cielo = "más limpio" if S[k4]["RMS en cielo vacío"] < w_sky_min else "comparable"
        d["param_txt"] = (rf"**¿Es un artefacto de los parámetros?** **Sólo en parte.** La conclusión sobre la **forma** no depende de ellos: para {lam_rango} la wavelet queda siempre con correlación ≤ {w_corr:.2f}, "
                          rf"mientras que el filtro gaussiano con $f_c\in[1/16,\,1/4]$ obtiene de {S[k4]['correlación de forma']:.2f} a {S[k16]['correlación de forma']:.2f}. "
                          rf"La del **cielo sí depende**: la wavelet queda entre {w_sky_min:.4f} y {w_sky:.4f} de RMS y el gaussiano entre {f_sky_min:.4f} (con $f_c=1/4$) y {S[k16]['RMS en cielo vacío']:.4f} (con $f_c=1/16$); "
                          rf"con $f_c=1/4$ el filtro gaussiano deja un cielo {cmp_cielo} que la wavelet, a cambio de conservar menos el pico ({S[k4]['retención del pico']:.2f}) y de localizar peor ({S[k4]['posición ≤1 px']:.0%}). "
                          "Hay, pues, un compromiso en la FFT: bajando $f_c$ mejora la forma pero deja más ruido en el cielo, y la ventaja de la wavelet en el fondo desaparece si se sube $f_c$.")
    else:
        d["param_txt"] = (rf"**¿Es un artefacto de los parámetros?** **No concluyente**: con estos datos los rangos de la wavelet (correlación ≤ {w_corr:.2f}, cielo ≤ {w_sky:.4f}) y del filtro gaussiano "
                          f"(correlación ≥ {f_corr_min:.2f}, cielo ≥ {f_sky_min:.4f}) se solapan; la conclusión depende de los parámetros (véase la tabla de sensibilidad).")

    if r > 1.15 and Wa["posición ≤1 px"] >= 0.85:
        d["ver_txt"] = (f"Para **detectar** fuentes (cielo limpio y posición correcta en el {Wa['posición ≤1 px']:.0%} de los casos) la wavelet es preferible; "
                        "para **medir forma y flujo** de las fuentes, el pasa-altas gaussiano.")
    elif r > 1.15:
        d["ver_txt"] = (f"La wavelet deja el cielo más limpio (×{r:.1f}) pero localiza peor que la FFT ({Wa['posición ≤1 px']:.0%} frente a {Ga['posición ≤1 px']:.0%} de máximos a ≤ 1 px), "
                        "así que su ventaja para **detectar** se limita a la limpieza del fondo; para **localizar y medir forma y flujo** de las fuentes, el pasa-altas gaussiano.")
    else:
        d["ver_txt"] = ("En esta imagen la wavelet no aporta una ventaja clara en la limpieza del fondo; para **localizar y medir forma y flujo** de las fuentes, el pasa-altas gaussiano.")
    dH = A["Wavelet db2 L=2 · umbral duro"]
    d["duro_txt"] = (" Aquí el cambio es pequeño, así que el sesgo del umbral suave **no** es la causa principal de la deformación."
                     if dH["retención del pico"] - Wa["retención del pico"] < 0.15
                     else " Esto confirma que el sesgo del umbral suave explica buena parte de la deformación.")
    d["tamano_txt"] = ("" if not fw or fw <= 4.5 else
                       rf"**Tamaño de las fuentes.** Las fuentes de referencia miden ≈ {fw:.1f} px de FWHM (área mediana {R['ej2_criterio_fuentes']['area_mediana']:.0f} px), frente a 1–3 px de una estrella a la escala de píxel natural del detector: "
                       rf"buena parte de su energía queda por debajo de $f_c=1/8$ ciclos/px, que ambos métodos eliminan (por eso la retención del pico es baja en los dos: {Wa['retención del pico']:.2f} y {Ga['retención del pico']:.2f}). "
                       "Con fuentes más finas (por ejemplo, promediando la imagen por bloques con `--reducir`) los números cambiarían; las comparaciones de esta pregunta valen para esta escala de fuentes.")
    return fill(r'''---
## Pregunta A — ¿qué método conserva con mayor fidelidad la forma y la posición de las fuentes sin alterar el fondo oscuro?

**Respuesta corta.** **Ninguno gana en todo**; depende de qué se entienda por "fidelidad". La **2D-FFT con pasa-altas gaussiano** conserva mejor la **forma** y la **posición** de las fuentes; la **wavelet** deja @@d>fondo_txt@@ pero deforma más las fuentes. Medido sobre @@ej2_n_fuentes@@ @@d>desc_fuentes@@:

| | Wavelet db2 L=2, umbral suave | FFT gaussiano ($f_c=1/8$) | FFT ideal |
|---|---|---|---|
| Posición del máximo a ≤ 1 px | @@ej2_tabla_A>Wavelet db2 L=2 · umbral suave (enunciado)>posición ≤1 px|.0%@@ | @@ej2_tabla_A>FFT pasa-altas gaussiano (enunciado)>posición ≤1 px|.0%@@ | @@ej2_tabla_A>FFT pasa-altas ideal (corte abrupto)>posición ≤1 px|.0%@@ |
| Correlación de forma (recorte 11×11) | @@ej2_tabla_A>Wavelet db2 L=2 · umbral suave (enunciado)>correlación de forma|.2f@@ | @@ej2_tabla_A>FFT pasa-altas gaussiano (enunciado)>correlación de forma|.2f@@ | @@ej2_tabla_A>FFT pasa-altas ideal (corte abrupto)>correlación de forma|.2f@@ |
| Retención del pico | @@ej2_tabla_A>Wavelet db2 L=2 · umbral suave (enunciado)>retención del pico|.2f@@ | @@ej2_tabla_A>FFT pasa-altas gaussiano (enunciado)>retención del pico|.2f@@ | @@ej2_tabla_A>FFT pasa-altas ideal (corte abrupto)>retención del pico|.2f@@ |
| RMS en cielo vacío | @@ej2_tabla_A>Wavelet db2 L=2 · umbral suave (enunciado)>RMS en cielo vacío|.4f@@ | @@ej2_tabla_A>FFT pasa-altas gaussiano (enunciado)>RMS en cielo vacío|.4f@@ | @@ej2_tabla_A>FFT pasa-altas ideal (corte abrupto)>RMS en cielo vacío|.4f@@ |

* **Posición:** @@d>pos_txt@@ (@@ej2_tabla_A>FFT pasa-altas gaussiano (enunciado)>posición ≤1 px|.0%@@ frente a @@ej2_tabla_A>Wavelet db2 L=2 · umbral suave (enunciado)>posición ≤1 px|.0%@@). Una posible causa de la menor precisión de la wavelet es que la DWT trabaja en una rejilla diádica (cada 4 px con L = 2).
* **Forma:** la FFT es claramente mejor (correlación @@ej2_tabla_A>FFT pasa-altas gaussiano (enunciado)>correlación de forma|.2f@@ frente a @@ej2_tabla_A>Wavelet db2 L=2 · umbral suave (enunciado)>correlación de forma|.2f@@). En los recortes y perfiles se ve que la wavelet @@d>pico_txt@@ y lo deja "pixelado".
* **Fondo oscuro:** @@d>fondo_gana@@: su RMS en cielo vacío es el @@d>ruido_wav|.0f@@ % de la dispersión robusta (MAD) de la imagen original ($\hat\sigma$ = @@ej2_sig_fondo|.4f@@; en un campo vacío es el ruido del fondo, en una imagen con nebulosa incluye su estructura) frente al @@d>ruido_fft|.0f@@ % con la FFT. El ruido del fondo tiene energía en todas las frecuencias, de modo que un pasa-altas sólo elimina su parte de baja frecuencia y deja casi todo lo demás; la umbralización suave, en cambio, pone a cero todos los coeficientes menores que $\lambda=$ @@ej2_lambda|.4f@@ (@@d>lam_k|.1f@@ $\hat\sigma_w$, con $\hat\sigma_w=$ @@ej2_sigma_w|.4f@@ estimada con MAD en $cD_1$), y sólo sobrevive el @@ej2_sobrev_pct|.2f@@ % de los coeficientes de detalle.

@@d>param_txt@@

**Por qué la wavelet deforma más las fuentes.** (i) El umbral suave **resta $\lambda$ a cada coeficiente que sobrevive**, lo que encoge sistemáticamente las fuentes (sesgo conocido de la umbralización suave): con umbral duro la retención del pico sube de @@ej2_tabla_A>Wavelet db2 L=2 · umbral suave (enunciado)>retención del pico|.2f@@ a @@ej2_tabla_A>Wavelet db2 L=2 · umbral duro>retención del pico|.2f@@ y la correlación de @@ej2_tabla_A>Wavelet db2 L=2 · umbral suave (enunciado)>correlación de forma|.2f@@ a @@ej2_tabla_A>Wavelet db2 L=2 · umbral duro>correlación de forma|.2f@@ (a costa de más ruido en el cielo: @@ej2_tabla_A>Wavelet db2 L=2 · umbral duro>RMS en cielo vacío|.4f@@).@@d>duro_txt@@ Ambos métodos eliminan además la baja frecuencia de la propia fuente ($cA_2=0$ en la wavelet, $H(0)=0$ en la FFT). (ii) La DWT **no es invariante a traslaciones** (Experimento 1b): la misma estrella queda con entre @@ej2_traslacion>Wavelet db2 L=2>pico/entrada · mín.|.0%@@ y @@ej2_traslacion>Wavelet db2 L=2>pico/entrada · máx.|.0%@@ de su pico según caiga respecto a la rejilla de 4 px, mientras que el filtro lineal da siempre @@ej2_traslacion>FFT gaussiano>pico/entrada · mín.|.0%@@.

@@d>tamano_txt@@

**Veredicto.** @@d>ver_txt@@ *Limitaciones:* @@d>lim_img@@, fuentes de referencia seleccionadas sobre la propia imagen original y una sola imagen.
''', R)


def disc2B(R):
    R = dict(R)
    B = R["ej2_tabla_B"]
    R["d"] = dict(sigma_r=1 / (2 * 3.141592653589793 * R["ej2_D0"]), per=1 / R["ej2_D0"])
    return fill(r'''---
## Pregunta B — ¿por qué la 2D-FFT produce *ringing* (Gibbs) y las wavelets no? Explicación matemática

**1 · El filtrado en frecuencia es una convolución.** Por el teorema de convolución,
$$f_{fft}=\mathcal F^{-1}\{F\cdot H\}=h*f,\qquad h=\mathcal F^{-1}\{H\}.$$
Cada fuente puntual se sustituye por una copia de la **respuesta al impulso** $h$. Si $H$ tiene una discontinuidad, $h$ no tiene soporte compacto y oscila.

**2 · Filtro ideal (corte abrupto): Gibbs.** Para $H=1-\mathbf 1_{\{D\le D_0\}}$ se tiene $h=\delta-h_{LP}$ con el núcleo "jinc"
$$h_{LP}(r)=D_0\,\frac{J_1(2\pi D_0 r)}{r}\;\sim\; r^{-3/2}\cos\!\left(2\pi D_0 r-\tfrac{3\pi}{4}\right),\quad r\to\infty,$$
que oscila con periodo $1/D_0\approx$ @@d>per|.1f@@ px y decae sólo como $r^{-3/2}$: el artefacto se extiende a **toda la imagen**. La causa es que una discontinuidad en el dominio de Fourier equivale a un decaimiento lento en el espacio. Es el fenómeno de Gibbs (por dualidad, el mismo que aparece al truncar la serie de Fourier de una discontinuidad en el espacio): la serie de Fourier truncada de una discontinuidad de salto 1 rebasa el valor en $\tfrac1\pi\mathrm{Si}(\pi)-\tfrac12\approx$ @@ej2_gibbs>limite_pct|.2f@@ % **sin importar cuántos términos se usen** (Experimento 2: @@ej2_gibbs>tabla>11|.3f@@ % con 11 armónicos y @@ej2_gibbs>tabla>1001|.3f@@ % con 1001). Una fuente brillante es un salto/impulso cuyo espectro es de banda ancha; al truncar o filtrar bruscamente ese espectro se pierde la coherencia de fases que lo mantenía localizado.

*Medido (Experimento 1):* filtro ideal → @@ej2_tabla_B>FFT ideal (corte abrupto)>cambios de signo|.0f@@ cambios de signo y un radio que contiene el 99 % de la energía de @@ej2_tabla_B>FFT ideal (corte abrupto)>radio 99 % energía [px]|.0f@@ px.

**3 · Matiz importante: el filtro gaussiano que pide el enunciado casi no produce Gibbs.** Para $H=1-e^{-D^2/2D_0^2}$, que es $C^\infty$ (sin discontinuidad), $h=\delta-2\pi D_0^2\,e^{-2\pi^2D_0^2r^2}$: un núcleo gaussiano de $\sigma_r=1/(2\pi D_0)\approx$ @@d>sigma_r|.1f@@ px, sin oscilaciones largas (@@ej2_tabla_B>FFT gaussiano>cambios de signo|.0f@@ cambio de signo y radio de 99 % de energía de @@ej2_tabla_B>FFT gaussiano>radio 99 % energía [px]|.0f@@ px). El "halo" negativo que sí se observa (≈ @@ej2_tabla_B>FFT gaussiano>lóbulo negativo mín. / pico|.2f@@ del pico) **no es Gibbs**: viene de que $H(0)=0$, es decir $\int h=0$, así que $h$ **debe** tener una parte negativa (el filtro quita la media local alrededor de cada fuente). La pregunta B se cumple, por tanto, para el pasa-altas ideal, y sólo en parte para el gaussiano.

**4 · Wavelets: artefactos acotados por el soporte.** La reconstrucción es $f=\sum_{j,k}c_{j,k}\,\psi_{j,k}$ con $\psi$ de soporte compacto (db2 tiene 4 coeficientes de filtro; a nivel $j$ su soporte es de $3(2^j-1)+1$ muestras, es decir, 10 px a $j=2$). Modificar o anular un coeficiente sólo puede alterar la imagen **dentro del soporte de su átomo**, así que el error no puede propagarse más allá de unos pocos píxeles; y la umbralización suave actúa coeficiente a coeficiente, sin un corte abrupto en frecuencia. *Medido:* radio del 99 % de la energía de la respuesta de la wavelet = @@ej2_tabla_B>Wavelet db2 L=2 (umbral suave)>radio 99 % energía [px]|.0f@@ px, frente a @@ej2_tabla_B>FFT ideal (corte abrupto)>radio 99 % energía [px]|.0f@@ px del filtro ideal.

**5 · Pero las wavelets *no* están libres de artefactos.** Dentro de su soporte, la respuesta de la wavelet tiene un lóbulo negativo de entre @@ej2_traslacion>Wavelet db2 L=2>lóbulo neg./pico · mejor|.2f@@ y @@ej2_traslacion>Wavelet db2 L=2>lóbulo neg./pico · peor|.2f@@ veces su pico (según la posición de la fuente), **mayor en términos relativos que el del filtro gaussiano** (@@ej2_traslacion>FFT gaussiano>lóbulo neg./pico · mejor|.2f@@). Causas: $cA_2=0$ también impone media cero; la DWT no es invariante a traslaciones; y db2 tiene sólo 2 momentos nulos y no es simétrica. Esos defectos se manifiestan como "bloques" (ver la compresión al 5 %). La conclusión correcta es: **las wavelets acotan *dónde* aparecen los artefactos (localidad), no su amplitud**.
''', R)


def disc2C(R):
    R = dict(R)
    T, T10 = R["ej2_mse_tab"]["5 %"], R["ej2_mse_tab"]["10 %"]
    comp = R["ej2_composicion_5pct"]
    F = T["FFT"]
    L2 = ["db2 L=2 (enunciado)", "db4 L=2", "bior2.2 L=2"]
    LH = ["db2 L=3", "db2 L=5", "db4 L=5", "bior2.2 L=5"]
    q = {k: T[k] / F for k in L2 + LH}
    nom = lambda k: "`" + k.replace(" (enunciado)", "") + "`"
    s2, sw2 = R["ej2_sig_fondo"] ** 2, R["ej2_sigma_w"] ** 2
    e5 = T["db2 L=5"]
    mejor = [k for k in LH if q[k] < 0.95]; peor = [k for k in LH if q[k] > 1.05]
    if len(mejor) == len(LH):
        veredicto_LH = "las wavelets ganan con claridad"
    elif not peor and mejor:
        veredicto_LH = "las wavelets empatan o ganan por poco (la FFT sigue siendo competitiva)"
    elif not peor and not mejor:
        veredicto_LH = "las wavelets empatan con la FFT"
    else:
        veredicto_LH = "el resultado es mixto"
    win10 = [k for k in L2 if T10[k] < T10["FFT"]]
    if len(win10) == len(L2):
        t10 = f"Con el 10 % de coeficientes todas las wavelets ganan, incluso a $L=2$ ({T10['db2 L=2 (enunciado)']:.2e} frente a {T10['FFT']:.2e})."
    elif not win10:
        t10 = (f"Con el 10 % de coeficientes, a $L=2$ **ninguna** de las tres gana todavía a la FFT ({T10['FFT']:.2e}): db2 {T10['db2 L=2 (enunciado)']:.2e}, db4 {T10['db4 L=2']:.2e}, bior2.2 {T10['bior2.2 L=2']:.2e}; "
               f"la brecha se reduce de ×{q['db2 L=2 (enunciado)']:.0f} a ×{T10['db2 L=2 (enunciado)'] / T10['FFT']:.2f} para db2 porque ya cabe toda la banda $cA_2$ y sobra presupuesto para detalles.")
    else:
        t10 = f"Con el 10 % de coeficientes, a $L=2$ ganan a la FFT ({T10['FFT']:.2e}) sólo " + ", ".join(nom(k) for k in win10) + "."
    # ¿cuánto del error de db2 L=5 es ruido a escala de píxel?
    if e5 > 0.2 * s2:
        ruido_txt = (f"El error de db2 a $L=5$ ({e5:.2e}) es del orden de la varianza del fondo ($\\hat\\sigma^2\\approx$ {s2:.1e}): el ruido de fondo no se concentra en pocos coeficientes en ninguna de estas bases, "
                     "lo que sugiere (sin demostrarlo aquí) que el error restante es en su mayor parte ruido.")
    else:
        ruido_txt = (f"El error de db2 a $L=5$ ({e5:.2e}) es sólo el {100 * e5 / s2:.2f} % de la varianza robusta de la imagen ($\\hat\\sigma^2\\approx$ {s2:.1e}) y {e5 / sw2:.1f} veces la varianza del ruido a escala de píxel "
                     f"($\\hat\\sigma_w^2\\approx$ {sw2:.1e}): la imagen es muy compresible y el 5 % de los coeficientes recoge casi toda la señal; que el error siga por encima del ruido a escala de píxel indica que queda algo de estructura sin representar (no demostrado aquí).")
    fw = fwhm_fuentes(R)
    if len(mejor) < len(LH) and fw and fw > 4.5:
        esparcidad_txt = (f"\n\n**Matiz para esta imagen.** La ventaja esperada de la localización es aquí pequeña o nula: las fuentes miden ≈ {fw:.1f} px de FWHM y el ruido a escala de píxel es muy pequeño (σ̂_w = {R['ej2_sigma_w']:.4f} frente a una dispersión global de {R['ej2_sig_fondo']:.4f}), de modo que la energía está en bajas frecuencias, "
                          "donde Fourier es muy eficiente. Es una explicación plausible, no verificada aquí; con fuentes más finas (`--reducir`) la comparación podría cambiar.")
    else:
        esparcidad_txt = ""
    lo, hi = min(q[k] for k in LH), max(q[k] for k in LH)
    mejor_k = nom(min(LH, key=lambda k: q[k]))
    if len(mejor) == len(LH):
        concl = (f"**una DWT con suficiente profundidad ($L\\ge3$) da un MSE ≈ {1 / hi:.2f}–{1 / lo:.2f} veces menor que la FFT al 5 %**, "
                 "por la localización de las fuentes y la concentración de energía en pocos coeficientes")
    elif mejor:
        concl = (f"**una DWT con $L\\ge3$ queda, al 5 %, entre {lo:.2f} y {hi:.2f} veces el MSE de la FFT** (la mejor, {mejor_k}): la ventaja de las wavelets sobre la FFT es pequeña "
                 "y depende de la familia y de $L$")
    else:
        concl = f"**con $L\\ge3$ la DWT no mejora a la FFT al 5 %** (cociente de MSE entre {lo:.2f} y {hi:.2f}): aquí la FFT es igual o mejor"
    R["d"] = dict(
        frase_L2="; ".join(f"{nom(k)}: {vs_fft(q[k])}" for k in L2), frase_LH="; ".join(f"{nom(k)}: {vs_fft(q[k])}" for k in LH), veredicto_LH=veredicto_LH,
        b2=vs_fft(q["db2 L=2 (enunciado)"]), bb=vs_fft(q["bior2.2 L=2"]), b3=vs_fft(q["db2 L=3"]), b5=vs_fft(q["db2 L=5"]), b5b=vs_fft(q["bior2.2 L=5"]),
        t10=t10, ruido_txt=ruido_txt, esparcidad_txt=esparcidad_txt,
        concl=concl,
        det=100 - comp["db2 L=2"]["retenidos_en_aprox_pct"],
        M=int(0.05 * R["ej2_shape"][0] * R["ej2_shape"][1]),
        cA=(R["ej2_shape"][0] // 4) * (R["ej2_shape"][1] // 4),
    )
    return fill(r'''---
## Pregunta C — compresión reteniendo sólo el 5 % de los coeficientes: ¿qué método garantiza menor MSE?

**Respuesta corta.** **Ninguno lo garantiza para cualquier configuración**; depende de la profundidad $L$ y de la familia. A $L=2$ (el presupuesto del 5 % es menor que la banda $cA_2$): @@d>frase_L2@@. Con $L\ge3$: @@d>frase_LH@@; en conjunto, @@d>veredicto_LH@@. MSE al 5 % (menor es mejor):

| FFT | db2 L=2 *(enunciado)* | db4 L=2 | bior2.2 L=2 | db2 L=3 | db2 L=5 | db4 L=5 | bior2.2 L=5 |
|---|---|---|---|---|---|---|---|
| @@ej2_mse_tab>5 %>FFT|.2e@@ | @@ej2_mse_tab>5 %>db2 L=2 (enunciado)|.2e@@ | @@ej2_mse_tab>5 %>db4 L=2|.2e@@ | @@ej2_mse_tab>5 %>bior2.2 L=2|.2e@@ | @@ej2_mse_tab>5 %>db2 L=3|.2e@@ | @@ej2_mse_tab>5 %>db2 L=5|.2e@@ | @@ej2_mse_tab>5 %>db4 L=5|.2e@@ | @@ej2_mse_tab>5 %>bior2.2 L=5|.2e@@ |

* Con db2 a $L=2$ (la configuración usada en la Parte 4): @@d>b2@@; con `bior2.2` a $L=2$: @@d>bb@@.
* Con $L=3$ (db2): @@d>b3@@; con $L=5$ (db2): @@d>b5@@; con $L=5$ y `bior2.2`: @@d>b5b@@. @@d>t10@@

**Por qué (el criterio).** Ambas transformadas son (casi) ortonormales, así que por Parseval $\mathrm{MSE}=\tfrac1n\sum(\text{coeficientes descartados})^2$, y la mejor aproximación con $M$ términos consiste en conservar los $M$ coeficientes de mayor módulo. Gana, por tanto, la base en la que la imagen es **más dispersa** (energía concentrada en menos coeficientes). Una imagen de campo profundo es un conjunto de **estructuras localizadas** (fuentes compactas) sobre un fondo suave más ruido: una base global como la de Fourier necesita muchas frecuencias coherentes para localizar cada fuente, mientras que una wavelet de soporte compacto la describe con pocos coeficientes por escala.@@d>esparcidad_txt@@

**Por qué a $L=2$ y 5 % no se aprecia la ventaja.** La banda de aproximación $cA_2$ tiene @@d>cA|.0f@@ coeficientes (el @@ej2_composicion_5pct>db2 L=2>banda_aprox_pct_de_N|.2f@@ % de $N$), más que todo el presupuesto de @@d>M|.0f@@ coeficientes (5 %). Como los coeficientes de aproximación son los de mayor módulo, el @@ej2_composicion_5pct>db2 L=2>retenidos_en_aprox_pct|.0f@@ % de lo retenido es de $cA_2$ y queda sólo ≈ @@d>det|.0f@@ % para los detalles que describen las fuentes; y ni siquiera la propia banda $cA_2$ se conserva entera (sólo el @@ej2_composicion_5pct>db2 L=2>aprox_retenida_pct|.0f@@ %). Resultado: reconstrucción con "bloques" oscuros y @@ej2_tabla_C>db2 L=2 (enunciado)>píxeles < 0 [%]|.1f@@ % de píxeles negativos (imposibles físicamente), frente a @@ej2_tabla_C>FFT>píxeles < 0 [%]|.1f@@ % con la FFT. Con $L=5$ la aproximación ocupa sólo el @@ej2_composicion_5pct>db2 L=5>banda_aprox_pct_de_N|.2f@@ % de $N$ y casi todo el presupuesto se dedica a detalles. @@d>ruido_txt@@

**Conclusión de C.** Para esta imagen, @@d>concl@@. No hay garantía general: depende de $L$, de la familia y del presupuesto. *Limitaciones:* una imagen, MSE medido sobre la imagen normalizada a [0, 1], y la selección por módulo no es la óptima para bases biortogonales (`bior2.2`).
''', R)
