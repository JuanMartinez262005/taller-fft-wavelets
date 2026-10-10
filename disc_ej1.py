# -*- coding: utf-8 -*-
from disc_util import fill, lim_img

N = 8192


def disc1(R):
    R = dict(R)
    ro = R["ej1_robustez"]
    R["d"] = dict(
        pF=100 * R["ej1_M_F"] / N, pW=100 * R["ej1_M_W"] / N,
        menos=100 * (1 - ro["M_W_mean"] / ro["M_F_mean"]),
        pierdeF=100 * (1 - ro["retR_F"]), pierdeW=100 * (1 - ro["retR_W"]),
        rQRS=ro["rmsQRS_F"] / ro["rmsQRS_W"],
        e21=R["ej1_niveles"]["D2"]["energia"] + R["ej1_niveles"]["D1"]["energia"],
    )
    return fill(r'''---
## Discusión · Ejercicio 1 — ¿qué método preserva mejor los eventos transitorios?

**Respuesta corta.** La **DWT (db4, L = 5)**. Con la misma energía retenida (95 %) y, por tanto, el mismo MSE, necesita bastante menos coeficientes y conserva mejor la amplitud de los complejos QRS, que son el evento transitorio de esta señal. Tiene matices que se detallan al final: la ventaja es moderada y no es uniforme.

**1 · Resultados.** Segmento de trabajo ($N=8192$):

| | FFT | DWT (db4, L = 5) |
|---|---|---|
| Coeficientes $M$ (% de $N$) | @@ej1_M_F@@ (@@d>pF|.2f@@ %) | @@ej1_M_W@@ (@@d>pW|.2f@@ %) |
| Razón de compresión $CR=N/M$ | @@ej1_CR_F|.1f@@ | @@ej1_CR_W|.1f@@ |
| MSE [mV²] | @@ej1_MSE_F|.5f@@ | @@ej1_MSE_W|.5f@@ |
| SNR [dB] | @@ej1_SNR_F|.1f@@ | @@ej1_SNR_W|.1f@@ |

Los MSE coinciden (diferencia relativa máxima entre métodos en los 13 segmentos: @@ej1_robustez>mse_gap_max|.1%@@) porque **ambas bases son ortonormales**: por Parseval, $\mathrm{MSE}=\tfrac1N\times$(energía descartada), y los dos métodos descartan ≈ 5 %. El MSE, por sí solo, **no distingue** los métodos en este experimento; lo que cambia es *cuántos* coeficientes hacen falta y *dónde* cae el error. En los 13 segmentos no solapados del registro, la DWT necesita @@ej1_robustez>M_W_mean|.0f@@ ± @@ej1_robustez>M_W_std|.0f@@ coeficientes frente a @@ej1_robustez>M_F_mean|.0f@@ ± @@ej1_robustez>M_F_std|.0f@@ de la FFT (@@d>menos|.0f@@ % menos; menor en @@ej1_robustez>wins_M|.0f@@/@@ej1_robustez>n_seg|.0f@@ segmentos), es decir, $CR$ medio de @@ej1_robustez>CR_W_mean|.1f@@ frente a @@ej1_robustez>CR_F_mean|.1f@@. Si se da a la FFT el mismo número de coeficientes que usa la DWT, su MSE es @@ej1_robustez>mseratio_igualM|.2f@@ veces mayor (la DWT gana en @@ej1_robustez>wins_igualM|.0f@@/@@ej1_robustez>n_seg|.0f@@).

**2 · Por qué: el principio de Heisenberg–Gabor.** Para cualquier función, $\Delta t\,\Delta\omega\ge \tfrac12$ (con $\Delta t,\Delta\omega$ las desviaciones típicas de $|f|^2$ y $|\hat f|^2$): no existe un átomo con localización perfecta en tiempo y en frecuencia a la vez, y cada base elige un punto de compromiso.

* **Fourier** elige $\Delta\omega\to0$ (frecuencia exacta) a costa de $\Delta t\to T$: sus átomos son sinusoides de soporte global. Un QRS dura ≈ 80 ms y es de banda ancha; para que una suma de sinusoides eternas sea cero fuera del QRS y grande dentro, hace falta que cientos de componentes **cancelen sus fases con precisión** en el resto del registro. Es una representación ineficiente (hay que conservar muchas líneas) y frágil: al descartar el 5 % de la energía la cancelación de fases deja de ser exacta y aparecen ondulaciones (*ringing*) alrededor del evento, visibles en el zoom de la figura de reconstrucción.
* **La DWT** usa átomos de soporte compacto cuya anchura crece con la escala: $\Delta t\propto 2^{j}$ y $\Delta\omega\propto 2^{-j}$ (el producto es constante; ver la figura de teselaciones). En las bandas de alta frecuencia (D3: 22.5–45 Hz, D4: 11.2–22.5 Hz) hay buena resolución temporal, justo donde vive la energía del QRS; en la banda lenta (A5: 0–5.6 Hz) hay buena resolución en frecuencia para la deriva de línea base y las ondas P/T.
* **Evidencia en los datos.** Los coeficientes de detalle conservados se alinean con los latidos: las ventanas QRS (±60 ms) ocupan el @@ej1_alineacion>frac_tiempo_QRS_pct|.0f@@ % del tiempo y contienen el @@ej1_alineacion>D3>en_QRS_pct|.0f@@ % de los coeficientes conservados de D3 y el @@ej1_alineacion>D4>en_QRS_pct|.0f@@ % de los de D4; D2 y D1 (45–180 Hz; juntas, @@d>e21|.2f@@ % de la energía) se descartan por completo. Los coeficientes de la aproximación A5 (@@ej1_niveles>A5>energia|.0f@@ % de la energía) se reparten por todo el registro (@@ej1_alineacion>A5>en_QRS_pct|.0f@@ % dentro de QRS). En consecuencia, la amplitud del pico R se conserva al @@ej1_robustez>retR_W|.1%@@ con la DWT frente a @@ej1_robustez>retR_F|.1%@@ con la FFT (mejor en @@ej1_robustez>wins_retR|.0f@@/@@ej1_robustez>n_seg|.0f@@ segmentos), y el RMS del error dentro de las ventanas QRS es @@ej1_robustez>rmsQRS_W|.3f@@ mV frente a @@ej1_robustez>rmsQRS_F|.3f@@ mV (menor en @@ej1_robustez>wins_QRS|.0f@@/@@ej1_robustez>n_seg|.0f@@).

**3 · Matices (lo que este experimento no permite afirmar).**

1. **Ninguno de los dos preserva bien el QRS con el criterio del 95 %**: ambos atenúan el pico R en un @@d>pierdeW|.0f@@–@@d>pierdeF|.0f@@ %. El criterio de energía está dominado por la línea base (A5 contiene el @@ej1_niveles>A5>energia|.0f@@ % de la energía), y el QRS, con poca energía total, es lo primero que se sacrifica. Para fidelidad diagnóstica convendría umbralizar por banda o subir el umbral de energía.
2. **En el segmento de trabajo la diferencia en el error del QRS es mínima** (@@ej1_seg0>rmsQRS_W|.3f@@ mV con DWT frente a @@ej1_seg0>rmsQRS_F|.3f@@ mV con FFT); su ventaja allí es sobre todo en $M$ y en la amplitud del pico R (@@ej1_seg0>retR_W|.2f@@ frente a @@ej1_seg0>retR_F|.2f@@). La ventaja en el error QRS aparece al promediar los 13 segmentos (cociente FFT/DWT de @@d>rQRS|.2f@@).
3. **Fuera de los QRS gana la FFT**: RMS @@ej1_robustez>rmsResto_F|.3f@@ mV frente a @@ej1_robustez>rmsResto_W|.3f@@ mV (la DWT sólo es mejor en @@ej1_robustez>wins_resto|.0f@@/@@ej1_robustez>n_seg|.0f@@ segmentos). La DWT concentra su presupuesto en los transitorios y sacrifica detalle de baja amplitud (ondas P/T, ruido); es un buen trato sólo si lo que importa son los eventos.
4. **A baja retención la FFT necesita menos coeficientes** (con el 80 % de la energía: @@ej1_curva>80.0>M_F@@ frente a @@ej1_curva>80.0>M_W@@). Una explicación plausible (no verificada aquí) es que con tan poco detalle basta describir la deriva lenta, casi sinusoidal, para la que Fourier es eficiente. La ventaja de la DWT aparece cuando hay que reconstruir transitorios (≥ 90 % de energía: @@ej1_curva>90.0>M_F@@ frente a @@ej1_curva>90.0>M_W@@ coeficientes; ver curva de compromiso).
5. **Límites del estudio:** un solo registro (MIT-BIH 208) y una sola familia (db4); $CR=N/M$ cuenta coeficientes, no bits (no incluye el coste de guardar las posiciones de los coeficientes retenidos —que podría penalizar más a la DWT, cuyos coeficientes conservados están más dispersos (figura (b))— ni codificación entrópica).

**Conclusión del Ejercicio 1.** Para señales con transitorios sobre un fondo lento (ECG), la DWT ofrece mejor compromiso tiempo–frecuencia: mismo MSE con ≈ @@d>menos|.0f@@ % menos coeficientes y mejor conservación del QRS. La FFT sigue siendo preferible para componentes estacionarias/lentas y cuando se exige poca energía retenida.
''', R)


def conclusiones(R):
    R = dict(R)
    R["d"] = dict(lim_img=lim_img(R))
    return fill(r'''---
## Conclusiones generales

* **Elegir la base es elegir qué es "pocos coeficientes".** Con bases ortonormales, el MSE de una compresión es la energía descartada (Parseval); gana la base en la que la señal es más *sparse*. Fourier lo es para componentes estacionarias; las wavelets, para estructuras localizadas (QRS, fuentes puntuales).
* **El MSE no basta.** En el Ejercicio 1 ambos métodos tienen el mismo MSE al 95 % de energía y, aun así, difieren en el número de coeficientes (@@ej1_robustez>M_F_mean|.0f@@ frente a @@ej1_robustez>M_W_mean|.0f@@) y en cómo tratan los transitorios.
* **"Las wavelets son mejores" depende de la tarea y de los parámetros.** En el Ejercicio 2 la wavelet db2 (L = 2, $cA_2=0$) deja el cielo mucho más limpio (RMS @@ej2_tabla_A>Wavelet db2 L=2 · umbral suave (enunciado)>RMS en cielo vacío|.4f@@ frente a @@ej2_tabla_A>FFT pasa-altas gaussiano (enunciado)>RMS en cielo vacío|.4f@@ con la FFT), pero deforma más las fuentes (correlación de forma @@ej2_tabla_A>Wavelet db2 L=2 · umbral suave (enunciado)>correlación de forma|.2f@@ frente a @@ej2_tabla_A>FFT pasa-altas gaussiano (enunciado)>correlación de forma|.2f@@); y en la compresión al 5 % sólo supera a la FFT con L ≥ 3 o con `bior2.2`.
* **Limitaciones comunes:** una sola señal y una sola imagen por ejercicio; @@d>lim_img@@; fuentes de referencia seleccionadas sobre la propia imagen original.
''', R)
