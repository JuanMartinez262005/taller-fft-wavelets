# -*- coding: utf-8 -*-
"""Comprueba si las afirmaciones cualitativas de las discusiones siguen siendo ciertas con los datos usados."""


def afirmaciones_fallidas(R):
    fallos = []
    ro = R["ej1_robustez"]
    A = R["ej2_tabla_A"]
    W = "Wavelet db2 L=2 · umbral suave (enunciado)"
    G = "FFT pasa-altas gaussiano (enunciado)"
    I = "FFT pasa-altas ideal (corte abrupto)"
    T5 = R["ej2_mse_tab"]["5 %"]
    B = R["ej2_tabla_B"]
    tr = R["ej2_traslacion"]
    cur = R["ej1_curva"]

    def chk(cond, texto):
        if not cond:
            fallos.append(texto)

    # Selección de fuentes (Ejercicio 2)
    c = R.get("ej2_criterio_fuentes")
    if c:
        chk(not c.get("relajado", c["modo"] != "estricto"),
            f"Ej. 2: «las fuentes de referencia son compactas y aisladas (área ≤ 12 px)». Con esta imagen el criterio estricto sólo dio {c['n_estricto']} fuentes y se usó el modo "
            f"«{c['modo']}» (área ≤ {c['area_max']} px, área mediana {c['area_mediana']:.0f} px). La Pregunta A (posición, forma, halo) está medida sobre fuentes más extensas que "
            "un punto: revisen la figura de selección y maticen en el texto que no son fuentes puntuales en sentido estricto.")
    # Ejercicio 1
    chk(ro["M_W_mean"] < ro["M_F_mean"], "Ej. 1: «la DWT necesita menos coeficientes que la FFT al 95 % de energía».")
    chk(ro["retR_W"] > ro["retR_F"], "Ej. 1: «la DWT conserva mejor la amplitud del pico R».")
    chk(ro["rmsQRS_W"] < ro["rmsQRS_F"], "Ej. 1: «el error en las ventanas QRS es menor con la DWT».")
    chk(ro["rmsResto_F"] < ro["rmsResto_W"], "Ej. 1: «fuera de los QRS gana la FFT».")
    chk(cur["80.0"]["M_F"] < cur["80.0"]["M_W"], "Ej. 1: «con el 80 % de energía la FFT necesita menos coeficientes».")
    chk(ro["mseratio_igualM"] > 1, "Ej. 1: «a igual M la FFT tiene mayor MSE que la DWT».")
    # Pregunta A
    chk(A[G]["correlación de forma"] > A[W]["correlación de forma"], "Pregunta A: «la FFT gaussiana conserva mejor la forma».")
    chk(A[G]["posición ≤1 px"] >= A[W]["posición ≤1 px"], "Pregunta A: «la FFT localiza algo mejor la posición».")
    chk(A[W]["RMS en cielo vacío"] < A[G]["RMS en cielo vacío"], "Pregunta A: «la wavelet deja el cielo más limpio».")
    chk(A["Wavelet db2 L=2 · umbral duro"]["retención del pico"] > A[W]["retención del pico"], "Pregunta A: «el umbral duro retiene más pico que el suave».")
    chk(tr["Wavelet db2 L=2"]["pico/entrada · máx."] - tr["Wavelet db2 L=2"]["pico/entrada · mín."] > 0.05, "Pregunta A: «la DWT no es invariante a traslaciones».")
    # Pregunta B
    chk(B["FFT ideal (corte abrupto)"]["radio 99 % energía [px]"] > B["Wavelet db2 L=2 (umbral suave)"]["radio 99 % energía [px]"] > 0, "Pregunta B: «el artefacto del filtro ideal se extiende más que el de la wavelet».")
    chk(B["FFT ideal (corte abrupto)"]["cambios de signo"] > B["FFT gaussiano"]["cambios de signo"], "Pregunta B: «el filtro ideal oscila más que el gaussiano».")
    chk(tr["Wavelet db2 L=2"]["lóbulo neg./pico · peor"] < tr["FFT gaussiano"]["lóbulo neg./pico · peor"], "Pregunta B: «el lóbulo negativo de la wavelet es mayor que el del filtro gaussiano».")
    # Pregunta C (los cocientes de MSE entre configuraciones se redactan a partir de los números; aquí sólo lo que el texto da por fijo)
    chk(T5["db2 L=2 (enunciado)"] >= T5["FFT"], "Pregunta C: «a L = 2 con db2, al 5 % no se aprecia la ventaja de la DWT» (aquí la DWT db2 L = 2 supera a la FFT).")
    chk(R["ej2_composicion_5pct"]["db2 L=2"]["retenidos_en_aprox_pct"] > 50, "Pregunta C: «a L = 2 la mayor parte del presupuesto se gasta en cA₂».")
    return fallos
