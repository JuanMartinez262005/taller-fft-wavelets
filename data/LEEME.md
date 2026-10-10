# Datos

| Archivo | Qué es | ¿Está en el repositorio? |
|---|---|---|
| `ecg_mitbih208.npy` | ECG del registro 208 de MIT-BIH (360 Hz, 5 min), copia de SciPy | Sí |
| `hubble_deep_field_rgb.npy` | Hubble eXtreme Deep Field (scikit-image). Solo sirve de imagen de respaldo para comprobar que el código corre | Sí |
| `ej1_comprimidos.npz` | Coeficientes retenidos del Ejercicio 1 (salida del cuaderno) | Sí |
| `epico.tif` | Imagen JWST NIRCam de IC 348 usada en el Ejercicio 2 (38 204 × 30 866 px) | **No**: GitHub no admite archivos de más de 100 MB |

Para reproducir el Ejercicio 2 con la imagen del taller, descárguela de nuevo, guárdela aquí como `data/epico.tif` y ejecute, desde la carpeta principal:

```
python construir.py --imagen data/epico.tif --reducir 3
```

Crédito de la imagen: https://science.nasa.gov/asset/webb/ic-348-nircam-image-2/