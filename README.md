# TeuB · De una secuencia sin etiqueta a un ranking de ligandos

Práctica del curso **Evolución y Diseño de Proteínas 2027-1**. A partir de la secuencia de TeuB (362 aa, *Rhizobium*), se identifica su familia, se predice su estructura, se recorta el péptido señal, se arma un panel de 16 ligandos con Foldseek y se ordenan por dos métodos independientes, **AutoDock Vina** (docking) y **Boltz-2** (co-plegamiento), midiendo el piso de ruido de cada uno.

**Conclusión:** los dos métodos colocan todos los ligandos en el mismo bolsillo, pero no coinciden en el orden (Spearman ρ = 0.39, p = 0.13). Boltz-2 empata a 12 de los 16 ligandos, y ninguno de los dos métodos rechaza a los señuelos. No se puede decir qué azúcar prefiere TeuB. La predicción previa (D-apiosa) no se confirma ni se descarta.

## Documentos

| Archivo | Qué es |
|---|---|
| **`ENTREGABLE_v2.pdf`** | **Entregable final** (3 páginas): ranking, pisos de ruido, comparación, cadena de inferencias, panel, péptido señal y veredicto, con una figura |
| `ENTREGABLE.pdf` / `ENTREGABLE.md` | Primera versión del entregable (mismo contenido, sin figura) |
| **`GUIA_PASO_A_PASO_v2.pdf`** / `.md` | Guía autosuficiente para reproducir la práctica: conceptos con diagramas, pasos web clic por clic y ruta alternativa por terminal |
| `GUIA_PASO_A_PASO.pdf` / `.md` | Primera versión de la guía |
| `CUADERNO.md` | Notas y resultados acumulados, etapa por etapa |
| `ESTADO.md` | Estado de cada etapa y decisiones tomadas |

## Entregable v2

Se genera a partir de dos archivos:

- `08_analisis/ENTREGABLE_v2.html`: la fuente, con el diseño de impresión en tamaño carta. Resultado: 3 páginas.
- `08_analisis/figura_entregable.py`: produce `figura_entregable.svg` y `figura_entregable.png`. La figura tiene tres paneles:
  - **A:** Vina ± DE entre semillas;
  - **B:** ipTM de Boltz-2 ± DE entre los 15 modelos, con una banda que marca los empates;
  - **C:** el puesto en Vina contra el puesto en Boltz-2.

Qué cambia respecto a la v1:

- **Diseño:** encabezado, secciones numeradas, señuelos en rojo, la predicción en ámbar y los empates sombreados.
- **Contenido nuevo:** un resumen al inicio, la figura y las cajas de ruido y de péptido señal.
- **Corrección:** el valor p del Spearman con ipTM es **0.13**; en la v1 decía 0.14 por redondeo.

Para regenerar el PDF y la figura después de editar:

```bash
~/miniforge3/bin/python 08_analisis/figura_entregable.py      # solo si cambian los datos
google-chrome --headless --no-pdf-header-footer --allow-file-access-from-files \
  --print-to-pdf=ENTREGABLE_v2.pdf 08_analisis/ENTREGABLE_v2.html
```

## Estructura del repositorio

```
00_datos_repartidos/   respaldo de la página de la práctica (el modelo tiene una Arg de menos en la posición 256)
01_blast/              teuB.fasta, swissprot_hits.tsv, nr_hits.tsv
02_estructura/         modelo AlphaFold2 (ColabFold, sin plantillas)
03_maduro/             modelo y secuencia madura 1-335 (sin péptido señal)
04_foldseek/           aciertos de Foldseek y ligandos por entrada del PDB
05_ligandos/           16 ligandos: SDF ideal, mol2 y pdbqt
06_docking/            receptor, run.sh y salidas de Vina (16 ligandos × 3 semillas)
07_coplegamiento/      ZIP de Boltz-2 válidos (+ _descartado_/ con las corridas malas)
08_analisis/           scripts de análisis, tablas, figura y fuente del entregable v2
guia/                  scripts de la guía paso a paso (p02_…–p07_…)
```

Las bases de datos locales (Swiss-Prot, 844 MB, y el PDB para Foldseek, 6.5 GB) no están en el repositorio. Cómo descargarlas se explica en la guía.

## Scripts de análisis

| Script | Qué hace |
|---|---|
| `08_analisis/revisar_boltz.py` | valida los ZIP de Boltz-2 (secuencia, plantillas, semillas, código) y mide la posición del ligando |
| `08_analisis/poses_vina.py` | comprueba que la mejor pose de Vina esté en el bolsillo |
| `08_analisis/comparar.py` | tabla comparativa, empates, Spearman y RMSD entre las poses de los dos métodos |
| `08_analisis/figura_entregable.py` | figura del entregable v2 |

Requisitos: Python con numpy, pandas, scipy, matplotlib y biopython (entorno de miniforge), además de `vina`, `obabel`, `foldseek` y `blastp` para rehacer los cálculos.
