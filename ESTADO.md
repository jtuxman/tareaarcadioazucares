# ESTADO DE LA PRÁCTICA — retomar desde aquí

Última actualización: 2026-09-27 noche. Predicción fijada: **D-apiosa (PAV)** (ver CUADERNO.md). Ojo: corridas Boltz con secuencia de 132 aa o código XIP descartadas en `07_coplegamiento/_descartado_*`.
Directorio raíz: `/home/jmanuel/cursos/proteinas/tareasugar`
Práctica: https://practica-ligandos-2027-1.netlify.app/ (copia local en el HTML de la raíz)

## CÓMO RETOMAR EN 30 SEGUNDOS

```bash
cd /home/jmanuel/cursos/proteinas/tareasugar
tmux attach -t tarea          # ver el docking corriendo  (Ctrl+B, D para salir sin matarlo)
ls 06_docking/out_*.pdbqt | wc -l   # progreso: debe llegar a 45
cat ESTADO.md CUADERNO.md
```

Si la sesión de tmux ya no existe pero faltan corridas:
```bash
cd /home/jmanuel/cursos/proteinas/tareasugar/06_docking
find . -name "out_*.pdbqt" -size -1k -delete      # limpia salidas truncadas
setsid tmux new-session -d -s tarea -c "$PWD" './run.sh; exec bash' < /dev/null
```
`run.sh` es **reanudable**: salta todo ligando/semilla que ya tenga su `out_*.pdbqt`.

---

## ESTADO POR ETAPA

| Etapa | Estado | Dónde está |
|---|---|---|
| 1 · BLAST Swiss-Prot | ✅ HECHA | `01_blast/swissprot_hits.tsv` |
| 1 · BLAST nr | ✅ HECHA (500 aciertos, 420 >40 % id, 259 Rhizobiaceae) | `01_blast/nr_hits.tsv` |
| 2 · Estructura AlphaFold2 | ✅ HECHA | `02_estructura/TeuB_full_146f4/` |
| 3 · Péptido señal + recorte | ✅ HECHA | `03_maduro/TeuB_maduro_1-335.pdb` |
| 4 · Foldseek + panel | ✅ HECHA | `04_foldseek/hits.tsv`, `tabla_ligandos.json` |
| 5 · Ligandos | ✅ HECHA | `05_ligandos/*.sdf .mol2 .pdbqt` (15) |
| 6 · Docking Vina | ✅ HECHA (45/45) — falta análisis de poses | `06_docking/`, `08_analisis/vina_scores.tsv` |
| 7 · Co-plegamiento Boltz-2 | 🔄 11 válidas (BGC RIP GAL XYP PAV INS AHR ALL WEB GZL 3VB); faltan ARA BDR FRU X9X HPA | `07_coplegamiento/`, `08_analisis/revisar_boltz.py`, `boltz_resumen.tsv` |
| 8-10 · Ranking y entregable | ⬜ PENDIENTE | `08_analisis/` |

---

## DECISIONES TOMADAS (hay que declararlas en el entregable)

1. **Todo local salvo las etapas 2 y 7.** La práctica está diseñada para navegador; aquí se usaron las mismas herramientas instaladas en la laptop, que es legítimo pero hay que decirlo.
   - `blastp` local contra swissprot descargada del NCBI (no la web del NCBI).
   - `foldseek` local contra la base **PDB completa** (6.5 GB), *no* PDB100. Ventaja: se ven todas las cadenas, así que no hace falta el paso manual de "buscar otras entradas de la misma proteína" — las entradas con ligando aparecen directamente.
   - `AutoDock Vina` local (el mismo motor que SwissDock usa por dentro), no el servidor.
   - Ligandos preparados con Open Babel desde los SDF `_ideal` del RCSB.
2. **Etapa 2: ColabFold en Colab con GPU**, no la versión WebGPU. `use_templates: false` ✅, `msa_mode: mmseqs2_uniref_env` ✅, `random_seed: 0`, 5 modelos, 3 reciclos.
   - Hubo una **corrida fallida previa** con la secuencia de ejemplo de la libreta (59 aa). Archivada en `02_estructura/_descartado_secuencia_de_ejemplo_59aa/`. No usar.
3. **Recorte del péptido señal: opción (b), recortar el PDB**, no volver a predecir. Coste a declarar: el resto de la proteína se plegó en presencia del péptido señal.
4. **NUMERACIÓN: cadena madura 1-335.** Residuo 1 = D de `DDTIA` = residuo 28 del precursor. Fijada, no cambiar.

## ⚠️ HALLAZGO IMPORTANTE: el modelo de respaldo de la práctica tiene un error

`00_datos_repartidos/TeuB_modelo_apo.pdb` (el que reparte el sitio web) tiene **334 residuos, no 335**.
En la posición 256 la secuencia real dice `...LNGWM**RR**WNDEK...` y el respaldo tiene `...LNGWM**R**WNDEK...`: **le falta una arginina**.
Consecuencia: todo lo posterior al residuo 256 va corrido en 1 respecto al modelo propio.
- El sitio que el guion da como `Tyr38 Ser43 His45 Arg174 Trp197 Asp222 Glu248 **Cys268**`
- en la numeración correcta (1-335) es `Y38 S43 H45 R174 W197 D222 E248 **C269**`
- Verificado: los 8 coinciden con la numeración propia. Esto va en el entregable.

---

## RESULTADOS YA OBTENIDOS

### Etapa 1 · BLAST (Swiss-Prot, E<1e-3)
- Secuencia verificada: **362 aa**.
- **17 aciertos**, todos proteínas periplásmicas de unión a sustrato (familia RbsB).
- Mejor acierto: **P36949** *Ribose import binding protein RbsB* (B. subtilis), **28.5 % id**, 61 % cobertura, E=7.7e-20.
- **6 ligandos distintos entre los 10 mejores**: ribosa, D-treitol, xilitol, ribosa/alosa, galactofuranosa, D-apiosa.
- Aciertos #1 (28.5 %) y #2 (31.5 %) tienen identidad casi igual y unen ligandos **distintos** (ribosa vs treitol).
- Identidades 22.6–31.5 % = **zona crepuscular**. Se puede afirmar familia y pliegue, NO el ligando.
- **qstart mínimo = 49, mediana = 88**: ningún homólogo alinea sobre los primeros 48 residuos.
- Aciertos #8–#11 son *D-apiose import binding protein*, y **#9 es de *Rhizobium rhizogenes*, #11 de *Rhizobium etli*** — mismo género que el organismo de origen.

### Etapa 2 · AlphaFold2 (ColabFold, job `TeuB_full_146f4`)
- Profundidad del alineamiento: **1841 secuencias**.
- Mejor modelo: `rank_001_alphafold2_ptm_model_3_seed_000`, **pLDDT 92.2, pTM 0.877**.
- Los 5 modelos dan pLDDT 91.9–92.2 y pTM 0.865–0.877: muy consistentes entre sí.
- pLDDT mínimo = **29.1 en el residuo 8** (dentro del péptido señal).
- **pLDDT medio residuos 1-27 = 36.3 · residuos 28-362 = 96.8 · diferencia de 60.5 puntos.**
- PAE: señal vs maduro = **29.6 Å** (el modelo no sabe dónde poner el péptido señal); dentro del maduro = 3.0 Å.
- Mejor partición en dos dominios: corte en el residuo 162, contraste inter–intra de solo **0.99 Å**
  → el modelo **sí** sabe cómo se orientan los dos lóbulos → **hendidura cerrada**, como anticipa el guion para la etapa 6.

### Etapa 3 · Péptido señal
Tres evidencias independientes coinciden en el corte 27/28:
- **Secuencia**: n(1-5) `MKRRT` carga +3 · h(6-22) `FLQTGSALIAAGAFGIP` Kyte-Doolittle +1.27 (pico +2.24 en 13-21) · c(23-27) `GILRA`. La hidropatía se desploma en 28.
- **pLDDT**: 36.3 vs 96.8. Sube y se estabiliza a partir del residuo 28-31.
- **BLAST**: ningún homólogo alinea antes del residuo 49.
- Nota para el entregable: BLAST empieza en 49, no en 28, porque los primeros ~20 residuos de la cadena madura tampoco alinean bien — son el extremo N variable del dominio, no parte del péptido señal.
- Cadena madura: **335 aa, empieza en `DDTIA`**.

### Etapa 4 · Foldseek
- **2811 aciertos en 7.6 s** contra la base PDB completa.
- Mejores: TM-score **0.76–0.92**, prob 1.00, E hasta 1e-23, pero **identidad de secuencia solo 17–25 %**.
  → Comparación directa con BLAST (22–31 % id, 17 aciertos): Foldseek encuentra dos órdenes de magnitud más parientes, en la zona crepuscular donde BLAST ya no resuelve.
- De **150 entradas únicas** analizadas vía la API GraphQL del RCSB: **67 apo (45 %)**, 83 con ligando.
- Ligandos más frecuentes: BGC 13 · RIP 9 · GAL 7 · HPA 6 · INS 4 · XYP 4 · PAV 4 · ARA 2 · FRU 2 · GLA 2 · ARB 2.
- **Fórmulas repetidas**: C5H10O5 con **11 ligandos distintos**, C6H12O6 con **10**. Esa es la dificultad, cuantificada.

### Panel de 15 elegido (confirmado por el usuario)

| cód | molécula | fórmula | n | criterio |
|---|---|---|---|---|
| BGC | β-D-glucopiranosa | C6H12O6 | 13 | el más frecuente |
| RIP | β-D-ribopiranosa | C5H10O5 | 9 | 2º; la etiqueta que daba BLAST |
| GAL | β-D-galactopiranosa | C6H12O6 | 7 | 3º; epímero de BGC en C4 |
| XYP | β-D-xilopiranosa | C5H10O5 | 4 | frecuente |
| PAV | **D-apiosa** | C5H10O5 | 4 | convergencia BLAST+Foldseek, y homólogos de *Rhizobium* |
| INS | *mio*-inositol | C6H12O6 | 4 | **señuelo 1**: ciclitol, misma fórmula que BGC |
| ARA | α-L-arabinopiranosa | C5H10O5 | 2 | anillo de 6 |
| AHR | α-L-arabinofuranosa | C5H10O5 | 1 | anillo de 5, **mismo azúcar que ARA** |
| BDR | β-D-ribofuranosa | C5H10O5 | 1 | anillo de 5 frente a RIP (6) |
| ALL | β-D-alopiranosa | C6H12O6 | 1 | aldosa |
| WEB | α-D-psicopiranosa | C6H12O6 | 1 | cetosa — **el par imposible con ALL** (anexo) |
| FRU | β-D-fructofuranosa | C6H12O6 | 2 | cetosa de anillo 5 |
| GZL | β-D-galactofuranosa | C6H12O6 | 1 | par piranosa/furanosa contra GAL |
| X9X | D-alitol | C6H14O6 | 1 | **señuelo 2**: poliol acíclico |
| HPA | hipoxantina | C5H4N4O | 6 | **señuelo 3**: base púrica, ni azúcar ni poliol |

En el panel: 6 ligandos comparten C5H10O5 y 7 comparten C6H12O6.
Los 15 pasaron las tres verificaciones del guion: fórmula = ficha, tienen hidrógenos, no son planos.

### Etapa 6 · Docking (parámetros fijos, idénticos para los 15)
```
receptor : 06_docking/receptor.pdbqt  (desde 03_maduro/TeuB_maduro_1-335.pdb, obabel -xr -p 7.4)
centro   : x=3.38  y=-3.17  z=0.79    <- CALCULADO en las coordenadas propias
size     : 26 x 26 x 26 Å
exhaustiveness 32, num_modes 9, cpu 8
semillas : 101, 202, 303  (3 réplicas por ligando -> piso de ruido medido, no estimado)
```
El centro es el centroide de los átomos de cadena lateral de `Y38 S43 H45 R174 W197 D222 E248 C269`.
**NO** son los 3.4/4.0/1.3 del guion: esos pertenecen al modelo de respaldo, que está en otro sistema de coordenadas.

Primeros resultados parciales: BGC ≈ -5.63/-5.67/-5.65 · RIP ≈ -6.00/-6.00
(la dispersión de BGC entre semillas es ~0.04 kcal/mol — ese es el orden del piso de ruido.)

---

## LO QUE FALTA

### A) Terminar el docking (automático, en tmux)
Cuando `ls 06_docking/out_*.pdbqt | wc -l` dé 45, analizar:
- mejor score por ligando y por semilla
- media ± dispersión de las 3 réplicas = **piso de ruido del método 1**
- si la mejor pose cayó **dentro del bolsillo o en la superficie**: medir la distancia del centroide de la pose al centro de la caja y los contactos con los 8 residuos del sitio. Marcar las superficiales.

### B) Etapa 7 · Co-plegamiento — REQUIERE ACCIÓN DEL USUARIO EN COLAB
La laptop **no tiene GPU NVIDIA** (Intel Iris Xe), así que Boltz-2 tiene que correr en Colab.
En la misma libreta ColabFold2_preview que ya se usó para la etapa 2:
- celda de instalación: `model = boltz2`
- `protein` = la secuencia madura de `03_maduro/teuB_maduro_solo_letras.txt` (335 aa, empieza DDTIA) — **la madura, no las 362**
- `ligand_ccd` = el código de 3 letras (uno por corrida), `ligand_smiles` vacío
- `jobname` = el código del ligando
- `msa_mode = mmseqs2_server`, `seeds = 1,2,3`
- num_recycles y num_diffusion_samples por omisión → 15 modelos por ligando
- **plantillas desactivadas** (si no, el servidor puede darle la estructura cristalográfica como molde y se anula el control de contaminación)
Repetir para los 15 códigos. Descargar cada ZIP a `~/Downloads`; se desempacan en `07_coplegamiento/`.
Leer de cada uno: `*_summary_confidences.json` → `iptm` y `has_clash`; `*_ranking_scores.csv` → las 15 filas (dispersión = **piso de ruido del método 2**); `*_model.cif` → mirar si el ligando quedó en la hendidura y si el anillo está deformado.

Si 15 ligandos × 15 modelos es demasiado tiempo de GPU, recortar a los 10 que se van a entregar y decirlo.

### C) BLAST contra nr (solo para contar)
Si `01_blast/nr_hits.tsv` está vacío:
```bash
cd /home/jmanuel/cursos/proteinas/tareasugar/01_blast
setsid nohup blastp -query teuB.fasta -db nr -remote -evalue 1e-5 -max_target_seqs 500 \
  -outfmt "6 sseqid stitle pident length qstart qend evalue bitscore qcovs" -out nr_hits.tsv &
```
Solo hace falta para: cuántos aciertos pasan de 40 % de identidad, y cuántos son de *Rhizobium* o géneros cercanos.

### D) Pendiente del usuario
- **Escribir la predicción antes de ver los números** (qué ligando cree que gana y por qué). Es un punto explícito del entregable y aún no la ha dado.

### E) Etapas 8-10 · Ranking y entregable
Tabla comparativa de 10 filas: ligando · puesto Vina · score Vina · puesto ipTM · ipTM · Δ de puestos.
Marcar empates con los dos pisos de ruido. Discutir: los de fórmula idéntica, dónde quedó cada señuelo, y si algún método acertó la predicción previa.
Entregable: 3 páginas máximo. Pesos: ranking 30 %, piso de ruido 20 %, comparación 20 %, construcción del panel 15 %, cadena de inferencias + péptido señal 10 %, veredicto 5 %.

---

## MAPA DE ARCHIVOS

```
00_datos_repartidos/   lo extraído del HTML: TeuB_modelo_apo.pdb (¡con el error de la R!), 10 .mol2, panel-azucares.smi
01_blast/              teuB.fasta, teuB_solo_letras.txt, swissprot_hits.tsv, nr_hits.tsv, base swissprot local
02_estructura/         TeuB_full_146f4/ (el bueno) · _descartado_secuencia_de_ejemplo_59aa/ (la corrida fallida)
03_maduro/             TeuB_maduro_1-335.pdb, teuB_maduro.fasta, teuB_maduro_solo_letras.txt
04_foldseek/           pdb_db* (6.5 GB), hits.tsv (2811), entradas_ligandos.json, tabla_ligandos.json
05_ligandos/           15 × {SDF ideal, .mol2, .pdbqt}
06_docking/            receptor.pdbqt, box.txt, run.sh, out_<COD>_s<SEED>.pdbqt, log_<COD>_s<SEED>.txt
07_coplegamiento/      (vacío, pendiente de Colab)
08_analisis/           (vacío)
CUADERNO.md            notas acumuladas para el entregable
ESTADO.md              este archivo
```

## HERRAMIENTAS INSTALADAS EN ESTA SESIÓN
`tmux` (apt) · `openbabel`, `foldseek`, `vina`, `numpy`, `pandas`, `matplotlib`, `biopython` (mamba, entorno base de miniforge3).
Ya estaban: `blastp`, `pymol`, `curl`, `jq`, `xclip`, `google-chrome`.
