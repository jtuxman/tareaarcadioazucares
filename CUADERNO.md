# Cuaderno · De una secuencia a un ranking de ligandos
Curso Evolución y Diseño de Proteínas 2027-1 · jmanuel@ibt.unam.mx

## Decisiones declaradas (para el entregable)
- Etapa 1: blastp 2.x **local** contra swissprot descargada del NCBI (no la web del NCBI). Contra nr: blastp -remote al servidor del NCBI.
- Etapa 2: AlphaFold2-WebGPU en Chrome, GPU Intel Iris Xe. Sin plantilla. MMseqs2 server. Semilla: _(pendiente)_
- Etapa 6: **AutoDock Vina local** (mismo motor que usa SwissDock), no la web. Permite medir el piso de ruido con réplicas reales.
- Etapa 7: Boltz-2 en Colab (la laptop no tiene GPU NVIDIA).
- NUMERACIÓN ELEGIDA: **cadena madura, 1-335** (residuo 1 = D de DDTIA = residuo 28 del precursor).

## Etapa 1 · BLAST
- Secuencia verificada: 362 aa.
- Swiss-Prot, E<1e-3: **17 aciertos**, todos proteínas periplásmicas de unión a sustrato (familia RbsB).
- Mejor acierto: P36949 Ribose import binding protein RbsB (B. subtilis), 28.5 % id, 61 % cobertura, E=7.7e-20.
- **6 ligandos distintos entre los 10 mejores**: ribosa, D-treitol, xilitol, ribosa/alosa, galactofuranosa, D-apiosa.
- Aciertos #1 (28.5 %) y #2 (31.5 %) tienen identidad casi igual y unen ligandos distintos (ribosa vs treitol).
- Zona crepuscular: identidades 22.6-31.5 %. Se puede afirmar familia y pliegue, NO el ligando.
- **qstart mínimo = 49, mediana = 88**: ningún homólogo alinea sobre los primeros 48 residuos.

## Etapa 3 · Péptido señal (evidencia de secuencia)
- n (1-5) MKRRT, carga neta +3
- h (6-22) FLQTGSALIAAGAFGIP, Kyte-Doolittle medio +1.27 (pico +2.24 en 13-21)
- c (23-27) GILRA, corte tras A27
- Caída de hidropatía justo en 28; maduro empieza en DDTIA, 335 aa. Coincide con la hipótesis del guion.

## Etapa 2 · AlphaFold2 (ColabFold, job TeuB_full_146f4, GPU de Colab)
- MSA: **1841 secuencias**. Sin plantilla. msa_mode mmseqs2_uniref_env. random_seed=0, 5 modelos, 3 reciclos.
- Mejor modelo rank_001 (model_3): **pLDDT 92.2, pTM 0.877**. Los 5 modelos: pLDDT 91.9-92.2, pTM 0.865-0.877.
- pLDDT mínimo 29.1 en el residuo 8.
- pLDDT 1-27 = **36.3** vs 28-362 = **96.8** (diferencia 60.5).
- PAE señal-vs-maduro = 29.6 A; dentro del maduro = 3.0 A.
- Dos dominios con corte en 162, contraste inter-intra 0.99 A => el modelo SÍ sabe orientarlos => **hendidura cerrada**.
- (Corrida previa fallida con la secuencia de ejemplo de la libreta, 59 aa. Archivada, no usada.)

## Etapa 3 · Recorte
- Corte 27/28. Cadena madura 335 aa, empieza DDTIA. Opción (b): recortar el PDB, no volver a predecir.
- **Numeración fijada: 1-335 sobre la cadena madura.**
- HALLAZGO: al modelo de respaldo repartido le falta una Arg en la posición 256 (334 residuos en vez de 335).
  Su Cys268 es la Cys269 en la numeración correcta. Los 8 residuos del sitio verificados: Y38 S43 H45 R174 W197 D222 E248 C269.

## Etapa 4 · Foldseek (base PDB completa, local)
- 2811 aciertos en 7.6 s. Mejores: TM 0.76-0.92, prob 1.00, E hasta 1e-23, **id de secuencia 17-25 %**.
- Contraste con BLAST: 22-31 % id y solo 17 aciertos. Foldseek resuelve donde BLAST ya no.
- 150 entradas únicas: 67 apo (45 %), 83 con ligando.
- Frecuencias: BGC 13, RIP 9, GAL 7, HPA 6, INS 4, XYP 4, PAV 4, ARA/FRU/GLA/ARB 2.
- Fórmulas repetidas: **C5H10O5 con 11 ligandos, C6H12O6 con 10**.
- PAV = D-apiosa: converge con los aciertos #8-#11 de BLAST (*D-apiose import binding protein*),
  dos de ellos de *Rhizobium rhizogenes* y *Rhizobium etli* — mismo género que el organismo de origen.

## Panel de 15
BGC RIP GAL XYP PAV INS ARA AHR BDR ALL WEB FRU GZL X9X HPA
3 pares piranosa/furanosa (ARA/AHR, RIP/BDR, GAL/GZL), 2 aldosa/cetosa (ALL/WEB, BGC/FRU),
3 señuelos de dureza creciente (INS ciclitol, X9X poliol, HPA base púrica).
Las 15 pasan: fórmula = ficha, tienen H, no son planas.

## Etapa 6 · Docking (parámetros)
receptor desde TeuB_maduro_1-335.pdb; centro 3.38 / -3.17 / 0.79 (centroide de cadenas laterales del sitio);
caja 26 A; exhaustiveness 32; semillas 101/202/303 para los 15 => piso de ruido medido.

## Predicción (punto explícito del entregable) — fijada 2026-09-27
**Ligando predicho: D-apiosa (PAV).** Segunda opción: ribosa (RIP).

Formulada con la evidencia de las etapas 1–4:
- El acierto #1 de BLAST (RbsB, ribosa) tiene 28.5 % de identidad, zona crepuscular; el #2, con identidad casi igual (31.5 %), une otra molécula (treitol). Una sola etiqueta de BLAST no basta.
- La D-apiosa es el único ligando presente a la vez en Swiss-Prot (aciertos #8–#11, *D-apiose import binding protein*, dos de ellos de *Rhizobium rhizogenes* y *R. etli*, el mismo género del organismo de origen) y en Foldseek (4 estructuras con PAV).
- Coherente con la ecología: *Rhizobium* vive asociado a raíces, y la apiosa es un azúcar característico de la pared celular vegetal.
- Los parientes casi idénticos en nr (97–100 % id, *Rhizobium*) anotados como "ribose transport system substrate-binding protein" NO se tomaron como evidencia: son anotaciones automáticas heredadas por homología, circulares.

**Declaración de honestidad:** los scores de Vina se vieron antes de escribir esta predicción, pero no se usaron para formularla.

## Etapa 1 · BLAST contra nr (remoto, E<1e-5, 500 aciertos máx.)
- 500 aciertos; **420 con >40 % de identidad**; **259 de Rhizobiaceae** (*Martinezella/Rhizobium*, *Mesorhizobium*, *Sinorhizobium*, *Agrobacterium*...).
- Mejores: 97–100 % id, todos de *Rhizobium*/*Martinezella*; anotaciones genéricas ("ABC transporter substrate-binding protein") o "ribose transport..." automáticas.

## Etapa 6 · Docking — resultados (16 ligandos × 3 semillas; 3VB añadido después con los mismos parámetros)
- Ruido entre semillas: rango máximo **0.035 kcal/mol** (BGC); DE típica 0.002–0.018. Rango del panel: **1.43 kcal/mol**.
- Con ese piso, las 16 medias son distinguibles entre sí — pero el piso de semillas solo mide la búsqueda, no el error de la función de puntaje (~2 kcal/mol frente a experimento), que es mayor que todo el rango del panel.
- **Todas las mejores poses caen en el bolsillo**: centroide a ≤1.5 Å del centro del sitio, 7–8/8 residuos del sitio a <4 Å. Ninguna superficial (la hendidura cerrada y la caja lo fuerzan).
- Ranking Vina: XYP −6.41 > ARA −6.14 > **HPA −6.10** > **INS −6.09** > GAL > RIP > ALL > GZL > **PAV −5.75 (9.º)** > FRU > BDR > BGC > AHR > WEB > **X9X −5.23** > 3VB −4.97.
- Señuelos: HPA 3.º e INS 4.º (¡por encima de casi todos los azúcares!); X9X 15.º. Vina no rechaza señuelos.

## Etapa 7 · Boltz-2 — resultados (16 ligandos × 15 modelos, sin plantillas, secuencia madura)
- ipTM 0.95–0.98 para TODOS. `ranking_score` medio: rango **0.034**; DE entre los 15 modelos 0.002–0.010 (media 0.0054) = **piso de ruido del método 2**.
- Criterio de empate: |Δ| < 2·√(sd₁²+sd₂²). **Los 12 primeros empatan con el líder (AHR)**; solo WEB, 3VB, FRU, X9X se separan por abajo. Solo 17 de 120 pares son distinguibles.
- Ligando en el sitio en los 15 modelos de todos: ≤2.2 Å del centro, ≥6/8 contactos.
- Señuelos: HPA 7.º, INS 10.º (dentro del empate), X9X último (16.º). Boltz solo rechaza (débilmente) el poliol acíclico.
- **HALLAZGO: Boltz-2 elimina el átomo llamado O1** en los 10 ligandos que lo tienen (RIP XYP ARA AHR BDR BGC GAL ALL GZL: el OH anomérico; FRU: el O del CH2OH). INS, X9X, HPA, 3VB, WEB, PAV salen completos. Para las aldosas el modelo no distingue α de β (el OH anomérico no está) y co-pliega un residuo glicosilo, no el azúcar libre. Esto hay que declararlo: la comparación Boltz entre azúcares no es entre las moléculas exactas.

## Etapas 8–10 · Comparación (`08_analisis/tabla_comparativa.tsv`)
| lig | puesto Vina | Vina | puesto Boltz | ranking_score | Δ |
|---|---|---|---|---|---|
| XYP | 1 | −6.41 | 5 | 0.963 | +4 |
| ARA | 2 | −6.14 | 3 | 0.963 | +1 |
| HPA* | 3 | −6.10 | 7 | 0.961 | +4 |
| INS* | 4 | −6.09 | 10 | 0.958 | +6 |
| GAL | 5 | −6.05 | 11 | 0.958 | +6 |
| RIP | 6 | −6.00 | 4 | 0.963 | −2 |
| ALL | 7 | −5.95 | 8 | 0.958 | +1 |
| GZL | 8 | −5.78 | 6 | 0.962 | −2 |
| PAV | 9 | −5.75 | 12 | 0.955 | +3 |
| FRU | 10 | −5.73 | 15 | 0.946 | +5 |
| BDR | 11 | −5.71 | 2 | 0.970 | −9 |
| BGC | 12 | −5.65 | 9 | 0.958 | −3 |
| AHR | 13 | −5.50 | 1 | 0.971 | −12 |
| WEB | 14 | −5.46 | 13 | 0.951 | −1 |
| X9X* | 15 | −5.23 | 16 | 0.937 | +1 |
| 3VB | 16 | −4.97 | 14 | 0.949 | −2 |
(* señuelo)
- Correlación de rangos: **Spearman ρ = 0.43 (p = 0.10), Kendall τ = 0.33 (p = 0.08)** → acuerdo débil, no significativo.
- Mayor desacuerdo: las **furanosas** (AHR 13.º→1.º, BDR 11.º→2.º): Vina las castiga, Boltz las prefiere.
- Poses: el mejor modelo Boltz (superpuesto por Cα sobre el receptor de Vina) y la pose 1 de Vina coinciden: centroides a 0.5–1.4 Å, RMSD de átomos pesados 1.1–2.0 Å. **Los dos métodos ponen el ligando en el mismo lugar; discrepan en el puntaje.**
- Predicción (PAV): 9.º en Vina, 12.º en Boltz (empatado con el líder). **Ningún método la confirma, y tampoco la descarta**: en Boltz está dentro del empate.
- Pares de fórmula idéntica: C5H10O5 (XYP, ARA, RIP, PAV, AHR, BDR) ocupan en Vina los puestos 1, 2, 6, 9, 13, 11 y en Boltz 5, 3, 4, 12, 1, 2 — el orden dentro de la misma fórmula no es reproducible entre métodos.
