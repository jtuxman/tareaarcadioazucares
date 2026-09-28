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
