---
title: "TeuB: de una secuencia sin etiqueta a un ranking de ligandos"
subtitle: "Evolución y Diseño de Proteínas 2027-1 · jmanuel@ibt.unam.mx"
---

**Numeración:** cadena madura 1-335. El residuo 1 es la D de `DDTIA`, que corresponde al residuo 28 del precursor. **Panel:** 15 ligandos elegidos a partir de Foldseek, más D-treitol (3VB), que viene del panel del tutorial; en total, 16. **Métodos:** AutoDock Vina 1.2 (local) y Boltz-2 (ColabFold2, sin plantillas).

## 1 · El ranking

Vina: la mejor pose de cada semilla, promediada sobre 3 semillas (101, 202, 303), ± desviación estándar. Boltz-2: ipTM promedio de los 15 modelos (semillas 1-3 × 5 muestras), ± desviación estándar. Los ligandos marcados con \* son señuelos. **Empates** según el criterio de la §2: en Vina no hay ninguno; en Boltz, los 10 del top-10 empatan con el primero, y también GAL y PAV.

| # | **Vina** | kcal/mol | | # | **Boltz-2** | ipTM (todos empatados con el 1.º) |
|---|---|---|---|---|---|---|
| 1 | XYP β-D-xilopiranosa | −6.408 ± 0.008 | | 1 | AHR α-L-arabinofuranosa | 0.975 ± 0.008 |
| 2 | ARA α-L-arabinopiranosa | −6.145 ± 0.012 | | 2 | BDR β-D-ribofuranosa | 0.973 ± 0.012 |
| 3 | HPA hipoxantina\* | −6.102 ± 0.004 | | 3 | GZL β-D-galactofuranosa | 0.963 ± 0.007 |
| 4 | INS *mio*-inositol\* | −6.091 ± 0.004 | | 4 | XYP β-D-xilopiranosa | 0.961 ± 0.003 |
| 5 | GAL β-D-galactopiranosa | −6.049 ± 0.011 | | 4 | RIP β-D-ribopiranosa | 0.961 ± 0.003 |
| 6 | RIP β-D-ribopiranosa | −6.003 ± 0.005 | | 4 | ARA α-L-arabinopiranosa | 0.961 ± 0.003 |
| 7 | ALL β-D-alopiranosa | −5.949 ± 0.004 | | 7 | HPA hipoxantina\* | 0.961 ± 0.007 |
| 8 | GZL β-D-galactofuranosa | −5.783 ± 0.008 | | 8 | ALL β-D-alopiranosa | 0.957 ± 0.010 |
| 9 | PAV D-apiosa | −5.751 ± 0.011 | | 9 | BGC β-D-glucopiranosa | 0.957 ± 0.009 |
| 10 | FRU β-D-fructofuranosa | −5.725 ± 0.002 | | 10 | INS *mio*-inositol\* | 0.957 ± 0.006 |

En Vina quedan fuera del top-10 BDR (−5.71), BGC (−5.65), AHR (−5.50), WEB (−5.46), X9X\* (−5.23) y 3VB (−4.97). En Boltz quedan fuera GAL (0.956) y PAV (0.954), que siguen dentro del empate, y WEB (0.949), 3VB (0.946), FRU (0.942) y X9X\* (0.929).

## 2 · Pisos de ruido y empates

**Criterio:** dos ligandos empatan cuando la diferencia entre sus medias es menor que 2·√(DE₁² + DE₂²), el doble del error de la diferencia.

- **Vina: medido con 3 réplicas de semilla por ligando.** Las desviaciones estándar van de 0.002 a 0.018 kcal/mol y el rango máximo entre semillas es de 0.035 kcal/mol (BGC). Con este piso, **los 120 pares son distinguibles y no hay empates**; el par más cercano es HPA/INS, con Δ = 0.011 frente a un umbral de 0.011. *Pero este piso solo mide la estocasticidad de la búsqueda, no el error de la función de puntaje.* Ese error ronda los 2 kcal/mol frente a datos experimentales y es mayor que el rango completo del panel (1.43 kcal/mol). El ranking de Vina es muy reproducible, pero eso no garantiza que sea correcto.
- **Boltz-2: medido con los 15 modelos por ligando.** Las desviaciones estándar del ipTM van de 0.003 a 0.012 y el rango del panel es de solo 0.046. **Los 12 primeros empatan con el líder (AHR).** Solo WEB, 3VB, FRU y X9X se separan, y lo hacen hacia abajo. Solo 15 de los 120 pares superan el piso. El ipTM se reporta con 2 decimales por modelo, lo que agrega ±0.005 de redondeo.

## 3 · Comparación de los dos métodos

| Lig. | P. Vina | Vina | P. Boltz | ipTM | **Δ** | | Lig. | P. Vina | Vina | P. Boltz | ipTM | **Δ** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| XYP | 1 | −6.41 | 4 | 0.961 | +3 | | PAV | 9 | −5.75 | 12 | 0.954 | +3 |
| ARA | 2 | −6.14 | 4 | 0.961 | +2 | | FRU | 10 | −5.73 | 15 | 0.942 | +5 |
| HPA\* | 3 | −6.10 | 7 | 0.961 | +4 | | BDR | 11 | −5.71 | 2 | 0.973 | **−9** |
| INS\* | 4 | −6.09 | 10 | 0.957 | +6 | | BGC | 12 | −5.65 | 9 | 0.957 | −3 |
| GAL | 5 | −6.05 | 11 | 0.956 | +6 | | AHR | 13 | −5.50 | 1 | 0.975 | **−12** |
| RIP | 6 | −6.00 | 4 | 0.961 | −2 | | WEB | 14 | −5.46 | 13 | 0.949 | −1 |
| ALL | 7 | −5.95 | 8 | 0.957 | +1 | | X9X\* | 15 | −5.23 | 16 | 0.929 | +1 |
| GZL | 8 | −5.78 | 3 | 0.963 | −5 | | 3VB | 16 | −4.97 | 14 | 0.946 | −2 |

- **El acuerdo es débil y no significativo:** Spearman ρ = 0.39 (p = 0.14) con ipTM, y ρ = 0.43 (p = 0.10) con `ranking_score`.
- **El mayor desacuerdo está en las furanosas.** AHR (Δ −12), BDR (Δ −9) y GZL (Δ −5) quedan abajo en Vina y hasta arriba en Boltz. Una explicación posible para Vina: cada furanosa tiene una torsión activa más que su piranosa (AHR 5 frente a ARA 4, BDR 5 frente a RIP 4, GZL 7 frente a GAL 6), y el término de entropía conformacional de Vina penaliza cada torsión. Boltz no tiene esa penalización.
- **Los dos métodos colocan el ligando en el mismo lugar.** Todas las poses están en el bolsillo: en Vina, el centroide queda a ≤ 1.5 Å del centro del sitio y la pose toca 7-8 de los 8 residuos (Y38 S43 H45 R174 W197 D222 E248 C269); en Boltz, a ≤ 2.2 Å y con ≥ 6 de 8 contactos, en los 15 modelos. Al superponer el mejor modelo de Boltz sobre el receptor de Vina (por Cα), la pose de Vina y la de Boltz difieren en **1.1-2.0 Å de RMSD de átomos pesados**. El desacuerdo está en el puntaje, no en la geometría.
- **Fórmulas idénticas.** Las 6 moléculas C₅H₁₀O₅ ocupan los puestos 1, 2, 6, 9, 11 y 13 en Vina, y 1, 2, 4, 4, 4 y 12 en Boltz. **El orden dentro de una misma fórmula no se reproduce entre métodos.** Los pares piranosa/furanosa del mismo azúcar se invierten: ARA/AHR queda 2.º/13.º en Vina y 4.º/1.º en Boltz; RIP/BDR, 6.º/11.º y 4.º/2.º.
- **Señuelos.** Ningún método rechaza a HPA ni a INS. En Vina quedan 3.º y 4.º, por delante de casi todos los azúcares. En Boltz quedan 7.º y 10.º, dentro del empate. Los dos métodos coinciden solo en poner al poliol acíclico X9X al fondo (15.º y 16.º).
- **⚠ Artefacto de Boltz-2.** En los 10 ligandos cuyo CCD tiene un átomo llamado **O1**, Boltz lo elimina: el OH anomérico de RIP, XYP, ARA, AHR, BDR, BGC, GAL, ALL y GZL, y el O del CH₂OH de FRU. INS, X9X, HPA, 3VB, WEB y PAV salen completos. Para las aldosas, **Boltz modela un residuo glicosilo y no el azúcar libre, así que no puede distinguir α de β**. Parte de la comparación con Boltz no se hace entre las mismas moléculas que dockeó Vina.

## 4 · Cadena de inferencias

| Etapa | Herramienta | Entrada | Qué salió | Qué concluyo · qué **no** |
|---|---|---|---|---|
| 1 BLAST | blastp 2.x local (Swiss-Prot) y remoto (nr) | 362 aa | Swiss-Prot: 17 aciertos (E < 1e-3), todos de la familia RbsB, 22.6-31.5 % de identidad; 6 ligandos distintos entre los 10 primeros. nr: 500 aciertos, 259 de Rhizobiaceae | Es una proteína periplásmica de unión a azúcar, con pliegue de tipo I. **No** asigna ligando: la identidad está en zona crepuscular y los aciertos #1 y #2 unen ligandos distintos. Las anotaciones "ribose" de nr son automáticas y circulares |
| 2 AlphaFold2 | ColabFold (GPU de Colab), sin plantillas, 5 modelos | 362 aa | pLDDT 92.2, pTM 0.877; MSA de 1841 secuencias; PAE entre lóbulos de 3.0 Å | Modelo confiable, con la hendidura cerrada. **No** muestra el ligando ni la conformación abierta |
| 3 Péptido señal | análisis n/h/c, pLDDT, alineamiento de BLAST | modelo y secuencia | Corte 27/28 | Cadena madura de 335 aa (ver §6) |
| 4 Foldseek | foldseek local contra el PDB completo | madura, 1-335 | 2811 aciertos; TM-score de 0.76 a 0.92 con solo 17-25 % de identidad; 150 entradas, 83 con ligando | Recupera parientes que BLAST no ve. **No** dice qué ligando: C₅H₁₀O₅ aparece con 11 códigos y C₆H₁₂O₆ con 10 |
| 5 Ligandos | CCD `_ideal.sdf` del RCSB y Open Babel | 16 códigos | mol2 y pdbqt; fórmula = ficha, con H, no planos | Panel listo. **No** incluye todos los tautómeros ni anómeros |
| 6 Docking | Vina 1.2, caja de 26 Å, exhaustiveness 32, 3 semillas | receptor maduro y 16 ligandos | −6.41 a −4.97 kcal/mol; poses en el sitio | Ranking reproducible. **No** es afinidad: el error del puntaje es mayor que el rango |
| 7 Co-plegamiento | Boltz-2, sin plantillas, 3 semillas × 5 muestras | secuencia madura y código CCD | ipTM medio de 0.93 a 0.98; 12 de 16 empatados | Todos los ligandos caben en el sitio. **No** discrimina entre azúcares. Artefacto O1 |

## 5 · Construcción del panel y predicción

De las 150 entradas únicas de Foldseek, 67 no tienen ligando (45 %). Conté los ligandos de las 83 restantes. El panel siguió cuatro criterios: **(a) frecuencia:** BGC (13: 2GBP, 3GBP, 2IPL…), RIP (9: 1DBP, 2FN8, 2IOY…), GAL (7: 2RJO, 8ABP…), XYP (4: 3C6Q, 5XSS…); **(b) convergencia con BLAST:** PAV (4: 3EJW, 3T95, 4PZ0, 6DSP), la única función que aparece en Swiss-Prot en el mismo género, *Rhizobium*; **(c) pares que pusieran a prueba la geometría con la misma fórmula:** piranosa/furanosa (ARA 1BAP/4RXT frente a AHR 5OCP, RIP frente a BDR 3KSM, GAL frente a GZL 2VK2) y aldosa/cetosa (ALL 5DTE frente a WEB 8WLB, BGC frente a FRU 2X7X/3BRQ); **(d) señuelos de dificultad creciente:** INS (4: 4IRX…, ciclitol con la misma fórmula que BGC), X9X (4WT7, poliol acíclico) y HPA (6: 1JFT…, base púrica). 3VB (4RSM) viene del panel del tutorial. El panel tiene 6 ligandos C₅H₁₀O₅ y 7 C₆H₁₂O₆.

**Predicción, formulada con las etapas 1-4: D-apiosa (PAV); segunda opción, ribosa.** El acierto #1 de BLAST (RbsB, ribosa, 28.5 %) está en zona crepuscular, y el #2 (31.5 %) une treitol. La apiosa es el único ligando que aparece a la vez en Swiss-Prot (aciertos #8-#11, dos de ellos de *R. rhizogenes* y *R. etli*) y en Foldseek, y es coherente con la ecología de *Rhizobium*, asociado a plantas; la apiosa es un azúcar de la pared celular vegetal. *Declaración:* vi los scores de Vina antes de escribir la predicción, pero no los usé para formularla. **Resultado:** PAV queda 9.º en Vina y 12.º en Boltz, dentro del empate. Ningún método la confirma ni la descarta.

## 6 · Péptido señal y numeración

El corte está **después de A27**, y hay tres evidencias independientes: **(i) secuencia:** región n `MKRRT` (carga +3), región h `FLQTGSALIAAGAFGIP` (Kyte-Doolittle +1.27, con máximo de +2.24) y región c `GILRA`, que termina en Ala −1, con una caída de hidropatía en el residuo 28; **(ii) pLDDT:** 36.3 en los residuos 1-27 frente a 96.8 en 28-362, y un PAE señal-maduro de 29.6 Å; **(iii) BLAST:** ningún homólogo alinea antes del residuo 49. **Recorté el PDB (opción b)** en lugar de volver a predecir, con un costo que hay que declarar: el resto de la proteína se plegó en presencia del péptido señal. **Numeración madura 1-335.** *Nota:* al modelo de respaldo que se reparte (`TeuB_modelo_apo.pdb`) le falta una Arg en la posición 256 (`WMRRW` → `WMRW`), así que tiene 334 residuos. Por eso su Cys268 es C269 en la numeración correcta, y el centro de la caja se recalculó sobre el modelo propio: (3.38, −3.17, 0.79).

## Veredicto

**No: con lo que medí, no puedo decir qué ligando prefiere TeuB.** Los dos métodos colocan todos los ligandos en el mismo bolsillo. Pero Boltz empata a 12 de 16 y Vina ordena con una precisión que su función de puntaje no respalda. Además, ninguno de los dos distingue un azúcar de un señuelo como la hipoxantina, y sus rankings no correlacionan (ρ = 0.39). Lo único que los dos sostienen es que **un poliol acíclico (X9X) es mal ligando**. La hipótesis mejor fundamentada, la D-apiosa, sigue apoyándose en la homología y en el contexto biológico, no en el cálculo.

<small>Datos, scripts y registros: https://github.com/jtuxman/tareaarcadioazucares · `08_analisis/comparar.py`, `tabla_comparativa.tsv`.</small>
