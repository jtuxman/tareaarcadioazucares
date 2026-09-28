#!/bin/bash
cd /home/jmanuel/cursos/proteinas/tareasugar/06_docking
CX=3.38; CY=-3.17; CZ=0.79; S=26
for c in BGC RIP GAL XYP PAV INS ARA AHR BDR ALL WEB FRU GZL X9X HPA; do
  for seed in 101 202 303; do
    o=out_${c}_s${seed}.pdbqt
    [ -s "$o" ] && continue
    vina --receptor receptor.pdbqt --ligand ../05_ligandos/${c}.pdbqt \
      --center_x $CX --center_y $CY --center_z $CZ --size_x $S --size_y $S --size_z $S \
      --exhaustiveness 32 --num_modes 9 --seed $seed --cpu 8 --out $o > log_${c}_s${seed}.txt 2>&1
    echo "$c seed=$seed $(grep -A3 '^-----' log_${c}_s${seed}.txt | sed -n '2p' | awk '{print $2}')"
  done
done
echo TODO_LISTO > .done
