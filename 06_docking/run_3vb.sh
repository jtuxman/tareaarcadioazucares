#!/bin/bash
cd /home/jmanuel/cursos/proteinas/tareasugar/06_docking
for seed in 101 202 303; do
  o=out_3VB_s${seed}.pdbqt; [ -s "$o" ] && continue
  vina --receptor receptor.pdbqt --ligand ../05_ligandos/3VB.pdbqt \
    --center_x 3.38 --center_y -3.17 --center_z 0.79 --size_x 26 --size_y 26 --size_z 26 \
    --exhaustiveness 32 --num_modes 9 --seed $seed --cpu 8 --out $o > log_3VB_s${seed}.txt 2>&1
done
