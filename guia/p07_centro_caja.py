"""Etapa 7 · Centro de la caja de docking = centroide de los atomos de cadena lateral de los
8 residuos del sitio, en las coordenadas del modelo PROPIO. Uso: python guia/p07_centro_caja.py"""
import numpy as np
SITIO = [38, 43, 45, 174, 197, 222, 248, 269]
P = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in open('03_maduro/TeuB_maduro_1-335.pdb')
              if l.startswith('ATOM') and int(l[22:26]) in SITIO and l[12:16].strip() not in ('N', 'CA', 'C', 'O')])
c = P.mean(0); ext = P.max(0) - P.min(0)
print(f"{len(P)} atomos de cadena lateral")
print(f"centro: x = {c[0]:.2f}  y = {c[1]:.2f}  z = {c[2]:.2f}")
print(f"extension del sitio: {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} A -> caja de {max(22, np.ceil(ext.max()) + 8):.0f} A por lado")
