"""Etapa 6 · Las tres verificaciones del guion para cada ligando: la formula del SDF coincide
con la ficha del CCD, tiene hidrogenos y no es plano. Cuenta ademas las torsiones activas.
Uso: python guia/p06_verificar_ligandos.py BGC RIP GAL ...   (sin argumentos: el panel completo)"""
import sys, re, json, collections
import numpy as np
PANEL = sys.argv[1:] or 'BGC RIP GAL XYP PAV INS ARA AHR BDR ALL WEB FRU GZL X9X HPA 3VB'.split()
info = json.load(open('04_foldseek/tabla_ligandos.json'))['info']
norm = lambda s: sorted(re.findall(r'[A-Z][a-z]?\d*', s))
print(f"{'cod':<5}{'ficha':<12}{'SDF':<12}{'nH':>4}{'tors':>6}{'grosor':>8}  veredicto")
for c in PANEL:
    L = open(f'05_ligandos/{c}_ideal.sdf').read().splitlines()
    na = int(L[3][:3]); at = L[4:4 + na]
    el = collections.Counter(l[31:34].strip() for l in at)
    xyz = [[float(l[0:10]), float(l[10:20]), float(l[20:30])] for l in at]
    es2d = all(abs(p[2]) < 1e-4 for p in xyz)          # un SDF "2D" trae todas las z = 0
    X = np.array(xyz) - np.mean(xyz, 0)
    grosor = 2 * np.sqrt(np.linalg.eigvalsh(X.T @ X / len(X))[0])   # espesor en el eje mas corto
    f = ''.join(f"{e}{el[e] if el[e] > 1 else ''}" for e in ('C', 'H', 'N', 'O') if el[e])
    ficha = (info.get(c, {}).get('formula') or '?').replace(' ', '')
    tors = open(f'05_ligandos/{c}.pdbqt').read().count('BRANCH') // 2
    ok = [('formula OK' if norm(f) == norm(ficha) else f'FORMULA? ficha={ficha}'),
          ('tiene H' if el['H'] else 'SIN H'), ('SDF 2D!' if es2d else ('3D' if grosor > 0.3 else 'plano (aromatico: esperado)'))]
    print(f"{c:<5}{ficha:<12}{f:<12}{el['H']:>4}{tors:>6}{grosor:>7.2f}A  {' | '.join(ok)}")
