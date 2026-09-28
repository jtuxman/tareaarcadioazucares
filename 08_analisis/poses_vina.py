"""Para la mejor pose (modo 1) de cada ligando/semilla de Vina: distancia del centroide
al centro del sitio y numero de residuos del sitio a <4 A. Marca las poses superficiales."""
import glob, re, numpy as np
from Bio.PDB import PDBParser
R = '/home/jmanuel/cursos/proteinas/tareasugar'
site = [38, 43, 45, 174, 197, 222, 248, 269]
A = PDBParser(QUIET=True).get_structure('r', f'{R}/03_maduro/TeuB_maduro_1-335.pdb')[0]
A = list(A)[0]
sc = np.array([a.coord for q in site for a in A[q] if a.get_id() not in ('N', 'C', 'O', 'CA')])
cen = sc.mean(0)
res = {q: np.array([a.coord for a in A[q]]) for q in site}
out = []
for f in sorted(glob.glob(f'{R}/06_docking/out_*_s*.pdbqt')):
    cod, seed = re.search(r'out_(\w+)_s(\d+)', f).groups()
    txt = open(f).read().split('ENDMDL')[0]
    sc_v = float(re.search(r'VINA RESULT:\s+(-?[\d.]+)', txt).group(1))
    L = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in txt.splitlines() if l.startswith(('ATOM', 'HETATM'))])
    d = np.linalg.norm(L.mean(0) - cen)
    c = [q for q in site if np.linalg.norm(L[:, None] - res[q][None], axis=2).min() < 4]
    out.append((cod, seed, sc_v, d, len(c)))
    print(f"{cod:4s} s{seed} {sc_v:7.3f}  dist={d:5.1f}A  sitio={len(c)}/8 {'SUPERFICIAL' if d > 6 or len(c) < 3 else ''}")
with open(f'{R}/08_analisis/vina_poses.tsv', 'w') as fh:
    fh.write('cod\tseed\tscore\tdist_A\tcontactos\n')
    for o in out: fh.write('\t'.join(map(str, o)) + '\n')
