"""Etapa 3 · Lee el modelo rank_001 de ColabFold: profundidad del MSA, pLDDT por region,
PAE senal-vs-madura y la mejor particion en dos dominios.
Uso: python guia/p03_plddt_pae.py [directorio_del_job]"""
import sys, glob, json, numpy as np
D = sys.argv[1] if len(sys.argv) > 1 else '02_estructura/TeuB_full_146f4'
CORTE = 27                                            # hipotesis: peptido senal = 1-27
a3m = glob.glob(f'{D}/*.a3m')[0]
print(f"MSA: {sum(l.startswith('>') for l in open(a3m))} secuencias")
pdb = sorted(glob.glob(f'{D}/*rank_001*.pdb'))[0]
bf = {}
for l in open(pdb):                                   # en ColabFold el pLDDT va en la columna B-factor
    if l.startswith('ATOM'): bf.setdefault(int(l[22:26]), []).append(float(l[60:66]))
p = {k: np.mean(v) for k, v in bf.items()}; N = max(p)
print(f"modelo: {pdb.split('/')[-1]}")
print(f"pLDDT medio {np.mean(list(p.values())):.1f} ; minimo {min(p.values()):.1f} en el residuo {min(p, key=p.get)}")
s = np.mean([p[i] for i in range(1, CORTE + 1)]); m = np.mean([p[i] for i in range(CORTE + 1, N + 1)])
print(f"pLDDT 1-{CORTE} = {s:.1f}   {CORTE+1}-{N} = {m:.1f}   diferencia = {m - s:.1f}")
print("\nperfil del extremo N:")
for i in range(1, 46): print(f"  {i:3d} {p[i]:5.1f} {'#' * int(p[i] / 3)}")
sc = glob.glob(f'{D}/*scores_rank_001*.json')[0]
pae = np.array(json.load(open(sc))['pae'])
S, M = slice(0, CORTE), slice(CORTE, N)
print(f"\nPAE senal-vs-madura = {pae[S, M].mean():.1f} A ; dentro de la madura = {pae[M, M].mean():.1f} A")
X = pae[CORTE:, CORTE:]; n = len(X); best = None
for c in range(40, n - 40):                            # corte que maximiza PAE inter - intra
    intra = (X[:c, :c].mean() + X[c:, c:].mean()) / 2; inter = (X[:c, c:].mean() + X[c:, :c].mean()) / 2
    if best is None or inter - intra > best[0]: best = (inter - intra, c)
print(f"mejor particion en 2 dominios: corte en el residuo {best[1] + CORTE} (precursor), contraste {best[0]:.2f} A")
print("contraste < 2 A -> el modelo sabe orientar los dos lobulos (hendidura cerrada)")
