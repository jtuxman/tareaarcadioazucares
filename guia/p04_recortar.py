"""Etapa 4 · Recorta el peptido senal del PDB (opcion b), renumera 1-335, escribe el FASTA
maduro y verifica los 8 residuos del sitio. Compara ademas con el modelo de respaldo.
Uso: python guia/p04_recortar.py"""
src = '02_estructura/TeuB_full_146f4/TeuB_full_146f4_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb'
out = '03_maduro/TeuB_maduro_1-335.pdb'; CORTE = 27
AA = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', ILE='I',
          LEU='L', LYS='K', MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V')
with open(out, 'w') as o:
    for l in open(src):
        if l.startswith(('ATOM', 'TER')) and int(l[22:26]) > CORTE:
            o.write(l[:22] + f"{int(l[22:26]) - CORTE:4d}" + l[26:])       # columnas 23-26 = numero de residuo
    o.write('END\n')
seq = open('01_blast/teuB_solo_letras.txt').read().strip()[CORTE:]
open('03_maduro/teuB_maduro_solo_letras.txt', 'w').write(seq + '\n')
open('03_maduro/teuB_maduro.fasta', 'w').write('>TeuB_maduro_1-335\n' + '\n'.join(seq[i:i+60] for i in range(0, len(seq), 60)) + '\n')
print(f"madura: {len(seq)} aa, empieza {seq[:5]}")
def leer(p):
    return {int(l[22:26]): AA[l[17:20]] for l in open(p) if l.startswith('ATOM')}
mio = leer(out)
print("sitio (numeracion madura):", ' '.join(f"{mio[n]}{n}" for n in (38, 43, 45, 174, 197, 222, 248, 269)))
resp = leer('00_datos_repartidos/TeuB_modelo_apo.pdb')
a = ''.join(resp[i] for i in sorted(resp)); b = ''.join(mio[i] for i in sorted(mio))
print(f"respaldo {len(a)} residuos, propio {len(b)}")
i = next((k for k in range(min(len(a), len(b))) if a[k] != b[k]), None)
if i is not None: print(f"primera diferencia en {i+1}: respaldo ...{a[i-6:i+6]}...  propio ...{b[i-6:i+6]}...")
