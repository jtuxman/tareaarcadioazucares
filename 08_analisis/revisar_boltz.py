"""Revisa todos los ZIP de Boltz-2 (teub_*.result.zip), valida los parametros,
copia los validos a 07_coplegamiento/ y mide si el ligando cae en el sitio."""
import zipfile, json, glob, csv, io, os, shutil
import numpy as np
from Bio.PDB import MMCIFParser
R = '/home/jmanuel/cursos/proteinas/tareasugar'
mat = open(f'{R}/03_maduro/teuB_maduro_solo_letras.txt').read().strip()
panel = "BGC RIP GAL XYP PAV INS ARA AHR BDR ALL WEB FRU GZL X9X HPA 3VB".split()
site = [38, 43, 45, 174, 197, 222, 248, 269]
files = glob.glob('/home/jmanuel/Downloads/teub_*.result.zip') + \
        glob.glob('/home/jmanuel/cursos/proteinas/**/teub_*.result.zip', recursive=True)
seen, ok, rows = set(), {}, []
for z in sorted(files, key=os.path.getmtime):
    b = os.path.basename(z)
    if b in seen: continue
    seen.add(b)
    zf = zipfile.ZipFile(z); names = zf.namelist()
    d = json.load(zf.open([n for n in names if n.endswith('_data.json')][0]))
    prot = [s['protein'] for s in d['sequences'] if 'protein' in s][0]
    lig = [c for s in d['sequences'] if 'ligand' in s for c in (s['ligand'].get('ccdCodes') or ['SMILES'])]
    s = json.load(zf.open([n for n in names if n.endswith('_summary_confidences.json') and n.count('/') == 2][0]))
    r = [float(x['ranking_score']) for x in csv.DictReader(io.TextIOWrapper(zf.open([n for n in names if n.endswith('ranking_scores.csv')][0])))]
    good = prot['sequence'] == mat and not prot.get('templates') and d['modelSeeds'] == [1, 2, 3] \
           and len(r) == 15 and len(lig) == 1 and lig[0] in panel
    why = [] if good else [t for t, c in [(f"{len(prot['sequence'])}aa", prot['sequence'] != mat),
           ('plantillas', bool(prot.get('templates'))), (f"codigo {lig}", not (len(lig) == 1 and lig[0] in panel)),
           ('semillas', d['modelSeeds'] != [1, 2, 3])] if c]
    print(f"{'OK ' if good else 'MAL'} {b:30s} {lig[0] if lig else '?':4s} iptm={s['iptm']:.2f} clash={s['has_clash']:.0f} "
          f"n={len(r)} rango={min(r):.3f}-{max(r):.3f} {' '.join(why)}")
    if not good: continue
    ok[lig[0]] = z
    dst = f'{R}/07_coplegamiento/{b}'
    if not os.path.exists(dst): shutil.copy2(z, dst)
    dist, cont = [], []
    for n in sorted(x for x in names if x.endswith('_model.cif') and 'seed-' in x):
        st = MMCIFParser(QUIET=True).get_structure('x', io.StringIO(zf.read(n).decode()))
        A, B = st[0]['A'], st[0]['B']
        L = np.array([a.coord for a in B.get_atoms()])
        sc = np.array([a.coord for q in site for a in A[q] if a.get_id() not in ('N', 'C', 'O', 'CA')])
        dist.append(np.linalg.norm(L.mean(0) - sc.mean(0)))
        cont.append(sum(np.linalg.norm(L[:, None] - np.array([a.coord for a in A[q]])[None], axis=2).min() < 4 for q in site))
    rows.append((lig[0], s['iptm'], np.mean(r), np.std(r), min(r), max(r), max(dist), min(cont)))
print("\nlig  iptm  rs_media  rs_sd   rs_min rs_max  dist_max  contactos_min")
for x in sorted(rows, key=lambda x: -x[2]):
    print(f"{x[0]:4s} {x[1]:.2f}  {x[2]:.4f}  {x[3]:.4f}  {x[4]:.3f}  {x[5]:.3f}  {x[6]:.1f}A     {x[7]}/8")
with open(f'{R}/08_analisis/boltz_resumen.tsv', 'w') as f:
    f.write("cod\tiptm\trs_media\trs_sd\trs_min\trs_max\tdist_max_A\tcontactos_min\n")
    for x in rows: f.write("\t".join(map(str, x)) + "\n")
print("\nVALIDOS:", len(ok), sorted(ok))
print("FALTAN (panel 15):", [p for p in panel[:15] if p not in ok])
