"""Etapa 5 · Toma los aciertos de Foldseek, se queda con las primeras 150 entradas
unicas del PDB, pregunta al RCSB (API GraphQL) que ligandos tiene cada una, descarta
iones/crioprotectores/tampones y cuenta. Imprime la tabla de ligandos y las formulas repetidas.
Uso: python guia/p05_rcsb_ligandos.py [salida.json]"""
import sys, json, urllib.request, collections
OUT = sys.argv[1] if len(sys.argv) > 1 else '04_foldseek/tabla_ligandos.json'
rows = [l.rstrip('\n').split('\t') for l in open('04_foldseek/hits.tsv')]
ids = []
for r in rows:                                   # target = '2ioy-assembly2_B' -> '2IOY'
    pid = r[1].split('-')[0].split('_')[0].upper()
    if pid not in ids: ids.append(pid)
    if len(ids) == 150: break
Q = '''{ entries(entry_ids: %s) { rcsb_id
  nonpolymer_entities { nonpolymer_comp { chem_comp { id name formula formula_weight } } } } }'''
E = {}
for i in range(0, len(ids), 50):
    req = urllib.request.Request('https://data.rcsb.org/graphql', headers={'Content-Type': 'application/json'},
                                 data=json.dumps({'query': Q % json.dumps(ids[i:i+50])}).encode())
    for e in json.load(urllib.request.urlopen(req, timeout=90))['data']['entries']:
        if e: E[e['rcsb_id']] = [ne['nonpolymer_comp']['chem_comp'] for ne in (e['nonpolymer_entities'] or [])]
# lo que NO es un ligando biologico: agua, iones, crioprotectores, tampones, detergentes, aminoacidos sueltos
BASURA = set('''HOH DOD SO4 PO4 GOL EDO PEG PG4 PGE 1PE 2PE MPD ACT ACY FMT CL BR IOD NA K MG CA ZN MN FE FE2 NI CO
CU CU1 CD HG CS RB SR BA LI AL F NO3 TRS EPE MES BTB CIT FLC TAR MLA SCN AZI NH4 IMD BME DTT DMS DMF URE ACE NAG BOG
LDA C8E TRT P6G 12P 15P XPE SIN OXL MAE EOH IPA ETX BU3 PDO MRD SGM TLA OCT HEZ PE4 PEU B3P CAC NH2 UNX UNL'''.split())
cnt, info, donde, apo = collections.Counter(), {}, collections.defaultdict(list), []
for pid, ligs in E.items():
    buenos = [l for l in ligs if l['id'] not in BASURA]
    if not buenos: apo.append(pid)
    for l in buenos: cnt[l['id']] += 1; info[l['id']] = l; donde[l['id']].append(pid)
print(f"entradas: {len(E)}   apo (sin ligando util): {len(apo)} = {100*len(apo)/len(E):.0f} %\n")
for c, n in cnt.most_common(30):
    i = info[c]; print(f"{c:<5}{n:>3}  {(i['formula'] or '?').replace(' ', ''):<14}{(i['name'] or '')[:45]:<46} {','.join(donde[c][:5])}")
porf = collections.defaultdict(list)
for c in cnt: porf[(info[c]['formula'] or '?').replace(' ', '')].append(c)
print("\nFORMULAS REPETIDAS:")
for f, cs in sorted(porf.items(), key=lambda x: -len(x[1])):
    if len(cs) > 1: print(f"  {f:<12} {len(cs)} ligandos: {' '.join(cs)}")
json.dump({'cnt': cnt, 'info': info, 'donde': donde, 'apo': apo}, open(OUT, 'w'), indent=1)
