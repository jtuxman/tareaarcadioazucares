"""Etapa 2 · Resume BLAST contra nr: identidades y organismos.
Uso: python guia/p02_nr_resumen.py"""
import re, collections
rows = [l.rstrip('\n').split('\t') for l in open('01_blast/nr_hits.tsv')]
ids = [float(r[2]) for r in rows]
print(f"aciertos: {len(rows)}   >40 % id: {sum(i > 40 for i in ids)}   >60 %: {sum(i > 60 for i in ids)}")
RHIZ = ('rhizob', 'agrobacterium', 'sinorhizobium', 'ensifer', 'mesorhizobium', 'martinezella', 'neorhizobium', 'rhizobiaceae')
org = [re.search(r'\[([^\]]+)\]', r[1]).group(1) if '[' in r[1] else '?' for r in rows]
print(f"Rhizobiaceae y parientes: {sum(any(k in o.lower() for k in RHIZ) for o in org)}")
for o, n in collections.Counter(org).most_common(10): print(f"  {n:4d}  {o}")
print("\ntop 5:")
for r in rows[:5]: print(f"  {float(r[2]):5.1f} %  {r[1][:100]}")
