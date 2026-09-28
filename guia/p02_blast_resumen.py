"""Etapa 2 · Resume la salida de BLAST contra Swiss-Prot: tabla de los 12 primeros,
qstart minimo/mediana y cuantos ligandos distintos aparecen en el top 10.
Uso (desde la raiz del proyecto):  python guia/p02_blast_resumen.py"""
import re
rows = [l.rstrip('\n').split('\t') for l in open('01_blast/swissprot_hits.tsv')]
# columnas: sseqid stitle pident length qstart qend evalue bitscore qcovs
print(f"Aciertos con E < 1e-3: {len(rows)}\n")
print(f"{'#':<3}{'accession':<12}{'nombre':<42}{'organismo':<30}{'%id':>6}{'cov':>5}{'qstart':>7}{'E':>10}")
for i, r in enumerate(rows[:12], 1):
    acc = r[0].split('|')[1]
    org = re.search(r'\[([^\]]+)\]$', r[1]); org = org.group(1) if org else ''
    name = re.sub(r'^RecName: Full=|;.*$|\s*\[.*$', '', r[1])
    print(f"{i:<3}{acc:<12}{name[:41]:<42}{org[:29]:<30}{float(r[2]):>6.1f}{r[8]:>5}{r[4]:>7}{float(r[6]):>10.1e}")
qs = sorted(int(r[4]) for r in rows)
print(f"\nqstart minimo = {qs[0]}  mediana = {qs[len(qs)//2]}  -> ningun homologo alinea antes del residuo {qs[0]}")
ids = [float(r[2]) for r in rows]
print(f"identidad: {min(ids):.1f}-{max(ids):.1f} %  (zona crepuscular si < ~35 %)")
