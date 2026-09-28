"""Tabla comparativa Vina vs Boltz-2, grupos de empate por piso de ruido,
correlacion de rangos y acuerdo de poses (Boltz superpuesto sobre el receptor de Vina)."""
import glob, re, io, zipfile, numpy as np, pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
from scipy.stats import spearmanr, kendalltau
from scipy.optimize import linear_sum_assignment
R = '/home/jmanuel/cursos/proteinas/tareasugar'
nombre = dict(BGC='β-D-glucopiranosa', RIP='β-D-ribopiranosa', GAL='β-D-galactopiranosa', XYP='β-D-xilopiranosa',
    PAV='D-apiosa', INS='mio-inositol (señuelo)', ARA='α-L-arabinopiranosa', AHR='α-L-arabinofuranosa',
    BDR='β-D-ribofuranosa', ALL='β-D-alopiranosa', WEB='α-D-psicopiranosa', FRU='β-D-fructofuranosa',
    GZL='β-D-galactofuranosa', X9X='D-alitol (señuelo)', HPA='hipoxantina (señuelo)', **{'3VB': 'D-treitol'})

# --- Vina
v = pd.read_csv(f'{R}/08_analisis/vina_poses.tsv', sep='\t')
V = v.groupby('cod').score.agg(['mean', 'std', 'min', 'max', 'count'])
# --- Boltz
b = pd.read_csv(f'{R}/08_analisis/boltz_resumen.tsv', sep='\t').set_index('cod')
cods = sorted(set(V.index) & set(b.index))
T = pd.DataFrame(index=cods)
T['nombre'] = [nombre[c] for c in cods]
T['vina'] = V.loc[cods, 'mean']; T['vina_sd'] = V.loc[cods, 'std']
T['boltz_rs'] = b.loc[cods, 'rs_media']; T['boltz_sd'] = b.loc[cods, 'rs_sd']; T['iptm'] = b.loc[cods, 'iptm']
T['p_vina'] = T.vina.rank(method='min').astype(int)            # mas negativo = mejor
T['p_boltz'] = (-T.boltz_rs).rank(method='min').astype(int)    # mas alto = mejor
T['delta'] = T.p_boltz - T.p_vina

def grupos(s, sd, asc):
    """Agrupa en empates: consecutivos cuya diferencia < 2*sqrt(sd1^2+sd2^2)."""
    o = s.sort_values(ascending=asc).index; g, k = {}, 1
    for i, c in enumerate(o):
        if i and abs(s[c] - s[o[i-1]]) >= 2 * np.hypot(sd[c], sd[o[i-1]]): k += 1
        g[c] = k
    return pd.Series(g)
T['g_vina'] = grupos(T.vina, T.vina_sd, True)
T['g_boltz'] = grupos(T.boltz_rs, T.boltz_sd, False)

# --- acuerdo de poses: superponer el mejor modelo Boltz sobre el receptor y comparar con modo 1 de Vina (semilla 101)
rec = list(PDBParser(QUIET=True).get_structure('r', f'{R}/03_maduro/TeuB_maduro_1-335.pdb')[0])[0]
site = [38, 43, 45, 174, 197, 222, 248, 269]
def heavy_vina(f):
    t = open(f).read().split('ENDMDL')[0]
    return [(l[77:79].strip()[0], np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
            for l in t.splitlines() if l.startswith(('ATOM', 'HETATM')) and l[77:79].strip() not in ('H', 'HD')]
# Nota: Boltz-2 elimina el atomo O1 (10 de 16 ligandos); la asignacion rectangular empareja
# los atomos que si existen en ambos y deja fuera el sobrante de Vina.
rms, dcen = {}, {}
import json
ZIPS = {}
for f in glob.glob(f'{R}/07_coplegamiento/teub_*.result.zip'):
    zf = zipfile.ZipFile(f); d = json.load(zf.open([n for n in zf.namelist() if n.endswith('_data.json')][0]))
    ZIPS[[s['ligand']['ccdCodes'][0] for s in d['sequences'] if 'ligand' in s][0]] = zf
for c in cods:
    zf = ZIPS[c]
    top = [n for n in zf.namelist() if n.count('/') == 2 and n.endswith('_model.cif')][0]
    st = MMCIFParser(QUIET=True).get_structure('b', io.StringIO(zf.read(top).decode()))[0]
    fixed = [rec[i]['CA'] for i in range(1, 336)]; moving = [st['A'][i]['CA'] for i in range(1, 336)]
    sup = Superimposer(); sup.set_atoms(fixed, moving)
    rot, tran = sup.rotran
    BL = [(a.element[0], a.coord @ rot + tran) for a in st['B'].get_atoms() if a.element != 'H']
    VL = heavy_vina(f'{R}/06_docking/out_{c}_s101.pdbqt')
    dcen[c] = np.linalg.norm(np.mean([x for _, x in BL], 0) - np.mean([x for _, x in VL], 0))
    # RMSD con asignacion optima por elemento (robusto a nombres de atomo distintos)
    sq = []
    for el in set(e for e, _ in VL):
        P = np.array([x for e, x in VL if e == el]); Q = np.array([x for e, x in BL if e == el])
        if len(Q) == 0: continue
        D = np.linalg.norm(P[:, None] - Q[None], axis=2) ** 2; i, j = linear_sum_assignment(D); sq += list(D[i, j])
    rms[c] = np.sqrt(np.mean(sq))
T['pose_rmsd'] = pd.Series(rms); T['pose_dcen'] = pd.Series(dcen)

T = T.sort_values('p_vina')
pd.set_option('display.width', 200)
print(T[['nombre', 'p_vina', 'vina', 'vina_sd', 'g_vina', 'p_boltz', 'boltz_rs', 'boltz_sd', 'iptm', 'g_boltz', 'delta', 'pose_rmsd', 'pose_dcen']].round(4).to_string())
rho, p = spearmanr(T.vina, -T.boltz_rs); tau, pt = kendalltau(T.vina, -T.boltz_rs)
print(f"\nN={len(T)}  Spearman rho={rho:.2f} (p={p:.2f})  Kendall tau={tau:.2f} (p={pt:.2f})")
print(f"Rango Vina: {T.vina.max()-T.vina.min():.2f} kcal/mol ; ruido max entre semillas (rango): {(V['max']-V['min']).max():.3f}")
print(f"Rango Boltz: {T.boltz_rs.max()-T.boltz_rs.min():.4f} ; sd media entre 15 modelos: {T.boltz_sd.mean():.4f}")
print(f"Grupos de empate: Vina {T.g_vina.max()} , Boltz {T.g_boltz.max()}")
T.round(4).to_csv(f'{R}/08_analisis/tabla_comparativa.tsv', sep='\t')
