"""Figura del entregable v2: (A) Vina, (B) Boltz-2 ipTM con banda de empate, (C) puesto vs puesto."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': 0.6, 'svg.fonttype': 'none'})
T = pd.read_csv('08_analisis/tabla_comparativa.tsv', sep='\t', index_col=0)
SEN = {'HPA', 'INS', 'X9X'}; AZUL, GRIS, ROJO, ORO = '#1f4e79', '#9aa5b1', '#b23a48', '#d4a017'
col = lambda c: ROJO if c in SEN else (ORO if c == 'PAV' else AZUL)
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.35), gridspec_kw={'width_ratios': [1, 1, 0.95]})
# A · Vina
S = T.sort_values('vina'); y = np.arange(len(S))[::-1]
ax[0].errorbar(S.vina, y, xerr=S.vina_sd, fmt='none', ecolor=GRIS, lw=0.8)
ax[0].scatter(S.vina, y, c=[col(c) for c in S.index], s=16, zorder=3)
ax[0].set_yticks(y, S.index); ax[0].set_xlabel('Vina (kcal/mol)   mejor →'); ax[0].set_title('A · Docking (3 semillas)', loc='left', fontweight='bold')
ax[0].invert_xaxis()
# B · Boltz
S = T.sort_values('iptm_m', ascending=False); y = np.arange(len(S))[::-1]
top = S.index[0]; thr = [S.iptm_m[top] - 2 * np.hypot(S.iptm_sd[top], S.iptm_sd[c]) for c in S.index]
tie = [S.iptm_m[c] > t for c, t in zip(S.index, thr)]
ax[1].axhspan(y[sum(tie) - 1] - 0.5, y[0] + 0.5, color='#e8eef5', zorder=0)
ax[1].text(S.iptm_m.min() + 0.001, y[sum(tie) // 2], f'empatados\ncon el 1.º\n({sum(tie)})', fontsize=6, color=AZUL, va='center')
ax[1].errorbar(S.iptm_m, y, xerr=S.iptm_sd, fmt='none', ecolor=GRIS, lw=0.8)
ax[1].scatter(S.iptm_m, y, c=[col(c) for c in S.index], s=16, zorder=3)
ax[1].set_yticks(y, S.index); ax[1].set_xlabel('ipTM medio (15 modelos)  mejor →'); ax[1].set_title('B · Co-plegamiento Boltz-2', loc='left', fontweight='bold')
# C · puestos
rho, p = spearmanr(T.vina, -T.iptm_m); print(f'rho={rho:.3f} p={p:.3f}')
ax[2].plot([1, 16], [1, 16], ls='--', c=GRIS, lw=0.7)
ax[2].scatter(T.p_vina, T.p_iptm, c=[col(c) for c in T.index], s=16, zorder=3)
for c in T.index:
    off = {'XYP': (-4, -7), 'ARA': (2.5, 3.5)}.get(c, (2.5, 1.5))
    ax[2].annotate(c, (T.p_vina[c], T.p_iptm[c]), xytext=off, textcoords='offset points', fontsize=5.3)
ax[2].set_xlim(0, 17.5); ax[2].set_ylim(0, 17.5); ax[2].invert_yaxis()
ax[2].set_xlabel('puesto Vina'); ax[2].set_ylabel('puesto Boltz-2')
ax[2].set_title('C · Puesto contra puesto', loc='left', fontweight='bold')
ax[2].text(0.8, 16.9, f'Spearman ρ = {rho:.2f}\np = {p:.2f}', ha='left', va='bottom', fontsize=6.5, color=AZUL)
ax[2].set_xticks([1, 4, 8, 12, 16]); ax[2].set_yticks([1, 4, 8, 12, 16])
for a in ax[:2]: a.tick_params(axis='y', length=0, labelsize=6)
fig.tight_layout(w_pad=1.2)
fig.savefig('08_analisis/figura_entregable.svg'); fig.savefig('08_analisis/figura_entregable.png', dpi=250)
