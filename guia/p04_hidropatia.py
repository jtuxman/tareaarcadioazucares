"""Etapa 4 · Evidencia de secuencia del peptido senal: regiones n/h/c y perfil de
hidropatia de Kyte-Doolittle (ventana de 9).  Uso: python guia/p04_hidropatia.py"""
seq = open('01_blast/teuB_solo_letras.txt').read().strip()
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
carga = lambda s: sum(c in 'KR' for c in s) - sum(c in 'DE' for c in s)
hid = lambda s: sum(KD[c] for c in s) / len(s)
print(f"n (1-5)   {seq[0:5]:<18} carga neta {carga(seq[0:5]):+d}")
print(f"h (6-22)  {seq[5:22]:<18} KD medio {hid(seq[5:22]):+.2f}")
print(f"c (23-27) {seq[22:27]:<18} corte despues del residuo 27")
print(f"madura: empieza {seq[27:32]}, {len(seq) - 27} aa\n")
W = 9
for i in range(0, 60 - W + 1):
    v = hid(seq[i:i + W]); print(f"{i+1:3d}-{i+W:<3d} {seq[i:i+W]} {v:+5.2f} {'#' * max(int(round((v + 4.5) * 4)), 0)}")
