# -*- coding: utf-8 -*-
"""Genera datos SINTÉTICOS de 1.100 'municipios' para la simulación (no son datos reales de Colombia).
Salida: datos.js (se incrusta luego en el HTML).  Uso: python generar_datos.py
"""
import json, numpy as np
from scipy.cluster.hierarchy import linkage

rng = np.random.default_rng(2026)
N = 1100

# (Versión anterior, solo para conservar el flujo del generador aleatorio y que las variables no cambien)
POLY_OLD = [(-71.3,12.45),(-71.9,11.7),(-72.9,11.4),(-74.2,11.25),(-74.85,11.0),(-75.55,10.4),(-75.7,9.5),
        (-76.4,8.7),(-76.9,8.1),(-77.35,8.65),(-77.9,7.2),(-77.5,6.0),(-77.3,4.9),(-77.1,3.9),(-78.4,2.6),
        (-78.8,1.7),(-77.7,0.8),(-76.3,0.4),(-75.3,-0.1),(-74.4,-0.9),(-73.2,-1.8),(-70.0,-4.2),(-69.4,-1.1),
        (-69.8,1.1),(-67.2,1.8),(-67.3,3.2),(-67.8,4.5),(-67.5,6.2),(-69.0,6.2),(-70.1,7.0),(-72.0,7.0),
        (-72.5,8.0),(-72.8,9.2),(-72.2,10.8),(-71.5,11.9)]

# Contorno de Colombia: silueta del archivo 'mapa colombia.jpg' (aportado por el usuario), vectorizada
# con trazado de contorno + simplificación y ajustada linealmente a la caja lon/lat aproximada del país.
# Sigue siendo un esquema (no cartografía oficial). Ver poly_colombia.json.
POLY = [tuple(p) for p in json.load(open('poly_colombia.json'))]

def in_poly(x, y, POLY=POLY):
    c = False; j = len(POLY) - 1
    for i in range(len(POLY)):
        xi, yi = POLY[i]; xj, yj = POLY[j]
        if (yi > y) != (yj > y) and x < (xj - xi) * (y - yi) / (yj - yi) + xi: c = not c
        j = i
    return c

# Variables: [%VA agro, %área agro, VA agro pc (mill COP), pob rural (hab), rendimiento (t/ha), conectividad (0-100), NBI %]
VARS = ["pct_va_agro","pct_area_agro","va_agro_pc","pob_rural","rendimiento","conectividad","nbi"]
# perfiles latentes (medias) y dispersión
P = {
 0: dict(m=[58,62,9.0,14000,14,66,30], s=[8,9,2.2,0.55,2.5,10,7]),   # vocación agrícola productiva
 1: dict(m=[62,46,2.4,9000,5,26,62],  s=[9,10,0.8,0.60,1.4,9,8]),    # agrícola baja productividad
 2: dict(m=[28,36,4.0,12000,9,56,38], s=[8,9,1.3,0.60,2.0,10,7]),    # economía mixta / transición
 3: dict(m=[6,10,1.6,25000,8,86,20],  s=[3,5,0.7,0.80,2.0,7,5]),     # diversificada, baja participación agro
}
# zonas (centros lon/lat) preferidas por perfil, solo para ubicar puntos de forma verosímil
ZONES = {
 0: [(-75.5,5.5),(-73.5,4.2),(-74.5,9.0),(-72.8,7.6)],
 1: [(-77.0,5.5),(-76.0,2.0),(-72.5,1.5),(-71.5,10.7),(-69.5,3.5)],
 2: [(-75.0,7.2),(-74.0,5.0),(-73.5,9.5),(-76.0,3.5)],
 3: [(-74.1,4.7),(-75.6,6.2),(-76.5,3.4),(-74.8,10.9),(-73.1,7.1)],
}
props = [.26,.30,.27,.17]
rows, lonlat, truth, blend = [], [], [], []
for i in range(N):
    k = rng.choice(4, p=props)
    w = np.zeros(4); w[k] = 1.0
    if rng.random() < 0.14:                       # ~14%: municipios "de frontera" (mezcla de dos perfiles)
        k2 = rng.choice([j for j in range(4) if j != k]); a = rng.uniform(.35,.55)
        w = np.zeros(4); w[k] = 1-a; w[k2] = a
    m = sum(w[j]*np.array(P[j]['m'],float) for j in range(4))
    s = sum(w[j]*np.array(P[j]['s'],float) for j in range(4))
    v = np.empty(7)
    for j in range(7):
        if VARS[j] == "pob_rural": v[j] = np.exp(rng.normal(np.log(m[j]), s[j]))   # lognormal (asimétrica)
        else: v[j] = rng.normal(m[j], s[j])
    v[0] = np.clip(v[0],0.5,90); v[1] = np.clip(v[1],1,95); v[2] = max(v[2],0.2)
    v[4] = max(v[4],1); v[5] = np.clip(v[5],5,100); v[6] = np.clip(v[6],3,95)
    # ubicación
    kk = k if w[k] >= .5 else k
    while True:
        zc = ZONES[kk][rng.integers(len(ZONES[kk]))]
        x = rng.normal(zc[0],0.9); y = rng.normal(zc[1],0.9)
        if in_poly(x,y,POLY_OLD): break
    rows.append(v); lonlat.append((x,y)); truth.append(int(k)); blend.append(round(float(w.max()),2))
# Ubicación definitiva sobre el nuevo contorno (generador aparte: NO altera las variables ya simuladas)
rng2 = np.random.default_rng(7)
lonlat = []
for k in truth:
    while True:
        zc = ZONES[k][rng2.integers(len(ZONES[k]))]
        x = rng2.normal(zc[0],0.9); y = rng2.normal(zc[1],0.9)
        if in_poly(x,y): break
    lonlat.append((x,y))
X = np.array(rows)
Z = (X - X.mean(0)) / X.std(0, ddof=1)

# PCA (sobre datos estandarizados) para la proyección 2D
C = np.cov(Z.T); ev, evec = np.linalg.eigh(C); o = np.argsort(ev)[::-1]; ev, evec = ev[o], evec[:,o]
scores = Z @ evec[:,:2]
# fijar signos para estabilidad
for c in range(2):
    if evec[np.argmax(abs(evec[:,c])), c] < 0: evec[:,c] *= -1
scores = Z @ evec[:,:2]

# Ward sobre Z: alturas de fusión al pasar a k grupos
L = linkage(Z, 'ward'); h = L[:,2]
fus = {k: float(h[N-1-k]) for k in range(1,13)}
mx = max(fus.values())
fusn = {k: round(v/mx,4) for k,v in fus.items()}

out = dict(
  vars=VARS, N=N,
  raw=[[round(float(a),3) for a in r] for r in X],
  lon=[round(float(a[0]),3) for a in lonlat], lat=[round(float(a[1]),3) for a in lonlat],
  pc=[[round(float(a),3) for a in r] for r in scores],
  evr=[round(float(a),4) for a in (ev/ev.sum())[:3]], load=[[round(float(a),3) for a in r] for r in evec[:,:2]],
  fusion=fusn, poly=POLY, blendmax=blend)
open("datos.js","w",encoding="utf-8").write("const DATA="+json.dumps(out,separators=(',',':'))+";")
import csv
with open("municipios_sinteticos.csv","w",newline="",encoding="utf-8") as fh:
    w = csv.writer(fh); w.writerow(["id"]+VARS+["lon","lat"])
    for i in range(N): w.writerow([f"M{i+1:04d}"]+[round(float(a),3) for a in X[i]]+[round(lonlat[i][0],3),round(lonlat[i][1],3)])
with open("poly_colombia.csv","w",newline="",encoding="utf-8") as fh:
    w = csv.writer(fh); w.writerow(["lon","lat"]); [w.writerow(p) for p in POLY]
print("N",N,"var expl PC1-2", out['evr'][:2]); print("fusion norm", fusn)
print("corr matriz:\n", np.round(np.corrcoef(X.T),2))

