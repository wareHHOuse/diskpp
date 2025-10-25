import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# --- Répertoire de travail ---
try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
except:
    pass

# --- Paramètres ---
results = {k: {'h': [], 'l2': [], 'dg': []} for k in range(1, 4)}  # k = 1,2,3

# Parcourir tous les fichiers implicit_l_*.txt
for filename in os.listdir('.'):
    if not filename.startswith("implicit_l_") or not filename.endswith(".txt"):
        continue

    # Extraire l et k du nom de fichier
    m = re.match(r"implicit_l_(\d+)_n_\d+_k_(\d+)_s_\d+\.txt", filename)
    if not m:
        continue
    l_number = int(m.group(1))
    k = int(m.group(2))
    if k not in results:
        continue

    h = l2_error = dg_error = None

    with open(filename, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            # h
            if "Characteristic h size" in line:
                match_h = re.search(r"Characteristic h size = ([0-9.eE+-]+)", line)
                if match_h:
                    h = float(match_h.group(1))

            # L2 error acoustique
            if line.strip() == "Acoustic region :":
                if i+1 < len(lines):
                    match_l2 = re.search(r"L2-norm error = ([0-9.eE+-]+)", lines[i+1])
                    if match_l2:
                        l2_error = float(match_l2.group(1))

            # dG error acoustique
            if line.strip().startswith("L2 errors of dG unknowns:"):
                for j in range(i+1, min(i+5, len(lines))):
                    if lines[j].strip().startswith("Acoustic region :"):
                        try:
                            dg_error = float(lines[j].split(":")[1].strip())
                        except:
                            dg_error = None
                        break

    if h is not None and l2_error is not None and dg_error is not None:
        results[k]['h'].append(h)
        results[k]['l2'].append(l2_error)
        results[k]['dg'].append(dg_error)
    else:
        print(f"⚠️ Données manquantes dans {filename}")

# --- Tracé ---
plt.figure(figsize=(8,8))
styles = {1:{'color':'darkgreen','marker':'^'}, 2:{'color':'orange','marker':'s'}, 3:{'color':'darkred','marker':'o'}}
for k in results:
    if results[k]['h']:
        h_sorted, l2_sorted, dg_sorted = zip(*sorted(zip(results[k]['h'], results[k]['l2'], results[k]['dg'])))
        plt.loglog(h_sorted, l2_sorted, marker=styles[k]['marker'], linestyle='-', color=styles[k]['color'], label=f'k={k} (L2)')
        plt.loglog(h_sorted, dg_sorted, marker=styles[k]['marker'], linestyle='--', color=styles[k]['color'], alpha=0.7, label=f'k={k} (dG)')

plt.xlabel("h")
plt.ylabel("L2 error / dG error")
plt.legend()
plt.grid(True, which='both')
plt.show()
