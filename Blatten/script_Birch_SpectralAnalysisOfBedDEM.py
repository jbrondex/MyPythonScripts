import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.path import Path as mpltPath
from scipy.ndimage import map_coordinates
from scipy.signal import detrend, welch

def amplitude_bosses(lam, f, P, bande=0.1):
    """RMS des bosses de longueur d'onde ~ lam (m), dans la bande lam ± bande*lam."""
    f0 = 1 / lam
    masque = (f >= f0 * (1 - bande)) & (f <= f0 * (1 + bande))
    df = f[1] - f[0]                       # f est uniformément espacée
    sigma = np.sqrt(np.sum(P[masque] * df))
    return sigma

# --- 1. Lecture du DEM et remise en grille 2D ---
Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/1_Mesh/MESH2025')
Filename_zbed = 'DEM_bed_Corrected.dat'
Beddem = pd.read_csv(Pathroot.joinpath(Filename_zbed), sep=r"\s+", header=None, names=["x", "y", "z"])

# pivot : lignes = y, colonnes = x (trie automatiquement les deux axes)
grid = Beddem.pivot(index="y", columns="x", values="z")
assert not grid.isna().all(axis=None), "Grille vide"
if grid.isna().any(axis=None):
    print("Attention : la grille contient des NaN (points manquants)")

xu, yu = grid.columns.to_numpy(), grid.index.to_numpy()
Z = grid.to_numpy()
dx, dy = np.median(np.diff(xu)), np.median(np.diff(yu))

# --- 2. Contour du glacier ---
Filename_GlacierContour = 'Contour_Birch_2025_Julien.dat'
contour = pd.read_csv(Pathroot.joinpath(Filename_GlacierContour), sep=r"\s+", header=None, names=["x", "y"])
polygon = mpltPath(contour[["x", "y"]].to_numpy())

# --- 3. Direction de l'écoulement ---
theta = np.deg2rad(135)                               # angle depuis l'axe x, sens trigo
u = np.array([np.cos(theta), np.sin(theta)])
v = np.array([-np.sin(theta), np.cos(theta)])

# --- 4. Profils parallèles à l'écoulement, restreints au glacier ---
ds = min(dx, dy)
nperseg = 32                                       # à adapter à la taille du glacier
xm, ym = contour["x"].mean(), contour["y"].mean()
half = np.hypot(np.ptp(contour["x"]), np.ptp(contour["y"])) / 2
s = np.arange(-half, half, ds)
offsets = np.arange(-half, half, ds)

spectres, longueurs = [], []
for t in offsets:
    px = xm + s * u[0] + t * v[0]
    py = ym + s * u[1] + t * v[1]

    inside = polygon.contains_points(np.column_stack([px, py]))
    if inside.sum() < nperseg:
        continue

    prof = map_coordinates(Z, [(py - yu[0]) / dy, (px - xu[0]) / dx],
                           order=1, mode="constant", cval=np.nan)
    valid = inside & ~np.isnan(prof)

    edges = np.flatnonzero(np.diff(np.r_[0, valid.astype(int), 0]))
    for i0, i1 in zip(edges[::2], edges[1::2]):
        if i1 - i0 < nperseg:
            continue
        f, Pxx = welch(detrend(prof[i0:i1]), fs=1/ds, window="hann",
                       nperseg=nperseg, scaling="density")
        spectres.append(Pxx)
        longueurs.append((i1 - i0) * ds)

print(f"{len(spectres)} tronçons utilisés, longueur médiane : {np.median(longueurs):.0f} m")
P = np.mean(spectres, axis=0)
f, P = f[1:], P[1:]
wavelength = 1 / f

for lam in [20, 30, 50, 100]:
    sigma = amplitude_bosses(lam, f, P)
    print(f"λ ≈ {lam:4d} m : RMS = {sigma:.2f} m, "
          f"hauteur bumps ≈ {np.sqrt(2)*sigma:.2f} m")

# --- 5. Résultats : DataFrame + export ---
spectre = pd.DataFrame({"wavelength_m": wavelength, "psd_m3": P})
# spectre.to_csv("spectre_longueurs_onde.dat", sep=" ", index=False)

fig, ax = plt.subplots()
ax.loglog(f, P)
ax.set_xlabel(r"$1/\lambda$ (m$^{-1}$)")
ax.set_ylabel(r"Power spectral density (m$^2$/Hz)")
ax.set_xlim(1e-3, 1e-1)
# fig.savefig("spectre_longueurs_onde.svg", bbox_inches="tight")

# Contrôle : DEM, contour et direction de l'écoulement
fig2, ax2 = plt.subplots()
ax2.pcolormesh(xu, yu, Z, shading="auto")
ax2.plot(contour["x"], contour["y"], "r-")
ax2.arrow(xm, ym, half/3*u[0], half/3*u[1], color="k", width=half/100)
ax2.set_aspect("equal")
# fig2.savefig("controle_domaine.png", dpi=150, bbox_inches="tight")

print("Longueur d'onde dominante :", spectre.loc[spectre["psd_m3"].idxmax(), "wavelength_m"], "m")
plt.show()