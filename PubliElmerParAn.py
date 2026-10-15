import numpy as np
import matplotlib.pyplot as plt

# Données
annees = [2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013,
           2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023,
           2024, 2025]

publis_par_an = [1, 0, 1, 3, 2, 4, 2, 7, 11, 11,
                 15, 9, 9, 15, 25, 22, 22, 22, 16, 12,
                 18, 8]

# Cumul
publis_cumulees = np.cumsum(publis_par_an)

# Figure
fig, ax = plt.subplots(figsize=(9,6))

ax.plot(
    annees,
    publis_cumulees,
    marker='o',
    markersize=10,
    linewidth=2.5
)

# Labels
ax.set_xlabel("Année", fontsize=22)
ax.set_ylabel("Nombre cumulé de publications", fontsize=22)
ax.set_title("Publications cumulées Elmer/Ice", fontsize=24)

# Ticks majeurs : une année sur deux
ax.set_xticks(annees[::2])

# Ticks mineurs : toutes les années
ax.set_xticks(annees, minor=True)

# Taille et rotation des ticks
ax.tick_params(axis='x', labelsize=18, rotation=45)
ax.tick_params(axis='y', labelsize=18)

# Grilles
ax.grid(True, which='major', linewidth=0.9, alpha=0.7)
ax.grid(True, which='minor', linewidth=0.4, alpha=0.3)

plt.tight_layout()
plt.show()