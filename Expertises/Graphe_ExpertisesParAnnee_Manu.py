import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_excel("/home/brondexj/Documents/Concours/Inventaire_demandes_Expertises_RisqueGL.xlsx")
annees = df.iloc[:, 0]
sollicitations = df.iloc[:, 1]

fig, ax = plt.subplots()

# ligne + étoiles
ax.plot(annees, sollicitations, linestyle='-', color='black')

mask = sollicitations > 0
ax.plot(
    annees[mask],
    sollicitations[mask],
    linestyle='None',
    marker='*',
    markersize=16,
    markerfacecolor=(95/255, 141/255, 211/255),
    markeredgecolor='black',
    markeredgewidth=0.8
)

# ticks majeurs (1 année sur 2)
ax.set_xticks(annees[::2])

# ticks mineurs (toutes les années)
ax.set_xticks(annees, minor=True)

ax.tick_params(axis='x', rotation=45, labelsize=18)
ax.tick_params(axis='y', labelsize=18)
# grille forte (années affichées)
ax.grid(True, which='major', linewidth=0.9, alpha=0.7)

# grille faible (années intermédiaires)
ax.grid(True, which='minor', linewidth=0.4, alpha=0.3)

ax.set_ylabel("Nombre de sites avec sollicitation ROG", fontsize=22)

plt.show()