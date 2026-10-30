# Modèle de Kuramoto : Synchronisation d'oscillateurs couplés

Ce projet reprend et fait fonctionner les codes présentés dans les **Annexes 5, 6 et 7** de l'exposé de TIPE (*Étude de la synchronisation des oscillateurs*, Paul Laurent-Verleye, 2025).

---

## 1. Structure des fichiers

| Fichier | Source | Description |
| :--- | :--- | :--- |
| `annexe5_animation.py` | Annexe 5 (page 28) | Animation dynamique des oscillateurs sur le cercle unité et de leur barycentre (paramètre d'ordre $\vec{p}$). |
| `annexe6_trouver_kc.py` | Annexe 6 (page 29) | Simulation de la transition de phase et tracé de $r_\infty$ en fonction de $K$ pour identifier le couplage critique $K_c$. |
| `annexe7_trouver_r_inf.py` | Annexe 7 (page 30) | Évolution temporelle du paramètre d'ordre $r(t)$ vers sa valeur limite stationnaire $r_\infty = \sqrt{1 - K_c/K}$. |
| `kuramoto_complet.py` | Synthèse & Optimisation | Version modulaire, vectorisée et propre (sans récursion profonde ni variables globales). Permet de reproduire toutes les figures (Figures 4 à 13). |

---

## 2. Corrections et améliorations apportées

### Dans `annexe5_animation.py` (Annexe 5) :
1. **Imports manquants** : Ajout de `numpy`, `matplotlib.pyplot`, `FuncAnimation`, `cauchy` (depuis `scipy.stats`), `uniform` (depuis `random`), ainsi que `atan2` et `sqrt`.
2. **Correction du tracé du barycentre** : Dans la version d'origine, `bari.set_data(xb, yb)` provoquait une erreur Matplotlib (`TypeError: float object is not iterable`), car `Line2D.set_data()` attend des séquences. Corrigé en `bari.set_data([xb], [yb])`.
3. **Persistance de l'animation** : `ani` est désormais retourné pour éviter qu'il ne soit supprimé par le garbage collector.
4. **Option d'export GIF** : Paramètre `save_path` pour exporter l'animation directement (ex. `animation_kuramoto.gif`).

### Dans `annexe6_trouver_kc.py` (Annexe 6) :
1. **Imports manquants** : Ajout de `numpy`, `matplotlib.pyplot`, `cauchy`, `uniform`, `atan2`, `sqrt`.
2. **Comptage des échantillons** : Dans la récursion initiale, le 10ᵉ échantillon était calculé mais non ajouté à la liste `P` avant de basculer sur le $K$ suivant. Tous les échantillons sont maintenant proprement collectés.
3. **Superposition théorie / expérience** : Le graphique affiche à la fois les points mesurés avec leur incertitude type ($\frac{\sigma}{\sqrt{N}}$) et la courbe théorique $r_\infty(K) = \sqrt{1 - \frac{K_c}{K}}$ pour $K > K_c = 2\gamma = 1.0$.
4. **Sauvegarde automatique** : Sauvegarde sous `figure_kc.png`.

### Dans `annexe7_trouver_r_inf.py` (Annexe 7) :
1. **Imports manquants** : Ajout de `numpy`, `matplotlib.pyplot`, `cauchy`, etc.
2. **Ajout de l'asymptote théorique** : Tracé de la droite horizontale en pointillés $r_\infty = \sqrt{1 - K_c/K} = \sqrt{1 - 1/3} \approx 0.8165$ pour comparaison directe avec la simulation (Figure 13 de la présentation).
3. **Sauvegarde automatique** : Sauvegarde sous `figure_r_inf.png`.

---

## 3. Exécution

### Installer les dépendances :
```bash
pip install numpy scipy matplotlib pillow
```

### Lancer les scripts des annexes :
```bash
# Annexe 5 : Animation
python3 annexe5_animation.py

# Annexe 6 : Détermination de Kc (courbe r vs K)
python3 annexe6_trouver_kc.py

# Annexe 7 : Relaxation temporelle r(t)
python3 annexe7_trouver_r_inf.py
```

### Script complet optimisé (`kuramoto_complet.py`) :
```bash
# Générer toutes les figures et animations en mode rapide (pour test)
python3 kuramoto_complet.py --fast

# Générer les figures avec la haute résolution complète (N = 3000)
python3 kuramoto_complet.py --mode all

# Générer uniquement l'évolution temporelle r(t)
python3 kuramoto_complet.py --mode r_inf

# Générer uniquement la transition de phase r vs K
python3 kuramoto_complet.py --mode kc

# Générer les instantanés du cercle unité (Figures 4, 5, 6)
python3 kuramoto_complet.py --mode snapshots
```

