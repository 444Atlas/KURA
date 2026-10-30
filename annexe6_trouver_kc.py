"""
Annexe 6 : Code pour trouver Kc (seuil critique de couplage dans le modèle de Kuramoto)
Fidèle au code de la planche 29 de la présentation, avec :
- Les imports nécessaires ajoutés
- La prise en compte propre de tous les échantillons
- Sauvegarde du graphique et affichage
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import cauchy
from random import uniform
from math import atan2, sqrt

P = []
lp = []
lk = []
perr = []

def test(N, Gamma, K, c):
    w0 = 1
    global phi, omega
    phi = np.array([uniform(-np.pi, np.pi) for _ in range(N)])
    omega = np.array(cauchy.rvs(loc=w0, scale=Gamma, size=N))
    for _ in range(3000):
        xb, yb = np.mean(np.cos(phi)), np.mean(np.sin(phi))
        psi = atan2(yb, xb)
        p = sqrt(xb**2 + yb**2)
        phi += (omega + K * p * np.sin(psi - phi)) * 0.01

    P.append(p)
    if c < 10:
        c += 1
        test(N, Gamma, K, c)
    elif K > 0.1:
        lp.append(np.mean(P))
        perr.append(np.std(P) / np.sqrt(len(P)))
        lk.append(K)
        P.clear()
        print(f"K = {K:.2f} -> <r> = {lp[-1]:.4f} +/- {perr[-1]:.4f}")
        test(N, Gamma, round(K - 0.1, 2), 1)
    else:
        display_graph()

def display_graph():
    fig_graph, ax_graph = plt.subplots(figsize=(8, 6))
    ax_graph.set_xlabel('K')
    ax_graph.set_ylabel('r')
    ax_graph.set_title('Comportement de r en fonction de K')
    
    # Points mesurés avec barres d'erreur
    ax_graph.errorbar(lk, lp, yerr=perr, fmt='.', label=r'$r_\infty$ mesuré avec incertitude')

    # Courbe théorique : r = 0 si K <= Kc, et r = sqrt(1 - Kc/K) si K > Kc
    # Avec Gamma = 0.5, Kc = 2 * Gamma = 1.0
    x1 = np.linspace(0.2, 1, 200)
    y1 = np.zeros_like(x1)
    x2 = np.linspace(1, 3, 900)
    y2 = np.sqrt(1 - 1 / x2)

    ax_graph.plot(x2, y2, color='blue', label=r'Théorie : $\sqrt{1 - \frac{K_c}{K}}$')
    ax_graph.plot(x1, y1, color='blue')
    ax_graph.grid(True, linestyle=':', alpha=0.6)
    ax_graph.legend()
    plt.savefig('figure_kc.png', dpi=300)
    print("Graphique sauvegardé sous 'figure_kc.png'")
    plt.show()

if __name__ == "__main__":
    test(N=3000, Gamma=0.5, K=3, c=1)

