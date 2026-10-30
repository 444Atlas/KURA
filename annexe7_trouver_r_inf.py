"""
Annexe 7 : Code pour trouver r_infini (évolution temporelle du paramètre d'ordre)
Fidèle au code de la planche 30 de la présentation, avec :
- Les imports nécessaires ajoutés
- Tracé de r(t) et comparaison avec la valeur théorique asymptotique r_infini
- Sauvegarde et affichage
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import cauchy
from random import uniform
from math import atan2, sqrt

lp = []
lt = []

def test(N, Gamma, K):
    w0 = 1
    global phi, omega
    phi = np.array([uniform(-np.pi, np.pi) for _ in range(N)])
    omega = np.array(cauchy.rvs(loc=w0, scale=Gamma, size=N))
    for s in range(3000):
        xb, yb = np.mean(np.cos(phi)), np.mean(np.sin(phi))
        psi = atan2(yb, xb)
        p = sqrt(xb**2 + yb**2)
        phi += (omega + K * p * np.sin(psi - phi)) * 0.01
        if s % 100 == 0:
            lp.append(p)
            lt.append(s // 100)
    display_graph(Gamma, K)

def display_graph(Gamma=0.5, K=3):
    fig_graph, ax_graph = plt.subplots(figsize=(8, 6))
    ax_graph.set_xlabel('t')
    ax_graph.set_ylabel('r')
    ax_graph.set_title('Comportement de r en fonction de t')
    ax_graph.plot(lt, lp, '.-', label='r(t) numérique')

    # Valeur théorique r_infini = sqrt(1 - Kc/K) avec Kc = 2*Gamma
    Kc = 2 * Gamma
    if K > Kc:
        r_inf = np.sqrt(1 - Kc / K)
        ax_graph.axhline(
            y=r_inf,
            color='coral',
            linestyle='--',
            label=rf'$r_\infty = \sqrt{{1 - K_c/K}} \approx {r_inf:.3f}$'
        )

    ax_graph.set_ylim(-0.05, 1.05)
    ax_graph.grid(True, linestyle=':', alpha=0.6)
    ax_graph.legend()
    plt.savefig('figure_r_inf.png', dpi=300)
    print("Graphique sauvegardé sous 'figure_r_inf.png'")
    plt.show()

if __name__ == "__main__":
    test(N=3000, Gamma=0.5, K=3)

