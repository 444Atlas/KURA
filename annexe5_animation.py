"""
Annexe 5 : Code pour l'animation de la synchronisation des oscillateurs (Kuramoto)
Fidèle au code de la planche 28 de la présentation, avec :
- Les imports nécessaires ajoutés
- La correction de `bari.set_data([xb], [yb])` (Line2D attend des séquences et non des scalaires)
- La possibilité optionnelle d'enregistrer en GIF / vidéo via `save_path`
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.stats import cauchy
from random import uniform
from math import atan2, sqrt

fig, ax = plt.subplots(figsize=(6, 6))
ax.set_aspect('equal')
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_axis_off()

circle = plt.Circle((0, 0), 1, fill=False, color='blue', linestyle='--')
ax.add_artist(circle)
points, = ax.plot([], [], 'ro', markersize=6)
bari, = ax.plot([], [], 'bo', markersize=8)

def test(N=30, Gamma=0.1, K=1.5, w0=1, save_path=None):
    dt = 0.01
    duree = 2000 * dt
    phi = np.array([uniform(-np.pi, np.pi) for _ in range(N)])
    omega = np.array(cauchy.rvs(loc=w0, scale=Gamma, size=N))

    def init():
        points.set_data([], [])
        bari.set_data([], [])
        return points, bari

    def update(frame):
        nonlocal phi, omega, K
        x = np.cos(phi)
        y = np.sin(phi)
        xb, yb = np.mean(np.cos(phi)), np.mean(np.sin(phi))
        psi = atan2(yb, xb)
        p = sqrt(xb**2 + yb**2)
        phi += (omega + K * p * np.sin(psi - phi)) * dt
        points.set_data(x, y)
        bari.set_data([xb], [yb])  # Correction : Line2D nécessite des listes [xb], [yb]
        return points, bari

    total_frames = int(duree / dt)
    ani = FuncAnimation(
        fig,
        update,
        frames=total_frames,
        init_func=init,
        interval=int(1000 * dt),
        blit=True
    )

    if save_path:
        print(f"Sauvegarde de l'animation dans {save_path} ({min(total_frames, 200)} frames)...")
        # Pour une sauvegarde rapide en aperçu
        ani_save = FuncAnimation(
            fig,
            update,
            frames=min(total_frames, 200),
            init_func=init,
            blit=True
        )
        ani_save.save(save_path, writer=PillowWriter(fps=30))
        print(f"Animation sauvegardée avec succès : {save_path}")
    else:
        plt.show()

    return ani

if __name__ == "__main__":
    test(N=3000, Gamma=0, K=0.3)
