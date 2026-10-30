#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Étude de la synchronisation des oscillateurs couplés - Modèle de Kuramoto
Auteur original : Paul LAURENT-VERLEYE (TIPE 2025)
Version restructurée, optimisée et documentée.

Ce module implémente :
1. La simulation numérique du modèle de Kuramoto avec méthode d'Euler.
2. La visualisation dynamique (animation) des oscillateurs sur le cercle unité.
3. La transition de phase et la détermination du couplage critique Kc = 2 * Gamma.
4. L'évolution temporelle du paramètre d'ordre r(t) vers son asymptote théorique r_infini.
5. La génération de figures comparatives à la théorie.
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.stats import cauchy


class KuramotoSimulation:
    """
    Simulateur vectorisé pour le modèle de Kuramoto :
    d(theta_i)/dt = omega_i + (K/N) * sum_j sin(theta_j - theta_i)
                  = omega_i + K * r(t) * sin(psi(t) - theta_i(t))
    """

    def __init__(self, n_oscillators=3000, gamma=0.5, omega_0=1.0, k=2.0, dt=0.01, seed=None):
        if seed is not None:
            np.random.seed(seed)
        self.n = n_oscillators
        self.gamma = gamma
        self.omega_0 = omega_0
        self.k = k
        self.dt = dt

        # Initialisation : phases uniformes sur [-pi, pi], fréquences propres lorentziennes (Cauchy)
        self.phi = np.random.uniform(-np.pi, np.pi, self.n)
        self.omega = cauchy.rvs(loc=self.omega_0, scale=self.gamma, size=self.n)

    def order_parameter(self, phi=None):
        """
        Calcule le paramètre d'ordre complexe :
        p = r * exp(i * psi) = (1/N) * sum_j exp(i * phi_j)
        Retourne (r, psi, xb, yb)
        """
        if phi is None:
            phi = self.phi
        xb = float(np.mean(np.cos(phi)))
        yb = float(np.mean(np.sin(phi)))
        psi = float(np.arctan2(yb, xb))
        r = float(np.hypot(xb, yb))
        return r, psi, xb, yb

    def step(self):
        """Avance d'un pas dt selon la méthode d'Euler."""
        r, psi, xb, yb = self.order_parameter()
        self.phi += (self.omega + self.k * r * np.sin(psi - self.phi)) * self.dt
        return r, psi, xb, yb

    def run(self, steps=3000, sample_interval=100):
        """Exécute la simulation pendant un nombre donné de pas."""
        t_list = []
        r_list = []
        for s in range(steps):
            r, _, _, _ = self.step()
            if s % sample_interval == 0:
                t_list.append(s * self.dt)
                r_list.append(r)
        return np.array(t_list), np.array(r_list)


def plot_snapshots(n=50, gamma=0.5, k=2.0, dt=0.01, seed=42, save_path="figures_snapshots.png"):
    """
    - État initial (désordonné, barycentre proche du centre)
    - État intermédiaire
    - État final (oscillateurs synchronisés en amas, barycentre vers la périphérie)
    """
    sim = KuramotoSimulation(n_oscillators=n, gamma=gamma, omega_0=1.0, k=k, dt=dt, seed=seed)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    circle_angles = np.linspace(0, 2 * np.pi, 200)

    stages = [
        ("État initial (t = 0)", 0),
        ("État intermédiaire (t = 5)", int(5.0 / dt)),
        ("État final (t = 25)", int(20.0 / dt))
    ]

    for idx, (title, steps_to_run) in enumerate(stages):
        for _ in range(steps_to_run):
            sim.step()
        r, psi, xb, yb = sim.order_parameter()
        ax = axes[idx]
        ax.set_aspect('equal')
        ax.set_xlim(-1.4, 1.4)
        ax.set_ylim(-1.4, 1.4)
        ax.plot(np.cos(circle_angles), np.sin(circle_angles), 'b--', alpha=0.5)
        ax.scatter(np.cos(sim.phi), np.sin(sim.phi), color='blue', s=25, label='Oscillateurs')
        ax.scatter([xb], [yb], color='red', s=60, zorder=5, label='Barycentre $\\vec{p}$')
        ax.set_title(f"{title}\n$r = {r:.2f}$", fontsize=11)
        ax.axis('off')
        if idx == 0:
            ax.legend(loc='lower left', frameon=True)

    plt.suptitle(f"Transition de phase de Kuramoto ($N={n}$, $\\gamma={gamma}$, $K={k}$)", fontsize=13)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"Instantanés sauvegardés sous '{save_path}'")
    plt.close(fig)


def simulate_r_vs_t(n=3000, gamma=0.5, k=3.0, dt=0.01, steps=3000, seed=42, save_path="figure_r_vs_t.png"):
    """
    Évolution de r(t) et comparaison avec r_infini = sqrt(1 - Kc / K)
    """
    print(f"Simulation de r(t) pour N={n}, Gamma={gamma}, K={k}...")
    sim = KuramotoSimulation(n_oscillators=n, gamma=gamma, omega_0=1.0, k=k, dt=dt, seed=seed)
    t_arr, r_arr = sim.run(steps=steps, sample_interval=100)

    kc = 2.0 * gamma
    r_inf_theo = np.sqrt(max(0.0, 1.0 - kc / k)) if k > kc else 0.0

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t_arr, r_arr, 'o-', color='#1f77b4', markersize=4, label='Mesure numérique $r(t)$')
    ax.axhline(
        y=r_inf_theo,
        color='#d62728',
        linestyle='--',
        linewidth=2,
        label=rf'Valeur théorique $r_\infty = \sqrt{{1 - \frac{{K_c}}{{K}}}} \approx {r_inf_theo:.3f}$'
    )

    ax.set_xlabel('Temps $t$', fontsize=12)
    ax.set_ylabel("Paramètre d'ordre $r$", fontsize=12)
    ax.set_title(f"Comportement de $r$ en fonction de $t$\n($N={n}$, $\\gamma={gamma}$, $K={k}$, $dt={dt}$)", fontsize=13)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"Graphique r(t) sauvegardé sous '{save_path}'")
    plt.close(fig)


def simulate_kc(n=3000, gamma=0.5, k_vals=None, n_repeats=10, dt=0.01, steps=3000, save_path="figure_kc.png"):
    """
    Transition de phase et détermination de Kc.
    Mesure r_final en fonction de K avec barres d'erreur et courbe théorique.
    """
    if k_vals is None:
        # Valeurs de K de 3.0 à 0.2 par pas de 0.1
        k_vals = np.round(np.arange(3.0, 0.1, -0.1), 2)

    kc_theo = 2.0 * gamma
    print(f"Simulation de r en fonction de K ({len(k_vals)} valeurs de K, {n_repeats} répétitions, N={n})...")

    lk = []
    lp = []
    perr = []

    for k in k_vals:
        p_samples = []
        for rep in range(n_repeats):
            sim = KuramotoSimulation(n_oscillators=n, gamma=gamma, omega_0=1.0, k=k, dt=dt)
            for _ in range(steps):
                sim.step()
            r, _, _, _ = sim.order_parameter()
            p_samples.append(r)

        mean_p = float(np.mean(p_samples))
        err_p = float(np.std(p_samples) / np.sqrt(len(p_samples)))
        lk.append(k)
        lp.append(mean_p)
        perr.append(err_p)
        print(f"  K = {k:4.2f} -> <r> = {mean_p:6.4f} +/- {err_p:6.4f}")

    # Tracé
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.errorbar(lk, lp, yerr=perr, fmt='o', color='#1f77b4', ecolor='#1f77b4',
                capsize=3, label=r'$r_\infty$ mesuré avec incertitude')

    # Courbe théorique
    x_sub = np.linspace(0.2, kc_theo, 200)
    y_sub = np.zeros_like(x_sub)
    x_sup = np.linspace(kc_theo, 3.0, 500)
    y_sup = np.sqrt(1.0 - kc_theo / x_sup)

    ax.plot(x_sub, y_sub, color='#2ca02c', linewidth=2)
    ax.plot(x_sup, y_sup, color='#2ca02c', linewidth=2,
            label=rf'Théorie : $r_\infty = \sqrt{{1 - \frac{{K_c}}{{K}}}}$ ($K_c={kc_theo}$)')
    ax.axvline(x=kc_theo, color='gray', linestyle=':', label=rf'$K_c = 2\gamma = {kc_theo}$')

    ax.set_xlabel('Coefficient de couplage $K$', fontsize=12)
    ax.set_ylabel(r"Paramètre d'ordre asymptotique $r_\infty$", fontsize=12)
    ax.set_title(f"Transition de phase : comportement de $r_\\infty$ en fonction de $K$\n"
                 f"($N={n}$, $\\gamma={gamma}$, $dt={dt}$)", fontsize=13)
    ax.set_xlim(0.1, 3.1)
    ax.set_ylim(-0.05, 0.95)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"Graphique Kc sauvegardé sous '{save_path}'")
    plt.close(fig)


def create_animation(n=50, gamma=0.1, k=1.5, dt=0.01, duration=20.0, save_path="animation.gif", fps=30):
    """
    Crée une animation fluide  des oscillateurs sur le cercle unité
    et du paramètre d'ordre représenté par le barycentre bleu.
    """
    print(f"Création de l'animation pour N={n}, Gamma={gamma}, K={k}...")
    sim = KuramotoSimulation(n_oscillators=n, gamma=gamma, omega_0=1.0, k=k, dt=dt)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect('equal')
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.axis('off')

    circle = plt.Circle((0, 0), 1, fill=False, color='#4682b4', linestyle='--', linewidth=1.5)
    ax.add_artist(circle)

    points, = ax.plot([], [], 'ro', markersize=6, alpha=0.8, label='Oscillateurs')
    bari, = ax.plot([], [], 'bo', markersize=9, label='Barycentre $\\vec{p}$')
    line_p, = ax.plot([], [], 'b-', linewidth=1.5, alpha=0.6)
    title_text = ax.text(0, 1.3, '', ha='center', fontsize=11)
    ax.legend(loc='lower left', frameon=False)

    n_frames = min(int(duration / dt), 300)

    def init():
        points.set_data([], [])
        bari.set_data([], [])
        line_p.set_data([], [])
        title_text.set_text('')
        return points, bari, line_p, title_text

    def update(frame):
        r, psi, xb, yb = sim.step()
        x = np.cos(sim.phi)
        y = np.sin(sim.phi)
        points.set_data(x, y)
        bari.set_data([xb], [yb])
        line_p.set_data([0, xb], [0, yb])
        title_text.set_text(f"$t = {frame * dt:.2f}$ s | $r(t) = {r:.2f}$")
        return points, bari, line_p, title_text

    ani = FuncAnimation(
        fig,
        update,
        frames=n_frames,
        init_func=init,
        interval=int(1000 * dt),
        blit=True
    )

    if save_path:
        print(f"Export en cours vers {save_path}...")
        ani.save(save_path, writer=PillowWriter(fps=fps))
        print(f"Animation sauvegardée sous '{save_path}'")
    else:
        plt.show()

    plt.close(fig)
    return ani


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulation du modèle de Kuramoto (TIPE)")
    parser.add_argument("--mode", choices=["all", "r_inf", "kc", "snapshots", "animation"], default="all",
                        help="Partie à exécuter (default: all)")
    parser.add_argument("--fast", action="store_true",
                        help="Mode rapide (N plus petit pour tester rapidement)")
    args = parser.parse_args()

    n_val = 600 if args.fast else 3000
    n_rep = 3 if args.fast else 5

    if args.mode in ["all", "snapshots"]:
        plot_snapshots(n=50, gamma=0.5, k=2.0)

    if args.mode in ["all", "r_inf"]:
        simulate_r_vs_t(n=n_val, gamma=0.5, k=3.0)

    if args.mode in ["all", "animation"]:
        create_animation(n=50, gamma=0.1, k=1.5, save_path="animation_kuramoto.gif")

    if args.mode in ["all", "kc"]:
        k_list = np.round(np.arange(3.0, 0.1, -0.2 if args.fast else -0.1), 2)
        simulate_kc(n=n_val, gamma=0.5, k_vals=k_list, n_repeats=n_rep)
