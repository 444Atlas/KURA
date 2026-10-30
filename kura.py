#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kura : Simulateur interactif du modèle de Kuramoto
Développé à partir du travail de TIPE sur la synchronisation des oscillateurs couplés.

Fonctionnalités :
- Animation fluide en temps réel des oscillateurs sur le cercle unité.
- Curseur (slider) interactif pour ajuster le couplage K entre 0.0 et 10.0 en direct.
- Curseur secondaire pour ajuster la dispersion gamma des fréquences propres.
- Tracé en direct de l'évolution du paramètre d'ordre r(t) avec asymptote théorique.
- Coloration des oscillateurs selon leur pulsation propre omega_i (spectre thermique).
- Boutons Pause/Reprendre, Réinitialisation des phases, et Tirage de nouvelles fréquences.
"""

from collections import deque
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.animation import FuncAnimation
from scipy.stats import cauchy


class KuramotoInteractiveApp:
    def __init__(self, n_oscillators=120, gamma=0.5, omega_0=1.0, k_init=2.0, dt=0.015):
        self.n = n_oscillators
        self.gamma = gamma
        self.omega_0 = omega_0
        self.k = k_init
        self.dt = dt
        self.is_running = True
        self.time = 0.0
        self.r = 0.0

        # Historique glissant pour le tracé de r(t)
        self.history_len = 300
        self.time_history = deque(maxlen=self.history_len)
        self.r_history = deque(maxlen=self.history_len)

        # Initialisation de la physique
        self._init_physics()

        # Configuration de l'interface graphique
        self._init_gui()

    def _init_physics(self):
        """Initialise les phases uniformément et les fréquences selon une loi de Cauchy."""
        self.phi = np.random.uniform(-np.pi, np.pi, self.n)
        # Écrêtage léger de Cauchy pour éviter des dérives infinies extrêmes
        raw_omega = cauchy.rvs(loc=self.omega_0, scale=self.gamma, size=self.n)
        self.omega = np.clip(raw_omega, self.omega_0 - 15 * self.gamma, self.omega_0 + 15 * self.gamma)
        self.time_history.clear()
        self.r_history.clear()
        self.time = 0.0

    def _reset_phases(self):
        """Réinitialise les phases sans modifier les fréquences propres."""
        self.phi = np.random.uniform(-np.pi, np.pi, self.n)
        self.time_history.clear()
        self.r_history.clear()
        self.time = 0.0

    def compute_order_parameter(self):
        """
        Calcule le paramètre d'ordre complexe :
        p = r * exp(i * psi) = (1/N) * sum_j exp(i * phi_j)
        """
        xb = float(np.mean(np.cos(self.phi)))
        yb = float(np.mean(np.sin(self.phi)))
        psi = float(np.arctan2(yb, xb))
        r = float(np.hypot(xb, yb))
        self.r = r
        return r, psi, xb, yb

    def step(self):
        """Un pas d'intégration d'Euler du modèle de Kuramoto."""
        r, psi, xb, yb = self.compute_order_parameter()
        # d(phi_i)/dt = omega_i + K * r * sin(psi - phi_i)
        dphi = self.omega + self.k * r * np.sin(psi - self.phi)
        self.phi = (self.phi + dphi * self.dt + np.pi) % (2 * np.pi) - np.pi
        self.time += self.dt
        self.time_history.append(self.time)
        self.r_history.append(r)
        return r, psi, xb, yb

    def _init_gui(self):
        # Palette sombre et moderne
        self.c_bg = '#0b0f19'         # Fond général
        self.c_card = '#151c2c'       # Fond des axes / cartes
        self.c_text = '#e2e8f0'       # Texte principal
        self.c_muted = '#64748b'      # Texte secondaire / axes
        self.c_accent = '#38bdf8'     # Bleu cyan
        self.c_sync = '#22c55e'       # Vert synchronisé
        self.c_desync = '#ef4444'     # Rouge désynchronisé
        self.c_gold = '#f59e0b'       # Or barycentre

        self.fig = plt.figure(figsize=(13, 7.5), facecolor=self.c_bg)
        self.fig.canvas.manager.set_window_title('KURA - Simulation Interactive de Kuramoto')

        # Disposition en grille :
        # Colonne gauche (0.05 -> 0.48) : Le cercle des oscillateurs
        # Colonne droite (0.56 -> 0.95) : Graphique r(t) + panneau d'informations
        # Bas (0.05 -> 0.95, hauteur 0.16) : Sliders & Boutons

        self.ax_circle = self.fig.add_axes([0.04, 0.24, 0.44, 0.68], facecolor=self.c_card)
        self.ax_circle.set_aspect('equal')
        self.ax_circle.set_xlim(-1.45, 1.45)
        self.ax_circle.set_ylim(-1.45, 1.45)
        self.ax_circle.axis('off')

        # Cercle unité
        angles = np.linspace(0, 2 * np.pi, 300)
        self.ax_circle.plot(np.cos(angles), np.sin(angles), color='#334155', linestyle='-', linewidth=2, zorder=1)
        # Axes de référence discrets
        self.ax_circle.axhline(0, color='#1e293b', linestyle=':', linewidth=1, zorder=0)
        self.ax_circle.axvline(0, color='#1e293b', linestyle=':', linewidth=1, zorder=0)

        # Tracé du cercle du barycentre (rayon r)
        self.r_circle, = self.ax_circle.plot([], [], color=self.c_gold, linestyle='--', linewidth=1.2, alpha=0.5, zorder=2)

        # Flèche / vecteur barycentre
        self.vec_line, = self.ax_circle.plot([], [], color=self.c_gold, linewidth=2.5, zorder=4)
        self.bary_pt, = self.ax_circle.plot([], [], 'o', color=self.c_gold, markersize=10, markeredgecolor='#ffffff', markeredgewidth=1.5, zorder=5)

        # Oscillateurs sur le cercle
        # Normalisation des couleurs selon omega_i
        norm_omega = (self.omega - self.omega_0) / (3.0 * self.gamma + 1e-5)
        norm_omega = np.clip(norm_omega, -1.0, 1.0)
        colors = plt.cm.coolwarm((norm_omega + 1.0) / 2.0)
        self.scat = self.ax_circle.scatter(
            np.cos(self.phi),
            np.sin(self.phi),
            c=colors,
            s=48,
            edgecolors='#0f172a',
            linewidths=0.6,
            zorder=3
        )

        self.ax_circle_title = self.ax_circle.text(
            0, 1.32, "Cercle des phases $\\theta_i$",
            ha='center', va='center', color=self.c_text, fontsize=13, fontweight='bold'
        )

        # Panneau droit : graphique de r(t)
        self.ax_graph = self.fig.add_axes([0.56, 0.44, 0.40, 0.48], facecolor=self.c_card)
        self.ax_graph.set_title("Évolution du paramètre d'ordre $r(t)$", color=self.c_text, fontsize=12, pad=10, fontweight='bold')
        self.ax_graph.set_xlabel("Temps écoulé (s)", color=self.c_muted, fontsize=10)
        self.ax_graph.set_ylabel("Cohérence $r$", color=self.c_muted, fontsize=10)
        self.ax_graph.set_ylim(-0.02, 1.05)
        self.ax_graph.tick_params(colors=self.c_muted, labelsize=9)
        for spine in self.ax_graph.spines.values():
            spine.set_color('#334155')
        self.ax_graph.grid(True, color='#1e293b', linestyle='--', alpha=0.7)

        self.line_r, = self.ax_graph.plot([], [], color=self.c_accent, linewidth=2, label="$r(t)$ instantané")
        self.line_r_inf = self.ax_graph.axhline(0, color=self.c_sync, linestyle=':', linewidth=1.8, label="$r_\\infty$ théorique")
        self.ax_graph.legend(loc='lower right', facecolor=self.c_card, edgecolor='#334155', labelcolor=self.c_text, fontsize=9)

        # Panneau d'informations textuel sous le graphique
        self.ax_info = self.fig.add_axes([0.56, 0.24, 0.40, 0.14], facecolor=self.c_card)
        self.ax_info.axis('off')
        for spine in self.ax_info.spines.values():
            spine.set_color('#334155')
        self.info_text = self.ax_info.text(
            0.04, 0.5, "",
            va='center', ha='left', color=self.c_text, fontsize=10, linespacing=1.6
        )

        # Contrôles en bas :
        # Sliders
        ax_slider_k = self.fig.add_axes([0.12, 0.13, 0.36, 0.035], facecolor=self.c_card)
        self.slider_k = Slider(
            ax=ax_slider_k,
            label='Couplage K ',
            valmin=0.0,
            valmax=10.0,
            valinit=self.k,
            valstep=0.05,
            color=self.c_accent,
            initcolor='none'
        )
        self._style_widget(self.slider_k, ax_slider_k)
        self.slider_k.on_changed(self._on_k_changed)

        ax_slider_gamma = self.fig.add_axes([0.12, 0.06, 0.36, 0.035], facecolor=self.c_card)
        self.slider_gamma = Slider(
            ax=ax_slider_gamma,
            label='Dispersion γ ',
            valmin=0.05,
            valmax=2.0,
            valinit=self.gamma,
            valstep=0.05,
            color=self.c_gold,
            initcolor='none'
        )
        self._style_widget(self.slider_gamma, ax_slider_gamma)
        self.slider_gamma.on_changed(self._on_gamma_changed)

        # Boutons
        ax_btn_pause = self.fig.add_axes([0.56, 0.07, 0.11, 0.07], facecolor=self.c_card)
        self.btn_pause = Button(ax_btn_pause, 'Pause', color=self.c_card, hovercolor='#334155')
        self._style_button(self.btn_pause, self.c_accent)
        self.btn_pause.on_clicked(self._toggle_pause)

        ax_btn_reset = self.fig.add_axes([0.69, 0.07, 0.12, 0.07], facecolor=self.c_card)
        self.btn_reset = Button(ax_btn_reset, 'Mélanger phases', color=self.c_card, hovercolor='#334155')
        self._style_button(self.btn_reset, self.c_text)
        self.btn_reset.on_clicked(lambda event: self._reset_phases())

        ax_btn_new = self.fig.add_axes([0.83, 0.07, 0.13, 0.07], facecolor=self.c_card)
        self.btn_new = Button(ax_btn_new, 'Nouvel échantillon', color=self.c_card, hovercolor='#334155')
        self._style_button(self.btn_new, self.c_text)
        self.btn_new.on_clicked(lambda event: self._new_frequencies())

    def _style_widget(self, slider, ax):
        ax.tick_params(colors=self.c_muted, labelsize=8)
        slider.label.set_color(self.c_text)
        slider.label.set_fontsize(10)
        slider.valtext.set_color(self.c_accent)
        slider.valtext.set_fontsize(10)
        for spine in ax.spines.values():
            spine.set_color('#334155')

    def _style_button(self, btn, text_color):
        btn.label.set_color(text_color)
        btn.label.set_fontsize(9)
        btn.label.set_fontweight('bold')
        for spine in btn.ax.spines.values():
            spine.set_color('#334155')

    def _on_k_changed(self, val):
        self.k = float(val)

    def _on_gamma_changed(self, val):
        self.gamma = float(val)
        # Recalculer les fréquences avec le nouveau gamma
        self._new_frequencies()

    def _toggle_pause(self, event=None):
        self.is_running = not self.is_running
        self.btn_pause.label.set_text('Reprendre' if not self.is_running else 'Pause')
        self.fig.canvas.draw_idle()

    def _new_frequencies(self):
        self._init_physics()
        # Mettre à jour les couleurs selon les nouvelles fréquences
        norm_omega = (self.omega - self.omega_0) / (3.0 * self.gamma + 1e-5)
        norm_omega = np.clip(norm_omega, -1.0, 1.0)
        colors = plt.cm.coolwarm((norm_omega + 1.0) / 2.0)
        self.scat.set_color(colors)

    def update(self, frame):
        if self.is_running:
            # Plusieurs micro-pas par frame d'animation pour fluidité physique
            for _ in range(2):
                r, psi, xb, yb = self.step()
        else:
            r, psi, xb, yb = self.compute_order_parameter()

        # Mise à jour des oscillateurs
        x = np.cos(self.phi)
        y = np.sin(self.phi)
        self.scat.set_offsets(np.column_stack((x, y)))

        # Cercle de rayon r
        angles = np.linspace(0, 2 * np.pi, 100)
        self.r_circle.set_data(r * np.cos(angles), r * np.sin(angles))

        # Vecteur et point barycentre
        self.vec_line.set_data([0, xb], [0, yb])
        self.bary_pt.set_data([xb], [yb])

        # Seuil critique analytique de Kuramoto : Kc = 2 * Gamma
        kc = 2.0 * self.gamma
        is_synchronized = self.k >= kc
        r_inf = np.sqrt(max(0.0, 1.0 - kc / self.k)) if is_synchronized and self.k > 0 else 0.0

        # Mise à jour graphique r(t)
        if len(self.time_history) > 1:
            t_data = list(self.time_history)
            r_data = list(self.r_history)
            self.line_r.set_data(t_data, r_data)
            self.ax_graph.set_xlim(max(0.0, t_data[-1] - 8.0), max(8.0, t_data[-1] + 0.5))

        if is_synchronized:
            self.line_r_inf.set_ydata([r_inf, r_inf])
            self.line_r_inf.set_visible(True)
        else:
            self.line_r_inf.set_visible(False)

        # Panneau texte dynamique
        self.r = r
        status_label = (
            "SYNCHRONISÉ" if is_synchronized and abs(self.r - r_inf) < 0.1
            else "TRANSITION DE PHASE" if is_synchronized and abs(self.r - r_inf) >= 0.1
            else "DÉSYNCHRONISÉ"
        )
        status_color = (
            self.c_sync if status_label == "SYNCHRONISÉ"
            else self.c_gold if status_label == "TRANSITION DE PHASE"
            else self.c_desync
        )

        info_str = (
            f"Régime : {status_label}\n"
            f"Paramètre d'ordre instantané : r(t) = {r:.3f}   |   Phase moyenne : ψ = {psi:.2f} rad\n"
            f"Couplage K = {self.k:.2f}   |   Seuil critique Kc = 2γ = {kc:.2f}   |   Asymptote r∞ = {r_inf:.3f}"
        )
        self.info_text.set_text(info_str)

        return (
            self.scat, self.r_circle, self.vec_line, self.bary_pt,
            self.line_r, self.line_r_inf, self.info_text
        )

    def run(self):
        """Lance l'animation interactive."""
        self.ani = FuncAnimation(
            self.fig,
            self.update,
            interval=25,
            blit=False,
            cache_frame_data=False
        )
        plt.show()


def main():
    app = KuramotoInteractiveApp(
        n_oscillators=120,
        gamma=0.5,
        omega_0=1.0,
        k_init=2.0,
        dt=0.015
    )
    app.run()


if __name__ == '__main__':
    main()

