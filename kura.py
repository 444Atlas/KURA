#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kura : Simulateur interactif du modèle de Kuramoto avec retard et connectome cérébral
Développé à partir du travail de TIPE sur la synchronisation des oscillateurs couplés.

Fonctionnalités :
- Animation fluide en temps réel des oscillateurs sur le cercle unité.
- Couleurs des oscillateurs sur le cercle reflétant leur fréquence (intrinsèque ou instantanée).
- Synchronisation temporelle stricte sur l'horloge réelle (1 s simulée = 1 s réelle).
- Curseur de pulsation w0 (de 0 à 10 rad/s) contrôlant fidèlement la vitesse réelle de rotation.
- Curseur de couplage adimensionné K / w0 (de 0 à 1).
- Curseur de retard adimensionné w0 * tau (de 0 à 10*pi).
- Sélection directe au clic/glisser sur l'image des zones de synchronisation pour fixer (w0*tau, K/w0).
- Matrice de connectivité structurelle du cerveau humain (AAL 90 régions HCP) issue de Hopf_Delay_Toolbox.
- Forçage automatique de N = 90 lors de la manipulation du retard.
- Diagramme des zones de synchronisation de Yeung & Strogatz (PRL 1999) verrouillé en déplacement,
  respectant l'échelle exacte de la publication et le color scheme sombre du graphe de r(t).
- Deux boutons de contrôle : Pause et Reset.
"""

import os
import time
from collections import deque
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.animation import FuncAnimation
import matplotlib.patches as mpatches
from scipy.stats import cauchy
import scipy.io as sio


class KuramotoInteractiveApp:
    def __init__(self, n_oscillators=90, gamma=0.0, omega_0=1.0, k_over_w0=0.8, w0_tau=0.0, dt=0.005):
        self.n = int(n_oscillators)
        self.gamma = float(gamma)
        self.omega_0 = float(omega_0)
        self.k_over_w0 = float(k_over_w0)
        self.w0_tau = float(w0_tau)
        self.dt = float(dt)

        # Grandeurs physiques déduites
        self.k = self.k_over_w0 * self.omega_0
        self.tau = (self.w0_tau / self.omega_0) if self.omega_0 > 1e-5 else 0.0

        self.is_running = True
        self.time = 0.0
        self.r = 0.0

        # Synchronisation stricte sur l'horloge temps réel
        self.last_real_time = None
        self.time_accumulator = 0.0

        # Chargement de la matrice de connectivité structurelle AAL (90 régions)
        self.C_aal = self._load_structural_connectivity()

        # Historique glissant pour le tracé de r(t)
        self.history_len = 300
        self.time_history = deque(maxlen=self.history_len)
        self.r_history = deque(maxlen=self.history_len)

        # Initialisation de la physique et du buffer circulaire de retard
        self._init_physics()

        # Configuration de l'interface graphique
        self._init_gui()

    def _load_structural_connectivity(self):
        """Charge et normalise la matrice de connectivité AAL 90 issue de Hopf_Delay_Toolbox."""
        script_dir = os.path.dirname(os.path.abspath(__file__))
        mat_path = os.path.join(script_dir, 'SC_90aal_32HCP.mat')
        if os.path.exists(mat_path):
            try:
                data = sio.loadmat(mat_path)
                raw_c = data['mat'].astype(np.float64)
                np.fill_diagonal(raw_c, 0.0)
                non_diag = raw_c[np.eye(raw_c.shape[0]) == 0]
                mean_val = float(np.mean(non_diag)) if np.mean(non_diag) > 0 else 1.0
                return raw_c / mean_val
            except Exception as e:
                print(f"Avertissement : Erreur lors du chargement de {mat_path}: {e}")
        # Fallback tous-à-tous normalisé si fichier absent
        fallback_c = np.ones((90, 90), dtype=np.float64)
        np.fill_diagonal(fallback_c, 0.0)
        return fallback_c

    def _init_physics(self):
        """Initialise les phases uniformément, les fréquences et le buffer circulaire d'historique."""
        self.phi = np.random.uniform(-np.pi, np.pi, self.n)

        # Dispersions intrinsèques des fréquences autour de omega_0
        if self.gamma <= 0.0:
            self.delta_omega = np.zeros(self.n)
            self.omega = np.full(self.n, self.omega_0)
        else:
            raw_delta = cauchy.rvs(loc=0.0, scale=self.gamma, size=self.n)
            self.delta_omega = np.clip(raw_delta, -15.0 * self.gamma, 15.0 * self.gamma)
            self.omega = self.omega_0 + self.delta_omega

        # Vitesses instantanées initiales
        self.dphi = np.copy(self.omega)

        # Matrice de connectivité adaptée au nombre d'oscillateurs
        if self.n == 90 and self.C_aal is not None:
            self.C = self.C_aal
        else:
            self.C = np.ones((self.n, self.n), dtype=np.float64)
            np.fill_diagonal(self.C, 0.0)

        # Buffer circulaire d'historique pour le retard tau (jusqu'à 10*pi s)
        self.max_buf = int(np.ceil(10.0 * np.pi / self.dt)) + 100
        self.history = np.tile(self.phi, (self.max_buf, 1))
        self.ptr = 0

        self.time_history.clear()
        self.r_history.clear()
        self.time = 0.0
        self.r = 0.0
        self.last_real_time = None
        self.time_accumulator = 0.0

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
        """
        Un pas d'intégration d'Euler du modèle de Kuramoto avec retard et matrice de connectivité :
        dphi_i / dt = omega_i + (K / N) * sum_j C_ij * sin(phi_j(t - tau) - phi_i(t))
        """
        r, psi, xb, yb = self.compute_order_parameter()
        delay_steps = int(round(self.tau / self.dt))

        if delay_steps == 0 and self.n != 90:
            # Traitement rapide champ moyen O(N) pour grand N sans retard
            dphi = self.omega + self.k * r * np.sin(psi - self.phi)
        else:
            # Traitement vectorisé avec matrice de connectivité et retard
            d_idx = (self.ptr - delay_steps) % self.max_buf
            delayed_phi = self.history[d_idx]
            diff = np.sin(delayed_phi[None, :] - self.phi[:, None])
            coupling = np.sum(self.C * diff, axis=1)
            dphi = self.omega + (self.k / self.n) * coupling

        self.dphi = dphi

        # Mise à jour des phases
        self.phi = (self.phi + dphi * self.dt + np.pi) % (2 * np.pi) - np.pi

        # Enregistrement dans le buffer circulaire
        self.ptr = (self.ptr + 1) % self.max_buf
        self.history[self.ptr] = self.phi

        self.time += self.dt
        self.time_history.append(self.time)
        self.r_history.append(r)
        return r, psi, xb, yb

    def _get_zone_name(self, x, y):
        """
        Identifie la zone de stabilité selon Yeung & Strogatz (PRL 1999) :
        x = omega_0 * tau, y = K / omega_0.
        """
        if y > 1.0:
            return "Zone Stable (Synchronie stable pour tout retard)"

        # Vérification langue noire (incohérence stable, synchronie instable)
        for m in range(1, 6):
            y_max_b = 1.0 / (2.0 * (2.0 * m - 1.0))
            if y <= y_max_b:
                x_left = (4.0 * m - 3.0) * np.pi / (2.0 - 2.0 * y)
                x_right = (4.0 * m - 1.0) * np.pi / (2.0 + 2.0 * y)
                if x_left <= x <= x_right:
                    return "Langue Noire (Synchronie instable, incohérence stable)"

        # Vérification langue grise (régime bistable)
        for m in range(1, 6):
            y_max_g = 1.0 / (2.0 * m - 1.0)
            if y <= y_max_g:
                x_left = (4.0 * m - 3.0) * np.pi / (2.0 - y)
                x_right = (4.0 * m - 1.0) * np.pi / (2.0 + y)
                if x_left <= x <= x_right:
                    return "Langue Grise (Bistable : synchronie ou incohérence selon CI)"

        return "Zone Stable (Synchronie stable)"

    def _get_point_size(self):
        if self.n <= 20:
            return 70
        elif self.n <= 60:
            return 50
        elif self.n <= 150:
            return 36
        elif self.n <= 350:
            return 22
        elif self.n <= 700:
            return 14
        else:
            return 8

    def _get_colors(self):
        """
        Calcule les couleurs des oscillateurs correspondant à leur fréquence :
        - Si gamma > 0 : fréquence intrinsèque omega_i (identité spectrale).
        - Si gamma == 0 : vitesse angulaire instantanée dphi/dt (permet de visualiser les clusters).
        """
        if self.gamma > 0.0:
            norm = (self.omega - self.omega_0) / (2.5 * self.gamma + 1e-5)
        else:
            spread = max(0.15, self.k * 0.35)
            norm = (self.dphi - self.omega_0) / spread
        norm = np.clip(norm, -1.0, 1.0)
        return plt.cm.coolwarm((norm + 1.0) / 2.0)

    def _create_scatter(self):
        colors = self._get_colors()
        size = self._get_point_size()
        return self.ax_circle.scatter(
            np.cos(self.phi),
            np.sin(self.phi),
            c=colors,
            s=size,
            edgecolors='#0f172a' if self.n <= 300 else 'none',
            linewidths=0.6,
            zorder=3
        )

    def _rebuild_scatter(self):
        if hasattr(self, 'scat') and self.scat is not None:
            self.scat.remove()
        self.scat = self._create_scatter()
        self.fig.canvas.draw_idle()

    def _init_gui(self):
        # Palette sombre et moderne
        self.c_bg = '#0b0f19'
        self.c_card = '#151c2c'
        self.c_text = '#e2e8f0'
        self.c_muted = '#64748b'
        self.c_accent = '#38bdf8'
        self.c_sync = '#22c55e'
        self.c_desync = '#ef4444'
        self.c_gold = '#f59e0b'
        self.c_purple = '#a855f7'
        self.c_rose = '#f43f5e'

        self.fig = plt.figure(figsize=(15.5, 8.8), facecolor=self.c_bg)
        self.fig.canvas.manager.set_window_title('KURA - Modèle de Kuramoto avec Retard et Connectome Cérébral')

        # 1. Panneau gauche : Cercle des phases theta_i
        self.ax_circle = self.fig.add_axes([0.04, 0.30, 0.28, 0.62], facecolor=self.c_card)
        self.ax_circle.set_aspect('equal')
        self.ax_circle.set_xlim(-1.45, 1.45)
        self.ax_circle.set_ylim(-1.45, 1.45)
        self.ax_circle.axis('off')

        # Cercle unité
        angles = np.linspace(0, 2 * np.pi, 300)
        self.ax_circle.plot(np.cos(angles), np.sin(angles), color='#334155', linestyle='-', linewidth=2, zorder=1)
        self.ax_circle.axhline(0, color='#1e293b', linestyle=':', linewidth=1, zorder=0)
        self.ax_circle.axvline(0, color='#1e293b', linestyle=':', linewidth=1, zorder=0)

        # Cercle barycentre et flèche
        self.r_circle, = self.ax_circle.plot([], [], color=self.c_gold, linestyle='--', linewidth=1.2, alpha=0.5, zorder=2)
        self.vec_line, = self.ax_circle.plot([], [], color=self.c_gold, linewidth=2.5, zorder=4)
        self.bary_pt, = self.ax_circle.plot([], [], 'o', color=self.c_gold, markersize=10, markeredgecolor='#ffffff', markeredgewidth=1.5, zorder=5)

        # Oscillateurs
        self.scat = self._create_scatter()
        self.ax_circle.text(
            0, 1.32, r"Cercle des phases $\theta_i$",
            ha='center', va='center', color=self.c_text, fontsize=12, fontweight='bold'
        )
        self.ax_circle.text(
            0, -1.30, "Couleur : fréquence (Bleu = lente, Rouge = rapide)",
            ha='center', va='center', color=self.c_muted, fontsize=8.5
        )

        # 2. Panneau central : graphique de r(t)
        self.ax_graph = self.fig.add_axes([0.36, 0.46, 0.29, 0.46], facecolor=self.c_card)
        self.ax_graph.set_title("Évolution du paramètre d'ordre $r(t)$", color=self.c_text, fontsize=11, pad=10, fontweight='bold')
        self.ax_graph.set_xlabel("Temps écoulé (s)", color=self.c_muted, fontsize=9)
        self.ax_graph.set_ylabel("Cohérence $r$", color=self.c_muted, fontsize=9)
        self.ax_graph.set_ylim(-0.02, 1.05)
        self.ax_graph.tick_params(colors=self.c_muted, labelsize=8)
        for spine in self.ax_graph.spines.values():
            spine.set_color('#334155')
        self.ax_graph.grid(True, color='#1e293b', linestyle='--', alpha=0.7)

        self.line_r, = self.ax_graph.plot([], [], color=self.c_accent, linewidth=2, label="$r(t)$ instantané")
        self.line_r_inf = self.ax_graph.axhline(0, color=self.c_sync, linestyle=':', linewidth=1.8, label=r"$r_\infty$ théorique (sans retard)")
        self.ax_graph.legend(loc='lower right', facecolor=self.c_card, edgecolor='#334155', labelcolor=self.c_text, fontsize=8)

        # 3. Panneau droit : Diagramme de synchronisation avec retard (Yeung & Strogatz PRL 1999)
        # Verrouillé en déplacement et interactif au clic/glisser
        self.ax_diag = self.fig.add_axes([0.69, 0.46, 0.28, 0.46], facecolor=self.c_card)
        self.ax_diag.set_navigate(False)
        self.ax_diag.set_autoscale_on(False)
        self.ax_diag.set_title("Zones de synchronisation (Yeung & Strogatz 1999)", color=self.c_text, fontsize=11, pad=10, fontweight='bold')

        # Tracé des langues avec le thème sombre (Zone blanche = fond, grise = ardoise, noire = noir)
        for m in range(1, 6):
            # Langue grise (bistabilité)
            y_max_g = 1.0 / (2.0 * m - 1.0)
            y_g = np.linspace(0.0, y_max_g, 200)
            x_left_g = (4.0 * m - 3.0) * np.pi / (2.0 - y_g)
            x_right_g = (4.0 * m - 1.0) * np.pi / (2.0 + y_g)
            self.ax_diag.fill_betweenx(y_g, x_left_g, x_right_g, color='#475569', edgecolor='#94a3b8', linewidth=1.1, zorder=2)

            # Langue noire (synchronie instable)
            y_max_b = 1.0 / (2.0 * (2.0 * m - 1.0))
            y_b = np.linspace(0.0, y_max_b, 200)
            x_left_b = (4.0 * m - 3.0) * np.pi / (2.0 - 2.0 * y_b)
            x_right_b = (4.0 * m - 1.0) * np.pi / (2.0 + 2.0 * y_b)
            self.ax_diag.fill_betweenx(y_b, x_left_b, x_right_b, color='#020617', edgecolor='#64748b', linewidth=1.1, zorder=3)

        self.ax_diag.set_xlim(0.0, 10.4 * np.pi)
        self.ax_diag.set_ylim(0.0, 1.05)

        # Échelle exacte de la photo : ticks intérieurs, labels identiques
        major_ticks_y = [0.0, 0.25, 0.5, 0.75, 1.0]
        major_labels_y = ['0', '1/4', '1/2', '3/4', '1']
        self.ax_diag.set_yticks(major_ticks_y)
        self.ax_diag.set_yticklabels(major_labels_y, fontsize=9.5, color=self.c_text)

        major_ticks_x = [0, 2 * np.pi, 4 * np.pi, 6 * np.pi, 8 * np.pi, 10 * np.pi]
        major_labels_x = ['0', r'$2\pi$', r'$4\pi$', r'$6\pi$', r'$8\pi$', r'$10\pi$']
        minor_ticks_x = [k * np.pi for k in range(1, 11, 2)]

        self.ax_diag.set_xticks(major_ticks_x)
        self.ax_diag.set_xticklabels(major_labels_x, fontsize=9.5, color=self.c_text)
        self.ax_diag.set_xticks(minor_ticks_x, minor=True)

        self.ax_diag.tick_params(axis='both', which='both', colors=self.c_muted, direction='in')
        self.ax_diag.tick_params(axis='both', which='major', length=6, width=1.2)
        self.ax_diag.tick_params(axis='x', which='minor', length=3.5, width=0.9)

        self.ax_diag.spines['left'].set_color('#475569')
        self.ax_diag.spines['left'].set_linewidth(1.4)
        self.ax_diag.spines['bottom'].set_color('#475569')
        self.ax_diag.spines['bottom'].set_linewidth(1.4)
        self.ax_diag.spines['top'].set_visible(False)
        self.ax_diag.spines['right'].set_visible(False)

        self.ax_diag.text(10.45 * np.pi, 0.0, r'$\omega_0 \tau$', fontsize=11, color=self.c_text, va='bottom', ha='left')
        self.ax_diag.text(0.0, 1.07, r'$K / \omega_0$', fontsize=11, color=self.c_text, va='bottom', ha='left')
        self.ax_diag.grid(True, color='#1e293b', linestyle='--', alpha=0.5, zorder=0)

        # Légende du diagramme
        legend_card = mpatches.Patch(facecolor=self.c_card, edgecolor='#64748b', label='Synchronie stable')
        legend_grey = mpatches.Patch(facecolor='#475569', edgecolor='#94a3b8', label='Bistable')
        legend_black = mpatches.Patch(facecolor='#020617', edgecolor='#64748b', label='Incohérence stable')
        marker_proxy = plt.Line2D([0], [0], marker='o', color='none', markerfacecolor=self.c_desync, markeredgecolor='#ffffff', markersize=6, label='Position courante (Cliquer pour régler)')

        self.ax_diag.legend(
            handles=[legend_card, legend_grey, legend_black, marker_proxy],
            loc='upper right',
            fontsize=7.5,
            framealpha=0.9,
            facecolor=self.c_card,
            edgecolor='#334155',
            labelcolor=self.c_text
        )

        # Point marqueur dynamique sur le diagramme
        init_x = self.w0_tau
        init_y = min(1.02, self.k_over_w0)
        self.diag_marker, = self.ax_diag.plot(
            [init_x], [init_y],
            'o',
            color=self.c_desync,
            markeredgecolor='#ffffff',
            markeredgewidth=1.4,
            markersize=8,
            zorder=10
        )

        # Connexion de l'interaction clic/glisser sur le diagramme pour fixer K/w0 et w0*tau
        self.fig.canvas.mpl_connect('button_press_event', self._on_diag_mouse)
        self.fig.canvas.mpl_connect('motion_notify_event', self._on_diag_mouse)

        # 4. Panneau d'informations textuel sous les graphiques
        self.ax_info = self.fig.add_axes([0.36, 0.30, 0.61, 0.11], facecolor=self.c_card)
        self.ax_info.axis('off')
        for spine in self.ax_info.spines.values():
            spine.set_color('#334155')
        self.info_text = self.ax_info.text(
            0.03, 0.5, "",
            va='center', ha='left', color=self.c_text, fontsize=9.5, linespacing=1.6
        )

        # 5. Contrôles en bas : 5 curseurs (3 à gauche, 2 au centre) et 2 boutons
        # Colonne 1 : Pulsation w0, Couplage adimensionné K / w0, Dispersion gamma
        ax_slider_w0 = self.fig.add_axes([0.10, 0.17, 0.21, 0.035], facecolor=self.c_card)
        self.slider_w0 = Slider(
            ax=ax_slider_w0,
            label='Pulsation ω₀ ',
            valmin=0.0,
            valmax=10.0,
            valinit=self.omega_0,
            valstep=0.05,
            color=self.c_rose,
            initcolor='none'
        )
        self._style_widget(self.slider_w0, ax_slider_w0)
        self.slider_w0.on_changed(self._on_w0_changed)

        ax_slider_k = self.fig.add_axes([0.10, 0.10, 0.21, 0.035], facecolor=self.c_card)
        self.slider_k_over_w0 = Slider(
            ax=ax_slider_k,
            label='Couplage K / ω₀ ',
            valmin=0.0,
            valmax=1.0,
            valinit=self.k_over_w0,
            valstep=0.01,
            color=self.c_accent,
            initcolor='none'
        )
        self._style_widget(self.slider_k_over_w0, ax_slider_k)
        self.slider_k_over_w0.on_changed(self._on_k_over_w0_changed)

        ax_slider_gamma = self.fig.add_axes([0.10, 0.03, 0.21, 0.035], facecolor=self.c_card)
        self.slider_gamma = Slider(
            ax=ax_slider_gamma,
            label='Dispersion γ ',
            valmin=0.0,
            valmax=2.0,
            valinit=self.gamma,
            valstep=0.02,
            color=self.c_gold,
            initcolor='none'
        )
        self._style_widget(self.slider_gamma, ax_slider_gamma)
        self.slider_gamma.on_changed(self._on_gamma_changed)

        # Colonne 2 : Oscillateurs N et Retard adimensionné w0 * tau (0 à 10 pi)
        ax_slider_n = self.fig.add_axes([0.44, 0.14, 0.21, 0.035], facecolor=self.c_card)
        self.slider_n = Slider(
            ax=ax_slider_n,
            label='Oscillateurs N ',
            valmin=1,
            valmax=1000,
            valinit=self.n,
            valstep=1,
            color=self.c_sync,
            initcolor='none'
        )
        self._style_widget(self.slider_n, ax_slider_n)
        self.slider_n.on_changed(self._on_n_changed)

        ax_slider_tau = self.fig.add_axes([0.44, 0.06, 0.21, 0.035], facecolor=self.c_card)
        self.slider_w0_tau = Slider(
            ax=ax_slider_tau,
            label='Retard ω₀·τ ',
            valmin=0.0,
            valmax=10.0 * np.pi,
            valinit=self.w0_tau,
            valstep=0.05,
            color=self.c_purple,
            initcolor='none'
        )
        self._style_widget(self.slider_w0_tau, ax_slider_tau)
        self.slider_w0_tau.on_changed(self._on_w0_tau_changed)

        # Boutons Pause et Reset
        ax_btn_pause = self.fig.add_axes([0.72, 0.08, 0.11, 0.08], facecolor=self.c_card)
        self.btn_pause = Button(ax_btn_pause, 'Pause', color=self.c_card, hovercolor='#334155')
        self._style_button(self.btn_pause, self.c_accent)
        self.btn_pause.on_clicked(self._toggle_pause)

        ax_btn_reset = self.fig.add_axes([0.85, 0.08, 0.11, 0.08], facecolor=self.c_card)
        self.btn_reset = Button(ax_btn_reset, 'Reset', color=self.c_card, hovercolor='#334155')
        self._style_button(self.btn_reset, self.c_gold)
        self.btn_reset.on_clicked(self._reset)

        self._update_diag_marker()

    def _style_widget(self, slider, ax):
        ax.tick_params(colors=self.c_muted, labelsize=8)
        slider.label.set_color(self.c_text)
        slider.label.set_fontsize(9.5)
        slider.valtext.set_color(self.c_accent)
        slider.valtext.set_fontsize(9.5)
        for spine in ax.spines.values():
            spine.set_color('#334155')

    def _style_button(self, btn, text_color):
        btn.label.set_color(text_color)
        btn.label.set_fontsize(10)
        btn.label.set_fontweight('bold')
        for spine in btn.ax.spines.values():
            spine.set_color('#334155')

    def _on_diag_mouse(self, event):
        """Permet de cliquer ou glisser sur l'image pour fixer K/w0 et w0*tau."""
        if event.inaxes == self.ax_diag and event.button == 1:
            if event.xdata is not None and event.ydata is not None:
                x_target = float(np.clip(event.xdata, 0.0, 10.0 * np.pi))
                y_target = float(np.clip(event.ydata, 0.0, 1.0))
                # Met à jour les curseurs qui synchronisent automatiquement la physique
                self.slider_w0_tau.set_val(x_target)
                self.slider_k_over_w0.set_val(y_target)

    def _update_diag_marker(self):
        """Met à jour le point de position sur le diagramme de Yeung & Strogatz."""
        x = self.w0_tau
        y = self.k_over_w0

        y_plot = min(1.02, y)
        marker_style = '^' if y > 1.0 else 'o'
        self.diag_marker.set_data([x], [y_plot])
        self.diag_marker.set_marker(marker_style)

        zone_name = self._get_zone_name(x, y)
        zone_short = "Synchronie" if "Stable" in zone_name and "Langue" not in zone_name else "Bistable" if "Grise" in zone_name else "Instable"
        self.ax_diag.set_title(f"Zones de synchronisation ({zone_short})", color=self.c_text, fontsize=11, pad=10, fontweight='bold')
        self.fig.canvas.draw_idle()

    def _on_w0_changed(self, val):
        self.omega_0 = float(val)
        # Met à jour K et tau physiques associés aux coordonnées adimensionnées
        self.k = self.k_over_w0 * self.omega_0
        self.tau = (self.w0_tau / self.omega_0) if self.omega_0 > 1e-5 else 0.0
        self.delay_steps = int(round(self.tau / self.dt))

        # Recalcul des fréquences intrinsèques sans réinitialiser les phases
        if self.gamma <= 0.0:
            self.omega = np.full(self.n, self.omega_0)
        else:
            self.omega = self.omega_0 + self.delta_omega
        self._update_diag_marker()

    def _on_k_over_w0_changed(self, val):
        self.k_over_w0 = float(val)
        self.k = self.k_over_w0 * self.omega_0
        self._update_diag_marker()

    def _on_gamma_changed(self, val):
        self.gamma = float(val)
        self._init_physics()
        self._rebuild_scatter()

    def _on_n_changed(self, val):
        new_n = int(round(val))
        if new_n != self.n:
            self.n = new_n
            self._init_physics()
            self._rebuild_scatter()

    def _on_w0_tau_changed(self, val):
        self.w0_tau = float(val)
        self.tau = (self.w0_tau / self.omega_0) if self.omega_0 > 1e-5 else 0.0
        self.delay_steps = int(round(self.tau / self.dt))
        # Forcer le nombre d'oscillateurs à 90 lors du déplacement du retard
        if self.n != 90:
            self.slider_n.set_val(90)
        self._update_diag_marker()

    def _toggle_pause(self, event=None):
        self.is_running = not self.is_running
        self.btn_pause.label.set_text('Reprendre' if not self.is_running else 'Pause')
        if self.is_running:
            self.last_real_time = None
            self.time_accumulator = 0.0
        self.fig.canvas.draw_idle()

    def _reset(self, event=None):
        self._init_physics()
        self._rebuild_scatter()
        self._update_diag_marker()

    def update(self, frame):
        if self.is_running:
            now = time.perf_counter()
            if self.last_real_time is None:
                self.last_real_time = now
                real_elapsed = 0.0
            else:
                real_elapsed = now - self.last_real_time
                self.last_real_time = now

            # Plafond de sauvegarde pour éviter un emballement
            real_elapsed = min(real_elapsed, 0.08)
            self.time_accumulator += real_elapsed

            stepped = False
            while self.time_accumulator >= self.dt:
                r, psi, xb, yb = self.step()
                self.time_accumulator -= self.dt
                stepped = True

            if not stepped:
                r, psi, xb, yb = self.compute_order_parameter()
        else:
            r, psi, xb, yb = self.compute_order_parameter()

        # Mise à jour des oscillateurs sur le cercle (positions et couleurs selon fréquence)
        x = np.cos(self.phi)
        y = np.sin(self.phi)
        self.scat.set_offsets(np.column_stack((x, y)))
        self.scat.set_facecolor(self._get_colors())

        # Cercle barycentrique de rayon r
        angles = np.linspace(0, 2 * np.pi, 100)
        self.r_circle.set_data(r * np.cos(angles), r * np.sin(angles))

        # Vecteur et point barycentre
        self.vec_line.set_data([0, xb], [0, yb])
        self.bary_pt.set_data([xb], [yb])

        # Mise à jour graphique r(t)
        if len(self.time_history) > 1:
            t_data = list(self.time_history)
            r_data = list(self.r_history)
            self.line_r.set_data(t_data, r_data)
            self.ax_graph.set_xlim(max(0.0, t_data[-1] - 8.0), max(8.0, t_data[-1] + 0.5))

        # Ligne théorique sans retard si applicable
        kc = 2.0 * self.gamma
        if self.tau == 0.0 and self.gamma > 0.0 and self.k > kc:
            r_inf = np.sqrt(max(0.0, 1.0 - (kc / self.k)))
            self.line_r_inf.set_ydata([r_inf, r_inf])
            self.line_r_inf.set_visible(True)
        else:
            self.line_r_inf.set_visible(False)

        # Panneau texte dynamique
        self.r = r
        if r >= 0.7:
            status_label = "SYNCHRONISÉ"
        elif r >= 0.35:
            status_label = "TRANSITION DE PHASE / BISTABLE"
        else:
            status_label = "DÉSYNCHRONISÉ / INCOHÉRENT"

        x_diag = self.w0_tau
        y_diag = self.k_over_w0
        zone_name = self._get_zone_name(x_diag, y_diag)
        network_info = "Connectome AAL (90 régions HCP)" if self.n == 90 else f"Tous-à-tous ({self.n} oscillateurs)"
        period_str = f"T = {2.0 * np.pi / self.omega_0:.2f} s" if self.omega_0 > 0.05 else "T = ∞ (statique)"

        info_str = (
            f"Régime : {status_label}   |   Réseau : {network_info}\n"
            f"Paramètres réduits : K/ω₀ = {self.k_over_w0:.2f}   |   ω₀·τ = {self.w0_tau / np.pi:.2f}π   |   Position : {zone_name}\n"
            f"Grandeurs physiques : ω₀ = {self.omega_0:.2f} rad/s ({period_str})   |   K = {self.k:.2f}   |   τ = {self.tau:.2f} s   |   Cohérence r(t) = {r:.3f}"
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
            interval=20,
            blit=False,
            cache_frame_data=False
        )
        plt.show()


def main():
    app = KuramotoInteractiveApp(
        n_oscillators=90,
        gamma=0.0,
        omega_0=1.0,
        k_over_w0=0.8,
        w0_tau=0.0,
        dt=0.005
    )
    app.run()


if __name__ == '__main__':
    main()
