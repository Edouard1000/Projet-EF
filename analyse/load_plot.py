import matplotlib.pyplot as plt
import numpy as np

# ==============================================
# == Données de Simulation (h=30) ==
# ==============================================
loads = np.array([
    1e-2, 100.00, 1000.00, 10000.00, 1e5, 1e6, 6250000.00, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12
])

ux_values = np.array([
    -1.090010e-05, -1.090069e-05, -1.090596e-05, -1.095868e-05, -1.031429e-05, -5.042000e-06,
    2.571302e-05, 4.768090e-05, 5.749099e-04, 5.847200e-03, 5.857010e-02, 5.857991e-01, 5.858089e+00
])

uy_values = np.array([
    -3.272980e-04, -3.272904e-04, -3.272221e-04, -3.265390e-04, -3.348873e-04, -4.031919e-04,
    -8.016348e-04, -1.086237e-03, -7.916688e-03, -7.622120e-02, -7.592663e-01, -7.589717e+00, -7.589423e+01
])

# Convertir les déplacements en millimètres
ux_values_mm = ux_values * 1000
uy_values_mm = uy_values * 1000

# ==============================================
# == Fonction pour Créer un Plot ==
# ==============================================
def create_plot(loads, ux_values_mm, uy_values_mm, title_suffix, filename, max_load=None):
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 8))

    # Filtrer les données si max_load est défini
    if max_load:
        mask = loads <= max_load
        loads = loads[mask]
        ux_values_mm = ux_values_mm[mask]
        uy_values_mm = uy_values_mm[mask]

    # Plot des déplacements
    line_uy, = ax.plot(loads, uy_values_mm, marker='o', markersize=8, linestyle='-', linewidth=2.5,
                       color='#1f77b4', label='Déplacement Vertical ($U_y$)')
    line_ux, = ax.plot(loads, ux_values_mm, marker='s', markersize=7, linestyle='--', linewidth=2.0,
                       color='#ff7f0e', label='Déplacement Horizontal ($U_x$)')

    # Configuration de l'échelle logarithmique pour l'axe X
    ax.set_xscale('log')

    # Labels et titre
    ax.set_xlabel('Charge Linéique Verticale Appliquée [N/m]', fontsize=14, labelpad=10)
    ax.set_ylabel("Déplacement du point d'application [mm]", fontsize=14, labelpad=10)
    ax.set_title(f"Réponse Élastique du Rail UIC60 (h=30mm) : Déplacement au Point (19.4, 171.0)\n"
                 f"en Fonction de la Charge Appliquée (Simulation EF Planar Strain, {title_suffix})",
                 fontsize=15, fontweight='bold', pad=20)

    # Grille et personnalisation des ticks
    ax.grid(True, which='major', linestyle='-', linewidth=0.6, color='gray', alpha=0.7)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.4, color='lightgray', alpha=0.5)
    ax.tick_params(axis='both', which='major', direction='in', length=6, width=1, labelsize=12)
    ax.tick_params(axis='both', which='minor', direction='in', length=4, width=0.8)

    # Lignes verticales pour les charges spécifiques
    reference_load = 6250000  # N/m
    ax.axvline(reference_load, color='darkred', linestyle='-.', linewidth=1.5, alpha=0.8,
               label=f'Train charge Réf. ({reference_load:.1e} N/m)')
    freight_train_load = 1.56e7  # N/m
    ax.axvline(freight_train_load, color='green', linestyle='--', linewidth=1.5, alpha=0.8,
               label=f'Train de marchandises ({freight_train_load:.2e} N/m)')
    overloaded_train_load = 3.13e7  # N/m
    ax.axvline(overloaded_train_load, color='purple', linestyle='--', linewidth=1.5, alpha=0.8,
               label=f'Train surchargé ({overloaded_train_load:.2e} N/m)')

    # Légende
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(handles=handles, labels=labels, fontsize=13, frameon=True,
                       framealpha=0.95, shadow=False, borderpad=0.8, loc='best')
    legend.set_title("Légende", prop={'size': 13, 'weight': 'semibold'})

    # Annotations pour les points
    for x, y in zip(loads, uy_values_mm):
        if x != 1e7:  # Exclure le point spécifique
            ax.annotate(f'{y:.2f}', xy=(x, y), xytext=(5, 5), textcoords='offset points',
                        fontsize=10, color='#1f77b4', alpha=0.8)
    for x, y in zip(loads, ux_values_mm):
        if x != 1e7:  # Exclure le point spécifique
            ax.annotate(f'{y:.2f}', xy=(x, y), xytext=(5, -10), textcoords='offset points',
                        fontsize=10, color='#ff7f0e', alpha=0.8)

    # Ajustements finaux
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['bottom'].set_linewidth(1.2)

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')  # Sauvegarder le plot
    plt.show()

# ==============================================
# == Premier Plot : Toutes les Charges ==
# ==============================================
create_plot(loads, ux_values_mm, uy_values_mm, "Toutes les Charges", "rail_load_vs_displacement_P29.png")

# ==============================================
# == Deuxième Plot : Charges jusqu'à 10^8 ==
# ==============================================
create_plot(loads, ux_values_mm, uy_values_mm, "Charges ≤ $10^8$", "rail_load_vs_displacement_P29_limited.png", max_load=1e8)