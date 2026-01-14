"""
Vergleich: Empirische Formeln vs. Shomate-basierte Berechnung
Zeigt den Unterschied in der Kruemmung der Phasengrenzen.

Ohne Aenderung am Hauptcode - nur Visualisierung.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================================
# SHOMATE-KOEFFIZIENTEN (kopiert aus thermodynamics.py)
# ============================================================================

R = 8.314462618  # J/(mol*K)

SHOMATE_DATA = {
    'Fe': {
        'dHf_298': 0.0,
        'ranges': [
            {'T_min': 298, 'T_max': 700,
             'A': 18.42868, 'B': 24.64301, 'C': -8.913720, 'D': 9.664706,
             'E': -0.012643, 'F': -6.573022, 'G': 42.51488, 'H': 0.0},
            {'T_min': 700, 'T_max': 1042,
             'A': -57767.65, 'B': 137919.7, 'C': -122773.2, 'D': 38682.42,
             'E': 3993.080, 'F': 24078.67, 'G': -87364.01, 'H': 0.0},
        ]
    },
    'FeO': {
        'dHf_298': -272.04,
        'ranges': [
            {'T_min': 298, 'T_max': 1650,
             'A': 45.75120, 'B': 18.78553, 'C': -5.952201, 'D': 0.852779,
             'E': -0.081265, 'F': -286.7429, 'G': 110.3120, 'H': -272.0441},
        ]
    },
    'Fe3O4': {
        'dHf_298': -1120.89,
        'ranges': [
            {'T_min': 298, 'T_max': 900,
             'A': 104.2096, 'B': 178.5108, 'C': 10.61510, 'D': 1.132534,
             'E': -0.994202, 'F': -1163.336, 'G': 212.0585, 'H': -1120.894},
            {'T_min': 900, 'T_max': 3000,
             'A': 200.8320, 'B': 1.586435e-7, 'C': -6.661682e-8, 'D': 9.452452e-9,
             'E': 3.186020e-8, 'F': -1174.135, 'G': 388.0790, 'H': -1120.894},
        ]
    },
    'H2': {
        'dHf_298': 0.0,
        'ranges': [
            {'T_min': 298, 'T_max': 1000,
             'A': 33.066178, 'B': -11.363417, 'C': 11.432816, 'D': -2.772874,
             'E': -0.158558, 'F': -9.980797, 'G': 172.707974, 'H': 0.0},
            {'T_min': 1000, 'T_max': 2500,
             'A': 18.563083, 'B': 12.257357, 'C': -2.859786, 'D': 0.268238,
             'E': 1.977990, 'F': -1.147438, 'G': 156.288133, 'H': 0.0},
        ]
    },
    'H2O': {
        'dHf_298': -241.83,
        'ranges': [
            {'T_min': 500, 'T_max': 1700,
             'A': 30.09200, 'B': 6.832514, 'C': 6.793435, 'D': -2.534480,
             'E': 0.082139, 'F': -250.8810, 'G': 223.3967, 'H': -241.8264},
        ]
    },
    'CO': {
        'dHf_298': -110.53,
        'ranges': [
            {'T_min': 298, 'T_max': 1300,
             'A': 25.56759, 'B': 6.096130, 'C': 4.054656, 'D': -2.671301,
             'E': 0.131021, 'F': -118.0089, 'G': 227.3665, 'H': -110.5271},
        ]
    },
    'CO2': {
        'dHf_298': -393.52,
        'ranges': [
            {'T_min': 298, 'T_max': 1200,
             'A': 24.99735, 'B': 55.18696, 'C': -33.69137, 'D': 7.948387,
             'E': -0.136638, 'F': -403.6075, 'G': 228.2431, 'H': -393.5224},
        ]
    },
}


def get_coeffs(species, T_K):
    """Holt Shomate-Koeffizienten fuer Temperatur."""
    data = SHOMATE_DATA[species]
    for r in data['ranges']:
        if r['T_min'] <= T_K <= r['T_max']:
            return r
    # Fallback
    if T_K < data['ranges'][0]['T_min']:
        return data['ranges'][0]
    return data['ranges'][-1]


def calc_H(species, T_K):
    """Enthalpie H(T) in kJ/mol."""
    c = get_coeffs(species, T_K)
    t = T_K / 1000.0
    return c['A']*t + c['B']*t**2/2 + c['C']*t**3/3 + c['D']*t**4/4 - c['E']/t + c['F']


def calc_S(species, T_K):
    """Entropie S(T) in J/(mol*K)."""
    c = get_coeffs(species, T_K)
    t = T_K / 1000.0
    return c['A']*np.log(t) + c['B']*t + c['C']*t**2/2 + c['D']*t**3/3 - c['E']/(2*t**2) + c['G']


def calc_G(species, T_K):
    """Gibbs-Energie G(T) in kJ/mol."""
    H = calc_H(species, T_K)
    S = calc_S(species, T_K)
    return H - T_K * S / 1000.0


# ============================================================================
# REAKTIONS-BERECHNUNGEN
# ============================================================================

def delta_G_reaction_shomate(reaction, T_K):
    """
    Berechnet Delta_G der Reaktion aus Shomate-Daten.
    Gibt Wert in J/mol zurueck.
    """
    if reaction == 'Fe3O4_Fe_H2':
        # 1/4 Fe3O4 + H2 = 3/4 Fe + H2O
        # Delta_G = 3/4*G(Fe) + G(H2O) - 1/4*G(Fe3O4) - G(H2)
        dG = 0.75*calc_G('Fe', T_K) + calc_G('H2O', T_K) - 0.25*calc_G('Fe3O4', T_K) - calc_G('H2', T_K)
    elif reaction == 'Fe3O4_Fe_CO':
        # 1/4 Fe3O4 + CO = 3/4 Fe + CO2
        dG = 0.75*calc_G('Fe', T_K) + calc_G('CO2', T_K) - 0.25*calc_G('Fe3O4', T_K) - calc_G('CO', T_K)
    elif reaction == 'FeO_Fe_H2':
        # FeO + H2 = Fe + H2O
        dG = calc_G('Fe', T_K) + calc_G('H2O', T_K) - calc_G('FeO', T_K) - calc_G('H2', T_K)
    elif reaction == 'FeO_Fe_CO':
        # FeO + CO = Fe + CO2
        dG = calc_G('Fe', T_K) + calc_G('CO2', T_K) - calc_G('FeO', T_K) - calc_G('CO', T_K)
    elif reaction == 'Fe3O4_FeO_H2':
        # 1/3 Fe3O4 + H2 = FeO + H2O  (normiert auf 1 mol H2)
        # Original: Fe3O4 + H2 = 3 FeO + H2O
        dG = calc_G('FeO', T_K) + (1/3)*calc_G('H2O', T_K) - (1/3)*calc_G('Fe3O4', T_K) - (1/3)*calc_G('H2', T_K)
        dG = dG * 3  # Fuer 1 mol H2
    elif reaction == 'Fe3O4_FeO_CO':
        # Fe3O4 + CO = 3 FeO + CO2
        dG = 3*calc_G('FeO', T_K) + calc_G('CO2', T_K) - calc_G('Fe3O4', T_K) - calc_G('CO', T_K)
    else:
        return np.nan

    return dG * 1000  # kJ -> J


def calc_K_shomate(reaction, T_K):
    """Gleichgewichtskonstante aus Shomate-Daten."""
    dG = delta_G_reaction_shomate(reaction, T_K)
    if np.isnan(dG):
        return np.nan
    return np.exp(-dG / (R * T_K))


def calc_pct_reducing_shomate(reaction, T_K):
    """Berechnet %H2 oder %CO bei Gleichgewicht aus Shomate."""
    K = calc_K_shomate(reaction, T_K)
    if np.isnan(K) or K <= 0:
        return np.nan
    # K = p_oxidized / p_reduced
    return 100.0 / (1.0 + K)


# ============================================================================
# EMPIRISCHE FORMELN (aus thermodynamics.py)
# ============================================================================

EMPIRICAL = {
    'Fe3O4_Fe_H2': {'A': -1071, 'B': 0.7460},
    'Fe3O4_Fe_CO': {'A': 0, 'B': -0.0052},
    'FeO_Fe_H2': {'A': -870, 'B': 0.5076},
    'FeO_Fe_CO': {'A': 919, 'B': -1.0895},
    'Fe3O4_FeO_H2': {'A': -3692, 'B': 3.8543},
    'Fe3O4_FeO_CO': {'A': -1881, 'B': 2.2305},
}


def calc_pct_reducing_empirical(reaction, T_K):
    """Berechnet %H2 oder %CO mit empirischer Formel."""
    if reaction not in EMPIRICAL:
        return np.nan
    A = EMPIRICAL[reaction]['A']
    B = EMPIRICAL[reaction]['B']
    log10_K = A / T_K + B
    K = 10 ** log10_K
    return 100.0 / (1.0 + K)


# ============================================================================
# VISUALISIERUNG
# ============================================================================

def create_comparison_plot():
    """Erstellt Vergleichsplot: Empirisch vs. Shomate."""

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle('Vergleich: Empirische Formeln (blau) vs. Shomate-Berechnung (rot)\n'
                 'Kruemmung entsteht durch temperaturabhaengige Delta_H und Delta_S',
                 fontsize=13, fontweight='bold')

    # Temperaturbereiche
    T_low_C = np.linspace(300, 570, 100)
    T_high_C = np.linspace(580, 1000, 100)
    T_low_K = T_low_C + 273.15
    T_high_K = T_high_C + 273.15

    reactions = [
        ('Fe3O4_Fe_H2', 'Fe3O4/Fe + H2 (< 570C)', T_low_C, T_low_K, axes[0, 0]),
        ('Fe3O4_Fe_CO', 'Fe3O4/Fe + CO (< 570C)', T_low_C, T_low_K, axes[0, 1]),
        ('FeO_Fe_H2', 'Fe/FeO + H2 (> 570C)', T_high_C, T_high_K, axes[0, 2]),
        ('FeO_Fe_CO', 'Fe/FeO + CO (> 570C)', T_high_C, T_high_K, axes[1, 0]),
        ('Fe3O4_FeO_H2', 'FeO/Fe3O4 + H2 (> 570C)', T_high_C, T_high_K, axes[1, 1]),
        ('Fe3O4_FeO_CO', 'FeO/Fe3O4 + CO (> 570C)', T_high_C, T_high_K, axes[1, 2]),
    ]

    for rxn, title, T_C, T_K, ax in reactions:
        # Empirische Berechnung
        pct_emp = [calc_pct_reducing_empirical(rxn, T) for T in T_K]
        god_emp = [1 - p/100 if not np.isnan(p) else np.nan for p in pct_emp]

        # Shomate-Berechnung
        pct_sho = [calc_pct_reducing_shomate(rxn, T) for T in T_K]
        god_sho = [1 - p/100 if not np.isnan(p) else np.nan for p in pct_sho]

        ax.plot(T_C, god_emp, 'b-', linewidth=2.5, label='Empirisch (A/T + B)')
        ax.plot(T_C, god_sho, 'r--', linewidth=2.5, label='Shomate (exakt)')

        ax.set_xlabel('Temperatur (C)')
        ax.set_ylabel('GOD')
        ax.set_title(title)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('shomate_vs_empirical.png', dpi=150, bbox_inches='tight')
    print("Gespeichert: shomate_vs_empirical.png")
    plt.close()


def create_full_diagram_comparison():
    """Erstellt vollstaendiges Baur-Glaessner Diagramm mit beiden Methoden."""

    fig, ax = plt.subplots(figsize=(14, 10))

    T_low = np.linspace(300, 570, 100)
    T_high = np.linspace(575, 1000, 150)

    # =========================================================================
    # H2-System
    # =========================================================================

    # Fe3O4/Fe (< 570C) - HIER IST DIE KRUEMMUNG!
    god_emp = [1 - calc_pct_reducing_empirical('Fe3O4_Fe_H2', T+273.15)/100 for T in T_low]
    god_sho = [1 - calc_pct_reducing_shomate('Fe3O4_Fe_H2', T+273.15)/100 for T in T_low]

    ax.plot(T_low, god_emp, 'b--', linewidth=2, label='Fe3O4/Fe H2 - Empirisch')
    ax.plot(T_low, god_sho, 'r--', linewidth=2, label='Fe3O4/Fe H2 - Shomate')

    # Fe/FeO (> 570C)
    god_emp = [1 - calc_pct_reducing_empirical('FeO_Fe_H2', T+273.15)/100 for T in T_high]
    god_sho = [1 - calc_pct_reducing_shomate('FeO_Fe_H2', T+273.15)/100 for T in T_high]

    ax.plot(T_high, god_emp, 'b--', linewidth=2, alpha=0.7)
    ax.plot(T_high, god_sho, 'r--', linewidth=2, alpha=0.7)

    # =========================================================================
    # CO-System
    # =========================================================================

    # Fe3O4/Fe (< 570C)
    god_emp = [1 - calc_pct_reducing_empirical('Fe3O4_Fe_CO', T+273.15)/100 for T in T_low]
    god_sho = [1 - calc_pct_reducing_shomate('Fe3O4_Fe_CO', T+273.15)/100 for T in T_low]

    ax.plot(T_low, god_emp, 'b-', linewidth=2, label='Fe3O4/Fe CO - Empirisch')
    ax.plot(T_low, god_sho, 'r-', linewidth=2, label='Fe3O4/Fe CO - Shomate')

    # Fe/FeO (> 570C)
    god_emp = [1 - calc_pct_reducing_empirical('FeO_Fe_CO', T+273.15)/100 for T in T_high]
    god_sho = [1 - calc_pct_reducing_shomate('FeO_Fe_CO', T+273.15)/100 for T in T_high]

    ax.plot(T_high, god_emp, 'b-', linewidth=2, alpha=0.7)
    ax.plot(T_high, god_sho, 'r-', linewidth=2, alpha=0.7)

    # =========================================================================
    # Formatierung
    # =========================================================================

    ax.axvline(x=570, color='gray', linestyle=':', linewidth=1.5, alpha=0.7)
    ax.text(575, 0.02, '570C', fontsize=10, color='gray')

    ax.set_xlabel('Temperatur (C)', fontsize=12)
    ax.set_ylabel('GOD', fontsize=12)
    ax.set_title('Baur-Glaessner Diagramm: Empirisch (blau) vs. Shomate (rot)\n'
                 'Die Shomate-Berechnung zeigt die korrekte Kruemmung',
                 fontsize=13, fontweight='bold')
    ax.set_xlim(300, 1000)
    ax.set_ylim(0, 0.6)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Phasenbeschriftungen
    ax.text(400, 0.08, 'Fe', fontsize=14, fontweight='bold')
    ax.text(800, 0.20, 'Fe', fontsize=14, fontweight='bold')
    ax.text(400, 0.45, 'Fe3O4', fontsize=12)

    plt.tight_layout()
    plt.savefig('baur_glassner_shomate_comparison.png', dpi=150, bbox_inches='tight')
    print("Gespeichert: baur_glassner_shomate_comparison.png")
    plt.close()


def print_numerical_comparison():
    """Druckt numerischen Vergleich."""

    print("\n" + "=" * 80)
    print("NUMERISCHER VERGLEICH: Empirisch vs. Shomate")
    print("=" * 80)

    print("\n--- Fe3O4/Fe + H2 (unter 570C) - HIER SIEHT MAN DIE KRUEMMUNG ---")
    print(f"{'T (C)':>8} | {'%H2 Emp':>10} | {'%H2 Shomate':>12} | {'Diff':>8} | {'GOD Emp':>8} | {'GOD Sho':>8}")
    print("-" * 75)

    for T_C in [300, 350, 400, 450, 500, 550, 570]:
        T_K = T_C + 273.15
        emp = calc_pct_reducing_empirical('Fe3O4_Fe_H2', T_K)
        sho = calc_pct_reducing_shomate('Fe3O4_Fe_H2', T_K)
        diff = emp - sho if not (np.isnan(emp) or np.isnan(sho)) else np.nan
        god_e = 1 - emp/100 if not np.isnan(emp) else np.nan
        god_s = 1 - sho/100 if not np.isnan(sho) else np.nan

        print(f"{T_C:>8} | {emp:>10.1f} | {sho:>12.1f} | {diff:>+8.1f} | {god_e:>8.3f} | {god_s:>8.3f}")

    print("\n--- Fe3O4/Fe + CO (unter 570C) ---")
    print(f"{'T (C)':>8} | {'%CO Emp':>10} | {'%CO Shomate':>12} | {'Diff':>8}")
    print("-" * 50)

    for T_C in [300, 350, 400, 450, 500, 550, 570]:
        T_K = T_C + 273.15
        emp = calc_pct_reducing_empirical('Fe3O4_Fe_CO', T_K)
        sho = calc_pct_reducing_shomate('Fe3O4_Fe_CO', T_K)
        diff = emp - sho if not (np.isnan(emp) or np.isnan(sho)) else np.nan

        print(f"{T_C:>8} | {emp:>10.1f} | {sho:>12.1f} | {diff:>+8.1f}")


# ============================================================================
# HAUPTPROGRAMM
# ============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 80)
    print("SHOMATE vs. EMPIRISCHE FORMELN - KRUEMMUNGSVERGLEICH")
    print("=" * 80)

    print_numerical_comparison()

    print("\n\nErstelle Diagramme...")
    create_comparison_plot()
    create_full_diagram_comparison()

    print("\n" + "=" * 80)
    print("ERKLAERUNG DER KRUEMMUNG")
    print("=" * 80)
    print("""
Die EMPIRISCHE Formel:  log10(K) = A/T + B
  -> Impliziert KONSTANTES Delta_H und Delta_S
  -> Ergibt eine fast gerade Linie im GOD-T Diagramm

Die SHOMATE-Berechnung:  Delta_G(T) = Sum[G_Produkte(T)] - Sum[G_Reaktanten(T)]
  -> H(T) und S(T) sind POLYNOME in T
  -> Delta_H(T) und Delta_S(T) aendern sich mit Temperatur
  -> Ergibt die GEKRUEMMTE Linie wie im Referenzbild

Der Unterschied ist besonders sichtbar bei:
  - Fe3O4/Fe + H2 Reaktion (unter 570C)
  - Grosse Temperaturspanne (300-570C)
""")
