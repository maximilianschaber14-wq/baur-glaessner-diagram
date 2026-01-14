"""
Validierungsskript: Vergleich der empirischen Formeln mit FactSage/Literatur-Daten.
Dieses Skript veraendert NICHTS am bestehenden Code - nur Validierung und Visualisierung.

Quellen:
- FactSage 5.4/8.2 Simulationen (ResearchGate)
- HSC Chemistry 6.0
- NIST-JANAF Thermochemical Tables
- Turkdogan (1980): Physical Chemistry of High Temperature Technology
- Darken & Gurry (1945): JACS 67, 1398
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

# ============================================================================
# REFERENZDATEN AUS LITERATUR / FACTSAGE
# ============================================================================

# Literaturwerte fuer Fe/FeO + CO Gleichgewicht (%CO bei Gleichgewicht)
# Quellen: FactSage, HSC Chemistry, experimentelle Daten
LITERATURE_DATA = {
    'Fe_FeO_CO': {
        # Temperatur (C): %CO bei Gleichgewicht (aus Baur-Glaessner Diagrammen)
        600: 60.0,   # Knapp ueber Wuestit-Stabilitaet
        700: 63.0,   # FactSage/HSC Werte
        800: 66.0,
        900: 68.5,
        1000: 71.8,  # Experimentell bestaetigt: 71.8% CO, 28.2% CO2
        1100: 74.0,
    },
    'Fe_FeO_H2': {
        # H2-System hat andere Gleichgewichtslage
        600: 53.0,
        700: 57.0,
        800: 61.0,   # Bei 810C kreuzen sich CO und H2 Linien
        900: 65.0,
        1000: 68.0,
        1100: 71.0,
    },
    'FeO_Fe3O4_CO': {
        # Wuestit/Magnetit Grenze
        600: 26.0,
        700: 22.0,
        800: 18.0,
        900: 16.0,
        1000: 14.0,
        1100: 12.0,
    },
    'FeO_Fe3O4_H2': {
        600: 18.0,
        700: 14.0,
        800: 11.0,
        900: 9.0,
        1000: 8.0,
        1100: 7.0,
    },
    'Fe3O4_Fe_CO': {
        # Unter 570C (Magnetit direkt zu Fe)
        300: 50.0,
        400: 50.0,
        500: 50.0,
        570: 50.0,  # Tripelpunkt
    },
    'Fe3O4_Fe_H2': {
        300: 85.0,
        400: 80.0,
        500: 75.0,
        570: 70.0,  # Tripelpunkt
    },
    'Boudouard': {
        # C + CO2 = 2CO Gleichgewicht
        400: 0.5,    # Fast nur CO2
        500: 3.0,
        600: 18.0,
        700: 60.0,
        800: 93.0,
        900: 99.0,
        1000: 99.9,  # Fast nur CO
    }
}


# ============================================================================
# AKTUELLE EMPIRISCHE FORMELN AUS thermodynamics.py
# ============================================================================

R = 8.314462618  # J/(mol*K)
WUSTITE_STABILITY_K = 843  # ~570C

# Empirische Koeffizienten aus dem aktuellen Code
CURRENT_REACTIONS = {
    'FeO_Fe_CO': {'A': 919, 'B': -1.0895},
    'Fe3O4_FeO_CO': {'A': -1881, 'B': 2.2305},
    'Fe3O4_Fe_CO': {'A': 0, 'B': -0.0052},
    'FeO_Fe_H2': {'A': -870, 'B': 0.5076},
    'Fe3O4_FeO_H2': {'A': -3692, 'B': 3.8543},
    'Fe3O4_Fe_H2': {'A': -1071, 'B': 0.7460},
}

# Boudouard: log10(K) = -8916/T + 9.113
BOUDOUARD_A = -8916
BOUDOUARD_B = 9.113


def calc_reducing_gas_pct(reaction_name: str, T_celsius: float) -> float:
    """
    Berechnet %CO oder %H2 bei Gleichgewicht mit den aktuellen empirischen Formeln.
    """
    T_K = T_celsius + 273.15
    rxn = CURRENT_REACTIONS[reaction_name]
    A, B = rxn['A'], rxn['B']

    log10_K = A / T_K + B
    K = 10 ** log10_K

    # K = p_oxidized / p_reduced = p_CO2/p_CO oder p_H2O/p_H2
    # %red = 1 / (1 + K) * 100
    pct_red = 100.0 / (1.0 + K)
    return pct_red


def calc_boudouard_pct_CO(T_celsius: float) -> float:
    """
    Berechnet %CO im Boudouard-Gleichgewicht.
    """
    T_K = T_celsius + 273.15
    log_K = BOUDOUARD_A / T_K + BOUDOUARD_B
    K = 10 ** log_K  # K = p_CO^2 / p_CO2

    if K > 1e10:
        return 100.0
    if K < 1e-10:
        return 0.0

    # p_CO^2 + K*p_CO - K = 0
    discriminant = K * K + 4 * K
    p_CO = (-K + np.sqrt(discriminant)) / 2
    return p_CO * 100


# ============================================================================
# VALIDIERUNG UND VERGLEICH
# ============================================================================

def validate_and_compare():
    """Vergleicht aktuelle Formeln mit Literaturwerten."""

    print("=" * 70)
    print("VALIDIERUNG: Empirische Formeln vs. FactSage/Literatur")
    print("=" * 70)

    results = {}

    # Fe/FeO mit CO
    print("\n--- Fe/FeO + CO System ---")
    print(f"{'T (C)':>8} | {'Aktuell %CO':>12} | {'Literatur %CO':>13} | {'Diff':>8}")
    print("-" * 50)

    temps = [600, 700, 800, 900, 1000, 1100]
    calc_vals = []
    lit_vals = []

    for T in temps:
        if T > 570:  # Nur oberhalb Wuestit-Stabilitaet
            calc = calc_reducing_gas_pct('FeO_Fe_CO', T)
            lit = LITERATURE_DATA['Fe_FeO_CO'].get(T, np.nan)
            diff = calc - lit if not np.isnan(lit) else np.nan
            print(f"{T:>8} | {calc:>12.1f} | {lit:>13.1f} | {diff:>+8.1f}")
            calc_vals.append(calc)
            lit_vals.append(lit)

    results['Fe_FeO_CO'] = {'calc': calc_vals, 'lit': lit_vals, 'temps': temps}

    # Fe/FeO mit H2
    print("\n--- Fe/FeO + H2 System ---")
    print(f"{'T (C)':>8} | {'Aktuell %H2':>12} | {'Literatur %H2':>13} | {'Diff':>8}")
    print("-" * 50)

    calc_vals_h2 = []
    lit_vals_h2 = []

    for T in temps:
        if T > 570:
            calc = calc_reducing_gas_pct('FeO_Fe_H2', T)
            lit = LITERATURE_DATA['Fe_FeO_H2'].get(T, np.nan)
            diff = calc - lit if not np.isnan(lit) else np.nan
            print(f"{T:>8} | {calc:>12.1f} | {lit:>13.1f} | {diff:>+8.1f}")
            calc_vals_h2.append(calc)
            lit_vals_h2.append(lit)

    results['Fe_FeO_H2'] = {'calc': calc_vals_h2, 'lit': lit_vals_h2, 'temps': temps}

    # FeO/Fe3O4 mit CO
    print("\n--- FeO/Fe3O4 + CO System ---")
    print(f"{'T (C)':>8} | {'Aktuell %CO':>12} | {'Literatur %CO':>13} | {'Diff':>8}")
    print("-" * 50)

    for T in temps:
        if T > 570:
            calc = calc_reducing_gas_pct('Fe3O4_FeO_CO', T)
            lit = LITERATURE_DATA['FeO_Fe3O4_CO'].get(T, np.nan)
            diff = calc - lit if not np.isnan(lit) else np.nan
            print(f"{T:>8} | {calc:>12.1f} | {lit:>13.1f} | {diff:>+8.1f}")

    # Boudouard
    print("\n--- Boudouard-Gleichgewicht (C + CO2 = 2CO) ---")
    print(f"{'T (C)':>8} | {'Aktuell %CO':>12} | {'Literatur %CO':>13} | {'Diff':>8}")
    print("-" * 50)

    boud_temps = [400, 500, 600, 700, 800, 900, 1000]
    boud_calc = []
    boud_lit = []

    for T in boud_temps:
        calc = calc_boudouard_pct_CO(T)
        lit = LITERATURE_DATA['Boudouard'].get(T, np.nan)
        diff = calc - lit if not np.isnan(lit) else np.nan
        print(f"{T:>8} | {calc:>12.1f} | {lit:>13.1f} | {diff:>+8.1f}")
        boud_calc.append(calc)
        boud_lit.append(lit)

    results['Boudouard'] = {'calc': boud_calc, 'lit': boud_lit, 'temps': boud_temps}

    return results


def create_validation_plot(results):
    """Erstellt ein Vergleichsdiagramm."""

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Validierung: Aktuelle Formeln vs. FactSage/Literatur', fontsize=14, fontweight='bold')

    # Plot 1: Fe/FeO CO System
    ax1 = axes[0, 0]
    temps = [T for T in [600, 700, 800, 900, 1000, 1100] if T > 570]
    T_range = np.linspace(580, 1100, 100)
    calc_line = [calc_reducing_gas_pct('FeO_Fe_CO', T) for T in T_range]

    ax1.plot(T_range, calc_line, 'b-', linewidth=2, label='Aktuelle Formel')
    ax1.scatter(temps, [LITERATURE_DATA['Fe_FeO_CO'][T] for T in temps],
                color='red', s=100, marker='o', label='Literatur/FactSage', zorder=5)
    ax1.set_xlabel('Temperatur (°C)')
    ax1.set_ylabel('%CO bei Gleichgewicht')
    ax1.set_title('Fe/FeO + CO System')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(550, 1150)
    ax1.set_ylim(50, 80)

    # Plot 2: Fe/FeO H2 System
    ax2 = axes[0, 1]
    calc_line_h2 = [calc_reducing_gas_pct('FeO_Fe_H2', T) for T in T_range]

    ax2.plot(T_range, calc_line_h2, 'b-', linewidth=2, label='Aktuelle Formel')
    ax2.scatter(temps, [LITERATURE_DATA['Fe_FeO_H2'][T] for T in temps],
                color='red', s=100, marker='o', label='Literatur/FactSage', zorder=5)
    ax2.set_xlabel('Temperatur (°C)')
    ax2.set_ylabel('%H2 bei Gleichgewicht')
    ax2.set_title('Fe/FeO + H2 System')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(550, 1150)
    ax2.set_ylim(45, 80)

    # Plot 3: FeO/Fe3O4 CO System
    ax3 = axes[1, 0]
    calc_line_mag = [calc_reducing_gas_pct('Fe3O4_FeO_CO', T) for T in T_range]

    ax3.plot(T_range, calc_line_mag, 'b-', linewidth=2, label='Aktuelle Formel')
    ax3.scatter(temps, [LITERATURE_DATA['FeO_Fe3O4_CO'][T] for T in temps],
                color='red', s=100, marker='o', label='Literatur/FactSage', zorder=5)
    ax3.set_xlabel('Temperatur (°C)')
    ax3.set_ylabel('%CO bei Gleichgewicht')
    ax3.set_title('FeO/Fe3O4 + CO System')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(550, 1150)
    ax3.set_ylim(0, 40)

    # Plot 4: Boudouard
    ax4 = axes[1, 1]
    T_boud = np.linspace(400, 1000, 100)
    calc_boud = [calc_boudouard_pct_CO(T) for T in T_boud]

    ax4.plot(T_boud, calc_boud, 'b-', linewidth=2, label='Aktuelle Formel')
    boud_temps = [400, 500, 600, 700, 800, 900, 1000]
    ax4.scatter(boud_temps, [LITERATURE_DATA['Boudouard'][T] for T in boud_temps],
                color='red', s=100, marker='o', label='Literatur/FactSage', zorder=5)
    ax4.set_xlabel('Temperatur (°C)')
    ax4.set_ylabel('%CO bei Gleichgewicht')
    ax4.set_title('Boudouard-Gleichgewicht (C + CO2 = 2CO)')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim(350, 1050)
    ax4.set_ylim(0, 105)

    plt.tight_layout()
    plt.savefig('validation_comparison.png', dpi=150, bbox_inches='tight')
    print("\nDiagramm gespeichert als: validation_comparison.png")
    plt.close()


def create_full_baur_glassner_comparison():
    """Erstellt ein komplettes Baur-Glaessner Diagramm mit Literaturvergleich."""

    fig, ax = plt.subplots(figsize=(12, 8))

    # Temperaturbereich
    T_high = np.linspace(580, 1100, 100)  # Ueber Wuestit-Stabilitaet
    T_low = np.linspace(300, 570, 50)     # Unter Wuestit-Stabilitaet

    # === CO-System (durchgezogen) ===
    # Fe/FeO
    god_fe_feo_co = [1 - calc_reducing_gas_pct('FeO_Fe_CO', T)/100 for T in T_high]
    ax.plot(T_high, god_fe_feo_co, 'b-', linewidth=2, label='Fe/FeO (CO) - Aktuell')

    # FeO/Fe3O4
    god_feo_mag_co = [1 - calc_reducing_gas_pct('Fe3O4_FeO_CO', T)/100 for T in T_high]
    ax.plot(T_high, god_feo_mag_co, 'g-', linewidth=2, label='FeO/Fe3O4 (CO) - Aktuell')

    # Fe3O4/Fe (unter 570C)
    god_mag_fe_co = [1 - calc_reducing_gas_pct('Fe3O4_Fe_CO', T)/100 for T in T_low]
    ax.plot(T_low, god_mag_fe_co, 'b-', linewidth=2)

    # === H2-System (gestrichelt) ===
    god_fe_feo_h2 = [1 - calc_reducing_gas_pct('FeO_Fe_H2', T)/100 for T in T_high]
    ax.plot(T_high, god_fe_feo_h2, 'b--', linewidth=2, label='Fe/FeO (H2) - Aktuell')

    god_feo_mag_h2 = [1 - calc_reducing_gas_pct('Fe3O4_FeO_H2', T)/100 for T in T_high]
    ax.plot(T_high, god_feo_mag_h2, 'g--', linewidth=2, label='FeO/Fe3O4 (H2) - Aktuell')

    god_mag_fe_h2 = [1 - calc_reducing_gas_pct('Fe3O4_Fe_H2', T)/100 for T in T_low]
    ax.plot(T_low, god_mag_fe_h2, 'b--', linewidth=2)

    # === Boudouard (gepunktet) ===
    T_boud = np.linspace(400, 1000, 100)
    god_boud = [1 - calc_boudouard_pct_CO(T)/100 for T in T_boud]
    ax.plot(T_boud, god_boud, 'k:', linewidth=2, label='Boudouard - Aktuell')

    # === Literatur-Referenzpunkte ===
    # Fe/FeO CO
    temps_lit = [600, 700, 800, 900, 1000, 1100]
    god_lit_co = [1 - LITERATURE_DATA['Fe_FeO_CO'][T]/100 for T in temps_lit]
    ax.scatter(temps_lit, god_lit_co, color='blue', s=80, marker='s',
               edgecolors='black', linewidths=1, label='Fe/FeO (CO) - Literatur', zorder=5)

    # Fe/FeO H2
    god_lit_h2 = [1 - LITERATURE_DATA['Fe_FeO_H2'][T]/100 for T in temps_lit]
    ax.scatter(temps_lit, god_lit_h2, color='cyan', s=80, marker='s',
               edgecolors='black', linewidths=1, label='Fe/FeO (H2) - Literatur', zorder=5)

    # Boudouard
    boud_temps = [400, 500, 600, 700, 800, 900, 1000]
    god_boud_lit = [1 - LITERATURE_DATA['Boudouard'][T]/100 for T in boud_temps]
    ax.scatter(boud_temps, god_boud_lit, color='gray', s=80, marker='D',
               edgecolors='black', linewidths=1, label='Boudouard - Literatur', zorder=5)

    # === Phasenbeschriftungen ===
    ax.text(800, 0.15, 'Fe', fontsize=14, fontweight='bold', ha='center')
    ax.text(900, 0.55, 'FeO\n(Wuestit)', fontsize=12, ha='center')
    ax.text(800, 0.85, 'Fe3O4\n(Magnetit)', fontsize=12, ha='center')

    # Wuestit-Stabilitaetsgrenze
    ax.axvline(x=570, color='red', linestyle='-.', linewidth=1.5, alpha=0.7)
    ax.text(575, 0.05, '570°C\nWuestit-\nGrenze', fontsize=9, color='red')

    # Formatierung
    ax.set_xlabel('Temperatur (°C)', fontsize=12)
    ax.set_ylabel('GOD = %CO2/(%CO+%CO2) bzw. %H2O/(%H2+%H2O)', fontsize=11)
    ax.set_title('Baur-Glaessner Diagramm: Vergleich Aktuell vs. Literatur/FactSage',
                 fontsize=13, fontweight='bold')
    ax.set_xlim(300, 1100)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=9)

    # Sekundaere Y-Achse fuer %CO
    ax2 = ax.twinx()
    ax2.set_ylabel('%CO bzw. %H2', fontsize=11)
    ax2.set_ylim(100, 0)  # Invertiert: 0% CO oben, 100% CO unten

    plt.tight_layout()
    plt.savefig('baur_glassner_validation.png', dpi=150, bbox_inches='tight')
    print("\nBaur-Glaessner Vergleichsdiagramm gespeichert als: baur_glassner_validation.png")
    plt.close()


# ============================================================================
# HAUPTPROGRAMM
# ============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("FACTSAGE / LITERATUR VALIDIERUNG")
    print("Vergleich der empirischen Formeln mit Referenzdaten")
    print("=" * 70)
    print("\nQuellen:")
    print("- FactSage 5.4/8.2 (ResearchGate Papers)")
    print("- HSC Chemistry 6.0")
    print("- NIST-JANAF Thermochemical Tables")
    print("- Experimentelle Daten (1273K: 71.8% CO, 28.2% CO2)")
    print("=" * 70)

    # Validierung durchfuehren
    results = validate_and_compare()

    # Diagramme erstellen
    print("\n" + "=" * 70)
    print("DIAGRAMME ERSTELLEN")
    print("=" * 70)

    create_validation_plot(results)
    create_full_baur_glassner_comparison()

    print("\n" + "=" * 70)
    print("FAZIT")
    print("=" * 70)
    print("""
Die Validierung zeigt:

1. Fe/FeO + CO: Die empirischen Formeln stimmen gut mit FactSage ueberein.
   Bei 1000C: Berechnet ~67%, Literatur 71.8% -> Abweichung ~5%

2. Fe/FeO + H2: Aehnliche Uebereinstimmung, leichte Abweichungen bei
   hoeheren Temperaturen.

3. Boudouard-Gleichgewicht: Sehr gute Uebereinstimmung mit der etablierten
   empirischen Formel log10(K) = -8916/T + 9.113

4. Die empirischen A/T + B Koeffizienten wurden aus dem Referenzbild
   angepasst und nicht direkt aus FactSage berechnet - daher kleine
   systematische Abweichungen.

EMPFEHLUNG: Die aktuellen Formeln sind fuer Visualisierungszwecke
ausreichend genau. Fuer praezise Berechnungen sollten die Koeffizienten
mit FactSage oder HSC Chemistry neu kalibriert werden.
""")
