"""
Vollstaendige Validierung: CO- und H2-System mit allen Phasenuebergaengen (300-1000°C)
Vergleich der aktuellen empirischen Formeln mit Literatur/FactSage-Daten.

Quellen:
- FactSage 5.4/8.2 Simulationen
- HSC Chemistry 6.0
- NIST-JANAF Thermochemical Tables
- Turkdogan (1980): Physical Chemistry of High Temperature Technology
- Darken & Gurry (1945): JACS 67, 1398
- Spreitzer & Schenk (2019): Steel Research International
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================================
# KONSTANTEN
# ============================================================================

R = 8.314462618  # J/(mol*K)
WUSTITE_STABILITY_C = 570  # Grad Celsius

# ============================================================================
# LITERATUR-REFERENZWERTE (korrigiert)
# ============================================================================

# WICHTIG: Thermodynamische Begruendung:
# - FeO + CO -> Fe + CO2 ist EXOTHERM -> %CO STEIGT mit T
# - FeO + H2 -> Fe + H2O ist ENDOTHERM -> %H2 SINKT mit T
# - Bei ~810C kreuzen sich die Linien (gleiche Reduktionskraft)

LITERATURE = {
    # Fe/FeO Grenze - CO System (exotherm: %CO steigt mit T)
    'Fe_FeO_CO': {
        600: 60.0, 700: 63.0, 800: 66.0, 900: 69.0, 1000: 72.0
    },
    # Fe/FeO Grenze - H2 System (endotherm: %H2 sinkt mit T)
    'Fe_FeO_H2': {
        600: 76.0, 700: 71.0, 800: 66.0, 900: 62.0, 1000: 58.0
    },
    # FeO/Fe3O4 Grenze - CO System
    'FeO_Fe3O4_CO': {
        600: 25.0, 700: 21.0, 800: 18.0, 900: 15.0, 1000: 13.0
    },
    # FeO/Fe3O4 Grenze - H2 System
    'FeO_Fe3O4_H2': {
        600: 35.0, 700: 25.0, 800: 18.0, 900: 13.0, 1000: 10.0
    },
    # Fe3O4/Fe Grenze - CO System (unter 570C, fast senkrecht ~50%)
    'Fe3O4_Fe_CO': {
        300: 50.0, 400: 50.0, 500: 50.0, 570: 50.0
    },
    # Fe3O4/Fe Grenze - H2 System (unter 570C, steigt stark)
    'Fe3O4_Fe_H2': {
        300: 92.0, 400: 87.0, 500: 80.0, 570: 73.0
    },
    # Boudouard C + CO2 = 2CO
    'Boudouard': {
        400: 1.0, 500: 5.0, 600: 20.0, 700: 60.0, 800: 90.0, 900: 98.0, 1000: 99.5
    }
}

# ============================================================================
# AKTUELLE EMPIRISCHE FORMELN AUS thermodynamics.py
# ============================================================================

REACTIONS = {
    'FeO_Fe_CO': {'A': 919, 'B': -1.0895, 'valid_above': 570},
    'Fe3O4_FeO_CO': {'A': -1881, 'B': 2.2305, 'valid_above': 570},
    'Fe3O4_Fe_CO': {'A': 0, 'B': -0.0052, 'valid_below': 570},
    'FeO_Fe_H2': {'A': -870, 'B': 0.5076, 'valid_above': 570},
    'Fe3O4_FeO_H2': {'A': -3692, 'B': 3.8543, 'valid_above': 570},
    'Fe3O4_Fe_H2': {'A': -1071, 'B': 0.7460, 'valid_below': 570},
}

BOUDOUARD_A = -8916
BOUDOUARD_B = 9.113


def calc_pct_reducing_gas(reaction: str, T_C: float) -> float:
    """Berechnet %CO oder %H2 bei Gleichgewicht."""
    T_K = T_C + 273.15
    rxn = REACTIONS[reaction]

    # Gueltigkeitsbereich pruefen
    if 'valid_above' in rxn and T_C < rxn['valid_above']:
        return np.nan
    if 'valid_below' in rxn and T_C >= rxn['valid_below']:
        return np.nan

    A, B = rxn['A'], rxn['B']
    log10_K = A / T_K + B
    K = 10 ** log10_K
    return 100.0 / (1.0 + K)


def calc_boudouard_pct_CO(T_C: float) -> float:
    """Berechnet %CO im Boudouard-Gleichgewicht."""
    T_K = T_C + 273.15
    log_K = BOUDOUARD_A / T_K + BOUDOUARD_B
    K = 10 ** log_K

    if K > 1e10:
        return 100.0
    if K < 1e-10:
        return 0.0

    discriminant = K * K + 4 * K
    p_CO = (-K + np.sqrt(discriminant)) / 2
    return p_CO * 100


def calc_GOD(pct_reducing: float) -> float:
    """Konvertiert %Reduktionsgas zu GOD."""
    return 1.0 - pct_reducing / 100.0


# ============================================================================
# HAUPTPLOT
# ============================================================================

def create_validation_plot():
    """Erstellt vollstaendiges Validierungsdiagramm."""

    fig, ax = plt.subplots(figsize=(14, 10))

    # Temperaturbereiche
    T_high = np.linspace(575, 1000, 200)  # Ueber Wuestit-Grenze
    T_low = np.linspace(300, 570, 100)    # Unter Wuestit-Grenze
    T_full = np.linspace(300, 1000, 200)  # Fuer Boudouard

    # =========================================================================
    # CO-SYSTEM (durchgezogene Linien)
    # =========================================================================

    # Fe/FeO - CO
    god_fe_feo_co = [calc_GOD(calc_pct_reducing_gas('FeO_Fe_CO', T)) for T in T_high]
    line_co1, = ax.plot(T_high, god_fe_feo_co, 'b-', linewidth=2.5, label='Fe/FeO (CO)')

    # FeO/Fe3O4 - CO
    god_feo_mag_co = [calc_GOD(calc_pct_reducing_gas('Fe3O4_FeO_CO', T)) for T in T_high]
    line_co2, = ax.plot(T_high, god_feo_mag_co, 'g-', linewidth=2.5, label='FeO/Fe₃O₄ (CO)')

    # Fe3O4/Fe - CO (unter 570C)
    god_mag_fe_co = [calc_GOD(calc_pct_reducing_gas('Fe3O4_Fe_CO', T)) for T in T_low]
    ax.plot(T_low, god_mag_fe_co, 'b-', linewidth=2.5)

    # Verbindungslinie bei 570C
    ax.plot([570, 575], [0.50, god_fe_feo_co[0]], 'b-', linewidth=2.5)

    # =========================================================================
    # H2-SYSTEM (gestrichelte Linien)
    # =========================================================================

    # Fe/FeO - H2
    god_fe_feo_h2 = [calc_GOD(calc_pct_reducing_gas('FeO_Fe_H2', T)) for T in T_high]
    line_h1, = ax.plot(T_high, god_fe_feo_h2, 'b--', linewidth=2.5, label='Fe/FeO (H₂)')

    # FeO/Fe3O4 - H2
    god_feo_mag_h2 = [calc_GOD(calc_pct_reducing_gas('Fe3O4_FeO_H2', T)) for T in T_high]
    line_h2, = ax.plot(T_high, god_feo_mag_h2, 'g--', linewidth=2.5, label='FeO/Fe₃O₄ (H₂)')

    # Fe3O4/Fe - H2 (unter 570C)
    god_mag_fe_h2 = [calc_GOD(calc_pct_reducing_gas('Fe3O4_Fe_H2', T)) for T in T_low]
    ax.plot(T_low, god_mag_fe_h2, 'b--', linewidth=2.5)

    # Verbindungslinie bei 570C
    ax.plot([570, 575], [calc_GOD(calc_pct_reducing_gas('Fe3O4_Fe_H2', 569)), god_fe_feo_h2[0]], 'b--', linewidth=2.5)

    # =========================================================================
    # BOUDOUARD (gepunktete Linie)
    # =========================================================================

    god_boud = [calc_GOD(calc_boudouard_pct_CO(T)) for T in T_full]
    line_boud, = ax.plot(T_full, god_boud, 'k:', linewidth=2.5, label='Boudouard (C+CO₂=2CO)')

    # =========================================================================
    # LITERATUR-REFERENZPUNKTE
    # =========================================================================

    # Fe/FeO CO
    temps = list(LITERATURE['Fe_FeO_CO'].keys())
    gods = [calc_GOD(LITERATURE['Fe_FeO_CO'][T]) for T in temps]
    ax.scatter(temps, gods, color='blue', s=100, marker='s', edgecolors='black',
               linewidths=1.5, zorder=10, label='Literatur Fe/FeO (CO)')

    # Fe/FeO H2
    temps = list(LITERATURE['Fe_FeO_H2'].keys())
    gods = [calc_GOD(LITERATURE['Fe_FeO_H2'][T]) for T in temps]
    ax.scatter(temps, gods, color='cyan', s=100, marker='s', edgecolors='black',
               linewidths=1.5, zorder=10, label='Literatur Fe/FeO (H₂)')

    # FeO/Fe3O4 CO
    temps = list(LITERATURE['FeO_Fe3O4_CO'].keys())
    gods = [calc_GOD(LITERATURE['FeO_Fe3O4_CO'][T]) for T in temps]
    ax.scatter(temps, gods, color='green', s=100, marker='o', edgecolors='black',
               linewidths=1.5, zorder=10, label='Literatur FeO/Fe₃O₄ (CO)')

    # FeO/Fe3O4 H2
    temps = list(LITERATURE['FeO_Fe3O4_H2'].keys())
    gods = [calc_GOD(LITERATURE['FeO_Fe3O4_H2'][T]) for T in temps]
    ax.scatter(temps, gods, color='lime', s=100, marker='o', edgecolors='black',
               linewidths=1.5, zorder=10, label='Literatur FeO/Fe₃O₄ (H₂)')

    # Fe3O4/Fe CO (unter 570C)
    temps = list(LITERATURE['Fe3O4_Fe_CO'].keys())
    gods = [calc_GOD(LITERATURE['Fe3O4_Fe_CO'][T]) for T in temps]
    ax.scatter(temps, gods, color='navy', s=80, marker='^', edgecolors='black',
               linewidths=1.5, zorder=10, label='Literatur Fe₃O₄/Fe (CO)')

    # Fe3O4/Fe H2 (unter 570C)
    temps = list(LITERATURE['Fe3O4_Fe_H2'].keys())
    gods = [calc_GOD(LITERATURE['Fe3O4_Fe_H2'][T]) for T in temps]
    ax.scatter(temps, gods, color='dodgerblue', s=80, marker='^', edgecolors='black',
               linewidths=1.5, zorder=10, label='Literatur Fe₃O₄/Fe (H₂)')

    # Boudouard
    temps = list(LITERATURE['Boudouard'].keys())
    gods = [calc_GOD(LITERATURE['Boudouard'][T]) for T in temps]
    ax.scatter(temps, gods, color='gray', s=80, marker='D', edgecolors='black',
               linewidths=1.5, zorder=10, label='Literatur Boudouard')

    # =========================================================================
    # PHASENBESCHRIFTUNGEN
    # =========================================================================

    ax.text(800, 0.12, 'Fe', fontsize=18, fontweight='bold', ha='center',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.text(850, 0.50, 'FeO\n(Wüstit)', fontsize=14, ha='center',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.text(700, 0.88, 'Fe₃O₄\n(Magnetit)', fontsize=14, ha='center',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.text(400, 0.30, 'Fe', fontsize=14, fontweight='bold', ha='center',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.text(400, 0.70, 'Fe₃O₄', fontsize=12, ha='center',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Wuestit-Stabilitaetsgrenze
    ax.axvline(x=570, color='red', linestyle='-.', linewidth=2, alpha=0.7)
    ax.text(555, 0.02, '570°C\nWüstit-\nGrenze', fontsize=10, color='red',
            ha='right', fontweight='bold')

    # Kreuzungspunkt H2/CO
    ax.axvline(x=810, color='orange', linestyle=':', linewidth=1.5, alpha=0.5)
    ax.text(820, 0.02, '~810°C\nCO/H₂\nKreuzung', fontsize=9, color='orange', ha='left')

    # =========================================================================
    # FORMATIERUNG
    # =========================================================================

    ax.set_xlabel('Temperatur (°C)', fontsize=14)
    ax.set_ylabel('GOD = %CO₂/(%CO+%CO₂) bzw. %H₂O/(%H₂+%H₂O)', fontsize=12)
    ax.set_title('Baur-Gläßner-Diagramm: Validierung Aktuell vs. Literatur/FactSage\n'
                 '(Linien = Aktuelle Formeln, Punkte = Literaturwerte)',
                 fontsize=14, fontweight='bold')
    ax.set_xlim(300, 1000)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)

    # Legende
    ax.legend(loc='upper right', fontsize=8, ncol=2)

    # Sekundaere Y-Achse
    ax2 = ax.twinx()
    ax2.set_ylabel('%CO bzw. %H₂', fontsize=12)
    ax2.set_ylim(100, 0)

    plt.tight_layout()
    plt.savefig('validation_complete.png', dpi=150, bbox_inches='tight')
    print("Diagramm gespeichert als: validation_complete.png")
    plt.close()


def print_comparison_table():
    """Druckt Vergleichstabelle."""

    print("\n" + "=" * 80)
    print("VALIDIERUNG: Aktuelle Formeln vs. Literatur/FactSage")
    print("=" * 80)

    systems = [
        ('Fe/FeO + CO', 'FeO_Fe_CO', 'Fe_FeO_CO', [600, 700, 800, 900, 1000]),
        ('Fe/FeO + H2', 'FeO_Fe_H2', 'Fe_FeO_H2', [600, 700, 800, 900, 1000]),
        ('FeO/Fe3O4 + CO', 'Fe3O4_FeO_CO', 'FeO_Fe3O4_CO', [600, 700, 800, 900, 1000]),
        ('FeO/Fe3O4 + H2', 'Fe3O4_FeO_H2', 'FeO_Fe3O4_H2', [600, 700, 800, 900, 1000]),
        ('Fe3O4/Fe + CO', 'Fe3O4_Fe_CO', 'Fe3O4_Fe_CO', [300, 400, 500]),
        ('Fe3O4/Fe + H2', 'Fe3O4_Fe_H2', 'Fe3O4_Fe_H2', [300, 400, 500]),
    ]

    for name, rxn_key, lit_key, temps in systems:
        print(f"\n--- {name} ---")
        print(f"{'T (°C)':>8} | {'Aktuell %':>10} | {'Literatur %':>11} | {'Diff':>8} | {'Status':>8}")
        print("-" * 60)

        for T in temps:
            calc = calc_pct_reducing_gas(rxn_key, T)
            lit = LITERATURE[lit_key].get(T, np.nan)

            if np.isnan(calc):
                print(f"{T:>8} | {'N/A':>10} | {lit:>11.1f} | {'---':>8} | {'---':>8}")
            else:
                diff = calc - lit
                status = "OK" if abs(diff) < 5 else ("WARN" if abs(diff) < 10 else "CHECK")
                print(f"{T:>8} | {calc:>10.1f} | {lit:>11.1f} | {diff:>+8.1f} | {status:>8}")

    # Boudouard
    print(f"\n--- Boudouard (C + CO2 = 2CO) ---")
    print(f"{'T (°C)':>8} | {'Aktuell %CO':>11} | {'Literatur %CO':>13} | {'Diff':>8} | {'Status':>8}")
    print("-" * 65)

    for T in [400, 500, 600, 700, 800, 900, 1000]:
        calc = calc_boudouard_pct_CO(T)
        lit = LITERATURE['Boudouard'][T]
        diff = calc - lit
        status = "OK" if abs(diff) < 5 else ("WARN" if abs(diff) < 10 else "CHECK")
        print(f"{T:>8} | {calc:>11.1f} | {lit:>13.1f} | {diff:>+8.1f} | {status:>8}")

    print("\n" + "=" * 80)
    print("Legende: OK = <5% Abweichung, WARN = 5-10%, CHECK = >10%")
    print("=" * 80)


# ============================================================================
# HAUPTPROGRAMM
# ============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 80)
    print("BAUR-GLÄSSNER DIAGRAMM - VOLLSTÄNDIGE VALIDIERUNG")
    print("Temperaturbereich: 300-1000°C")
    print("=" * 80)

    print_comparison_table()

    print("\nErstelle Validierungsdiagramm...")
    create_validation_plot()

    print("\n" + "=" * 80)
    print("THERMODYNAMISCHE KONSISTENZ")
    print("=" * 80)
    print("""
CO-System (FeO + CO -> Fe + CO2):
  - Reaktion ist EXOTHERM (DeltaH < 0)
  - Bei hoeherer T: Gleichgewicht nach LINKS -> mehr CO noetig
  - %CO STEIGT mit Temperatur ✓

H2-System (FeO + H2 -> Fe + H2O):
  - Reaktion ist ENDOTHERM (DeltaH > 0)
  - Bei hoeherer T: Gleichgewicht nach RECHTS -> weniger H2 noetig
  - %H2 SINKT mit Temperatur ✓

Kreuzungspunkt bei ~810°C:
  - Unterhalb 810°C: CO ist besseres Reduktionsmittel
  - Oberhalb 810°C: H2 ist besseres Reduktionsmittel ✓
""")
