"""
Test-Script fuer die thermodynamischen Berechnungen.
Vergleicht berechnete Werte mit Referenzwerten.
"""

from thermodynamics import BaurGlassnerThermo
import numpy as np

print("=" * 60)
print("TEST: Baur-Glaessner Thermodynamik")
print("=" * 60)

thermo = BaurGlassnerThermo("NBS_Tables_Library.xlsx")

# Debug-Info ausgeben
thermo.print_debug_info()

print("\n" + "=" * 60)
print("VERGLEICH MIT REFERENZWERTEN")
print("=" * 60)

print("\n=== Erwartete Werte (Referenz) ===")
print("Fe/FeO (CO) bei 800C:  GOD ~ 0.28-0.35")
print("Fe/FeO (CO) bei 1000C: GOD ~ 0.35-0.40")
print("Fe/FeO (H2) bei 800C:  GOD ~ 0.33")
print("Fe3O4/Fe (CO) bei 400C: GOD ~ 0.25")
print("Boudouard bei 700C:    GOD ~ 0.35-0.50")
print("Boudouard bei 600C:    GOD ~ 0.10-0.20")

print("\n=== Berechnete Werte ===")

# Fe/FeO Tests
for T_C in [600, 700, 800, 900, 1000]:
    T_K = T_C + 273.15
    god_co = thermo.GOD('FeO_Fe_CO', T_K)
    god_h2 = thermo.GOD('FeO_Fe_H2', T_K)

    co_str = f"{god_co:.3f}" if not np.isnan(god_co) else "N/A"
    h2_str = f"{god_h2:.3f}" if not np.isnan(god_h2) else "N/A"

    print(f"Fe/FeO bei {T_C}C: CO={co_str}, H2={h2_str}")

print()

# Fe3O4/Fe Tests (unter 570C)
for T_C in [400, 500]:
    T_K = T_C + 273.15
    god_co = thermo.GOD('Fe3O4_Fe_CO', T_K)
    god_h2 = thermo.GOD('Fe3O4_Fe_H2', T_K)

    co_str = f"{god_co:.3f}" if not np.isnan(god_co) else "N/A"
    h2_str = f"{god_h2:.3f}" if not np.isnan(god_h2) else "N/A"

    print(f"Fe3O4/Fe bei {T_C}C: CO={co_str}, H2={h2_str}")

print()

# Boudouard Tests
print("Boudouard-Gleichgewicht:")
for T_C in [500, 600, 700, 800, 900, 1000]:
    T_K = T_C + 273.15
    god = thermo.boudouard_GOD(T_K)
    print(f"  {T_C}C: GOD = {god:.3f}")

print("\n" + "=" * 60)
print("VALIDIERUNG")
print("=" * 60)

# Check: Bei 570C sollten sich Fe/FeO und Fe3O4/Fe treffen
T_570 = 843  # K
god_fe_feo_co = thermo.GOD('FeO_Fe_CO', T_570 + 1)
god_fe3o4_fe_co = thermo.GOD('Fe3O4_Fe_CO', T_570 - 1)

print(f"\n570C Tripelpunkt (CO-System):")
print(f"  Fe/FeO bei 571C:    GOD = {god_fe_feo_co:.3f}")
print(f"  Fe3O4/Fe bei 569C:  GOD = {god_fe3o4_fe_co:.3f}")
print(f"  Differenz: {abs(god_fe_feo_co - god_fe3o4_fe_co):.4f}")

# Check: H2 vs CO Kreuzungspunkt
print("\nKreuzung H2/CO (Fe/FeO-Linie):")
for T_C in [700, 750, 800, 850, 900]:
    T_K = T_C + 273.15
    god_co = thermo.GOD('FeO_Fe_CO', T_K)
    god_h2 = thermo.GOD('FeO_Fe_H2', T_K)
    diff = god_co - god_h2
    print(f"  {T_C}C: CO={god_co:.3f}, H2={god_h2:.3f}, Diff={diff:+.3f}")
