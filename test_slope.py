"""Test Fe2O3 vs Fe3O4 slope comparison"""
from mass_balance_h2 import *

# Test mit Fe2O3 (97%)
ore1 = OreCompositionH2(Fe2O3=97.0, Fe3O4=0.0, SiO2=2.0, CaO=1.0)
dri = DRIComposition(Fe_met=92.0, FeO=3.0, C=0.0, Gangue=5.0)
gas = ReductionGasComposition(H2=95.0, H2O=5.0)

ore_amount1 = calculate_ore_amount_from_fe_h2(ore1, dri)
r1 = calculate_h2_mass_balance(ore1, ore_amount1, gas, 800.0, dri)

print("=== Fe2O3 (Hematite) ===")
print(f"  O/Fe Erz:    {r1.o_fe_ore:.4f}")
print(f"  H2 consumed: {r1.h2_amount_nm3:.1f} Nm3/t")
print(f"  Slope k:     {r1.slope_k:.4f}")
print(f"  Intercept d: {r1.intercept_d:.4f}")
print(f"  GOD out:     {r1.god_output*100:.1f}%")

# Test mit Fe3O4 (97%)
ore2 = OreCompositionH2(Fe2O3=0.0, Fe3O4=97.0, SiO2=2.0, CaO=1.0)
ore_amount2 = calculate_ore_amount_from_fe_h2(ore2, dri)
r2 = calculate_h2_mass_balance(ore2, ore_amount2, gas, 800.0, dri)

print()
print("=== Fe3O4 (Magnetite) ===")
print(f"  O/Fe Erz:    {r2.o_fe_ore:.4f}")
print(f"  H2 consumed: {r2.h2_amount_nm3:.1f} Nm3/t")
print(f"  Slope k:     {r2.slope_k:.4f}")
print(f"  Intercept d: {r2.intercept_d:.4f}")
print(f"  GOD out:     {r2.god_output*100:.1f}%")

print()
print("=== Differenz ===")
print(f"  Δ Slope:     {r2.slope_k - r1.slope_k:.4f}")
print(f"  Δ H2:        {r2.h2_amount_nm3 - r1.h2_amount_nm3:.1f} Nm3/t")
