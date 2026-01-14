"""
Massenbilanz-Berechnungen fuer das H2-Rist-Diagramm.

Berechnet die Betriebslinie einer H2-Direktreduktionsanlage basierend auf:
- Erzzusammensetzung (Fe2O3, SiO2, CaO)
- Reduktionsgas-Zusammensetzung (H2, H2O)
- Produktzusammensetzung (DRI/HBI: Fe, ggf. Gangue)

Die Betriebslinie im H2-Rist-Diagramm:
    O/Fe = k * (O/H2) + d

wobei:
    k = H2/Fe (Steigung, mol H2 pro mol Fe)
    d = 0 (kein externer Sauerstoff bei reiner H2-Reduktion)
"""

from dataclasses import dataclass
from typing import Tuple, Optional


# Molare Massen [kg/kmol]
MOLAR_MASS = {
    'Fe': 56.0,
    'Fe2O3': 160.0,
    'FeO': 72.0,
    'Fe3O4': 232.0,
    'Si': 28.0,
    'SiO2': 60.0,
    'Ca': 40.0,
    'CaO': 56.0,
    'O': 16.0,
    'H2': 2.0,
    'H2O': 18.0,
}

# Molares Volumen ideales Gas bei Normalbedingungen [Nm³/kmol]
MOLAR_VOLUME_NM3 = 22.414


@dataclass
class OreCompositionH2:
    """Erzzusammensetzung in Massen-% fuer H2-Reduktion.
    
    Kann sowohl Fe2O3 (Haematit) als auch Fe3O4 (Magnetit) enthalten.
    """
    Fe2O3: float = 94.4  # % Haematit
    Fe3O4: float = 0.0   # % Magnetit
    SiO2: float = 2.5    # %
    CaO: float = 3.1     # %

    def validate(self) -> bool:
        """Prueft ob die Summe ~100% ergibt."""
        total = self.Fe2O3 + self.Fe3O4 + self.SiO2 + self.CaO
        return 99.0 <= total <= 101.0


@dataclass
class ReductionGasComposition:
    """Reduktionsgas-Zusammensetzung in Vol-%.

    Typische Werte fuer H2-Direktreduktion:
    - Frischgas: 95-100% H2, 0-5% H2O
    - Verbrauchtes Gas: 50-70% H2, 30-50% H2O
    """
    H2: float = 95.0     # %
    H2O: float = 5.0     # %

    def validate(self) -> bool:
        """Prueft ob die Summe ~100% ergibt."""
        total = self.H2 + self.H2O
        return 99.0 <= total <= 101.0

    @property
    def god(self) -> float:
        """Gas Oxidation Degree = H2O / (H2 + H2O)."""
        total = self.H2 + self.H2O
        if total <= 0:
            return 0.0
        return self.H2O / total


@dataclass
class DRIComposition:
    """DRI/HBI Produkt-Zusammensetzung in Massen-%.

    Typische Werte fuer DRI:
    - Metallisierung 92-95%
    - Kohlenstoffgehalt 0-4% (je nach Prozess)
    """
    Fe_met: float = 92.0   # % metallisches Eisen
    FeO: float = 3.0       # % Rest-FeO (nicht reduziert)
    C: float = 0.0         # % Kohlenstoff (bei reiner H2-Reduktion = 0)
    Gangue: float = 5.0    # % Gangart (SiO2, CaO, etc.)

    def validate(self) -> bool:
        """Prueft ob die Summe ~100% ergibt."""
        total = self.Fe_met + self.FeO + self.C + self.Gangue
        return 99.0 <= total <= 101.0

    @property
    def metallization(self) -> float:
        """Metallisierungsgrad in %."""
        # Fe_total = Fe_met + Fe in FeO
        fe_in_feo = self.FeO * MOLAR_MASS['Fe'] / MOLAR_MASS['FeO']
        fe_total = self.Fe_met + fe_in_feo
        if fe_total <= 0:
            return 0.0
        return (self.Fe_met / fe_total) * 100.0


@dataclass
class H2MassBalanceResults:
    """Ergebnisse der H2-Massenbilanz-Berechnung."""
    # Stoffmengen [kmol/tDRI]
    n_Fe: float          # Eisen im DRI (metallisch + in FeO)
    n_H2_consumed: float # Verbrauchter Wasserstoff
    n_O_removed: float   # Entfernter Sauerstoff
    n_H2O_produced: float  # Produziertes Wasser
    
    # Gas Ein-/Ausgang [kmol/tDRI]
    n_gas_input: float   # Gesamt-Gas Eingang
    n_H2_input: float    # H2 im Eingangsgas
    n_H2O_input: float   # H2O im Eingangsgas
    n_H2_output: float   # H2 im Ausgangsgas
    n_H2O_output: float  # H2O im Ausgangsgas

    # Betriebslinien-Parameter
    slope_k: float       # Steigung k = H2/Fe
    intercept_d: float   # Ordinatenabschnitt d (normalerweise 0)

    # Ein- und Ausgangs-Gaszusammensetzung
    god_input: float     # GOD am Eingang (frisches Gas)
    god_output: float    # GOD am Ausgang (verbrauchtes Gas)
    o_h2_input: float    # O/H2 am Eingang
    o_h2_output: float   # O/H2 am Ausgang
    
    # Gasnutzung
    gas_utilization: float  # Gasnutzungsgrad [%]

    # Reduktionsanteile
    indirect_reduction_pct: float  # Anteil indirekte Reduktion [%]
    direct_reduction_pct: float    # Anteil direkte Reduktion [%] (bei H2: normalerweise 100% indirekt)

    # Massen
    ore_amount: float    # Erzmenge [kg/tDRI]
    dri_amount: float    # DRI-Menge [kg] (Bezugsbasis)
    h2_amount_nm3: float # H2-Verbrauch [Nm³/tDRI]
    gas_input_nm3: float # Gesamt-Gas Eingang [Nm³/tDRI]
    o_fe_ore: float      # O/Fe des Erzes (abhängig von Fe2O3/Fe3O4 Mischung)


def calculate_ore_amount_from_fe_h2(
    ore: OreCompositionH2,
    dri: DRIComposition,
    dri_amount: float = 1000.0,
) -> float:
    """
    Berechnet die benoetigte Erzmenge basierend auf Fe-Bilanz.

    Args:
        ore: Erzzusammensetzung
        dri: DRI-Zusammensetzung
        dri_amount: DRI-Menge [kg/tDRI]

    Returns:
        Erzmenge [kg/tDRI]
    """
    # Fe im DRI [kg]
    m_Fe_met = dri_amount * dri.Fe_met / 100.0
    m_FeO = dri_amount * dri.FeO / 100.0
    m_Fe_in_FeO = m_FeO * MOLAR_MASS['Fe'] / MOLAR_MASS['FeO']
    m_Fe_total = m_Fe_met + m_Fe_in_FeO

    # Fe-Gehalt im Erz [kg Fe / kg Erz]
    # Fe2O3 -> 2 Fe, Fe-Anteil = 2*56 / 160 = 0.7
    # Fe3O4 -> 3 Fe, Fe-Anteil = 3*56 / 232 = 0.724
    fe_in_fe2o3 = 2 * MOLAR_MASS['Fe'] / MOLAR_MASS['Fe2O3']  # = 0.7
    fe_in_fe3o4 = 3 * MOLAR_MASS['Fe'] / 232.0  # M_Fe3O4 = 232
    fe_content_ore = (ore.Fe2O3 / 100.0 * fe_in_fe2o3 + 
                      ore.Fe3O4 / 100.0 * fe_in_fe3o4)

    # Erzmenge = Fe im DRI / Fe-Gehalt im Erz
    if fe_content_ore > 0:
        ore_amount = m_Fe_total / fe_content_ore
    else:
        ore_amount = 0.0

    return ore_amount


def calculate_h2_mass_balance(
    ore: OreCompositionH2,
    ore_amount: float,
    input_gas: ReductionGasComposition,
    gas_amount_nm3: float,
    dri: DRIComposition,
    dri_amount: float = 1000.0,
) -> H2MassBalanceResults:
    """
    Berechnet die Massenbilanz und Betriebslinien-Parameter fuer H2-Reduktion.

    Args:
        ore: Erzzusammensetzung [Massen-%]
        ore_amount: Erzmenge [kg/tDRI]
        input_gas: Eingangs-Gaszusammensetzung [Vol-%]
        gas_amount_nm3: Gesamte Reduktionsgasmenge [Nm³/tDRI]
        dri: DRI-Zusammensetzung [Massen-%]
        dri_amount: DRI-Menge [kg/tDRI] (Bezugsbasis, default 1000)

    Returns:
        H2MassBalanceResults mit allen berechneten Werten
    """
    # === 1. Stoffmengen im DRI [kmol/tDRI] ===
    m_Fe_met = dri_amount * dri.Fe_met / 100.0  # kg Fe metallisch
    m_FeO = dri_amount * dri.FeO / 100.0        # kg FeO

    n_Fe_met = m_Fe_met / MOLAR_MASS['Fe']      # kmol Fe metallisch
    n_FeO = m_FeO / MOLAR_MASS['FeO']           # kmol FeO
    n_Fe_in_FeO = n_FeO                         # kmol Fe in FeO (1:1)

    n_Fe_total = n_Fe_met + n_Fe_in_FeO         # kmol Fe gesamt

    # Sauerstoff im DRI (aus Rest-FeO)
    n_O_in_dri = n_FeO  # kmol O (1 O pro FeO)

    # === 2. Stoffmengen im Erz [kmol/tDRI] ===
    m_Fe2O3_ore = ore_amount * ore.Fe2O3 / 100.0  # kg Fe2O3
    m_Fe3O4_ore = ore_amount * ore.Fe3O4 / 100.0  # kg Fe3O4

    n_Fe2O3_ore = m_Fe2O3_ore / MOLAR_MASS['Fe2O3']  # kmol Fe2O3
    n_Fe3O4_ore = m_Fe3O4_ore / 232.0  # kmol Fe3O4 (M = 232)

    # Sauerstoff im Erz (Fe2O3: 3 O, Fe3O4: 4 O)
    n_O_ore = 3 * n_Fe2O3_ore + 4 * n_Fe3O4_ore  # kmol O
    
    # Fe im Erz (Fe2O3: 2 Fe, Fe3O4: 3 Fe)
    n_Fe_ore = 2 * n_Fe2O3_ore + 3 * n_Fe3O4_ore  # kmol Fe
    
    # O/Fe des Erzes berechnen (dynamisch basierend auf Erzzusammensetzung)
    # Fe2O3: O/Fe = 3/2 = 1.5, Fe3O4: O/Fe = 4/3 = 1.333
    o_fe_ore = n_O_ore / n_Fe_ore if n_Fe_ore > 0 else 1.5

    # === 3. Sauerstoff-Bilanz ===
    # Entfernter Sauerstoff = O im Erz - O im DRI
    n_O_removed = n_O_ore - n_O_in_dri

    # === 4. Wasserstoff-Bilanz ===
    # Jedes entfernte O wird zu H2O: H2 + 1/2 O2 -> H2O
    # Also: n_H2_consumed = n_O_removed (1:1 Verhaeltnis)
    n_H2_consumed = n_O_removed
    n_H2O_produced = n_O_removed

    # === 5. Gaszusammensetzung mit tatsächlicher Gasmenge ===
    god_input = input_gas.god
    o_h2_input = 1.0 + god_input
    
    # Gas-Stoffmengen aus eingegebener Gasmenge [kmol/tDRI]
    n_gas_input = gas_amount_nm3 / MOLAR_VOLUME_NM3
    n_H2_input = n_gas_input * input_gas.H2 / 100.0
    n_H2O_input = n_gas_input * input_gas.H2O / 100.0
    
    # Ausgangsgas: H2 wird verbraucht, H2O wird produziert
    n_H2_output = n_H2_input - n_H2_consumed
    n_H2O_output = n_H2O_input + n_H2O_produced
    
    # Sicherheitsprüfung: Genug H2 vorhanden?
    if n_H2_output < 0:
        # Mehr H2 verbraucht als eingegeben - nicht genug Gas!
        n_H2_output = 0.0
    
    # Ausgangs-GOD aus tatsächlicher Gaszusammensetzung
    n_gas_output = n_H2_output + n_H2O_output
    god_output = n_H2O_output / n_gas_output if n_gas_output > 0 else 1.0
    o_h2_output = 1.0 + god_output
    
    # Gasnutzungsgrad = verbrauchtes H2 / eingesetztes H2
    gas_utilization = (n_H2_consumed / n_H2_input * 100.0) if n_H2_input > 0 else 0.0

    # === 6. Betriebslinien-Parameter ===
    # Im H2-Rist-Diagramm verbindet die Betriebslinie:
    #   - Punkt 1: (O/H2_input, O/Fe_product) = (o_h2_input, ~0)
    #   - Punkt 2: (O/H2_output, O/Fe_ore) = (o_h2_output, o_fe_ore)
    #
    # Die Steigung ist: slope = Δ(O/Fe) / Δ(O/H2)
    # Das Rist-Diagramm zeigt die Massenbilanz des Gegenstrom-Reaktors:
    # - Frisches Gas (niedrig GOD) trifft auf reduziertes Produkt (O/Fe ≈ 0)
    # - Verbrauchtes Gas (hoch GOD) verlässt den Reaktor mit dem Erz (O/Fe = o_fe_ore)
    
    # O/Fe im Produkt (aus Rest-FeO)
    o_fe_product = n_O_in_dri / n_Fe_total if n_Fe_total > 0 else 0
    # o_fe_ore wurde oben aus Erzzusammensetzung berechnet
    
    # Betriebslinie: O/Fe = k * (O/H2) + d
    # Steigung k = Δ(O/Fe) / Δ(O/H2)
    delta_o_fe = o_fe_ore - o_fe_product
    delta_o_h2 = o_h2_output - o_h2_input
    
    if delta_o_h2 > 0.001:  # Vermeidung Division durch 0
        slope_k = delta_o_fe / delta_o_h2
    else:
        slope_k = o_fe_ore  # Fallback (dynamisch: 1.5 für Fe2O3, 1.333 für Fe3O4)
    
    # Intercept d: Die Linie geht durch (o_h2_input, o_fe_product)
    # o_fe_product = k * o_h2_input + d  =>  d = o_fe_product - k * o_h2_input
    intercept_d = o_fe_product - slope_k * o_h2_input

    # === 7. Reduktionsanteile ===
    # Bei reiner H2-Reduktion ist alles indirekte Reduktion (Gas)
    # Es gibt keine direkte Reduktion wie beim Hochofen (Koks)
    indirect_reduction_pct = 100.0
    direct_reduction_pct = 0.0

    # Falls es Rest-FeO gibt, ist die Reduktion nicht 100%
    metallization = dri.metallization
    actual_reduction_pct = (o_fe_ore - o_fe_product) / o_fe_ore * 100.0 if o_fe_ore > 0 else 100.0

    # === 8. H2-Verbrauch in Nm³ ===
    h2_amount_nm3 = n_H2_consumed * MOLAR_VOLUME_NM3  # Nm³/tDRI

    return H2MassBalanceResults(
        n_Fe=n_Fe_total,
        n_H2_consumed=n_H2_consumed,
        n_O_removed=n_O_removed,
        n_H2O_produced=n_H2O_produced,
        n_gas_input=n_gas_input,
        n_H2_input=n_H2_input,
        n_H2O_input=n_H2O_input,
        n_H2_output=n_H2_output,
        n_H2O_output=n_H2O_output,
        slope_k=slope_k,
        intercept_d=intercept_d,
        god_input=god_input,
        god_output=god_output,
        o_h2_input=o_h2_input,
        o_h2_output=o_h2_output,
        gas_utilization=gas_utilization,
        indirect_reduction_pct=indirect_reduction_pct,
        direct_reduction_pct=direct_reduction_pct,
        ore_amount=ore_amount,
        dri_amount=dri_amount,
        h2_amount_nm3=h2_amount_nm3,
        gas_input_nm3=gas_amount_nm3,
        o_fe_ore=o_fe_ore,
    )


def calculate_h2_operating_line_points(
    results: H2MassBalanceResults,
    o_fe_start: float = 0.0,
    o_fe_end: Optional[float] = None,
) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """
    Berechnet die Endpunkte der Betriebslinie.

    Die Betriebslinie: O/Fe = k * (O/H2) + d
    Umgestellt: O/H2 = (O/Fe - d) / k

    Args:
        results: Ergebnisse der Massenbilanz
        o_fe_start: Start O/Fe-Wert (default: 0 = reines Fe)
        o_fe_end: End O/Fe-Wert (default: None = aus Erzzusammensetzung)

    Returns:
        ((o_h2_start, o_fe_start), (o_h2_end, o_fe_end))
    """
    k = results.slope_k
    d = results.intercept_d
    
    # Wenn o_fe_end nicht angegeben, aus Erz-O/Fe nehmen
    if o_fe_end is None:
        o_fe_end = results.o_fe_ore

    if k <= 0:
        return ((1.0, o_fe_start), (2.0, o_fe_end))

    # O/H2 = (O/Fe - d) / k
    o_h2_start = (o_fe_start - d) / k
    o_h2_end = (o_fe_end - d) / k

    return ((o_h2_start, o_fe_start), (o_h2_end, o_fe_end))


if __name__ == "__main__":
    print("=== H2 Massenbilanz Test ===\n")

    # Standardwerte
    ore = OreCompositionH2(Fe2O3=94.4, SiO2=2.5, CaO=3.1)
    input_gas = ReductionGasComposition(H2=95.0, H2O=5.0)
    dri = DRIComposition(Fe_met=92.0, FeO=3.0, C=0.0, Gangue=5.0)

    # Erzmenge aus Fe-Bilanz berechnen
    ore_amount = calculate_ore_amount_from_fe_h2(ore, dri)
    print(f"Berechnete Erzmenge: {ore_amount:.2f} kg/tDRI")

    # Massenbilanz berechnen (mit typischer Gasmenge von 800 Nm³/tDRI)
    results = calculate_h2_mass_balance(
        ore=ore,
        ore_amount=ore_amount,
        input_gas=input_gas,
        gas_amount_nm3=800.0,
        dri=dri,
    )

    print(f"\n=== Stoffmengen [kmol/tDRI] ===")
    print(f"  Fe (gesamt):       {results.n_Fe:.3f}")
    print(f"  H2 (verbraucht):   {results.n_H2_consumed:.3f}")
    print(f"  O (entfernt):      {results.n_O_removed:.3f}")
    print(f"  H2O (produziert):  {results.n_H2O_produced:.3f}")

    print(f"\n=== Betriebslinie ===")
    print(f"  Steigung k (H2/Fe):        {results.slope_k:.3f}")
    print(f"  Ordinatenabschnitt d:      {results.intercept_d:.3f}")
    print(f"  Geradengleichung: O/Fe = {results.slope_k:.3f} * (O/H2) + ({results.intercept_d:.3f})")

    print(f"\n=== Gaszusammensetzung ===")
    print(f"  Eingang: GOD = {results.god_input*100:.1f}%, O/H2 = {results.o_h2_input:.3f}")
    print(f"  Ausgang: GOD = {results.god_output*100:.1f}%, O/H2 = {results.o_h2_output:.3f}")

    print(f"\n=== Verbrauch ===")
    print(f"  H2: {results.h2_amount_nm3:.1f} Nm³/tDRI")

    # Betriebslinie Punkte
    start, end = calculate_h2_operating_line_points(results)
    print(f"\n=== Betriebslinie Punkte ===")
    print(f"  Start (Fe):    O/H2 = {start[0]:.3f}, O/Fe = {start[1]:.3f}")
    print(f"  Ende (Fe2O3):  O/H2 = {end[0]:.3f}, O/Fe = {end[1]:.3f}")
