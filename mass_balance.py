"""
Massenbilanz-Berechnungen fuer das Rist-Diagramm.

Berechnet die Betriebslinie eines Hochofens basierend auf:
- Erzzusammensetzung (Fe2O3, SiO2, CaO)
- Kokszusammensetzung (C, SiO2)
- Heisswindzusammensetzung (O2, N2)
- Roheisenzusammensetzung (Fe, Si, C)
- Schlackenbasizitaet (B2 = CaO/SiO2)

Die Betriebslinie im Rist-Diagramm:
    O/Fe = k * (O/C) + d

wobei:
    k = C/Fe (Steigung, mol C pro mol Fe)
    d = -O_extern/Fe (Ordinatenabschnitt)
"""

from dataclasses import dataclass
from typing import Tuple, Optional


# Molare Massen [kg/kmol]
MOLAR_MASS = {
    'Fe': 56.0,
    'Fe2O3': 160.0,
    'Si': 28.0,
    'SiO2': 60.0,
    'Ca': 40.0,
    'CaO': 56.0,
    'CaCO3': 100.0,
    'O': 16.0,
    'O2': 32.0,
    'N2': 28.0,
    'C': 12.0,
    'CO': 28.0,
    'CO2': 44.0,
}

# Molares Volumen ideales Gas bei Normalbedingungen [Nm³/kmol]
MOLAR_VOLUME_NM3 = 22.414


@dataclass
class OreComposition:
    """Erzzusammensetzung in Massen-%.
    
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
class CokeComposition:
    """Kokszusammensetzung in Massen-%."""
    C: float = 87.0      # %
    SiO2: float = 13.0   # %

    def validate(self) -> bool:
        """Prueft ob die Summe ~100% ergibt."""
        total = self.C + self.SiO2
        return 99.0 <= total <= 101.0


@dataclass
class HotBlastComposition:
    """Heisswind-Zusammensetzung in Vol-%."""
    O2: float = 21.0     # %
    N2: float = 79.0     # %

    def validate(self) -> bool:
        """Prueft ob die Summe ~100% ergibt."""
        total = self.O2 + self.N2
        return 99.0 <= total <= 101.0


@dataclass
class HotMetalComposition:
    """Roheisen-Zusammensetzung in Massen-%."""
    Fe: float = 95.7     # %
    Si: float = 0.3      # %
    C: float = 4.0       # %

    def validate(self) -> bool:
        """Prueft ob die Summe ~100% ergibt."""
        total = self.Fe + self.Si + self.C
        return 99.0 <= total <= 101.0


@dataclass
class FluxComposition:
    """Zuschlag-Zusammensetzung in Massen-%.

    Typischer Zuschlag ist Kalkstein (CaCO3) oder gebrannter Kalk (CaO).
    """
    CaCO3: float = 100.0  # % Kalkstein
    CaO: float = 0.0      # % gebrannter Kalk (alternativ)

    def validate(self) -> bool:
        """Prueft ob die Summe ~100% ergibt."""
        total = self.CaCO3 + self.CaO
        return 99.0 <= total <= 101.0


@dataclass
class MassBalanceInputs:
    """Eingabedaten fuer die Massenbilanz."""
    ore: OreComposition
    ore_amount: float  # kg/tHM

    coke: CokeComposition
    coke_amount: float  # kg/tHM

    hot_blast: HotBlastComposition
    hot_blast_amount: float  # Nm³/tHM

    hot_metal: HotMetalComposition
    hot_metal_amount: float = 1000.0  # kg/tHM (Bezugsbasis)

    flux: FluxComposition = None  # Zuschlaege
    flux_amount: float = 0.0  # kg/tHM

    slag_basicity_B2: float = 1.0  # CaO/SiO2


@dataclass
class MassBalanceResults:
    """Ergebnisse der Massenbilanz-Berechnung."""
    # Stoffmengen [kmol/tHM]
    n_Fe: float          # Eisen im Roheisen
    n_C_total: float     # Gesamter Kohlenstoff im Koks
    n_C_reduction: float # Kohlenstoff fuer Reduktion (ohne C im Roheisen)
    n_O_wind: float      # Sauerstoff aus Heisswind
    n_O_Si: float        # Sauerstoff fuer Si-Reduktion
    n_O_ore: float       # Sauerstoff aus Erz (aus Fe2O3)
    n_O_extern: float    # Externer Sauerstoff (Wind + Si)

    # Betriebslinien-Parameter
    slope_k: float       # Steigung k = C/Fe
    intercept_d: float   # Ordinatenabschnitt d = -O_extern/Fe

    # Reduktionsanteile
    indirect_reduction_pct: float  # Anteil indirekte Reduktion [%]
    direct_reduction_pct: float    # Anteil direkte Reduktion [%]

    # Gichtgas
    god_gichtgas: float  # GOD am Gichtgas (bei O/Fe = 1.5)
    o_c_gichtgas: float  # O/C am Gichtgas

    # Massen
    ore_amount: float    # Erzmenge [kg/tHM]
    flux_amount: float   # Zuschlagsmenge [kg/tHM]
    slag_amount: float   # Schlackenmenge [kg/tHM]

    # Schlacke Details
    slag_cao: float      # CaO in Schlacke [kg/tHM]
    slag_sio2: float     # SiO2 in Schlacke [kg/tHM]
    slag_basicity_B2: float  # Tatsaechliche Basizitaet
    
    # Erz O/Fe
    o_fe_ore: float      # O/Fe des Erzes (abhängig von Fe2O3/Fe3O4 Mischung)
    
    # Gichtgas Details [NEU]
    n_CO_gichtgas: float = 0.0      # kmol CO im Gichtgas
    n_CO2_gichtgas: float = 0.0     # kmol CO2 im Gichtgas
    n_N2_gichtgas: float = 0.0      # kmol N2 im Gichtgas
    gichtgas_volume_nm3: float = 0.0  # Gichtgasmenge [Nm³/tHM]
    gichtgas_co_pct: float = 0.0    # CO-Anteil im Gichtgas [Vol-%]
    gichtgas_co2_pct: float = 0.0   # CO2-Anteil im Gichtgas [Vol-%]
    gichtgas_n2_pct: float = 0.0    # N2-Anteil im Gichtgas [Vol-%]


def calculate_flux_amount(
    ore: OreComposition,
    ore_amount: float,
    coke: CokeComposition,
    coke_amount: float,
    hot_metal: HotMetalComposition,
    slag_basicity_B2: float,
    flux: 'FluxComposition' = None,
    hot_metal_amount: float = 1000.0,
) -> float:
    """
    Berechnet die benoetigte Zuschlagsmenge basierend auf Basizitaet.

    Die Basizitaet B2 = CaO/SiO2 in der Schlacke.

    Bilanz:
    - SiO2 in Schlacke = SiO2(Erz) + SiO2(Koks) - Si(Roheisen) * M_SiO2/M_Si
    - CaO in Schlacke = CaO(Erz) + CaO(Zuschlag)
    - B2 = CaO_Schlacke / SiO2_Schlacke

    Daraus: CaO_Zuschlag = B2 * SiO2_Schlacke - CaO_Erz

    Args:
        ore: Erzzusammensetzung
        ore_amount: Erzmenge [kg/tHM]
        coke: Kokszusammensetzung
        coke_amount: Koksmenge [kg/tHM]
        hot_metal: Roheisenzusammensetzung
        slag_basicity_B2: Ziel-Basizitaet
        flux: Zuschlagszusammensetzung (default: 100% CaCO3)
        hot_metal_amount: Roheisenmenge [kg/tHM]

    Returns:
        Zuschlagsmenge [kg/tHM]
    """
    if flux is None:
        flux = FluxComposition()

    # SiO2 in Schlacke [kg]
    m_SiO2_ore = ore_amount * ore.SiO2 / 100.0
    m_SiO2_coke = coke_amount * coke.SiO2 / 100.0
    m_Si_HM = hot_metal_amount * hot_metal.Si / 100.0
    # Si -> SiO2: m_SiO2 = m_Si * M_SiO2 / M_Si
    m_SiO2_reduced = m_Si_HM * MOLAR_MASS['SiO2'] / MOLAR_MASS['Si']
    m_SiO2_slag = m_SiO2_ore + m_SiO2_coke - m_SiO2_reduced

    # CaO aus Erz [kg]
    m_CaO_ore = ore_amount * ore.CaO / 100.0

    # Benoetigtes CaO fuer Ziel-Basizitaet [kg]
    m_CaO_required = slag_basicity_B2 * m_SiO2_slag

    # CaO das aus Zuschlag kommen muss [kg]
    m_CaO_from_flux = m_CaO_required - m_CaO_ore

    if m_CaO_from_flux <= 0:
        return 0.0  # Kein Zuschlag noetig, Erz hat genug CaO

    # Zuschlagsmenge berechnen
    # CaCO3 -> CaO + CO2: M_CaO/M_CaCO3 = 56/100 = 0.56
    # Also: m_CaCO3 = m_CaO / 0.56
    cao_from_caco3 = flux.CaCO3 / 100.0 * (MOLAR_MASS['CaO'] / MOLAR_MASS['CaCO3'])
    cao_from_cao = flux.CaO / 100.0
    cao_per_flux = cao_from_caco3 + cao_from_cao

    if cao_per_flux <= 0:
        return 0.0

    flux_amount = m_CaO_from_flux / cao_per_flux
    return flux_amount


def calculate_mass_balance(
    ore: OreComposition,
    ore_amount: float,
    coke: CokeComposition,
    coke_amount: float,
    hot_blast: HotBlastComposition,
    hot_blast_amount: float,
    hot_metal: HotMetalComposition,
    slag_basicity_B2: float = 1.0,
    flux: 'FluxComposition' = None,
    flux_amount: float = None,  # None = automatisch berechnen
    hot_metal_amount: float = 1000.0,
) -> MassBalanceResults:
    """
    Berechnet die Massenbilanz und Betriebslinien-Parameter.

    Args:
        ore: Erzzusammensetzung [Massen-%]
        ore_amount: Erzmenge [kg/tHM]
        coke: Kokszusammensetzung [Massen-%]
        coke_amount: Koksmenge [kg/tHM]
        hot_blast: Heisswind-Zusammensetzung [Vol-%]
        hot_blast_amount: Heisswind-Menge [Nm³/tHM]
        hot_metal: Roheisen-Zusammensetzung [Massen-%]
        slag_basicity_B2: Schlacken-Basizitaet CaO/SiO2
        flux: Zuschlagszusammensetzung (default: 100% CaCO3)
        flux_amount: Zuschlagsmenge [kg/tHM], None = automatisch aus B2 berechnen
        hot_metal_amount: Roheisen-Menge [kg/tHM] (Bezugsbasis, default 1000)

    Returns:
        MassBalanceResults mit allen berechneten Werten
    """
    if flux is None:
        flux = FluxComposition()

    # === 1. Stoffmengen im Roheisen [kmol/tHM] ===
    m_Fe_HM = hot_metal_amount * hot_metal.Fe / 100.0  # kg Fe im Roheisen
    m_Si_HM = hot_metal_amount * hot_metal.Si / 100.0  # kg Si im Roheisen
    m_C_HM = hot_metal_amount * hot_metal.C / 100.0    # kg C im Roheisen

    n_Fe = m_Fe_HM / MOLAR_MASS['Fe']  # kmol Fe
    n_Si_HM = m_Si_HM / MOLAR_MASS['Si']  # kmol Si im Roheisen
    n_C_HM = m_C_HM / MOLAR_MASS['C']  # kmol C im Roheisen

    # === 2. Stoffmengen im Koks [kmol/tHM] ===
    m_C_coke = coke_amount * coke.C / 100.0  # kg C im Koks
    m_SiO2_coke = coke_amount * coke.SiO2 / 100.0  # kg SiO2 im Koks

    n_C_total = m_C_coke / MOLAR_MASS['C']  # kmol C gesamt
    n_SiO2_coke = m_SiO2_coke / MOLAR_MASS['SiO2']  # kmol SiO2 aus Koks

    # Kohlenstoff fuer Reduktion = Gesamt-C minus C im Roheisen
    n_C_reduction = n_C_total - n_C_HM

    # === 3. Stoffmengen im Heisswind [kmol/tHM] ===
    # Umrechnung Nm³ -> kmol: n = V / V_m
    n_wind_total = hot_blast_amount / MOLAR_VOLUME_NM3  # kmol Gesamtgas
    n_O2_wind = n_wind_total * hot_blast.O2 / 100.0  # kmol O2
    n_N2_wind = n_wind_total * hot_blast.N2 / 100.0  # kmol N2

    # Sauerstoff aus Wind (als O-Atome)
    n_O_wind = 2 * n_O2_wind  # kmol O (2 O pro O2)

    # === 4. Stoffmengen im Erz [kmol/tHM] ===
    m_Fe2O3_ore = ore_amount * ore.Fe2O3 / 100.0  # kg Fe2O3
    m_Fe3O4_ore = ore_amount * ore.Fe3O4 / 100.0  # kg Fe3O4
    m_SiO2_ore = ore_amount * ore.SiO2 / 100.0    # kg SiO2
    m_CaO_ore = ore_amount * ore.CaO / 100.0      # kg CaO

    n_Fe2O3_ore = m_Fe2O3_ore / MOLAR_MASS['Fe2O3']  # kmol Fe2O3
    n_Fe3O4_ore = m_Fe3O4_ore / 232.0  # kmol Fe3O4 (M = 3*56 + 4*16 = 232)
    n_SiO2_ore = m_SiO2_ore / MOLAR_MASS['SiO2']    # kmol SiO2
    n_CaO_ore = m_CaO_ore / MOLAR_MASS['CaO']       # kmol CaO

    # Sauerstoff aus Erz (Fe2O3: 3 O, Fe3O4: 4 O)
    n_O_ore = 3 * n_Fe2O3_ore + 4 * n_Fe3O4_ore  # kmol O
    
    # Fe im Erz (Fe2O3: 2 Fe, Fe3O4: 3 Fe)
    n_Fe_ore = 2 * n_Fe2O3_ore + 3 * n_Fe3O4_ore  # kmol Fe
    
    # O/Fe des Erzes berechnen (dynamisch basierend auf Erzzusammensetzung)
    # Fe2O3: O/Fe = 3/2 = 1.5, Fe3O4: O/Fe = 4/3 = 1.333
    o_fe_ore = n_O_ore / n_Fe_ore if n_Fe_ore > 0 else 1.5

    # === 4b. Zuschlaege (CaCO3/CaO) ===
    # SiO2 in Schlacke berechnen (fuer Basizitaet)
    m_Si_HM = hot_metal_amount * hot_metal.Si / 100.0
    m_SiO2_reduced = m_Si_HM * MOLAR_MASS['SiO2'] / MOLAR_MASS['Si']
    m_SiO2_slag = m_SiO2_ore + m_SiO2_coke - m_SiO2_reduced

    # Zuschlagsmenge berechnen wenn nicht angegeben
    if flux_amount is None:
        flux_amount = calculate_flux_amount(
            ore, ore_amount, coke, coke_amount,
            hot_metal, slag_basicity_B2, flux, hot_metal_amount
        )

    # CaO aus Zuschlag
    # CaCO3 -> CaO + CO2
    m_CaCO3_flux = flux_amount * flux.CaCO3 / 100.0
    m_CaO_flux_direct = flux_amount * flux.CaO / 100.0
    m_CaO_from_CaCO3 = m_CaCO3_flux * MOLAR_MASS['CaO'] / MOLAR_MASS['CaCO3']
    m_CaO_flux = m_CaO_from_CaCO3 + m_CaO_flux_direct

    # CO2 aus Kalzinierung (CaCO3 -> CaO + CO2) - geht ins Gichtgas
    n_CaCO3_flux = m_CaCO3_flux / MOLAR_MASS['CaCO3']
    # Dieses CO2 beeinflusst das GOD, aber nicht die Betriebslinie im Rist-Diagramm
    # da es nicht aus der Reduktion stammt

    # Gesamtes CaO in Schlacke
    m_CaO_slag = m_CaO_ore + m_CaO_flux

    # Tatsaechliche Basizitaet
    actual_basicity = m_CaO_slag / m_SiO2_slag if m_SiO2_slag > 0 else 0.0

    # === 5. Si-Reduktion und externer Sauerstoff ===
    # Si im Roheisen kommt aus SiO2-Reduktion
    # SiO2 + 2C -> Si + 2CO
    # Der Sauerstoff aus SiO2 wird zu CO, zaehlt als "externer" O
    n_O_Si = 2 * n_Si_HM  # kmol O fuer Si-Reduktion (2 O pro Si)

    # Externer Sauerstoff = Wind-O + Si-Reduktions-O
    n_O_extern = n_O_wind + n_O_Si

    # === 6. Betriebslinien-Parameter ===
    # Geradengleichung: O/Fe = k * (O/C) + d
    # k = C/Fe (Steigung)
    # d = -O_extern/Fe (Ordinatenabschnitt bei O/C = 0)

    slope_k = n_C_reduction / n_Fe
    intercept_d = -n_O_extern / n_Fe

    # === 7. Reduktionsanteile ===
    # Bei O/C = 1 (100% CO): O/Fe_1 = k * 1 + d = k + d
    # Dieser Wert gibt den Anteil indirekter Reduktion
    # O/Fe = 1.5 ist Haematit (Fe2O3), O/Fe = 0 ist reines Fe
    # Indirekte Reduktion: Anteil der Reduktion durch Gas (O/C >= 1)

    o_fe_at_oc1 = slope_k * 1.0 + intercept_d  # O/Fe bei O/C = 1

    # Gesamte Reduktion: von O/Fe = o_fe_ore (Erz) zu O/Fe = 0 (Fe)
    # Delta O/Fe = o_fe_ore (dynamisch: 1.5 für Fe2O3, 1.333 für Fe3O4)
    # Indirekte Reduktion ist der Teil von O/Fe bei O/C=1 bis O/Fe=o_fe_ore
    # Das entspricht dem Teil der Reduktion durch CO/CO2 Gas

    # Wenn o_fe_at_oc1 < 0, dann beginnt die "Gasphase" erst bei O/C > 1
    # Wenn o_fe_at_oc1 > 0, ist bei O/C=1 bereits etwas reduziert

    # Der Schnittpunkt mit O/Fe = o_fe_ore gibt den O/C-Wert am Gichtgas
    # o_fe_ore = k * o_c + d  =>  o_c = (o_fe_ore - d) / k
    o_c_gichtgas = (o_fe_ore - intercept_d) / slope_k

    # GOD = (O/C - 1) fuer O/C >= 1
    god_gichtgas = (o_c_gichtgas - 1.0) * 100.0 if o_c_gichtgas >= 1.0 else 0.0

    # Anteil Reduktion:
    # - Direkte Reduktion: Anteil von O/Fe=0 bis o_fe_at_oc1 (mit festem C)
    # - Indirekte Reduktion: Anteil von o_fe_at_oc1 bis O/Fe=o_fe_ore (mit Gas CO/CO2)
    #
    # Direkte Reduktion = o_fe_at_oc1 / o_fe_ore * 100%
    # Indirekte Reduktion = (o_fe_ore - o_fe_at_oc1) / o_fe_ore * 100%
    if o_fe_at_oc1 > 0:
        direct_reduction_pct = (o_fe_at_oc1 / o_fe_ore) * 100.0
    else:
        direct_reduction_pct = 0.0

    direct_reduction_pct = min(100.0, max(0.0, direct_reduction_pct))
    indirect_reduction_pct = 100.0 - direct_reduction_pct

    # === 8. Schlackenmenge ===
    # Schlacke besteht aus CaO und SiO2
    # m_SiO2_slag und m_CaO_slag wurden bereits in Schritt 4b berechnet
    slag_amount = m_SiO2_slag + m_CaO_slag

    # === 9. Gichtgas-Berechnung ===
    # Das Gichtgas besteht aus CO, CO2 und N2
    
    # N2 kommt unverändert aus dem Heißwind
    n_N2_gichtgas = n_N2_wind
    
    # Gesamter Kohlenstoff im Gichtgas = C für Reduktion (wird zu CO/CO2)
    n_C_gas = n_C_reduction
    
    # O/C Verhältnis am Gichtgas bestimmt das CO/CO2 Verhältnis
    # O/C = (n_CO + 2*n_CO2) / (n_CO + n_CO2)
    # Aus o_c_gichtgas können wir CO und CO2 berechnen:
    # n_CO + n_CO2 = n_C_gas
    # n_CO + 2*n_CO2 = o_c_gichtgas * n_C_gas
    # => n_CO2 = (o_c_gichtgas - 1) * n_C_gas
    # => n_CO = n_C_gas - n_CO2 = (2 - o_c_gichtgas) * n_C_gas
    
    n_CO2_gichtgas = (o_c_gichtgas - 1.0) * n_C_gas if o_c_gichtgas > 1.0 else 0.0
    n_CO_gichtgas = n_C_gas - n_CO2_gichtgas
    
    # CO2 aus Kalzinierung addieren (CaCO3 -> CaO + CO2)
    n_CO2_gichtgas += n_CaCO3_flux
    
    # Gesamte Stoffmenge im Gichtgas
    n_gichtgas_total = n_CO_gichtgas + n_CO2_gichtgas + n_N2_gichtgas
    
    # Gichtgas-Zusammensetzung [Vol-%]
    if n_gichtgas_total > 0:
        gichtgas_co_pct = n_CO_gichtgas / n_gichtgas_total * 100.0
        gichtgas_co2_pct = n_CO2_gichtgas / n_gichtgas_total * 100.0
        gichtgas_n2_pct = n_N2_gichtgas / n_gichtgas_total * 100.0
    else:
        gichtgas_co_pct = 0.0
        gichtgas_co2_pct = 0.0
        gichtgas_n2_pct = 0.0
    
    # Gichtgas-Volumen [Nm³/tHM]
    gichtgas_volume_nm3 = n_gichtgas_total * MOLAR_VOLUME_NM3

    return MassBalanceResults(
        n_Fe=n_Fe,
        n_C_total=n_C_total,
        n_C_reduction=n_C_reduction,
        n_O_wind=n_O_wind,
        n_O_Si=n_O_Si,
        n_O_ore=n_O_ore,
        n_O_extern=n_O_extern,
        slope_k=slope_k,
        intercept_d=intercept_d,
        indirect_reduction_pct=indirect_reduction_pct,
        direct_reduction_pct=direct_reduction_pct,
        god_gichtgas=god_gichtgas,
        o_c_gichtgas=o_c_gichtgas,
        ore_amount=ore_amount,
        flux_amount=flux_amount,
        slag_amount=slag_amount,
        slag_cao=m_CaO_slag,
        slag_sio2=m_SiO2_slag,
        slag_basicity_B2=actual_basicity,
        o_fe_ore=o_fe_ore,
        n_CO_gichtgas=n_CO_gichtgas,
        n_CO2_gichtgas=n_CO2_gichtgas,
        n_N2_gichtgas=n_N2_gichtgas,
        gichtgas_volume_nm3=gichtgas_volume_nm3,
        gichtgas_co_pct=gichtgas_co_pct,
        gichtgas_co2_pct=gichtgas_co2_pct,
        gichtgas_n2_pct=gichtgas_n2_pct,
    )


def calculate_operating_line_points(
    results: MassBalanceResults,
    o_fe_start: float = 0.0,
    o_fe_end: Optional[float] = None,
) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """
    Berechnet die Endpunkte der Betriebslinie.

    Die Betriebslinie: O/Fe = k * (O/C) + d
    Umgestellt: O/C = (O/Fe - d) / k

    Args:
        results: Ergebnisse der Massenbilanz
        o_fe_start: Start O/Fe-Wert (default: 0 = reines Fe)
        o_fe_end: End O/Fe-Wert (default: None = aus Erzzusammensetzung)

    Returns:
        ((o_c_start, o_fe_start), (o_c_end, o_fe_end))
    """
    k = results.slope_k
    d = results.intercept_d
    
    # Wenn o_fe_end nicht angegeben, aus Erz-O/Fe nehmen
    if o_fe_end is None:
        o_fe_end = results.o_fe_ore

    # O/C = (O/Fe - d) / k
    o_c_start = (o_fe_start - d) / k
    o_c_end = (o_fe_end - d) / k

    return ((o_c_start, o_fe_start), (o_c_end, o_fe_end))


def calculate_ore_amount_from_fe(
    ore: OreComposition,
    hot_metal: HotMetalComposition,
    hot_metal_amount: float = 1000.0,
) -> float:
    """
    Berechnet die benoetigte Erzmenge basierend auf Fe-Bilanz.

    Args:
        ore: Erzzusammensetzung
        hot_metal: Roheisenzusammensetzung
        hot_metal_amount: Roheisenmenge [kg/tHM]

    Returns:
        Erzmenge [kg/tHM]
    """
    # Fe im Roheisen [kg]
    m_Fe_HM = hot_metal_amount * hot_metal.Fe / 100.0

    # Fe-Gehalt im Erz [kg Fe / kg Erz]
    # Fe2O3 -> 2 Fe, Fe-Anteil = 2*56 / 160 = 0.7
    # Fe3O4 -> 3 Fe, Fe-Anteil = 3*56 / 232 = 0.724
    fe_in_fe2o3 = 2 * MOLAR_MASS['Fe'] / MOLAR_MASS['Fe2O3']
    fe_in_fe3o4 = 3 * MOLAR_MASS['Fe'] / 232.0  # M_Fe3O4 = 232
    fe_content_ore = (ore.Fe2O3 / 100.0 * fe_in_fe2o3 + 
                      ore.Fe3O4 / 100.0 * fe_in_fe3o4)

    # Erzmenge = Fe im Roheisen / Fe-Gehalt im Erz
    if fe_content_ore > 0:
        ore_amount = m_Fe_HM / fe_content_ore
    else:
        ore_amount = 0.0

    return ore_amount


if __name__ == "__main__":
    print("=== Massenbilanz Test ===\n")

    # Standardwerte aus Excel
    ore = OreComposition(Fe2O3=94.4, SiO2=2.5, CaO=3.1)
    coke = CokeComposition(C=87.0, SiO2=13.0)
    hot_blast = HotBlastComposition(O2=21.0, N2=79.0)
    hot_metal = HotMetalComposition(Fe=95.7, Si=0.3, C=4.0)

    # Erzmenge aus Fe-Bilanz berechnen
    ore_amount = calculate_ore_amount_from_fe(ore, hot_metal)
    print(f"Berechnete Erzmenge: {ore_amount:.2f} kg/tHM")

    # Massenbilanz berechnen
    results = calculate_mass_balance(
        ore=ore,
        ore_amount=ore_amount,
        coke=coke,
        coke_amount=500.0,
        hot_blast=hot_blast,
        hot_blast_amount=1200.0,
        hot_metal=hot_metal,
    )

    print(f"\n=== Stoffmengen [kmol/tHM] ===")
    print(f"  Fe (Roheisen):     {results.n_Fe:.3f}")
    print(f"  C (gesamt):        {results.n_C_total:.3f}")
    print(f"  C (Reduktion):     {results.n_C_reduction:.3f}")
    print(f"  O (Wind):          {results.n_O_wind:.3f}")
    print(f"  O (Si-Reduktion):  {results.n_O_Si:.3f}")
    print(f"  O (extern):        {results.n_O_extern:.3f}")

    print(f"\n=== Betriebslinie ===")
    print(f"  Steigung k (C/Fe):         {results.slope_k:.3f}")
    print(f"  Ordinatenabschnitt d:      {results.intercept_d:.3f}")
    print(f"  Geradengleichung: O/Fe = {results.slope_k:.3f} * (O/C) + ({results.intercept_d:.3f})")

    print(f"\n=== Reduktionsanteile ===")
    print(f"  Indirekte Reduktion: {results.indirect_reduction_pct:.1f}%")
    print(f"  Direkte Reduktion:   {results.direct_reduction_pct:.1f}%")

    print(f"\n=== Gichtgas ===")
    print(f"  O/C am Gichtgas:  {results.o_c_gichtgas:.3f}")
    print(f"  GOD am Gichtgas:  {results.god_gichtgas:.1f}%")

    # Betriebslinie Punkte
    start, end = calculate_operating_line_points(results)
    print(f"\n=== Betriebslinie Punkte ===")
    print(f"  Start (Fe):    O/C = {start[0]:.3f}, O/Fe = {start[1]:.3f}")
    print(f"  Ende (Fe2O3):  O/C = {end[0]:.3f}, O/Fe = {end[1]:.3f}")
