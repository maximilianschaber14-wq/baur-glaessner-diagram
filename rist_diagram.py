"""
Rist-Diagramm Berechnungen fuer Eisenoxid-Reduktionsprozesse.

Das Rist-Diagramm stellt die Massenbilanz von Reduktionsprozessen dar.
Achsen:
- Y-Achse: O/Fe Verhaeltnis (Sauerstoff pro Mol Eisen)
- X-Achse: O/C Verhaeltnis (Sauerstoff pro Mol Kohlenstoff)

O/Fe Werte:
- Fe (Eisen): 0
- FeO (Wuestit): ~1.05 (nicht-stoechiometrisch wegen Eisenleerstellen)
- Fe3O4 (Magnetit): 4/3 = 1.333
- Fe2O3 (Haematit): 3/2 = 1.5

O/C Werte:
- O/C = 0: Direktreduktion (reiner Kohlenstoff, kein Gas)
- O/C = 1: 100% CO
- O/C = 2: 100% CO2

Der Wuestit-Punkt (W) und Magnetit-Punkt (M) werden aus dem Baur-Glaessner-Diagramm
abgeleitet und definieren die thermodynamischen Grenzen.

Quellen:
- Rist, A. (1963-1967): Blast Furnace Operating Diagram
- Spreitzer & Schenk (2019): Steel Research International
"""

from typing import Tuple, Optional, List
from thermodynamics import BaurGlassnerThermo, WUSTITE_STABILITY_K


# O/Fe Verhaeltnisse fuer die Eisenoxide
O_FE_RATIOS = {
    'Fe': 0.0,
    'FeO': 1.05,      # Nicht-stoechiometrisch (Fe_{1-x}O)
    'Fe3O4': 4/3,     # = 1.333
    'Fe2O3': 3/2,     # = 1.5
}


def god_to_o_c(god: float) -> float:
    """
    Konvertiert GOD (Gas Oxidation Degree) zu O/C Verhaeltnis.

    GOD = CO2 / (CO + CO2)

    Fuer CO: O/C = 1
    Fuer CO2: O/C = 2

    Bei GOD = 0: 100% CO -> O/C = 1
    Bei GOD = 1: 100% CO2 -> O/C = 2

    O/C = (n_CO * 1 + n_CO2 * 2) / (n_CO + n_CO2)
        = (1-GOD)*1 + GOD*2
        = 1 + GOD

    Args:
        god: Gas Oxidation Degree (0-1)

    Returns:
        O/C Verhaeltnis (1-2)
    """
    return 1.0 + god


def o_c_to_god(o_c: float) -> float:
    """
    Konvertiert O/C Verhaeltnis zu GOD.

    Nur gueltig fuer O/C zwischen 1 und 2 (Gasphase).
    Fuer O/C < 1 ist GOD nicht definiert (Direktreduktion).

    Args:
        o_c: O/C Verhaeltnis (1-2)

    Returns:
        GOD (0-1), oder None wenn O/C < 1
    """
    if o_c < 1.0:
        return None  # Direktreduktionsbereich
    return o_c - 1.0


class RistDiagram:
    """
    Klasse zur Berechnung des Rist-Diagramms.

    Das Rist-Diagramm zeigt:
    - Horizontale Linien fuer die Eisenoxide bei ihren O/Fe Werten
    - Den Wuestit-Punkt (W) und Magnetit-Punkt (M) aus dem Baur-Glaessner-Diagramm
    - Die Betriebslinie eines idealen Gegenstromreaktors
    - Den "thermodynamisch verbotenen" Bereich
    """

    def __init__(self, thermo: BaurGlassnerThermo):
        """
        Initialisiert das Rist-Diagramm.

        Args:
            thermo: BaurGlassnerThermo-Instanz fuer thermodynamische Berechnungen
        """
        self.thermo = thermo

    def get_wustite_point(self, T_kelvin: float) -> Tuple[float, float]:
        """
        Berechnet den Wuestit-Punkt (W) fuer eine gegebene Temperatur.

        Der Wuestit-Punkt liegt auf der Fe/FeO Gleichgewichtslinie
        bei O/Fe = 1.05 (Wuestit).

        Args:
            T_kelvin: Temperatur in Kelvin

        Returns:
            (O/C, O/Fe) Koordinaten des Wuestit-Punkts
        """
        if T_kelvin < WUSTITE_STABILITY_K:
            # Unterhalb 570C: Verwende Fe/Fe3O4 Gleichgewicht
            god = self.thermo.GOD('Fe3O4_Fe_CO', T_kelvin)
            o_c = god_to_o_c(god) if god == god else 1.5  # NaN check
            return (o_c, O_FE_RATIOS['FeO'])
        else:
            # Oberhalb 570C: Verwende Fe/FeO Gleichgewicht
            god = self.thermo.GOD('FeO_Fe_CO', T_kelvin)
            o_c = god_to_o_c(god)
            return (o_c, O_FE_RATIOS['FeO'])

    def get_magnetite_point(self, T_kelvin: float) -> Tuple[float, float]:
        """
        Berechnet den Magnetit-Punkt (M) fuer eine gegebene Temperatur.

        Der Magnetit-Punkt liegt auf der FeO/Fe3O4 Gleichgewichtslinie
        bei O/Fe = 1.333 (Magnetit).

        Args:
            T_kelvin: Temperatur in Kelvin

        Returns:
            (O/C, O/Fe) Koordinaten des Magnetit-Punkts
        """
        if T_kelvin < WUSTITE_STABILITY_K:
            # Unterhalb 570C: FeO ist instabil, M-Punkt existiert nicht sinnvoll
            god = self.thermo.GOD('Fe3O4_Fe_CO', T_kelvin)
            o_c = god_to_o_c(god) if god == god else 1.5
            return (o_c, O_FE_RATIOS['Fe3O4'])
        else:
            # Oberhalb 570C: Verwende FeO/Fe3O4 Gleichgewicht
            god = self.thermo.GOD('Fe3O4_FeO_CO', T_kelvin)
            o_c = god_to_o_c(god)
            return (o_c, O_FE_RATIOS['Fe3O4'])

    def get_ideal_operating_line(
        self,
        T_kelvin: float,
        input_o_fe: float = 1.5,  # Haematit als Eingangsmaterial
        output_o_fe: float = 0.0,  # Reines Eisen als Ausgang
    ) -> Tuple[Tuple[float, float], Tuple[float, float], float]:
        """
        Berechnet die Betriebslinie eines Gegenstromreaktors.

        Die Betriebslinie im Rist-Diagramm:
        - Geht durch den Wuestit-Punkt (W) - das ist die thermodynamische Grenze
        - Die Steigung wird durch den W-Punkt bestimmt
        - Der C/Fe-Verbrauch (Reduktionsmittelverbrauch) ergibt sich aus der Steigung

        Physikalische Bedeutung im Gegenstromreaktor (Hochofen):
        - UNTEN: Gas tritt als CO ein (O/C = 1), Eisen (Fe) kommt raus (O/Fe = 0)
        - OBEN: Erz (Fe2O3) kommt rein (O/Fe = 1.5), verbrauchtes Gas (CO/CO2) geht raus

        Die Betriebslinie zeigt die Massenbilanz:
        - Sie muss durch den W-Punkt gehen (thermodynamischer Engpass bei FeO/Fe)
        - Die Steigung m = dO/Fe / dO/C
        - C/Fe = 1/m (Kohlenstoffverbrauch pro mol Eisen)

        Berechnung:
        - Unten: Gas-Eingang bei O/C = 1 (100% CO), Produkt O/Fe = 0 (Fe)
        - Die Linie geht durch W und durch den Punkt (1, 0) am unteren Ende
        - Dann wird der Endpunkt oben berechnet wo O/Fe = 1.5

        Args:
            T_kelvin: Temperatur in Kelvin
            input_o_fe: O/Fe des Eingangsmaterials (Standard: Haematit = 1.5)
            output_o_fe: O/Fe des Ausgangsmaterials (Standard: Fe = 0)

        Returns:
            ((start_o_c, start_o_fe), (end_o_c, end_o_fe), c_per_fe)
        """
        # Wuestit-Punkt bestimmen (aus Baur-Glaessner-Diagramm)
        w_point = self.get_wustite_point(T_kelvin)
        w_o_c, w_o_fe = w_point

        # Die Betriebslinie im Rist-Diagramm:
        #
        # Die Betriebslinie zeigt die Massenbilanz im Gegenstromreaktor.
        # Sie MUSS durch den W-Punkt gehen (thermodynamischer Engpass).
        #
        # Die "ideale" oder "minimale" Betriebslinie ist diejenige, die
        # den geringsten Kohlenstoffverbrauch hat. Diese geht durch:
        # - Den Ursprung (0, 0): Hier waere das Produkt Fe mit 100% C (Direktreduktion)
        # - Den W-Punkt: Thermodynamischer Engpass
        #
        # Die Steigung dieser Linie ist: m = w_o_fe / w_o_c
        #
        # Diese Linie verlaengert bis O/Fe = 1.5 ergibt den Endpunkt.
        #
        # ACHTUNG: Im Bild das der User gezeigt hat, ist die Betriebslinie
        # NICHT durch den Ursprung, sondern hat einen Startpunkt bei O/C > 0.
        # Das entspricht einem realen Hochofenprozess mit einer Mischung aus
        # Direktreduktion (C) und indirekter Reduktion (CO).

        # Fuer die ideale (minimale) Betriebslinie:
        # Steigung durch Ursprung und W
        slope = w_o_fe / w_o_c

        # Start bei O/Fe = 0:
        # Die Linie geht durch den Ursprung (0, 0) mit Steigung slope
        # Also: O/Fe = slope * O/C
        # Bei O/Fe = 0: start_o_c = 0
        start_o_c = 0.0
        start_o_fe = output_o_fe

        # Ende bei O/Fe = 1.5:
        # 1.5 = slope * end_o_c
        # end_o_c = 1.5 / slope = 1.5 * w_o_c / w_o_fe
        end_o_c = input_o_fe / slope
        end_o_fe = input_o_fe

        # C/Fe Verbrauch = horizontale Distanz / vertikale Distanz
        c_per_fe = (end_o_c - start_o_c) / (end_o_fe - start_o_fe)

        start_point = (start_o_c, start_o_fe)
        end_point = (end_o_c, end_o_fe)

        return (start_point, end_point, c_per_fe)

    def get_forbidden_zone_w(self, T_kelvin: float) -> List[Tuple[float, float]]:
        """
        Berechnet die Eckpunkte des thermodynamisch verbotenen Bereichs
        beim Wuestit-Punkt (W).

        Der verbotene Bereich liegt rechts vom Wuestit-Punkt im Bereich
        zwischen der Fe-Linie (O/Fe=0) und der Wuestit-Linie (O/Fe=1.05).

        Args:
            T_kelvin: Temperatur in Kelvin

        Returns:
            Liste von (O/C, O/Fe) Punkten
        """
        w_point = self.get_wustite_point(T_kelvin)
        w_o_c, w_o_fe = w_point

        # Der verbotene Bereich beim W-Punkt ist ein Viereck:
        # - W-Punkt
        # - Rechts oben bei Wuestit-Linie (O/C=2, O/Fe=1.05)
        # - Rechts unten bei Fe-Linie (O/C=2, O/Fe=0)
        # - Unter dem W-Punkt (w_o_c, O/Fe=0)
        vertices = [
            w_point,                    # W-Punkt
            (2.0, O_FE_RATIOS['FeO']),  # Rechts oben (bei Wuestit-Linie)
            (2.0, O_FE_RATIOS['Fe']),   # Rechts unten (bei Fe-Linie)
            (w_o_c, O_FE_RATIOS['Fe']), # Unter dem W-Punkt
        ]

        return vertices

    def get_forbidden_zone_m(self, T_kelvin: float) -> List[Tuple[float, float]]:
        """
        Berechnet die Eckpunkte des thermodynamisch verbotenen Bereichs
        beim Magnetit-Punkt (M).

        Der verbotene Bereich liegt rechts vom Magnetit-Punkt im Bereich
        zwischen der Wuestit-Linie (O/Fe=1.05) und der Magnetit-Linie (O/Fe=1.33).

        Args:
            T_kelvin: Temperatur in Kelvin

        Returns:
            Liste von (O/C, O/Fe) Punkten
        """
        m_point = self.get_magnetite_point(T_kelvin)
        m_o_c, m_o_fe = m_point

        # Der verbotene Bereich beim M-Punkt:
        # - M-Punkt
        # - Rechts oben bei Magnetit-Linie (O/C=2, O/Fe=1.33)
        # - Rechts unten bei Wuestit-Linie (O/C=2, O/Fe=1.05)
        # - Unter dem M-Punkt auf der Wuestit-Linie (m_o_c, O/Fe=1.05)
        vertices = [
            m_point,                        # M-Punkt
            (2.0, O_FE_RATIOS['Fe3O4']),    # Rechts oben (bei Magnetit-Linie)
            (2.0, O_FE_RATIOS['FeO']),      # Rechts bei Wuestit-Linie
            (m_o_c, O_FE_RATIOS['FeO']),    # Unter dem M-Punkt
        ]

        return vertices

    def get_forbidden_zone_combined(self, T_kelvin: float) -> List[Tuple[float, float]]:
        """
        Berechnet den kombinierten thermodynamisch verbotenen Bereich.

        Dieser umfasst:
        - Den Bereich rechts vom W-Punkt (Fe -> FeO Reduktion)
        - Den Bereich rechts vom M-Punkt (FeO -> Fe3O4 Reduktion)

        Args:
            T_kelvin: Temperatur in Kelvin

        Returns:
            Liste von (O/C, O/Fe) Punkten fuer den gesamten verbotenen Bereich
        """
        w_point = self.get_wustite_point(T_kelvin)
        m_point = self.get_magnetite_point(T_kelvin)
        w_o_c, w_o_fe = w_point
        m_o_c, m_o_fe = m_point

        # Kombinierter verbotener Bereich:
        # Beginnt beim W-Punkt, geht hoch zum M-Punkt, dann nach rechts
        vertices = [
            (w_o_c, O_FE_RATIOS['Fe']),     # Unter W-Punkt auf Fe-Linie
            w_point,                         # W-Punkt
            (m_o_c, O_FE_RATIOS['FeO']),    # Unter M-Punkt auf FeO-Linie
            m_point,                         # M-Punkt
            (2.0, O_FE_RATIOS['Fe3O4']),    # Rechter Rand bei Fe3O4
            (2.0, O_FE_RATIOS['Fe']),       # Rechter Rand bei Fe
        ]

        return vertices

    def calculate_carbon_consumption(
        self,
        T_kelvin: float,
        input_oxide: str = 'Fe2O3',
        output_state: str = 'Fe',
    ) -> float:
        """
        Berechnet den theoretischen Kohlenstoffverbrauch (mol C / mol Fe).

        Args:
            T_kelvin: Temperatur in Kelvin
            input_oxide: Eingangsmaterial ('Fe2O3', 'Fe3O4', 'FeO')
            output_state: Ausgangszustand ('Fe', 'FeO', 'Fe3O4')

        Returns:
            Kohlenstoffverbrauch in mol C pro mol Fe
        """
        _, _, c_per_fe = self.get_ideal_operating_line(
            T_kelvin,
            input_o_fe=O_FE_RATIOS[input_oxide],
            output_o_fe=O_FE_RATIOS[output_state],
        )
        return c_per_fe


if __name__ == "__main__":
    # Test
    print("=== Rist-Diagramm Test ===\n")

    thermo = BaurGlassnerThermo()
    rist = RistDiagram(thermo)

    for T_C in [700, 800, 900, 1000]:
        T_K = T_C + 273.15
        print(f"\nTemperatur: {T_C} C")

        w_point = rist.get_wustite_point(T_K)
        m_point = rist.get_magnetite_point(T_K)

        print(f"  Wuestit-Punkt (W): O/C = {w_point[0]:.3f}, O/Fe = {w_point[1]:.3f}")
        print(f"  Magnetit-Punkt (M): O/C = {m_point[0]:.3f}, O/Fe = {m_point[1]:.3f}")

        # GOD zurueckrechnen (nur fuer O/C >= 1)
        god_w = o_c_to_god(w_point[0])
        god_m = o_c_to_god(m_point[0])
        if god_w is not None:
            print(f"  GOD am W-Punkt: {god_w:.3f} ({(1-god_w)*100:.1f}% CO)")
        if god_m is not None:
            print(f"  GOD am M-Punkt: {god_m:.3f} ({(1-god_m)*100:.1f}% CO)")

        # Betriebslinie
        start, end, c_fe = rist.get_ideal_operating_line(T_K)
        print(f"  Ideale Betriebslinie: ({start[0]:.2f}, {start[1]:.2f}) -> ({end[0]:.2f}, {end[1]:.2f})")
        print(f"  C-Verbrauch: {c_fe:.3f} mol C / mol Fe")
