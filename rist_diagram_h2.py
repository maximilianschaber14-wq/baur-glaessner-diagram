"""
Rist-Diagramm Berechnungen fuer Wasserstoff-Direktreduktion (H2-DR).

Das H2-Rist-Diagramm stellt die Massenbilanz von H2-Reduktionsprozessen dar.
Achsen:
- Y-Achse: O/Fe Verhaeltnis (Sauerstoff pro Mol Eisen)
- X-Achse: O/H2 Verhaeltnis (Sauerstoff pro Mol Wasserstoff)
  - O/H2 = 1: 100% H2
  - O/H2 = 2: 100% H2O

Die W- und M-Punkte werden aus dem Baur-Glaessner-Diagramm (H2-System) uebernommen.
"""

from typing import Tuple, List
import numpy as np

from thermodynamics import BaurGlassnerThermo, WUSTITE_STABILITY_K


# O/Fe Verhaeltnisse fuer die Eisenoxide (gleich wie CO-System)
O_FE_RATIOS = {
    'Fe': 0.0,
    'FeO': 1.05,      # Nicht-stoechiometrisch (Fe_{1-x}O)
    'Fe3O4': 4/3,     # = 1.333
    'Fe2O3': 3/2,     # = 1.5
}


def god_to_o_h2(god: float) -> float:
    """
    Konvertiert GOD (Gas Oxidation Degree) zu O/H2 Verhaeltnis.

    GOD = H2O / (H2 + H2O)

    Fuer H2: O/H2 = 1
    Fuer H2O: O/H2 = 2 (weil H2O ein O und 2 H hat, also 1 O pro H2-Aequivalent)

    Bei GOD = 0: 100% H2 -> O/H2 = 1
    Bei GOD = 1: 100% H2O -> O/H2 = 2

    O/H2 = (n_H2 * 1 + n_H2O * 2) / (n_H2 + n_H2O)
         = (1-GOD)*1 + GOD*2
         = 1 + GOD

    Args:
        god: Gas Oxidation Degree (0-1)

    Returns:
        O/H2 Verhaeltnis (1-2)
    """
    return 1.0 + god


def o_h2_to_god(o_h2: float) -> float:
    """
    Konvertiert O/H2 Verhaeltnis zu GOD.

    Nur gueltig fuer O/H2 zwischen 1 und 2 (Gasphase).
    Fuer O/H2 < 1 ist GOD nicht definiert.

    Args:
        o_h2: O/H2 Verhaeltnis (1-2)

    Returns:
        GOD (0-1), oder None wenn O/H2 < 1
    """
    if o_h2 < 1.0:
        return None
    return o_h2 - 1.0


class RistDiagramH2:
    """
    Klasse zur Berechnung des H2-Rist-Diagramms.

    Das H2-Rist-Diagramm zeigt:
    - Horizontale Linien fuer die Eisenoxide bei ihren O/Fe Werten
    - Den Wuestit-Punkt (W) und Magnetit-Punkt (M) aus dem Baur-Glaessner-Diagramm (H2-System)
    - Die Betriebslinie basierend auf der Massenbilanz
    """

    def __init__(self, thermo: BaurGlassnerThermo):
        """
        Initialisiert das H2-Rist-Diagramm.

        Args:
            thermo: BaurGlassnerThermo-Instanz fuer thermodynamische Berechnungen
        """
        self.thermo = thermo

    def get_wustite_point(self, T_kelvin: float) -> Tuple[float, float]:
        """
        Berechnet den Wuestit-Punkt (W) fuer eine gegebene Temperatur (H2-System).

        Der Wuestit-Punkt liegt auf der Fe/FeO Gleichgewichtslinie
        bei O/Fe = 1.05 (Wuestit).

        Args:
            T_kelvin: Temperatur in Kelvin

        Returns:
            (O/H2, O/Fe) Koordinaten des Wuestit-Punkts
        """
        if T_kelvin < WUSTITE_STABILITY_K:
            # Unterhalb 570C: Verwende Fe/Fe3O4 Gleichgewicht
            god = self.thermo.GOD('Fe3O4_Fe_H2', T_kelvin)
            o_h2 = god_to_o_h2(god) if god == god else 1.5  # NaN check
            return (o_h2, O_FE_RATIOS['FeO'])
        else:
            # Oberhalb 570C: Verwende Fe/FeO Gleichgewicht
            god = self.thermo.GOD('FeO_Fe_H2', T_kelvin)
            o_h2 = god_to_o_h2(god)
            return (o_h2, O_FE_RATIOS['FeO'])

    def get_magnetite_point(self, T_kelvin: float) -> Tuple[float, float]:
        """
        Berechnet den Magnetit-Punkt (M) fuer eine gegebene Temperatur (H2-System).

        Der Magnetit-Punkt liegt auf der FeO/Fe3O4 Gleichgewichtslinie
        bei O/Fe = 1.333 (Magnetit).

        Args:
            T_kelvin: Temperatur in Kelvin

        Returns:
            (O/H2, O/Fe) Koordinaten des Magnetit-Punkts
        """
        if T_kelvin < WUSTITE_STABILITY_K:
            # Unterhalb 570C: FeO ist instabil, M-Punkt existiert nicht sinnvoll
            god = self.thermo.GOD('Fe3O4_Fe_H2', T_kelvin)
            o_h2 = god_to_o_h2(god) if god == god else 1.5
            return (o_h2, O_FE_RATIOS['Fe3O4'])
        else:
            # Oberhalb 570C: Verwende FeO/Fe3O4 Gleichgewicht
            god = self.thermo.GOD('Fe3O4_FeO_H2', T_kelvin)
            o_h2 = god_to_o_h2(god)
            return (o_h2, O_FE_RATIOS['Fe3O4'])

    def get_forbidden_zone_w(self, T_kelvin: float) -> List[Tuple[float, float]]:
        """
        Berechnet die Eckpunkte des thermodynamisch verbotenen Bereichs
        beim Wuestit-Punkt (W) im H2-System.

        Der verbotene Bereich liegt rechts vom Wuestit-Punkt im Bereich
        zwischen der Fe-Linie (O/Fe=0) und der Wuestit-Linie (O/Fe=1.05).

        Args:
            T_kelvin: Temperatur in Kelvin

        Returns:
            Liste von (O/H2, O/Fe) Punkten
        """
        w_point = self.get_wustite_point(T_kelvin)
        o_h2_w = w_point[0]

        # Verbotener Bereich: Rechts vom W-Punkt bis O/H2 = 2
        return [
            (o_h2_w, O_FE_RATIOS['Fe']),      # W-Punkt auf Fe-Linie
            (2.0, O_FE_RATIOS['Fe']),          # Rechte untere Ecke
            (2.0, O_FE_RATIOS['FeO']),         # Rechte obere Ecke
            (o_h2_w, O_FE_RATIOS['FeO']),     # W-Punkt
        ]

    def get_forbidden_zone_m(self, T_kelvin: float) -> List[Tuple[float, float]]:
        """
        Berechnet die Eckpunkte des thermodynamisch verbotenen Bereichs
        beim Magnetit-Punkt (M) im H2-System.

        Der verbotene Bereich liegt rechts vom Magnetit-Punkt im Bereich
        zwischen der Wuestit-Linie (O/Fe=1.05) und der Magnetit-Linie (O/Fe=1.33).

        Args:
            T_kelvin: Temperatur in Kelvin

        Returns:
            Liste von (O/H2, O/Fe) Punkten
        """
        m_point = self.get_magnetite_point(T_kelvin)
        o_h2_m = m_point[0]

        return [
            (o_h2_m, O_FE_RATIOS['FeO']),     # M-Punkt auf FeO-Linie
            (2.0, O_FE_RATIOS['FeO']),         # Rechte untere Ecke
            (2.0, O_FE_RATIOS['Fe3O4']),       # Rechte obere Ecke
            (o_h2_m, O_FE_RATIOS['Fe3O4']),   # M-Punkt
        ]

    def get_forbidden_zone_combined(self, T_kelvin: float) -> List[Tuple[float, float]]:
        """
        Berechnet den kombinierten thermodynamisch verbotenen Bereich.

        Dieser umfasst:
        - Den Bereich rechts vom W-Punkt (Fe -> FeO Reduktion)
        - Den Bereich rechts vom M-Punkt (FeO -> Fe3O4 Reduktion)

        Args:
            T_kelvin: Temperatur in Kelvin

        Returns:
            Liste von (O/H2, O/Fe) Punkten fuer den gesamten verbotenen Bereich
        """
        w_point = self.get_wustite_point(T_kelvin)
        m_point = self.get_magnetite_point(T_kelvin)

        o_h2_w = w_point[0]
        o_h2_m = m_point[0]

        # Kombiniertes Polygon
        if T_kelvin >= WUSTITE_STABILITY_K:
            # Oberhalb 570C: W und M haben unterschiedliche X-Werte
            return [
                (o_h2_w, O_FE_RATIOS['Fe']),      # W auf Fe-Linie
                (2.0, O_FE_RATIOS['Fe']),          # Untere rechte Ecke
                (2.0, O_FE_RATIOS['Fe3O4']),       # Obere rechte Ecke
                (o_h2_m, O_FE_RATIOS['Fe3O4']),   # M-Punkt
                (o_h2_m, O_FE_RATIOS['FeO']),     # M auf FeO-Linie
                (o_h2_w, O_FE_RATIOS['FeO']),     # W-Punkt
            ]
        else:
            # Unterhalb 570C: Nur eine vertikale Linie
            return [
                (o_h2_w, O_FE_RATIOS['Fe']),
                (2.0, O_FE_RATIOS['Fe']),
                (2.0, O_FE_RATIOS['Fe3O4']),
                (o_h2_w, O_FE_RATIOS['Fe3O4']),
            ]

    def calculate_h2_consumption(
        self,
        T_kelvin: float,
        input_oxide: str = 'Fe2O3',
        output_state: str = 'Fe',
    ) -> float:
        """
        Berechnet den theoretischen Wasserstoffverbrauch (mol H2 / mol Fe).

        Args:
            T_kelvin: Temperatur in Kelvin
            input_oxide: Eingangsmaterial ('Fe2O3', 'Fe3O4', 'FeO')
            output_state: Ausgangszustand ('Fe', 'FeO', 'Fe3O4')

        Returns:
            Wasserstoffverbrauch in mol H2 pro mol Fe
        """
        o_fe_in = O_FE_RATIOS.get(input_oxide, O_FE_RATIOS['Fe2O3'])
        o_fe_out = O_FE_RATIOS.get(output_state, O_FE_RATIOS['Fe'])

        # Sauerstoff der entfernt werden muss
        delta_o_fe = o_fe_in - o_fe_out

        # Bei idealer Reduktion (100% H2O am Ausgang) waere O/H2 = 2
        # Bei realistischer Reduktion am W-Punkt
        w_point = self.get_wustite_point(T_kelvin)
        o_h2_equilibrium = w_point[0]

        # H2 Verbrauch = delta_O / (O/H2)
        if o_h2_equilibrium > 0:
            # Am W-Punkt: H2/Fe = delta_O_Fe / (O/H2 - 1)
            # Weil nur der oxidierte Anteil des H2 (als H2O) Sauerstoff aufnimmt
            god_w = o_h2_to_god(o_h2_equilibrium)
            if god_w and god_w > 0:
                h2_per_fe = delta_o_fe / god_w
                return h2_per_fe

        return delta_o_fe  # Fallback: 1:1 Verhaeltnis


if __name__ == "__main__":
    # Test
    print("=== H2 Rist-Diagramm Test ===\n")

    thermo = BaurGlassnerThermo()
    rist = RistDiagramH2(thermo)

    for T_C in [700, 800, 900, 1000]:
        T_K = T_C + 273.15
        print(f"\nTemperatur: {T_C} C")

        # W- und M-Punkte
        w = rist.get_wustite_point(T_K)
        m = rist.get_magnetite_point(T_K)
        print(f"  W-Punkt (H2): O/H2 = {w[0]:.3f}, O/Fe = {w[1]:.3f}")
        print(f"  M-Punkt (H2): O/H2 = {m[0]:.3f}, O/Fe = {m[1]:.3f}")

        # GOD-Werte
        god_w = o_h2_to_god(w[0])
        god_m = o_h2_to_god(m[0])
        print(f"  GOD_W = {god_w:.3f}, GOD_M = {god_m:.3f}")

        # H2-Verbrauch
        h2_fe = rist.calculate_h2_consumption(T_K)
        print(f"  H2-Verbrauch: {h2_fe:.3f} mol H2 / mol Fe")
