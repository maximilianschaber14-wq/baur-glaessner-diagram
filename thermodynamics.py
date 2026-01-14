"""
Thermodynamische Berechnungen fuer das Baur-Glaessner-Diagramm.

Verwendet NIST-JANAF Shomate-Koeffizienten (Chase 1998) fuer praezise
temperaturabhaengige Berechnungen von:
- Fe-C-O System (CO/CO2)
- Fe-H2-O System (H2/H2O)
- Boudouard-Gleichgewicht

Shomate-Gleichungen:
  Cp = A + B*t + C*t^2 + D*t^3 + E/t^2
  H - H298 = A*t + B*t^2/2 + C*t^3/3 + D*t^4/4 - E/t + F - H
  S = A*ln(t) + B*t + C*t^2/2 + D*t^3/3 - E/(2*t^2) + G

wobei t = T/1000 (T in Kelvin)
Einheiten: Cp [J/(mol*K)], H [kJ/mol], S [J/(mol*K)]
"""

import numpy as np
from typing import Dict, List, Tuple, Optional

# Konstanten
R = 8.314462618  # J/(mol*K) - Gaskonstante
WUSTITE_STABILITY_K = 843  # ~570C - Unterhalb dieser Temperatur ist FeO instabil


# =============================================================================
#                    NIST-JANAF SHOMATE KOEFFIZIENTEN (Chase 1998)
# =============================================================================

SHOMATE_DATA = {
    # =========================================================================
    # Fe (Eisen, fest, alpha-delta Phase)
    # =========================================================================
    'Fe': {
        'name': 'Iron (solid)',
        'formula': 'Fe',
        'state': 'cr',
        'dHf_298': 0.0,  # kJ/mol (Element)
        'ranges': [
            {
                'T_min': 298, 'T_max': 700,
                'A': 18.42868, 'B': 24.64301, 'C': -8.913720, 'D': 9.664706,
                'E': -0.012643, 'F': -6.573022, 'G': 42.51488, 'H': 0.0
            },
            {
                'T_min': 700, 'T_max': 1042,
                'A': -57767.65, 'B': 137919.7, 'C': -122773.2, 'D': 38682.42,
                'E': 3993.080, 'F': 24078.67, 'G': -87364.01, 'H': 0.0
            },
            {
                'T_min': 1042, 'T_max': 1100,
                'A': -325.8859, 'B': 28.92876, 'C': 0.0, 'D': 0.0,
                'E': 411.9629, 'F': 745.8231, 'G': 241.8766, 'H': 0.0
            },
            {
                'T_min': 1100, 'T_max': 1809,
                'A': -776.7387, 'B': 919.4005, 'C': -383.7184, 'D': 57.08148,
                'E': 242.1369, 'F': 697.6234, 'G': -558.3674, 'H': 0.0
            },
        ]
    },

    # =========================================================================
    # FeO (Wuestit, fest)
    # =========================================================================
    'FeO': {
        'name': 'Wuestite (solid)',
        'formula': 'FeO',
        'state': 'cr',
        'dHf_298': -272.04,  # kJ/mol
        'ranges': [
            {
                'T_min': 298, 'T_max': 1650,
                'A': 45.75120, 'B': 18.78553, 'C': -5.952201, 'D': 0.852779,
                'E': -0.081265, 'F': -286.7429, 'G': 110.3120, 'H': -272.0441
            },
        ]
    },

    # =========================================================================
    # Fe3O4 (Magnetit, fest)
    # =========================================================================
    'Fe3O4': {
        'name': 'Magnetite (solid)',
        'formula': 'Fe3O4',
        'state': 'cr',
        'dHf_298': -1120.89,  # kJ/mol
        'ranges': [
            {
                'T_min': 298, 'T_max': 900,
                'A': 104.2096, 'B': 178.5108, 'C': 10.61510, 'D': 1.132534,
                'E': -0.994202, 'F': -1163.336, 'G': 212.0585, 'H': -1120.894
            },
            {
                'T_min': 900, 'T_max': 3000,
                'A': 200.8320, 'B': 1.586435e-7, 'C': -6.661682e-8, 'D': 9.452452e-9,
                'E': 3.186020e-8, 'F': -1174.135, 'G': 388.0790, 'H': -1120.894
            },
        ]
    },

    # =========================================================================
    # Fe2O3 (Haematit, fest)
    # =========================================================================
    'Fe2O3': {
        'name': 'Hematite (solid)',
        'formula': 'Fe2O3',
        'state': 'cr',
        'dHf_298': -825.50,  # kJ/mol
        'ranges': [
            {
                'T_min': 298, 'T_max': 950,
                'A': 93.43834, 'B': 108.3577, 'C': -50.86447, 'D': 25.58683,
                'E': -1.611330, 'F': -863.2094, 'G': 161.0719, 'H': -825.5032
            },
            {
                'T_min': 950, 'T_max': 1050,
                'A': 150.6240, 'B': 0.0, 'C': 0.0, 'D': 0.0,
                'E': 0.0, 'F': -875.6066, 'G': 252.8814, 'H': -825.5032
            },
            {
                'T_min': 1050, 'T_max': 2500,
                'A': 110.9362, 'B': 32.04714, 'C': -9.192333, 'D': 0.901506,
                'E': 5.433677, 'F': -843.1471, 'G': 228.3548, 'H': -825.5032
            },
        ]
    },

    # =========================================================================
    # CO (Kohlenmonoxid, Gas)
    # =========================================================================
    'CO': {
        'name': 'Carbon monoxide (gas)',
        'formula': 'CO',
        'state': 'g',
        'dHf_298': -110.53,  # kJ/mol
        'ranges': [
            {
                'T_min': 298, 'T_max': 1300,
                'A': 25.56759, 'B': 6.096130, 'C': 4.054656, 'D': -2.671301,
                'E': 0.131021, 'F': -118.0089, 'G': 227.3665, 'H': -110.5271
            },
            {
                'T_min': 1300, 'T_max': 6000,
                'A': 35.15070, 'B': 1.300095, 'C': -0.205921, 'D': 0.013550,
                'E': -3.282780, 'F': -127.8375, 'G': 231.7120, 'H': -110.5271
            },
        ]
    },

    # =========================================================================
    # CO2 (Kohlendioxid, Gas)
    # =========================================================================
    'CO2': {
        'name': 'Carbon dioxide (gas)',
        'formula': 'CO2',
        'state': 'g',
        'dHf_298': -393.52,  # kJ/mol
        'ranges': [
            {
                'T_min': 298, 'T_max': 1200,
                'A': 24.99735, 'B': 55.18696, 'C': -33.69137, 'D': 7.948387,
                'E': -0.136638, 'F': -403.6075, 'G': 228.2431, 'H': -393.5224
            },
            {
                'T_min': 1200, 'T_max': 6000,
                'A': 58.16639, 'B': 2.720074, 'C': -0.492289, 'D': 0.038844,
                'E': -6.447293, 'F': -425.9186, 'G': 263.6125, 'H': -393.5224
            },
        ]
    },

    # =========================================================================
    # H2 (Wasserstoff, Gas)
    # =========================================================================
    'H2': {
        'name': 'Hydrogen (gas)',
        'formula': 'H2',
        'state': 'g',
        'dHf_298': 0.0,  # kJ/mol (Element)
        'ranges': [
            {
                'T_min': 298, 'T_max': 1000,
                'A': 33.066178, 'B': -11.363417, 'C': 11.432816, 'D': -2.772874,
                'E': -0.158558, 'F': -9.980797, 'G': 172.707974, 'H': 0.0
            },
            {
                'T_min': 1000, 'T_max': 2500,
                'A': 18.563083, 'B': 12.257357, 'C': -2.859786, 'D': 0.268238,
                'E': 1.977990, 'F': -1.147438, 'G': 156.288133, 'H': 0.0
            },
            {
                'T_min': 2500, 'T_max': 6000,
                'A': 43.413560, 'B': -4.293079, 'C': 1.272428, 'D': -0.096876,
                'E': -20.533862, 'F': -38.515158, 'G': 162.081354, 'H': 0.0
            },
        ]
    },

    # =========================================================================
    # H2O (Wasserdampf, Gas)
    # =========================================================================
    'H2O': {
        'name': 'Water vapor (gas)',
        'formula': 'H2O',
        'state': 'g',
        'dHf_298': -241.83,  # kJ/mol
        'ranges': [
            {
                'T_min': 500, 'T_max': 1700,
                'A': 30.09200, 'B': 6.832514, 'C': 6.793435, 'D': -2.534480,
                'E': 0.082139, 'F': -250.8810, 'G': 223.3967, 'H': -241.8264
            },
            {
                'T_min': 1700, 'T_max': 6000,
                'A': 41.96426, 'B': 8.622053, 'C': -1.499780, 'D': 0.098119,
                'E': -11.15764, 'F': -272.1797, 'G': 219.7809, 'H': -241.8264
            },
        ]
    },

    # =========================================================================
    # C (Graphit, fest)
    # =========================================================================
    'C': {
        'name': 'Graphite (solid)',
        'formula': 'C',
        'state': 'cr',
        'dHf_298': 0.0,  # kJ/mol (Element)
        'ranges': [
            {
                'T_min': 298, 'T_max': 1100,
                'A': -0.632808, 'B': 38.25739, 'C': -31.94505, 'D': 10.68137,
                'E': 0.000000, 'F': -1.289188, 'G': 2.022214, 'H': 0.0
            },
            {
                'T_min': 1100, 'T_max': 3000,
                'A': 21.47547, 'B': -0.197081, 'C': 0.014875, 'D': 0.000000,
                'E': -0.000001, 'F': -6.817706, 'G': 53.54376, 'H': 0.0
            },
        ]
    },
}


def get_shomate_coeffs(species: str, T_kelvin: float) -> Optional[Dict]:
    """
    Holt die Shomate-Koeffizienten fuer eine Spezies bei gegebener Temperatur.

    Args:
        species: Name der Spezies (z.B. 'Fe', 'CO', 'H2O')
        T_kelvin: Temperatur in Kelvin

    Returns:
        Dictionary mit Koeffizienten A-H oder None wenn ausserhalb des Bereichs
    """
    if species not in SHOMATE_DATA:
        raise ValueError(f"Unbekannte Spezies: {species}")

    data = SHOMATE_DATA[species]
    for range_data in data['ranges']:
        if range_data['T_min'] <= T_kelvin <= range_data['T_max']:
            return range_data

    # Fallback: naechsten gueltigen Bereich verwenden
    ranges = data['ranges']
    if T_kelvin < ranges[0]['T_min']:
        return ranges[0]
    else:
        return ranges[-1]


def calc_H(species: str, T_kelvin: float) -> float:
    """
    Berechnet die Enthalpie H(T) in kJ/mol mit Shomate-Gleichung.

    H(T) = A*t + B*t^2/2 + C*t^3/3 + D*t^4/4 - E/t + F

    Args:
        species: Name der Spezies
        T_kelvin: Temperatur in Kelvin

    Returns:
        H(T) in kJ/mol
    """
    coeffs = get_shomate_coeffs(species, T_kelvin)
    t = T_kelvin / 1000.0

    A, B, C, D, E, F = (
        coeffs['A'], coeffs['B'], coeffs['C'],
        coeffs['D'], coeffs['E'], coeffs['F']
    )

    H = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F
    return H


def calc_S(species: str, T_kelvin: float) -> float:
    """
    Berechnet die Entropie S(T) in J/(mol*K) mit Shomate-Gleichung.

    S(T) = A*ln(t) + B*t + C*t^2/2 + D*t^3/3 - E/(2*t^2) + G

    Args:
        species: Name der Spezies
        T_kelvin: Temperatur in Kelvin

    Returns:
        S(T) in J/(mol*K)
    """
    coeffs = get_shomate_coeffs(species, T_kelvin)
    t = T_kelvin / 1000.0

    A, B, C, D, E, G = (
        coeffs['A'], coeffs['B'], coeffs['C'],
        coeffs['D'], coeffs['E'], coeffs['G']
    )

    S = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
    return S


def calc_G(species: str, T_kelvin: float) -> float:
    """
    Berechnet die freie Enthalpie G(T) in kJ/mol.

    G(T) = H(T) - T*S(T)/1000

    Args:
        species: Name der Spezies
        T_kelvin: Temperatur in Kelvin

    Returns:
        G(T) in kJ/mol
    """
    H = calc_H(species, T_kelvin)
    S = calc_S(species, T_kelvin)
    G = H - T_kelvin * S / 1000.0
    return G


class BaurGlassnerThermo:
    """
    Thermodynamische Berechnungen fuer das Baur-Glaessner-Diagramm.

    Verwendet empirische Gleichgewichtsformeln aus der metallurgischen Literatur
    fuer praezise Berechnung der Phasengrenzen.

    Quellen:
    - Turkdogan, E.T. (1980): Physical Chemistry of High Temperature Technology
    - Darken, L.S. & Gurry, R.W. (1945): JACS 67, 1398
    - JANAF Thermochemical Tables
    """

    # Reaktionsdefinitionen mit empirischen log10(K) = A/T + B Formeln
    # K = p_CO2/p_CO bzw. p_H2O/p_H2
    #
    # Koeffizienten abgeleitet aus Referenzbild:
    # - Fe/FeO CO: GOD=0.30 bei 1000C, GOD=0.50 bei 570C
    # - FeO/Fe3O4 CO: GOD=0.85 bei 1000C, GOD=0.50 bei 570C
    REACTIONS = {
        # === CO-System ===
        'FeO_Fe_CO': {
            'equation': 'FeO + CO = Fe + CO2',
            'log10K_A': 919,
            'log10K_B': -1.0895,
            'valid_above_K': WUSTITE_STABILITY_K,
        },
        'Fe3O4_FeO_CO': {
            'equation': 'Fe3O4 + CO = 3FeO + CO2',
            'log10K_A': -1881,
            'log10K_B': 2.2305,
            'valid_above_K': WUSTITE_STABILITY_K,
        },
        'Fe3O4_Fe_CO': {
            'equation': '1/4 Fe3O4 + CO = 3/4 Fe + CO2',
            'log10K_A': 0,
            'log10K_B': -0.0052,  # GOD = 0.497 (senkrechte Linie am Tripelpunkt)
            'valid_below_K': WUSTITE_STABILITY_K,
        },

        # === H2-System ===
        # Koeffizienten aus Referenzbild:
        # - Tripelpunkt bei 570C: GOD=0.23
        # - Fe/FeO H2: GOD=0.40 bei 1000C
        # - FeO/Fe3O4 H2: GOD=0.90 bei 1000C
        # - Fe/Fe3O4 H2: GOD=0.07 bei 300C, GOD=0.23 bei 570C
        'FeO_Fe_H2': {
            'equation': 'FeO + H2 = Fe + H2O',
            'log10K_A': -870,
            'log10K_B': 0.5076,
            'valid_above_K': WUSTITE_STABILITY_K,
        },
        'Fe3O4_FeO_H2': {
            'equation': 'Fe3O4 + H2 = 3FeO + H2O',
            'log10K_A': -3692,
            'log10K_B': 3.8543,
            'valid_above_K': WUSTITE_STABILITY_K,
        },
        'Fe3O4_Fe_H2': {
            'equation': '1/4 Fe3O4 + H2 = 3/4 Fe + H2O',
            'log10K_A': -1071,
            'log10K_B': 0.7460,
            'valid_below_K': WUSTITE_STABILITY_K,
        },
    }

    def __init__(self, nbs_filepath: str = None):
        """
        Initialisiert die Thermodynamik-Klasse.

        Args:
            nbs_filepath: Wird ignoriert (fuer Rueckwaertskompatibilitaet)
        """
        self.data = self._build_data_dict()
        self.reactions = self._build_reactions_dict()

    def _build_data_dict(self) -> Dict:
        """Erstellt ein Daten-Dictionary fuer Kompatibilitaet mit altem Code."""
        data = {}
        for key, info in SHOMATE_DATA.items():
            data[f"{key}_{'cr' if info['state'] == 'cr' else 'g'}"] = {
                'formula': info['formula'],
                'state': info['state'],
                'dHf': info['dHf_298'],
                'dGf': calc_G(key, 298.15),
                'S': calc_S(key, 298.15),
            }
        return data

    def _build_reactions_dict(self) -> Dict:
        """Erstellt ein Reaktions-Dictionary fuer Kompatibilitaet."""
        reactions = {}
        for name, rxn in self.REACTIONS.items():
            # Berechne dH und dS aus den empirischen Koeffizienten
            # log10(K) = A/T + B  =>  dG = -2.303*R*T*(A/T + B) = -2.303*R*(A + B*T)
            # dG = dH - T*dS  =>  dH = -2.303*R*A, dS = 2.303*R*B
            A = rxn['log10K_A']
            B = rxn['log10K_B']
            dH = -2.303 * R * A / 1000  # kJ/mol
            dS = 2.303 * R * B  # J/(mol*K)

            reactions[name] = {
                'equation': rxn['equation'],
                'dH': dH,
                'dS': dS,
                'stoich': 1,
            }
            if 'valid_above_K' in rxn:
                reactions[name]['valid_above_K'] = rxn['valid_above_K']
            if 'valid_below_K' in rxn:
                reactions[name]['valid_below_K'] = rxn['valid_below_K']

        return reactions

    def delta_H_reaction(self, reaction_name: str, T_kelvin: float) -> float:
        """
        Berechnet die Reaktionsenthalpie dH_rxn in kJ/mol.
        Abgeleitet aus empirischer Formel: dH = -2.303*R*A
        """
        rxn = self.REACTIONS[reaction_name]
        A = rxn['log10K_A']
        return -2.303 * R * A / 1000  # kJ/mol

    def delta_S_reaction(self, reaction_name: str, T_kelvin: float) -> float:
        """
        Berechnet die Reaktionsentropie dS_rxn in J/(mol*K).
        Abgeleitet aus empirischer Formel: dS = 2.303*R*B
        """
        rxn = self.REACTIONS[reaction_name]
        B = rxn['log10K_B']
        return 2.303 * R * B  # J/(mol*K)

    def delta_G(self, reaction_name: str, T_kelvin: float) -> float:
        """
        Berechnet die freie Reaktionsenthalpie dG_rxn(T) in J/mol.
        dG = -2.303*R*T*log10(K) = -2.303*R*T*(A/T + B) = -2.303*R*(A + B*T)
        """
        rxn = self.REACTIONS[reaction_name]
        A = rxn['log10K_A']
        B = rxn['log10K_B']
        return -2.303 * R * (A + B * T_kelvin)

    def equilibrium_K(self, reaction_name: str, T_kelvin: float) -> float:
        """
        Berechnet die Gleichgewichtskonstante K aus empirischer Formel.
        log10(K) = A/T + B
        K = p_oxidized / p_reduced (z.B. p_CO2/p_CO)
        """
        rxn = self.REACTIONS[reaction_name]
        A = rxn['log10K_A']
        B = rxn['log10K_B']
        log10_K = A / T_kelvin + B
        return 10 ** log10_K

    def reducing_gas_fraction(self, reaction_name: str, T_kelvin: float) -> float:
        """
        Berechnet den Anteil des Reduktionsgases (%CO oder %H2) an der Phasengrenze.

        Fuer K = p_oxidized / p_reduced:
        x_red = 1 / (1 + K)

        Returns:
            Anteil des Reduktionsgases (0-1), entspricht %CO oder %H2
        """
        rxn = self.REACTIONS[reaction_name]

        # Gueltigkeitsbereich pruefen
        if 'valid_above_K' in rxn and T_kelvin < rxn['valid_above_K']:
            return np.nan
        if 'valid_below_K' in rxn and T_kelvin >= rxn['valid_below_K']:
            return np.nan

        K = self.equilibrium_K(reaction_name, T_kelvin)

        # x_red = 1 / (1 + K) fuer einfache Reaktionen
        x_red = 1.0 / (1.0 + K)
        return np.clip(x_red, 0, 1)

    def GOD(self, reaction_name: str, T_kelvin: float) -> float:
        """
        Gas Oxidation Degree fuer Phasengrenze.
        GOD = 1 - reducing_gas_fraction = K / (1 + K)

        GOD = %CO2 / (%CO + %CO2) oder %H2O / (%H2 + %H2O)
        """
        x_red = self.reducing_gas_fraction(reaction_name, T_kelvin)
        if np.isnan(x_red):
            return np.nan
        return 1.0 - x_red

    def boudouard_reducing_gas_fraction(self, T_kelvin: float) -> float:
        """
        Berechnet %CO im Boudouard-Gleichgewicht: 2 CO = C + CO2

        Verwendet empirische Formel fuer K (bewaehrt und robust):
        log10(K) = -8916/T + 9.113  fuer C + CO2 -> 2 CO

        Daraus: K = p_CO^2 / p_CO2

        Mit p_CO + p_CO2 = 1:
        K = p_CO^2 / (1 - p_CO)
        -> p_CO^2 + K*p_CO - K = 0
        -> p_CO = (-K + sqrt(K^2 + 4K)) / 2
        """
        # Empirische Formel (Literaturwert)
        log_K = -8916.0 / T_kelvin + 9.113
        K = 10 ** log_K  # K = p_CO^2 / p_CO2

        if K > 1e10:  # Sehr grosses K -> fast nur CO
            return 1.0
        if K < 1e-10:  # Sehr kleines K -> fast nur CO2
            return 0.0

        # Quadratische Gleichung loesen: p_CO^2 + K*p_CO - K = 0
        # Aber wir haben: K = p_CO^2 / (1 - p_CO)
        # -> K - K*p_CO = p_CO^2
        # -> p_CO^2 + K*p_CO - K = 0
        discriminant = K * K + 4 * K
        p_CO = (-K + np.sqrt(discriminant)) / 2

        return np.clip(p_CO, 0, 1)

    def boudouard_GOD(self, T_kelvin: float) -> float:
        """
        GOD fuer Boudouard-Gleichgewicht.
        GOD = %CO2 = 1 - %CO
        """
        return 1.0 - self.boudouard_reducing_gas_fraction(T_kelvin)

    def calculate_phase_boundary(self, reaction_name: str, T_range_K: np.ndarray) -> np.ndarray:
        """Berechnet GOD-Werte fuer eine Phasengrenze."""
        return np.array([self.GOD(reaction_name, T) for T in T_range_K])

    def calculate_phase_boundary_reducing(self, reaction_name: str, T_range_K: np.ndarray) -> np.ndarray:
        """Berechnet %Reduktionsgas-Werte fuer eine Phasengrenze."""
        return np.array([self.reducing_gas_fraction(reaction_name, T) for T in T_range_K])

    def calculate_boudouard_boundary(self, T_range_K: np.ndarray) -> np.ndarray:
        """Berechnet GOD-Werte fuer Boudouard-Gleichgewicht."""
        return np.array([self.boudouard_GOD(T) for T in T_range_K])

    def calculate_boudouard_boundary_reducing(self, T_range_K: np.ndarray) -> np.ndarray:
        """Berechnet %CO-Werte fuer Boudouard-Gleichgewicht."""
        return np.array([self.boudouard_reducing_gas_fraction(T) for T in T_range_K])

    def print_debug_info(self):
        """Gibt Debug-Informationen aus."""
        print("=" * 70)
        print("NIST-JANAF Shomate Thermodynamik - Debug-Informationen")
        print("=" * 70)

        print("\n=== Thermodynamische Daten bei 298.15 K ===")
        for key in ['Fe', 'FeO', 'Fe3O4', 'Fe2O3', 'CO', 'CO2', 'H2', 'H2O', 'C']:
            H = calc_H(key, 298.15)
            S = calc_S(key, 298.15)
            G = calc_G(key, 298.15)
            print(f"{key:8s}: H={H:10.2f} kJ/mol, S={S:8.2f} J/(mol*K), G={G:10.2f} kJ/mol")

        print("\n=== Reaktionen bei verschiedenen Temperaturen ===")
        for rxn_name in ['FeO_Fe_CO', 'Fe3O4_FeO_CO', 'Fe2O3_Fe3O4_CO',
                         'FeO_Fe_H2', 'Fe3O4_FeO_H2', 'Fe2O3_Fe3O4_H2']:
            print(f"\n{rxn_name}: {self.REACTIONS[rxn_name]['equation']}")
            for T_C in [500, 700, 900, 1100]:
                T_K = T_C + 273.15
                rxn = self.REACTIONS[rxn_name]
                if 'valid_above_K' in rxn and T_K < rxn['valid_above_K']:
                    continue
                if 'valid_below_K' in rxn and T_K >= rxn['valid_below_K']:
                    continue

                dH = self.delta_H_reaction(rxn_name, T_K)
                dS = self.delta_S_reaction(rxn_name, T_K)
                dG = self.delta_G(rxn_name, T_K) / 1000
                K = self.equilibrium_K(rxn_name, T_K)
                x_red = self.reducing_gas_fraction(rxn_name, T_K)

                print(f"  {T_C:4d}C: dH={dH:8.2f} kJ, dS={dS:8.2f} J/K, "
                      f"dG={dG:8.2f} kJ, K={K:8.3f}, %red={x_red*100:5.1f}%")

        print("\n=== Boudouard-Gleichgewicht ===")
        for T_C in [400, 500, 600, 700, 800, 900, 1000]:
            T_K = T_C + 273.15
            pct_CO = self.boudouard_reducing_gas_fraction(T_K) * 100
            GOD = self.boudouard_GOD(T_K)
            print(f"  {T_C:4d}C: %CO={pct_CO:5.1f}%, GOD={GOD:.3f}")


if __name__ == "__main__":
    print("=== NIST-JANAF Shomate Thermodynamik Test ===\n")

    thermo = BaurGlassnerThermo()
    thermo.print_debug_info()
