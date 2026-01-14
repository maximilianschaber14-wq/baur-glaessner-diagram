"""
Datenloader für NBS Thermochemical Tables.

Lädt thermodynamische Daten aus der Excel-Datei NBS_Tables_Library.xlsx.
"""

import pandas as pd
from typing import Optional


def load_nbs_data(filepath: str = "NBS_Tables_Library.xlsx") -> dict:
    """
    Lädt thermodynamische Daten aus NBS Excel-Tabelle.

    Args:
        filepath: Pfad zur Excel-Datei

    Returns:
        dict mit Substanzen und ihren thermodynamischen Eigenschaften:
        {
            'Fe_cr': {'dHf': 0, 'dGf': 0, 'S': 27.28, 'formula': 'Fe', 'state': 'cr'},
            ...
        }
    """
    # Excel einlesen (Header in Zeile 2, also index 1)
    df = pd.read_excel(filepath, sheet_name='NBS Tables', header=1)

    substances = {}

    # Zu suchende Substanzen: (Formula in Excel, State, interner Key)
    targets = [
        ('Fe', 'cr', 'Fe_cr'),
        ('FeO0.947O', 'cr', 'FeO_cr'),      # Wüstit
        ('Fe2O3', 'cr', 'Fe2O3_cr'),         # Hämatit
        ('Fe3O4', 'cr', 'Fe3O4_cr'),         # Magnetit
        ('CO', 'g', 'CO_g'),
        ('CO2', 'g', 'CO2_g'),
        ('H2', 'g', 'H2_g'),
        ('H2O', 'g', 'H2O_g'),
        ('C', 'cr', 'C_gr'),                 # Graphit
    ]

    for formula, state, key in targets:
        if 'FeO0.947' in formula:
            # Spezielle Suche für Wüstit
            mask = df['Formula'].str.contains('FeO0.947', na=False) & (df['State'] == state)
        else:
            mask = (df['Formula'] == formula) & (df['State'] == state)

        rows = df[mask]

        if len(rows) > 0:
            row = rows.iloc[0]

            # Werte extrahieren aus den korrekten Spalten-Indizes
            # Index 7: ΔfH°, Index 8: ΔfG°, Index 10: S°/Cp
            dHf = parse_value(row.iloc[7])
            dGf = parse_value(row.iloc[8])
            S = parse_entropy(row.iloc[10])

            substances[key] = {
                'dHf': dHf,     # kJ/mol
                'dGf': dGf,     # kJ/mol
                'S': S,         # J/(mol·K)
                'formula': row['Formula'],
                'state': state,
            }

    return substances


def parse_value(val) -> Optional[float]:
    """Konvertiert Zellenwert zu float, None wenn '-' oder leer."""
    if pd.isna(val) or val == '-' or val == '':
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def parse_entropy(val) -> Optional[float]:
    """
    Extrahiert S° aus der S°/Cp Spalte.
    Format kann sein: "27.28" oder "27.28    25.10" (S und Cp zusammen)
    """
    if pd.isna(val) or val == '-' or val == '':
        return None
    try:
        # Nimm ersten Wert (S°), falls mehrere vorhanden
        s = str(val).split()[0]
        return float(s)
    except (ValueError, TypeError, IndexError):
        return None


if __name__ == '__main__':
    # Test: Daten laden und anzeigen
    data = load_nbs_data()
    print("Geladene NBS-Daten:")
    print("-" * 60)
    for key, vals in data.items():
        print(f"{key:12} | dHf={vals['dHf']:>10} | dGf={vals['dGf']:>10} | S={vals['S']:>8}")
