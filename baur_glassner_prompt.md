# Baur-Glässner-Diagramm – Interaktive Streamlit-App

## Projektübersicht
Entwickle eine interaktive Web-App mit Python und Streamlit, die das Baur-Glässner-Diagramm für die Systeme **Fe-C-O** und **Fe-H₂-O** parallel darstellt.

---

## Diagramm-Spezifikation

### Achsen
- **X-Achse**: Gas Oxidation Degree (GOD) von 0.0 bis 1.0
  - Für CO-System: GOD = CO₂/(CO+CO₂)
  - Für H₂-System: GOD = H₂O/(H₂+H₂O)
  - Links (GOD=0): reduzierende Atmosphäre
  - Rechts (GOD=1): oxidierende Atmosphäre
- **Y-Achse**: Temperatur von 300°C bis 1000°C

### Phasenfelder
- **Fe** (metallisches Eisen): links
- **Fe₃O₄** (Magnetit): Mitte-rechts, unterhalb ~570°C
- **Fe₍₁₋ₓ₎O** (Wüstit): großes Feld, oberhalb ~570°C
- **Fe₂O₃** (Hämatit): ganz rechts (schmaler Streifen)

### Linientypen
1. **System Fe-C-O**: Durchgezogene Linien (solid)
2. **System Fe-H₂-O**: Gestrichelte Linien (dashed)
3. **Boudouard-Linie**: Gepunktete Linie (dotted)

---

## Thermodynamische Daten

### Grundformel
Für jede Reaktion: ΔG° = A + B·T [J/mol]

Gleichgewichtskonstante: K = exp(-ΔG°/(R·T))

Umrechnung zu GOD: 
- Für K = p_ox/p_red: GOD = K/(1+K)
- Für K = (p_ox/p_red)^n: GOD = K^(1/n)/(1 + K^(1/n))

### Thermodynamische Koeffizienten (ΔG° = A + B·T in J/mol)

```json
{
  "reactions": {
    "CO_system": {
      "Fe2O3_Fe3O4": {
        "equation": "3 Fe2O3 + CO = 2 Fe3O4 + CO2",
        "A": -52131,
        "B": -41.0,
        "stoichiometry": 1,
        "notes": "Hämatit zu Magnetit, über gesamten T-Bereich"
      },
      "Fe3O4_FeO": {
        "equation": "Fe3O4 + CO = 3 FeO + CO2",
        "A": 35380,
        "B": -40.16,
        "stoichiometry": 1,
        "valid_above_K": 843,
        "notes": "Magnetit zu Wüstit, nur oberhalb 570°C (843 K)"
      },
      "FeO_Fe": {
        "equation": "FeO + CO = Fe + CO2",
        "A": -22800,
        "B": 24.26,
        "stoichiometry": 1,
        "valid_above_K": 843,
        "notes": "Wüstit zu Eisen, nur oberhalb 570°C (843 K)"
      },
      "Fe3O4_Fe_direct": {
        "equation": "Fe3O4 + 4 CO = 3 Fe + 4 CO2",
        "A": -33020,
        "B": 32.62,
        "stoichiometry": 4,
        "valid_below_K": 843,
        "notes": "Direktreduktion unterhalb 570°C"
      }
    },
    "H2_system": {
      "Fe2O3_Fe3O4": {
        "equation": "3 Fe2O3 + H2 = 2 Fe3O4 + H2O",
        "A": -8700,
        "B": -58.7,
        "stoichiometry": 1,
        "notes": "Hämatit zu Magnetit"
      },
      "Fe3O4_FeO": {
        "equation": "Fe3O4 + H2 = 3 FeO + H2O",
        "A": 78100,
        "B": -57.6,
        "stoichiometry": 1,
        "valid_above_K": 843,
        "notes": "Magnetit zu Wüstit"
      },
      "FeO_Fe": {
        "equation": "FeO + H2 = Fe + H2O",
        "A": 23600,
        "B": 7.0,
        "stoichiometry": 1,
        "valid_above_K": 843,
        "notes": "Wüstit zu Eisen"
      },
      "Fe3O4_Fe_direct": {
        "equation": "Fe3O4 + 4 H2 = 3 Fe + 4 H2O",
        "A": 149000,
        "B": -115.0,
        "stoichiometry": 4,
        "valid_below_K": 843,
        "notes": "Direktreduktion unterhalb 570°C"
      }
    },
    "boudouard": {
      "equation": "C + CO2 = 2 CO (oder umgekehrt: 2 CO = C + CO2)",
      "log10_Keq_formula": "log10(K) = -8916/T + 9.113",
      "notes": "Gültig 500-2200 K. K=1 bei ca. 975 K (702°C). Unterhalb: Kohlenstoffabscheidung möglich."
    }
  },
  "constants": {
    "R": 8.314,
    "wustite_stability_K": 843,
    "wustite_stability_C": 570
  },
  "data_sources": [
    "NIST-JANAF Thermochemical Tables",
    "HSC Chemistry Database",
    "Kubaschewski: Metallurgical Thermochemistry",
    "Wikipedia: Boudouard reaction"
  ]
}
```

### Alternative Formulierung für Boudouard
Für die Boudouard-Reaktion (2 CO ⇌ C + CO₂):
```
log₁₀(K) = -8916/T + 9.113  (T in Kelvin)
```
Daraus GOD berechnen:
```python
K_boudouard = 10**(-8916/T + 9.113)
# Für 2CO = C + CO2: K = p_CO2 / p_CO^2 bei Gesamtdruck 1 bar
# Wenn p_CO + p_CO2 = 1:
# K = p_CO2 / (1-p_CO2)^2
# Diese quadratische Gleichung lösen für p_CO2 (=GOD)
```

---

## Implementierung

### Projektstruktur
```
baur-glassner-app/
├── app.py                 # Streamlit Hauptanwendung
├── thermodynamics.py      # Gleichgewichtsberechnungen
├── plotting.py            # Plotly-Diagramm
├── data/
│   └── thermo_data.json   # Thermodynamische Konstanten
├── requirements.txt
└── README.md
```

### thermodynamics.py - Kernfunktionen

```python
import numpy as np

R = 8.314  # J/(mol·K)
WUSTITE_STABILITY_K = 843  # ~570°C

def delta_G(A: float, B: float, T_kelvin: float) -> float:
    """Berechnet ΔG° = A + B·T [J/mol]"""
    return A + B * T_kelvin

def equilibrium_constant(A: float, B: float, T_kelvin: float) -> float:
    """Berechnet K aus ΔG° = -RT·ln(K)"""
    dG = delta_G(A, B, T_kelvin)
    return np.exp(-dG / (R * T_kelvin))

def GOD_from_K(K: float, stoichiometry: int = 1) -> float:
    """
    Berechnet GOD aus K.
    Für K = (p_ox/p_red)^n: 
    GOD = K^(1/n) / (1 + K^(1/n))
    """
    K_eff = K ** (1.0 / stoichiometry)
    return K_eff / (1.0 + K_eff)

def boudouard_GOD(T_kelvin: float) -> float:
    """
    Berechnet GOD für Boudouard-Gleichgewicht.
    2 CO ⇌ C + CO2, K = p_CO2 / p_CO^2
    Mit p_CO + p_CO2 = 1 (Gesamtdruck normiert)
    """
    log_K = -8916.0 / T_kelvin + 9.113
    K = 10 ** log_K
    # Löse: K = GOD / (1-GOD)^2
    # K(1-GOD)^2 = GOD
    # K - 2K·GOD + K·GOD^2 = GOD
    # K·GOD^2 - (2K+1)·GOD + K = 0
    a = K
    b = -(2*K + 1)
    c = K
    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        return np.nan
    GOD = (-b - np.sqrt(discriminant)) / (2*a)
    return np.clip(GOD, 0, 1)

def calculate_phase_boundary(system: str, reaction: str, T_range: np.ndarray) -> np.ndarray:
    """
    Berechnet GOD-Werte für eine Phasengrenze.
    system: 'CO' oder 'H2'
    reaction: 'Fe2O3_Fe3O4', 'Fe3O4_FeO', 'FeO_Fe', 'Fe3O4_Fe_direct'
    """
    # Thermodynamische Daten (siehe JSON oben)
    data = get_thermo_data(system, reaction)
    
    GOD_values = []
    for T in T_range:
        # Prüfe Gültigkeitsbereich
        if 'valid_above_K' in data and T < data['valid_above_K']:
            GOD_values.append(np.nan)
            continue
        if 'valid_below_K' in data and T >= data['valid_below_K']:
            GOD_values.append(np.nan)
            continue
            
        K = equilibrium_constant(data['A'], data['B'], T)
        GOD = GOD_from_K(K, data['stoichiometry'])
        GOD_values.append(GOD)
    
    return np.array(GOD_values)
```

### plotting.py - Diagrammerstellung

```python
import plotly.graph_objects as go
import numpy as np

def create_baur_glassner_diagram(
    show_CO_system: bool = True,
    show_H2_system: bool = True,
    show_boudouard: bool = True,
    T_min_C: float = 300,
    T_max_C: float = 1000
) -> go.Figure:
    """Erstellt das interaktive Baur-Glässner-Diagramm."""
    
    fig = go.Figure()
    
    T_range_C = np.linspace(T_min_C, T_max_C, 500)
    T_range_K = T_range_C + 273.15
    
    # CO-System (durchgezogene Linien)
    if show_CO_system:
        # Fe/FeO Grenze (T > 570°C)
        # Fe/Fe3O4 Grenze (T < 570°C)
        # Fe3O4/FeO Grenze
        # FeO/Fe2O3 Grenze
        pass  # Implementiere mit calculate_phase_boundary()
    
    # H2-System (gestrichelte Linien)
    if show_H2_system:
        pass  # Analog zu CO-System
    
    # Boudouard-Linie (gepunktet)
    if show_boudouard:
        GOD_boudouard = [boudouard_GOD(T) for T in T_range_K]
        fig.add_trace(go.Scatter(
            x=GOD_boudouard, y=T_range_C,
            mode='lines',
            line=dict(dash='dot', color='black', width=1.5),
            name='Boudouard (C + CO₂ ⇌ 2CO)'
        ))
    
    # Phasenfelder beschriften
    # Legende und Layout
    fig.update_layout(
        xaxis_title="Gas oxidation degree GOD [-]",
        yaxis_title="Temperature [°C]",
        xaxis=dict(range=[0, 1]),
        yaxis=dict(range=[T_min_C, T_max_C]),
        legend=dict(x=0.6, y=0.15),
        template="plotly_white"
    )
    
    return fig
```

### app.py - Streamlit App

```python
import streamlit as st
from plotting import create_baur_glassner_diagram

st.set_page_config(page_title="Baur-Glässner Diagram", layout="wide")
st.title("Baur-Glässner Diagram")
st.markdown("Interactive phase diagram for Fe-O-C and Fe-O-H systems")

# Sidebar
with st.sidebar:
    st.header("Options")
    show_CO = st.checkbox("Show Fe-C-O system (solid lines)", value=True)
    show_H2 = st.checkbox("Show Fe-H₂-O system (dashed lines)", value=True)
    show_boudouard = st.checkbox("Show Boudouard equilibrium (dotted)", value=True)
    
    st.header("Temperature Range")
    T_min = st.slider("Min Temperature [°C]", 200, 500, 300)
    T_max = st.slider("Max Temperature [°C]", 800, 1200, 1000)

# Hauptdiagramm
fig = create_baur_glassner_diagram(
    show_CO_system=show_CO,
    show_H2_system=show_H2,
    show_boudouard=show_boudouard,
    T_min_C=T_min,
    T_max_C=T_max
)

st.plotly_chart(fig, use_container_width=True)

# Erklärung
with st.expander("About this diagram"):
    st.markdown("""
    ### Baur-Glässner Diagram
    
    Shows the stability regions of iron and its oxides as a function of:
    - **Temperature** (y-axis)
    - **Gas oxidation degree** (x-axis): CO₂/(CO+CO₂) or H₂O/(H₂+H₂O)
    
    **Phase regions:**
    - **Fe**: Metallic iron
    - **FeO (Wüstite)**: Only stable above ~570°C
    - **Fe₃O₄ (Magnetite)**: Stable over wide range
    - **Fe₂O₃ (Hematite)**: Most oxidized form
    
    **Boudouard line**: Below this line, carbon deposition can occur.
    """)
```

### requirements.txt
```
streamlit>=1.28.0
plotly>=5.18.0
numpy>=1.24.0
```

---

## Validierung

Vergleiche die berechneten Kurven mit dem Referenzbild:
1. **570°C-Punkt**: Tripelpunkt Fe/FeO/Fe₃O₄ muss korrekt sein
2. **Kurvenverläufe**: H₂-System liegt bei hohen Temperaturen links vom CO-System (H₂ ist dort reduktiver)
3. **Boudouard-Linie**: Schneidet Y-Achse bei ca. 700°C, verläuft nach rechts unten
4. **GOD-Werte**: Bei 1000°C sollte Fe/FeO-Grenze für CO bei ca. GOD=0.25-0.30 liegen

---

## Diagramm-Styling (entsprechend Referenzbild)

- **Fe-C-O System**: Schwarze durchgezogene Linien
- **Fe-H₂-O System**: Schwarze gestrichelte Linien  
- **Boudouard**: Schwarze gepunktete Linie
- **Achsen**: Deutsche Beschriftung optional
- **Legende**: Im Plot unten rechts
- **Reaktionsgleichungen**: Als Annotationen an den Linien (optional)

---

## Hinweise zur Implementierung

1. **Starte mit den Berechnungen**: Implementiere zuerst `thermodynamics.py` und validiere die Gleichgewichtskonstanten gegen bekannte Werte
2. **Teste einzelne Kurven**: Plotte erst eine Kurve (z.B. FeO/Fe für CO), bevor alle kombiniert werden
3. **Achte auf den 570°C-Punkt**: Hier treffen sich mehrere Linien - das ist ein guter Validierungspunkt
4. **Boudouard ist kritisch**: Die Formel log₁₀(K) = -8916/T + 9.113 ist validiert (Wikipedia-Quelle)

Viel Erfolg!
