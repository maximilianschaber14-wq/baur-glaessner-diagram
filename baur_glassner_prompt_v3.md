# Baur-Glässner-Diagramm – Interaktive Streamlit-App

## Projektübersicht
Erstelle eine interaktive Web-App mit Python und Streamlit, die das Baur-Glässner-Diagramm für die Systeme **Fe-C-O** und **Fe-H₂-O** parallel darstellt.

**Datenquelle:** `NBS_Tables_Library.xlsx` (NBS Thermochemical Tables)

---

## Datenquelle: NBS Thermochemical Tables (Excel)

Die Datei `NBS_Tables_Library.xlsx` enthält die offiziellen NBS (National Bureau of Standards) thermochemischen Daten bei 298.15 K.

### Excel-Struktur

- **Sheet:** `NBS Tables`
- **Header:** Zeile 2 (Index 1)
- **Spalten:** Formula, Solvent, Name, State Description, State, Molar Mass, ΔfH0° (0K), ΔfH°, ΔfG°, H°-H0°, S°/Cp

### Relevante Substanzen

| Formula | Name | State | ΔfH° [kJ/mol] | ΔfG° [kJ/mol] | S° [J/(mol·K)] |
|---------|------|-------|---------------|---------------|----------------|
| Fe | Eisen | cr | 0 | 0 | 27.28 |
| FeO0.947O | Wüstit | cr | -266.27 | -245.12 | 57.49 |
| Fe2O3 | Hämatit | cr | -824.2 | -742.2 | 87.40 |
| Fe3O4 | Magnetit | cr | -1118.4 | -1015.4 | 146.4 |
| CO | Kohlenmonoxid | g | -110.525 | -137.168 | 197.674 |
| CO2 | Kohlendioxid | g | -393.509 | -394.359 | 213.74 |
| H2 | Wasserstoff | g | 0 | 0 | 130.684 |
| H2O | Wasser | g | -241.818 | -228.572 | 188.825 |
| C | Graphit | cr | 0 | 0 | 5.740 |

---

## Thermodynamische Berechnungen

### Grundprinzip

Für eine Reaktion bei Temperatur T gilt:

```
ΔG°(T) = ΔH°(298) - T · ΔS°(298)
```

Dabei:
- **ΔH°(Reaktion)** = Σ ΔfH°(Produkte) - Σ ΔfH°(Reaktanten)
- **ΔS°(Reaktion)** = Σ S°(Produkte) - Σ S°(Reaktanten)

### Gleichgewichtskonstante und GOD

```
K = exp(-ΔG° / (R·T))
GOD = K^(1/n) / (1 + K^(1/n))
```

wobei n der stöchiometrische Koeffizient ist.

### Reaktionen

**CO-System:**
1. FeO + CO → Fe + CO₂ (T > 570°C)
2. Fe₃O₄ + CO → 3 FeO + CO₂ (T > 570°C)
3. Fe₃O₄ + 4 CO → 3 Fe + 4 CO₂ (T < 570°C)
4. 3 Fe₂O₃ + CO → 2 Fe₃O₄ + CO₂

**H₂-System:**
1. FeO + H₂ → Fe + H₂O (T > 570°C)
2. Fe₃O₄ + H₂ → 3 FeO + H₂O (T > 570°C)
3. Fe₃O₄ + 4 H₂ → 3 Fe + 4 H₂O (T < 570°C)
4. 3 Fe₂O₃ + H₂ → 2 Fe₃O₄ + H₂O

**Boudouard:**
- 2 CO ⇌ C + CO₂
- Empirisch: log₁₀(K) = -8916/T + 9.113

---

## Projektstruktur

```
baur-glassner-app/
├── app.py                      # Streamlit Hauptanwendung
├── thermodynamics.py           # Berechnungen
├── data_loader.py              # Excel-Parser
├── plotting.py                 # Diagramm
├── NBS_Tables_Library.xlsx     # Datenquelle
├── requirements.txt
└── README.md
```

---

## Implementierung

### data_loader.py

```python
import pandas as pd

def load_nbs_data(filepath: str = "NBS_Tables_Library.xlsx") -> dict:
    """
    Lädt thermodynamische Daten aus NBS Excel-Tabelle.
    
    Returns:
        dict mit Substanzen und ihren thermodynamischen Eigenschaften
    """
    # Excel einlesen (Header in Zeile 2)
    df = pd.read_excel(filepath, sheet_name='NBS Tables', header=1)
    
    # Spaltennamen vereinfachen
    df.columns = ['Formula', 'Solvent', 'Name', 'State_Desc', 'State', 
                  'Molar_Mass', 'dHf0_0K', 'dHf', 'dGf', 'H_H0', 'S_Cp', 
                  'Extra1', 'Extra2']
    
    substances = {}
    
    # Zu suchende Substanzen: (Formula in Excel, State, interner Key)
    targets = [
        ('Fe', 'cr', 'Fe_cr'),
        ('FeO0.947O', 'cr', 'FeO_cr'),    # Wüstit
        ('Fe2O3', 'cr', 'Fe2O3_cr'),       # Hämatit
        ('Fe3O4', 'cr', 'Fe3O4_cr'),       # Magnetit
        ('CO', 'g', 'CO_g'),
        ('CO2', 'g', 'CO2_g'),
        ('H2', 'g', 'H2_g'),
        ('H2O', 'g', 'H2O_g'),
        ('C', 'cr', 'C_gr'),               # Graphit
    ]
    
    for formula, state, key in targets:
        mask = (df['Formula'] == formula) & (df['State'] == state)
        rows = df[mask]
        
        if len(rows) > 0:
            row = rows.iloc[0]
            
            # Werte extrahieren (behandle '-' als None)
            dHf = parse_value(row['dHf'])
            dGf = parse_value(row['dGf'])
            S = parse_entropy(row['S_Cp'])
            
            substances[key] = {
                'dHf': dHf,   # kJ/mol
                'dGf': dGf,   # kJ/mol
                'S': S,       # J/(mol·K)
                'formula': formula,
                'state': state,
            }
    
    return substances


def parse_value(val) -> float:
    """Konvertiert Zellenwert zu float, None wenn '-' oder leer."""
    if pd.isna(val) or val == '-' or val == '':
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def parse_entropy(val) -> float:
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
```

### thermodynamics.py

```python
import numpy as np
from data_loader import load_nbs_data

R = 8.314  # J/(mol·K)
WUSTITE_STABILITY_K = 843  # 570°C - unterhalb ist FeO instabil


class BaurGlassnerThermo:
    """Thermodynamische Berechnungen für das Baur-Glässner-Diagramm."""
    
    def __init__(self, nbs_filepath: str = "NBS_Tables_Library.xlsx"):
        self.data = load_nbs_data(nbs_filepath)
        self.reactions = self._calculate_all_reactions()
    
    def _calculate_all_reactions(self) -> dict:
        """Berechnet ΔH° und ΔS° für alle Reaktionen aus NBS-Daten."""
        d = self.data
        reactions = {}
        
        # === CO-System ===
        
        # FeO + CO → Fe + CO2
        reactions['FeO_Fe_CO'] = {
            'equation': 'FeO + CO → Fe + CO₂',
            'dH': self._calc_dH([d['Fe_cr'], d['CO2_g']], 
                                [d['FeO_cr'], d['CO_g']]),
            'dS': self._calc_dS([d['Fe_cr'], d['CO2_g']], 
                                [d['FeO_cr'], d['CO_g']]),
            'stoich': 1,
            'valid_above_K': WUSTITE_STABILITY_K,
        }
        
        # Fe3O4 + CO → 3 FeO + CO2
        reactions['Fe3O4_FeO_CO'] = {
            'equation': 'Fe₃O₄ + CO → 3 FeO + CO₂',
            'dH': self._calc_dH([d['FeO_cr'], d['FeO_cr'], d['FeO_cr'], d['CO2_g']], 
                                [d['Fe3O4_cr'], d['CO_g']]),
            'dS': self._calc_dS([d['FeO_cr'], d['FeO_cr'], d['FeO_cr'], d['CO2_g']], 
                                [d['Fe3O4_cr'], d['CO_g']]),
            'stoich': 1,
            'valid_above_K': WUSTITE_STABILITY_K,
        }
        
        # Fe3O4 + 4 CO → 3 Fe + 4 CO2 (unterhalb 570°C)
        reactions['Fe3O4_Fe_CO'] = {
            'equation': 'Fe₃O₄ + 4 CO → 3 Fe + 4 CO₂',
            'dH': self._calc_dH([d['Fe_cr']]*3 + [d['CO2_g']]*4, 
                                [d['Fe3O4_cr']] + [d['CO_g']]*4),
            'dS': self._calc_dS([d['Fe_cr']]*3 + [d['CO2_g']]*4, 
                                [d['Fe3O4_cr']] + [d['CO_g']]*4),
            'stoich': 4,
            'valid_below_K': WUSTITE_STABILITY_K,
        }
        
        # 3 Fe2O3 + CO → 2 Fe3O4 + CO2
        reactions['Fe2O3_Fe3O4_CO'] = {
            'equation': '3 Fe₂O₃ + CO → 2 Fe₃O₄ + CO₂',
            'dH': self._calc_dH([d['Fe3O4_cr']]*2 + [d['CO2_g']], 
                                [d['Fe2O3_cr']]*3 + [d['CO_g']]),
            'dS': self._calc_dS([d['Fe3O4_cr']]*2 + [d['CO2_g']], 
                                [d['Fe2O3_cr']]*3 + [d['CO_g']]),
            'stoich': 1,
        }
        
        # === H2-System ===
        
        # FeO + H2 → Fe + H2O
        reactions['FeO_Fe_H2'] = {
            'equation': 'FeO + H₂ → Fe + H₂O',
            'dH': self._calc_dH([d['Fe_cr'], d['H2O_g']], 
                                [d['FeO_cr'], d['H2_g']]),
            'dS': self._calc_dS([d['Fe_cr'], d['H2O_g']], 
                                [d['FeO_cr'], d['H2_g']]),
            'stoich': 1,
            'valid_above_K': WUSTITE_STABILITY_K,
        }
        
        # Fe3O4 + H2 → 3 FeO + H2O
        reactions['Fe3O4_FeO_H2'] = {
            'equation': 'Fe₃O₄ + H₂ → 3 FeO + H₂O',
            'dH': self._calc_dH([d['FeO_cr']]*3 + [d['H2O_g']], 
                                [d['Fe3O4_cr'], d['H2_g']]),
            'dS': self._calc_dS([d['FeO_cr']]*3 + [d['H2O_g']], 
                                [d['Fe3O4_cr'], d['H2_g']]),
            'stoich': 1,
            'valid_above_K': WUSTITE_STABILITY_K,
        }
        
        # Fe3O4 + 4 H2 → 3 Fe + 4 H2O (unterhalb 570°C)
        reactions['Fe3O4_Fe_H2'] = {
            'equation': 'Fe₃O₄ + 4 H₂ → 3 Fe + 4 H₂O',
            'dH': self._calc_dH([d['Fe_cr']]*3 + [d['H2O_g']]*4, 
                                [d['Fe3O4_cr']] + [d['H2_g']]*4),
            'dS': self._calc_dS([d['Fe_cr']]*3 + [d['H2O_g']]*4, 
                                [d['Fe3O4_cr']] + [d['H2_g']]*4),
            'stoich': 4,
            'valid_below_K': WUSTITE_STABILITY_K,
        }
        
        # 3 Fe2O3 + H2 → 2 Fe3O4 + H2O
        reactions['Fe2O3_Fe3O4_H2'] = {
            'equation': '3 Fe₂O₃ + H₂ → 2 Fe₃O₄ + H₂O',
            'dH': self._calc_dH([d['Fe3O4_cr']]*2 + [d['H2O_g']], 
                                [d['Fe2O3_cr']]*3 + [d['H2_g']]),
            'dS': self._calc_dS([d['Fe3O4_cr']]*2 + [d['H2O_g']], 
                                [d['Fe2O3_cr']]*3 + [d['H2_g']]),
            'stoich': 1,
        }
        
        return reactions
    
    def _calc_dH(self, products: list, reactants: list) -> float:
        """ΔH° = Σ ΔfH°(Produkte) - Σ ΔfH°(Reaktanten) in kJ/mol"""
        prod_sum = sum(p['dHf'] for p in products)
        react_sum = sum(r['dHf'] for r in reactants)
        return prod_sum - react_sum
    
    def _calc_dS(self, products: list, reactants: list) -> float:
        """ΔS° = Σ S°(Produkte) - Σ S°(Reaktanten) in J/(mol·K)"""
        prod_sum = sum(p['S'] for p in products)
        react_sum = sum(r['S'] for r in reactants)
        return prod_sum - react_sum
    
    def delta_G(self, reaction_name: str, T_kelvin: float) -> float:
        """
        ΔG°(T) = ΔH° - T·ΔS° in J/mol
        (ΔH° wird von kJ zu J konvertiert)
        """
        rxn = self.reactions[reaction_name]
        return rxn['dH'] * 1000 - T_kelvin * rxn['dS']
    
    def equilibrium_K(self, reaction_name: str, T_kelvin: float) -> float:
        """Gleichgewichtskonstante K = exp(-ΔG°/RT)"""
        dG = self.delta_G(reaction_name, T_kelvin)
        return np.exp(-dG / (R * T_kelvin))
    
    def GOD(self, reaction_name: str, T_kelvin: float) -> float:
        """
        Gas Oxidation Degree aus Gleichgewichtskonstante.
        GOD = K^(1/n) / (1 + K^(1/n))
        """
        rxn = self.reactions[reaction_name]
        
        # Gültigkeitsbereich prüfen
        if 'valid_above_K' in rxn and T_kelvin < rxn['valid_above_K']:
            return np.nan
        if 'valid_below_K' in rxn and T_kelvin >= rxn['valid_below_K']:
            return np.nan
        
        K = self.equilibrium_K(reaction_name, T_kelvin)
        stoich = rxn['stoich']
        K_eff = K ** (1.0 / stoich)
        return np.clip(K_eff / (1.0 + K_eff), 0, 1)
    
    def boudouard_GOD(self, T_kelvin: float) -> float:
        """
        Boudouard-Gleichgewicht: 2 CO ⇌ C + CO₂
        Empirische Formel: log₁₀(K) = -8916/T + 9.113
        K = p_CO₂ / p_CO²
        """
        log_K = -8916.0 / T_kelvin + 9.113
        K = 10 ** log_K
        
        # Löse K = GOD / (1-GOD)² für GOD
        # K·(1-GOD)² = GOD
        # K - 2K·GOD + K·GOD² = GOD
        # K·GOD² - (2K+1)·GOD + K = 0
        a, b, c = K, -(2*K + 1), K
        discriminant = b**2 - 4*a*c
        
        if discriminant < 0:
            return np.nan
        
        GOD = (-b - np.sqrt(discriminant)) / (2*a)
        return np.clip(GOD, 0, 1)
    
    def calculate_phase_boundary(self, reaction_name: str, T_range_K: np.ndarray) -> np.ndarray:
        """Berechnet GOD-Werte für eine Phasengrenze."""
        return np.array([self.GOD(reaction_name, T) for T in T_range_K])
    
    def calculate_boudouard_boundary(self, T_range_K: np.ndarray) -> np.ndarray:
        """Berechnet GOD-Werte für Boudouard-Gleichgewicht."""
        return np.array([self.boudouard_GOD(T) for T in T_range_K])
```

### plotting.py

```python
import plotly.graph_objects as go
import numpy as np
from thermodynamics import BaurGlassnerThermo


def create_baur_glassner_diagram(
    thermo: BaurGlassnerThermo,
    show_CO_system: bool = True,
    show_H2_system: bool = True,
    show_boudouard: bool = True,
    show_labels: bool = True,
    T_min_C: float = 300,
    T_max_C: float = 1000,
) -> go.Figure:
    """Erstellt das interaktive Baur-Glässner-Diagramm."""
    
    fig = go.Figure()
    
    # Temperaturbereich
    T_range_C = np.linspace(T_min_C, T_max_C, 500)
    T_range_K = T_range_C + 273.15
    
    # === CO-System (durchgezogene Linien) ===
    if show_CO_system:
        co_lines = [
            ('FeO_Fe_CO', 'Fe/FeO'),
            ('Fe3O4_FeO_CO', 'Fe₃O₄/FeO'),
            ('Fe3O4_Fe_CO', 'Fe₃O₄/Fe'),
            ('Fe2O3_Fe3O4_CO', 'Fe₂O₃/Fe₃O₄'),
        ]
        
        for i, (rxn, label) in enumerate(co_lines):
            GOD = thermo.calculate_phase_boundary(rxn, T_range_K)
            fig.add_trace(go.Scatter(
                x=GOD, y=T_range_C,
                mode='lines',
                line=dict(color='black', width=2),
                name=f'{label} (CO)' if i == 0 else None,
                legendgroup='CO',
                showlegend=(i == 0),
                hovertemplate=f'{label}<br>GOD: %{{x:.3f}}<br>T: %{{y:.0f}}°C<extra></extra>'
            ))
    
    # === H2-System (gestrichelte Linien) ===
    if show_H2_system:
        h2_lines = [
            ('FeO_Fe_H2', 'Fe/FeO'),
            ('Fe3O4_FeO_H2', 'Fe₃O₄/FeO'),
            ('Fe3O4_Fe_H2', 'Fe₃O₄/Fe'),
            ('Fe2O3_Fe3O4_H2', 'Fe₂O₃/Fe₃O₄'),
        ]
        
        for i, (rxn, label) in enumerate(h2_lines):
            GOD = thermo.calculate_phase_boundary(rxn, T_range_K)
            fig.add_trace(go.Scatter(
                x=GOD, y=T_range_C,
                mode='lines',
                line=dict(color='black', width=2, dash='dash'),
                name=f'{label} (H₂)' if i == 0 else None,
                legendgroup='H2',
                showlegend=(i == 0),
                hovertemplate=f'{label}<br>GOD: %{{x:.3f}}<br>T: %{{y:.0f}}°C<extra></extra>'
            ))
    
    # === Boudouard-Linie (gepunktet) ===
    if show_boudouard:
        GOD_boud = thermo.calculate_boudouard_boundary(T_range_K)
        fig.add_trace(go.Scatter(
            x=GOD_boud, y=T_range_C,
            mode='lines',
            line=dict(color='black', width=1.5, dash='dot'),
            name='Boudouard',
            hovertemplate='Boudouard<br>GOD: %{x:.3f}<br>T: %{y:.0f}°C<extra></extra>'
        ))
    
    # === Phasenbeschriftungen ===
    if show_labels:
        labels = [
            (0.12, 750, 'Fe', 'black'),
            (0.50, 480, 'Fe₃O₄', 'darkblue'),
            (0.60, 820, 'Fe₍₁₋ₓ₎O', 'darkblue'),
            (0.92, 650, 'Fe₂O₃', 'darkred'),
        ]
        for x, y, text, color in labels:
            if T_min_C <= y <= T_max_C:
                fig.add_annotation(
                    x=x, y=y, text=text,
                    showarrow=False,
                    font=dict(size=16, color=color, family='Arial Black'),
                )
    
    # === Layout ===
    fig.update_layout(
        title=dict(
            text='Baur-Glässner Diagram',
            font=dict(size=20),
        ),
        xaxis=dict(
            title='Gas oxidation degree GOD [-]',
            range=[0, 1],
            dtick=0.1,
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
        ),
        yaxis=dict(
            title='Temperature [°C]',
            range=[T_min_C, T_max_C],
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
        ),
        legend=dict(
            x=0.02,
            y=0.98,
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor='gray',
            borderwidth=1,
        ),
        template='plotly_white',
        height=700,
        margin=dict(l=80, r=40, t=60, b=60),
    )
    
    return fig
```

### app.py

```python
import streamlit as st
import numpy as np
from thermodynamics import BaurGlassnerThermo
from plotting import create_baur_glassner_diagram

# Seitenkonfiguration
st.set_page_config(
    page_title="Baur-Glässner Diagram",
    page_icon="🔥",
    layout="wide"
)

# Titel
st.title("🔥 Baur-Glässner Diagram")
st.markdown("*Interactive phase equilibrium diagram for iron oxide reduction*")

# Thermodynamik laden
@st.cache_resource
def load_thermo():
    return BaurGlassnerThermo("NBS_Tables_Library.xlsx")

thermo = load_thermo()

# === Sidebar ===
with st.sidebar:
    st.header("⚙️ Display Options")
    
    show_CO = st.checkbox("Fe-C-O system (solid lines)", value=True)
    show_H2 = st.checkbox("Fe-H₂-O system (dashed lines)", value=True)
    show_boudouard = st.checkbox("Boudouard equilibrium (dotted)", value=True)
    show_labels = st.checkbox("Show phase labels", value=True)
    
    st.header("🌡️ Temperature Range")
    T_min = st.slider("Min Temperature [°C]", 200, 500, 300)
    T_max = st.slider("Max Temperature [°C]", 700, 1200, 1000)
    
    st.header("📊 Data Source")
    st.info("NBS Thermochemical Tables\n(National Bureau of Standards)")

# === Hauptdiagramm ===
fig = create_baur_glassner_diagram(
    thermo=thermo,
    show_CO_system=show_CO,
    show_H2_system=show_H2,
    show_boudouard=show_boudouard,
    show_labels=show_labels,
    T_min_C=T_min,
    T_max_C=T_max,
)

st.plotly_chart(fig, use_container_width=True)

# === Legende ===
col1, col2, col3 = st.columns(3)

with col1:
    st.markdown("""
    **Line Types:**
    - ━━━ Solid: CO/CO₂ system
    - ┅┅┅ Dashed: H₂/H₂O system
    - ····· Dotted: Boudouard equilibrium
    """)

with col2:
    st.markdown("""
    **X-Axis (GOD):**
    - CO system: CO₂/(CO+CO₂)
    - H₂ system: H₂O/(H₂+H₂O)
    """)

with col3:
    st.markdown("""
    **Key Temperature:**
    - 570°C: FeO stability limit
    - Below: FeO → Fe + Fe₃O₄
    """)

# === Info-Box ===
with st.expander("ℹ️ About this diagram"):
    st.markdown("""
    ### Baur-Glässner Diagram
    
    The Baur-Glässner diagram shows the thermodynamic stability regions of iron 
    and its oxides as a function of temperature and gas composition.
    
    **Phase Regions:**
    - **Fe**: Metallic iron (most reduced)
    - **FeO (Wüstite)**: Only stable above ~570°C
    - **Fe₃O₄ (Magnetite)**: Stable over wide temperature range
    - **Fe₂O₃ (Hematite)**: Most oxidized form
    
    **Boudouard Line:**
    Below this line, carbon can deposit from CO-rich gas mixtures.
    This is important for blast furnace operation.
    
    **Applications:**
    - Blast furnace process optimization
    - Direct reduction of iron ore
    - Gas-based steelmaking (H₂-DRI)
    
    **Data Source:**
    Thermodynamic data from NBS (National Bureau of Standards) 
    Thermochemical Tables at 298.15 K.
    """)

# === Thermodynamische Daten anzeigen ===
with st.expander("📋 Thermodynamic Data (NBS Tables)"):
    st.markdown("**Substances used for calculations:**")
    
    data_table = []
    for key, vals in thermo.data.items():
        data_table.append({
            'Key': key,
            'Formula': vals['formula'],
            'State': vals['state'],
            'ΔfH° [kJ/mol]': vals['dHf'],
            'ΔfG° [kJ/mol]': vals['dGf'],
            'S° [J/(mol·K)]': vals['S'],
        })
    
    st.dataframe(data_table, use_container_width=True)
```

### requirements.txt

```
streamlit>=1.28.0
plotly>=5.18.0
numpy>=1.24.0
pandas>=2.0.0
openpyxl>=3.1.0
```

### README.md

```markdown
# Baur-Glässner Diagram App

Interactive Streamlit application for visualizing the Baur-Glässner phase 
equilibrium diagram for iron oxide reduction.

## Features

- Fe-C-O system (CO/CO₂ equilibria)
- Fe-H₂-O system (H₂/H₂O equilibria)
- Boudouard equilibrium line
- Interactive temperature range selection
- Thermodynamic data from NBS Tables

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
streamlit run app.py
```

## Data Source

Thermodynamic data from:
- NBS (National Bureau of Standards) Thermochemical Tables
- File: `NBS_Tables_Library.xlsx`

## Thermodynamic Basis

Calculations use:
- ΔG°(T) = ΔH°(298) - T·ΔS°(298)
- K = exp(-ΔG°/RT)
- GOD = K^(1/n) / (1 + K^(1/n))

## References

- Baur, E., Glässner, A. (1903)
- NBS Thermochemical Tables
```

---

## Dateien für Claude Code

Kopiere in deinen Projektordner:

1. **`baur_glassner_prompt_v3.md`** (dieser Prompt)
2. **`NBS_Tables_Library.xlsx`** (deine Excel-Datei)

**Sage Claude Code:**

> Lies `baur_glassner_prompt_v3.md` und erstelle die Baur-Glässner App. 
> Die thermodynamischen Daten kommen aus `NBS_Tables_Library.xlsx`.
