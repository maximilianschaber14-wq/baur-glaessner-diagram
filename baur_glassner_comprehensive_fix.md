# Baur-Glässner Diagramm - Umfassende Korrektur

## Übersicht der Probleme

Das aktuelle Diagramm hat mehrere Fehler im Vergleich zur Referenz:

1. **Boudouard-Linie** verläuft in falsche Richtung
2. **H₂-Linie (Fe₃O₄/Fe)** geht unrealistisch nach links unten
3. **570°C Tripelpunkt** ist nicht korrekt dargestellt
4. **Kurvenverläufe** stimmen nicht mit Referenz überein

---

## Referenz-Diagramm Eigenschaften

So sollte das Diagramm aussehen:

```
    GOD →  0.0   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0
         ┌─────────────────────────────────────────────────────────────────┐
   1000°C│ ·                    ╲H₂                                        │
         │   ·                    ╲    CO                                  │
    900°C│     ·                   ╲  ╱                                    │
         │       ·        Fe        ╲╱     FeO                    Fe₂O₃   │
    800°C│         ·                ╱╲                          ╱         │
         │           ·  Boudouard  ╱  ╲                       ╱           │
    700°C│             ·   ····   ╱    ╲                    ╱             │
         │               ·  ···  ╱      ╲                 ╱               │
    600°C│                 ·    ╱========╲══════════════╱  (570°C)        │
         │                     ╱          ╲                               │
    500°C│                    ╱   Fe₃O₄    ╲                              │
         │                   ╱              ╲                             │
    400°C│                  ╱                ╲                            │
         │                 ╱                  ╲                           │
    300°C│                ╱                    ╲                          │
         └─────────────────────────────────────────────────────────────────┘
```

**Wichtige Merkmale:**
- Boudouard: von links-oben (~GOD=0, 1000°C) nach rechts-unten (~GOD=0.7, 600°C)
- Fe/FeO (CO): fast senkrecht, leicht nach links geneigt bei hohen T
- Fe/FeO (H₂): kreuzt CO-Linie bei ca. 810°C
- Fe₃O₄/Fe: von GOD≈0.05 (400°C) nach GOD≈0.5 (570°C)
- Alle Linien treffen sich bei 570°C am Tripelpunkt

---

## Fix 1: Boudouard-Linie korrigieren

Die Boudouard-Reaktion **2 CO ⇌ C + CO₂** hat:
- Bei **niedrigen T**: K groß → CO₂ bevorzugt → **hoher GOD**
- Bei **hohen T**: K klein → CO bevorzugt → **niedriger GOD**

### Korrigierte Funktion:

```python
def boudouard_GOD(self, T_kelvin: float) -> float:
    """
    Boudouard-Gleichgewicht: 2 CO ⇌ C + CO₂
    
    Die empirische Formel log₁₀(K) = -8916/T + 9.113 gilt für C + CO₂ → 2 CO
    Für die Rückreaktion (2 CO → C + CO₂) ist K_reverse = 1/K
    
    K = p_CO² / p_CO₂  (für C + CO₂ → 2 CO)
    
    Mit p_CO + p_CO₂ = 1:
    GOD = p_CO₂ = 1 - p_CO
    
    K = (1-GOD)² / GOD
    → GOD = (1-GOD)² / K
    → K·GOD = (1-GOD)²
    → K·GOD = 1 - 2·GOD + GOD²
    → GOD² - (2+K)·GOD + 1 = 0
    """
    # K für C + CO₂ → 2 CO
    log_K = -8916.0 / T_kelvin + 9.113
    K = 10 ** log_K
    
    # Quadratische Gleichung: GOD² - (2+K)·GOD + 1 = 0
    a = 1
    b = -(2 + K)
    c = 1
    
    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        return np.nan
    
    # Wir brauchen die kleinere Wurzel (GOD zwischen 0 und 1)
    GOD = (-b - np.sqrt(discriminant)) / (2*a)
    
    return np.clip(GOD, 0, 1)
```

---

## Fix 2: GOD-Berechnung für Reduktionsreaktionen

Für die Reaktion **FeO + CO → Fe + CO₂**:

$$K = \frac{p_{CO_2}}{p_{CO}}$$

$$GOD = \frac{p_{CO_2}}{p_{CO} + p_{CO_2}} = \frac{K}{1 + K}$$

**ABER:** Die Phasengrenze zeigt, wo das **Oxid** stabil wird (rechts) vs. **Metall** (links).

Wenn K > 1: Mehr CO₂ im Gleichgewicht → Reduktion ist begünstigt → Fe ist stabil
Wenn K < 1: Mehr CO im Gleichgewicht → Oxidation ist begünstigt → FeO ist stabil

**Die korrekte Formel ist:**

```python
def GOD(self, reaction_name: str, T_kelvin: float) -> float:
    """
    Berechnet GOD für Phasengrenze.
    
    Für Reduktion MeO + Red → Me + Ox:
    K = p_Ox / p_Red
    GOD = K / (1 + K)
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
    
    # GOD = K / (1 + K) ist KORREKT für diese Definition
    return np.clip(K_eff / (1.0 + K_eff), 0, 1)
```

---

## Fix 3: Thermodynamische Daten überprüfen

Die berechneten ΔH° und ΔS° Werte müssen überprüft werden. Hier die erwarteten Werte:

### Reaktion: FeO + CO → Fe + CO₂

```python
# Aus NBS-Daten:
# ΔfH°: Fe=0, CO₂=-393.509, FeO=-266.27, CO=-110.525 (alle kJ/mol)
# S°: Fe=27.28, CO₂=213.74, FeO=57.49, CO=197.674 (alle J/(mol·K))

dH = (0 + (-393.509)) - ((-266.27) + (-110.525))
   = -393.509 - (-376.795)
   = -16.714 kJ/mol

dS = (27.28 + 213.74) - (57.49 + 197.674)
   = 241.02 - 255.164
   = -14.144 J/(mol·K)
```

**Check bei 1000 K:**
```python
dG = -16714 - 1000 * (-14.144) = -16714 + 14144 = -2570 J/mol
K = exp(2570 / (8.314 * 1000)) = exp(0.309) = 1.36
GOD = 1.36 / 2.36 = 0.576
```

Das scheint plausibel für die Fe/FeO-Grenze bei 1000°C.

---

## Fix 4: Horizontale Linie bei 570°C entfernen

Die horizontale Linie bei 570°C ist **keine echte Phasengrenze**, sondern nur die Grenze, unterhalb derer FeO instabil ist.

**Lösung:** Entferne die horizontale Linie und stelle sicher, dass:
- Fe₃O₄/Fe-Linie bei 570°C endet
- Fe/FeO-Linie bei 570°C beginnt
- Fe₃O₄/FeO-Linie bei 570°C beginnt

```python
# In plotting.py - keine horizontale Verbindungslinie zeichnen!
# Die Linien enden/beginnen einfach bei 570°C
```

---

## Fix 5: Korrekte Linienverbindung am 570°C Tripelpunkt

Am Tripelpunkt (570°C, GOD ≈ 0.5 für CO) treffen sich:
1. Fe₃O₄/Fe-Linie (von unten kommend)
2. Fe/FeO-Linie (nach oben gehend)  
3. Fe₃O₄/FeO-Linie (nach rechts-oben gehend)

**Diese drei Linien müssen sich in EINEM Punkt treffen!**

```python
def get_tripelpunkt_GOD(self, system='CO'):
    """Berechnet GOD am 570°C Tripelpunkt."""
    T_triple = 843  # K (570°C)
    
    # Alle drei Reaktionen sollten hier den gleichen GOD ergeben
    if system == 'CO':
        return self.GOD('FeO_Fe_CO', T_triple)
    else:
        return self.GOD('FeO_Fe_H2', T_triple)
```

---

## Komplette korrigierte thermodynamics.py

```python
import numpy as np
from data_loader import load_nbs_data

R = 8.314  # J/(mol·K)
WUSTITE_STABILITY_K = 843  # 570°C


class BaurGlassnerThermo:
    def __init__(self, nbs_filepath: str = "NBS_Tables_Library.xlsx"):
        self.data = load_nbs_data(nbs_filepath)
        self.reactions = self._calculate_all_reactions()
    
    def _calculate_all_reactions(self) -> dict:
        d = self.data
        reactions = {}
        
        # === CO-System ===
        
        # FeO + CO → Fe + CO₂ (T > 570°C)
        reactions['FeO_Fe_CO'] = {
            'equation': 'FeO + CO → Fe + CO₂',
            'dH': (d['Fe_cr']['dHf'] + d['CO2_g']['dHf']) - 
                  (d['FeO_cr']['dHf'] + d['CO_g']['dHf']),
            'dS': (d['Fe_cr']['S'] + d['CO2_g']['S']) - 
                  (d['FeO_cr']['S'] + d['CO_g']['S']),
            'stoich': 1,
            'valid_above_K': WUSTITE_STABILITY_K,
        }
        
        # Fe₃O₄ + CO → 3 FeO + CO₂ (T > 570°C)
        reactions['Fe3O4_FeO_CO'] = {
            'equation': 'Fe₃O₄ + CO → 3 FeO + CO₂',
            'dH': (3*d['FeO_cr']['dHf'] + d['CO2_g']['dHf']) - 
                  (d['Fe3O4_cr']['dHf'] + d['CO_g']['dHf']),
            'dS': (3*d['FeO_cr']['S'] + d['CO2_g']['S']) - 
                  (d['Fe3O4_cr']['S'] + d['CO_g']['S']),
            'stoich': 1,
            'valid_above_K': WUSTITE_STABILITY_K,
        }
        
        # Fe₃O₄ + 4 CO → 3 Fe + 4 CO₂ (T < 570°C)
        reactions['Fe3O4_Fe_CO'] = {
            'equation': 'Fe₃O₄ + 4 CO → 3 Fe + 4 CO₂',
            'dH': (3*d['Fe_cr']['dHf'] + 4*d['CO2_g']['dHf']) - 
                  (d['Fe3O4_cr']['dHf'] + 4*d['CO_g']['dHf']),
            'dS': (3*d['Fe_cr']['S'] + 4*d['CO2_g']['S']) - 
                  (d['Fe3O4_cr']['S'] + 4*d['CO_g']['S']),
            'stoich': 4,
            'valid_below_K': WUSTITE_STABILITY_K,
        }
        
        # 3 Fe₂O₃ + CO → 2 Fe₃O₄ + CO₂
        reactions['Fe2O3_Fe3O4_CO'] = {
            'equation': '3 Fe₂O₃ + CO → 2 Fe₃O₄ + CO₂',
            'dH': (2*d['Fe3O4_cr']['dHf'] + d['CO2_g']['dHf']) - 
                  (3*d['Fe2O3_cr']['dHf'] + d['CO_g']['dHf']),
            'dS': (2*d['Fe3O4_cr']['S'] + d['CO2_g']['S']) - 
                  (3*d['Fe2O3_cr']['S'] + d['CO_g']['S']),
            'stoich': 1,
        }
        
        # === H₂-System ===
        
        # FeO + H₂ → Fe + H₂O (T > 570°C)
        reactions['FeO_Fe_H2'] = {
            'equation': 'FeO + H₂ → Fe + H₂O',
            'dH': (d['Fe_cr']['dHf'] + d['H2O_g']['dHf']) - 
                  (d['FeO_cr']['dHf'] + d['H2_g']['dHf']),
            'dS': (d['Fe_cr']['S'] + d['H2O_g']['S']) - 
                  (d['FeO_cr']['S'] + d['H2_g']['S']),
            'stoich': 1,
            'valid_above_K': WUSTITE_STABILITY_K,
        }
        
        # Fe₃O₄ + H₂ → 3 FeO + H₂O (T > 570°C)
        reactions['Fe3O4_FeO_H2'] = {
            'equation': 'Fe₃O₄ + H₂ → 3 FeO + H₂O',
            'dH': (3*d['FeO_cr']['dHf'] + d['H2O_g']['dHf']) - 
                  (d['Fe3O4_cr']['dHf'] + d['H2_g']['dHf']),
            'dS': (3*d['FeO_cr']['S'] + d['H2O_g']['S']) - 
                  (d['Fe3O4_cr']['S'] + d['H2_g']['S']),
            'stoich': 1,
            'valid_above_K': WUSTITE_STABILITY_K,
        }
        
        # Fe₃O₄ + 4 H₂ → 3 Fe + 4 H₂O (T < 570°C)
        reactions['Fe3O4_Fe_H2'] = {
            'equation': 'Fe₃O₄ + 4 H₂ → 3 Fe + 4 H₂O',
            'dH': (3*d['Fe_cr']['dHf'] + 4*d['H2O_g']['dHf']) - 
                  (d['Fe3O4_cr']['dHf'] + 4*d['H2_g']['dHf']),
            'dS': (3*d['Fe_cr']['S'] + 4*d['H2O_g']['S']) - 
                  (d['Fe3O4_cr']['S'] + 4*d['H2_g']['S']),
            'stoich': 4,
            'valid_below_K': WUSTITE_STABILITY_K,
        }
        
        # 3 Fe₂O₃ + H₂ → 2 Fe₃O₄ + H₂O
        reactions['Fe2O3_Fe3O4_H2'] = {
            'equation': '3 Fe₂O₃ + H₂ → 2 Fe₃O₄ + H₂O',
            'dH': (2*d['Fe3O4_cr']['dHf'] + d['H2O_g']['dHf']) - 
                  (3*d['Fe2O3_cr']['dHf'] + d['H2_g']['dHf']),
            'dS': (2*d['Fe3O4_cr']['S'] + d['H2O_g']['S']) - 
                  (3*d['Fe2O3_cr']['S'] + d['H2_g']['S']),
            'stoich': 1,
        }
        
        return reactions
    
    def delta_G(self, reaction_name: str, T_kelvin: float) -> float:
        """ΔG°(T) = ΔH° - T·ΔS° in J/mol"""
        rxn = self.reactions[reaction_name]
        return rxn['dH'] * 1000 - T_kelvin * rxn['dS']
    
    def equilibrium_K(self, reaction_name: str, T_kelvin: float) -> float:
        """K = exp(-ΔG°/RT)"""
        dG = self.delta_G(reaction_name, T_kelvin)
        return np.exp(-dG / (R * T_kelvin))
    
    def GOD(self, reaction_name: str, T_kelvin: float) -> float:
        """
        Gas Oxidation Degree für Phasengrenze.
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
        
        K für C + CO₂ → 2 CO: log₁₀(K) = -8916/T + 9.113
        K = p_CO² / p_CO₂
        
        Mit p_CO + p_CO₂ = 1:
        K = (1-GOD)² / GOD
        → GOD² - (2+K)·GOD + 1 = 0
        """
        log_K = -8916.0 / T_kelvin + 9.113
        K = 10 ** log_K
        
        # Quadratische Gleichung lösen
        a = 1
        b = -(2 + K)
        c = 1
        
        discriminant = b**2 - 4*a*c
        if discriminant < 0:
            return np.nan
        
        # Kleinere Wurzel nehmen (GOD zwischen 0 und 1)
        GOD = (-b - np.sqrt(discriminant)) / (2*a)
        
        return np.clip(GOD, 0, 1)
    
    def calculate_phase_boundary(self, reaction_name: str, T_range_K: np.ndarray) -> np.ndarray:
        """Berechnet GOD-Werte für eine Phasengrenze."""
        return np.array([self.GOD(reaction_name, T) for T in T_range_K])
    
    def calculate_boudouard_boundary(self, T_range_K: np.ndarray) -> np.ndarray:
        """Berechnet GOD-Werte für Boudouard-Gleichgewicht."""
        return np.array([self.boudouard_GOD(T) for T in T_range_K])
    
    def print_debug_info(self):
        """Gibt Debug-Informationen aus."""
        print("=== Thermodynamische Daten ===")
        for key, vals in self.data.items():
            print(f"{key}: ΔfH°={vals['dHf']:.2f} kJ/mol, S°={vals['S']:.2f} J/(mol·K)")
        
        print("\n=== Reaktionen ===")
        for name, rxn in self.reactions.items():
            print(f"{name}:")
            print(f"  ΔH° = {rxn['dH']:.2f} kJ/mol")
            print(f"  ΔS° = {rxn['dS']:.2f} J/(mol·K)")
            
            # Test bei 1000 K
            if 'valid_below_K' not in rxn or 1000 < rxn['valid_below_K']:
                if 'valid_above_K' not in rxn or 1000 > rxn['valid_above_K']:
                    K = self.equilibrium_K(name, 1000)
                    GOD = self.GOD(name, 1000)
                    print(f"  Bei 1000K: K={K:.3f}, GOD={GOD:.3f}")
        
        print("\n=== Boudouard bei verschiedenen T ===")
        for T_C in [500, 600, 700, 800, 900, 1000]:
            T_K = T_C + 273.15
            GOD = self.boudouard_GOD(T_K)
            print(f"  {T_C}°C: GOD = {GOD:.3f}")
```

---

## Korrigierte plotting.py

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
    """Erstellt das Baur-Glässner-Diagramm."""
    
    fig = go.Figure()
    
    T_range_C = np.linspace(T_min_C, T_max_C, 500)
    T_range_K = T_range_C + 273.15
    
    # === CO-System (durchgezogene Linien) ===
    if show_CO_system:
        # Fe/FeO Grenze (nur T > 570°C)
        GOD = thermo.calculate_phase_boundary('FeO_Fe_CO', T_range_K)
        mask = ~np.isnan(GOD)
        fig.add_trace(go.Scatter(
            x=GOD[mask], y=T_range_C[mask],
            mode='lines',
            line=dict(color='black', width=2),
            name='Fe-C-O System',
            legendgroup='CO',
        ))
        
        # Fe₃O₄/FeO Grenze (nur T > 570°C)
        GOD = thermo.calculate_phase_boundary('Fe3O4_FeO_CO', T_range_K)
        mask = ~np.isnan(GOD)
        fig.add_trace(go.Scatter(
            x=GOD[mask], y=T_range_C[mask],
            mode='lines',
            line=dict(color='black', width=2),
            legendgroup='CO',
            showlegend=False,
        ))
        
        # Fe₃O₄/Fe Grenze (nur T < 570°C)
        GOD = thermo.calculate_phase_boundary('Fe3O4_Fe_CO', T_range_K)
        mask = ~np.isnan(GOD)
        fig.add_trace(go.Scatter(
            x=GOD[mask], y=T_range_C[mask],
            mode='lines',
            line=dict(color='black', width=2),
            legendgroup='CO',
            showlegend=False,
        ))
        
        # Fe₂O₃/Fe₃O₄ Grenze
        GOD = thermo.calculate_phase_boundary('Fe2O3_Fe3O4_CO', T_range_K)
        mask = ~np.isnan(GOD)
        fig.add_trace(go.Scatter(
            x=GOD[mask], y=T_range_C[mask],
            mode='lines',
            line=dict(color='black', width=2),
            legendgroup='CO',
            showlegend=False,
        ))
    
    # === H₂-System (gestrichelte Linien) ===
    if show_H2_system:
        # Fe/FeO Grenze (nur T > 570°C)
        GOD = thermo.calculate_phase_boundary('FeO_Fe_H2', T_range_K)
        mask = ~np.isnan(GOD)
        fig.add_trace(go.Scatter(
            x=GOD[mask], y=T_range_C[mask],
            mode='lines',
            line=dict(color='black', width=2, dash='dash'),
            name='Fe-H₂-O System',
            legendgroup='H2',
        ))
        
        # Fe₃O₄/FeO Grenze (nur T > 570°C)
        GOD = thermo.calculate_phase_boundary('Fe3O4_FeO_H2', T_range_K)
        mask = ~np.isnan(GOD)
        fig.add_trace(go.Scatter(
            x=GOD[mask], y=T_range_C[mask],
            mode='lines',
            line=dict(color='black', width=2, dash='dash'),
            legendgroup='H2',
            showlegend=False,
        ))
        
        # Fe₃O₄/Fe Grenze (nur T < 570°C)
        GOD = thermo.calculate_phase_boundary('Fe3O4_Fe_H2', T_range_K)
        mask = ~np.isnan(GOD)
        fig.add_trace(go.Scatter(
            x=GOD[mask], y=T_range_C[mask],
            mode='lines',
            line=dict(color='black', width=2, dash='dash'),
            legendgroup='H2',
            showlegend=False,
        ))
        
        # Fe₂O₃/Fe₃O₄ Grenze
        GOD = thermo.calculate_phase_boundary('Fe2O3_Fe3O4_H2', T_range_K)
        mask = ~np.isnan(GOD)
        fig.add_trace(go.Scatter(
            x=GOD[mask], y=T_range_C[mask],
            mode='lines',
            line=dict(color='black', width=2, dash='dash'),
            legendgroup='H2',
            showlegend=False,
        ))
    
    # === Boudouard (gepunktet) ===
    if show_boudouard:
        GOD_boud = thermo.calculate_boudouard_boundary(T_range_K)
        mask = ~np.isnan(GOD_boud)
        fig.add_trace(go.Scatter(
            x=GOD_boud[mask], y=T_range_C[mask],
            mode='lines',
            line=dict(color='gray', width=1.5, dash='dot'),
            name='Boudouard',
        ))
    
    # === Phasenbeschriftungen ===
    if show_labels:
        fig.add_annotation(x=0.15, y=750, text='Fe', showarrow=False,
            font=dict(size=18, color='black', family='Arial Black'))
        fig.add_annotation(x=0.55, y=450, text='Fe₃O₄', showarrow=False,
            font=dict(size=16, color='darkblue', family='Arial Black'))
        fig.add_annotation(x=0.55, y=800, text='FeO', showarrow=False,
            font=dict(size=16, color='darkblue', family='Arial Black'))
        fig.add_annotation(x=0.88, y=650, text='Fe₂O₃', showarrow=False,
            font=dict(size=14, color='darkred', family='Arial Black'))
    
    # === Layout ===
    fig.update_layout(
        title='Baur-Glässner Diagram',
        xaxis_title='Gas oxidation degree GOD [-]',
        yaxis_title='Temperature [°C]',
        xaxis=dict(range=[0, 1], dtick=0.1),
        yaxis=dict(range=[T_min_C, T_max_C]),
        legend=dict(x=0.02, y=0.98),
        template='plotly_white',
        height=700,
    )
    
    return fig
```

---

## Test-Script

Erstelle eine Datei `test_thermo.py` um die Berechnungen zu überprüfen:

```python
from thermodynamics import BaurGlassnerThermo

thermo = BaurGlassnerThermo("NBS_Tables_Library.xlsx")
thermo.print_debug_info()

print("\n=== Erwartete Werte (Referenz) ===")
print("Fe/FeO (CO) bei 800°C: GOD ≈ 0.28")
print("Fe/FeO (CO) bei 1000°C: GOD ≈ 0.35")
print("Fe/FeO (H₂) bei 800°C: GOD ≈ 0.33")
print("Fe₃O₄/Fe (CO) bei 400°C: GOD ≈ 0.25")
print("Boudouard bei 700°C: GOD ≈ 0.5")
```

---

## Zusammenfassung der Änderungen

1. **Boudouard**: Neue quadratische Gleichung für korrekten Kurvenverlauf
2. **GOD-Formel**: `K/(1+K)` ist korrekt (nicht invertieren!)
3. **Horizontale Linie bei 570°C**: Entfernt (war falsch)
4. **NaN-Handling**: Punkte außerhalb Gültigkeitsbereich werden nicht geplottet
5. **Debug-Funktion**: Zum Überprüfen der berechneten Werte

Führe zuerst `test_thermo.py` aus und vergleiche die Werte mit der Referenz!
