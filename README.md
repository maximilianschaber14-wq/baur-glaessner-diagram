# Baur-Glässner Diagram Generator

Interactive Streamlit application for visualizing thermodynamic phase diagrams and mass balances in iron oxide reduction processes.

## Features

### Baur-Glässner Diagram
- Phase equilibrium diagram for Fe-C-O and Fe-H₂-O systems
- Temperature-dependent stability regions of iron oxides
- Boudouard equilibrium visualization
- Custom point calculator

### CO Rist Diagram (Blast Furnace)
- Mass balance visualization for blast furnace processes
- Operating line calculation
- Direct/indirect reduction fractions
- Top gas (Gichtgas) composition and volume

### H₂ Rist Diagram (Direct Reduction)
- Mass balance for hydrogen-based direct reduction (H₂-DR)
- Gas utilization calculation
- Support for Fe₂O₃ (hematite) and Fe₃O₄ (magnetite) input

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
streamlit run app.py
```

## Thermodynamic Data

All calculations use **NIST-JANAF Thermochemical Tables** (Chase 1998) with Shomate polynomial coefficients for:
- Iron phases: Fe, FeO, Fe₃O₄, Fe₂O₃
- Gas species: CO, CO₂, H₂, H₂O, C (graphite)

## Author

Maximilian Schaber

## Disclaimer

All information is provided without guarantee. No responsibility is taken for the accuracy, completeness, or timeliness of the content.
