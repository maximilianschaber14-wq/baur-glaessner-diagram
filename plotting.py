"""
Plotting-Modul fuer das Baur-Glaessner-Diagramm.

Erstellt interaktive Plotly-Diagramme fuer die Fe-C-O und Fe-H2-O Systeme.
Basierend auf dem klassischen Baur-Glaessner-Diagramm mit GOD-Darstellung.
"""

import numpy as np
import plotly.graph_objects as go
from thermodynamics import BaurGlassnerThermo, WUSTITE_STABILITY_K


# Farbschema
COLORS = {
    'CO_system': '#000000',      # Schwarz fuer CO-System (durchgezogen)
    'H2_system': '#000000',      # Schwarz fuer H2-System (gestrichelt)
    'boudouard': '#000000',      # Schwarz fuer Boudouard (gepunktet)
    'wustite_line': '#7f7f7f',   # Grau fuer Wuestit-Stabilitaet
    'Fe': '#4a4a4a',
    'FeO': '#4a4a4a',
    'Fe3O4': '#4a4a4a',
    'Fe2O3': '#4a4a4a',
}


from typing import List, Tuple, Optional


def create_baur_glassner_diagram(
    thermo: BaurGlassnerThermo,
    show_CO_system: bool = True,
    show_H2_system: bool = True,
    show_boudouard: bool = True,
    show_labels: bool = True,
    show_wustite_line: bool = True,
    T_min_C: float = 300,
    T_max_C: float = 1000,
    custom_points: Optional[List[Tuple[float, float, int]]] = None,
) -> go.Figure:
    """
    Erstellt das interaktive Baur-Glaessner-Diagramm im GOD-Format.

    Args:
        thermo: BaurGlassnerThermo-Instanz
        show_CO_system: CO/CO2-Gleichgewichte anzeigen
        show_H2_system: H2/H2O-Gleichgewichte anzeigen
        show_boudouard: Boudouard-Linie anzeigen
        show_labels: Phasenbeschriftungen anzeigen
        show_wustite_line: Wuestit-Stabilitaetslinie bei 570C anzeigen
        T_min_C: Minimale Temperatur in Celsius
        T_max_C: Maximale Temperatur in Celsius
        custom_points: Liste von (GOD, Temperatur_C, Nummer) Tupeln

    Returns:
        Plotly Figure-Objekt
    """
    fig = go.Figure()

    # Temperaturbereich
    T_range_C = np.linspace(T_min_C, T_max_C, 500)
    T_range_K = T_range_C + 273.15

    # Wuestit-Stabilitaetstemperatur
    T_wustite_C = WUSTITE_STABILITY_K - 273.15  # ~570C

    # Berechne GOD am Tripelpunkt (570C) fuer Fe/FeO CO-System
    T_triple_K = WUSTITE_STABILITY_K
    GOD_triple = thermo.GOD('FeO_Fe_CO', T_triple_K + 1)  # knapp oberhalb

    # === CO-System (durchgezogene Linien) ===
    if show_CO_system:
        # Fe/FeO Grenze (oberhalb 570C)
        T_above = T_range_C[T_range_C > T_wustite_C]
        T_above_K = T_above + 273.15
        if len(T_above) > 0:
            god_fe_feo = np.array([thermo.GOD('FeO_Fe_CO', T) for T in T_above_K])
            valid = ~np.isnan(god_fe_feo)
            if np.any(valid):
                fig.add_trace(go.Scatter(
                    x=god_fe_feo[valid],
                    y=T_above[valid],
                    mode='lines',
                    line=dict(color=COLORS['CO_system'], width=2),
                    name='System CO/CO₂',
                    legendgroup='CO',
                    showlegend=True,
                    hovertemplate='FeO + CO = Fe + CO2<br>GOD: %{x:.2f}<br>T: %{y:.0f}C<extra></extra>'
                ))

        # FeO/Fe3O4 Grenze (oberhalb 570C)
        if len(T_above) > 0:
            god_feo_fe3o4 = np.array([thermo.GOD('Fe3O4_FeO_CO', T) for T in T_above_K])
            valid = ~np.isnan(god_feo_fe3o4)
            if np.any(valid):
                fig.add_trace(go.Scatter(
                    x=god_feo_fe3o4[valid],
                    y=T_above[valid],
                    mode='lines',
                    line=dict(color=COLORS['CO_system'], width=2),
                    name='System CO/CO₂',
                    legendgroup='CO',
                    showlegend=False,
                    hovertemplate='Fe₃O₄ + CO = 3FeO + CO₂<br>GOD: %{x:.2f}<br>T: %{y:.0f}C<extra></extra>'
                ))

        # Fe/Fe3O4 Grenze (unterhalb 570C) - SENKRECHTE LINIE bei GOD des Tripelpunkts
        T_below = T_range_C[T_range_C <= T_wustite_C]
        if len(T_below) > 0:
            # Senkrechte Linie bei GOD_triple
            fig.add_trace(go.Scatter(
                x=[GOD_triple, GOD_triple],
                y=[T_min_C, T_wustite_C],
                mode='lines',
                line=dict(color=COLORS['CO_system'], width=2),
                name='System CO/CO₂',
                legendgroup='CO',
                showlegend=False,
                hovertemplate='¼Fe₃O₄ + CO = ¾Fe + CO₂<br>GOD: %{x:.2f}<br>T: %{y:.0f}C<extra></extra>'
            ))

    # === H2-System (gestrichelte Linien) ===
    if show_H2_system:
        # Fe/FeO Grenze (oberhalb 570C)
        T_above = T_range_C[T_range_C > T_wustite_C]
        T_above_K = T_above + 273.15
        if len(T_above) > 0:
            god_fe_feo_h2 = np.array([thermo.GOD('FeO_Fe_H2', T) for T in T_above_K])
            valid = ~np.isnan(god_fe_feo_h2)
            if np.any(valid):
                fig.add_trace(go.Scatter(
                    x=god_fe_feo_h2[valid],
                    y=T_above[valid],
                    mode='lines',
                    line=dict(color=COLORS['H2_system'], width=2, dash='dash'),
                    name='System H₂/H₂O',
                    legendgroup='H2',
                    showlegend=True,
                    hovertemplate='FeO + H2 = Fe + H2O<br>GOD: %{x:.2f}<br>T: %{y:.0f}C<extra></extra>'
                ))

        # FeO/Fe3O4 Grenze (oberhalb 570C)
        if len(T_above) > 0:
            god_feo_fe3o4_h2 = np.array([thermo.GOD('Fe3O4_FeO_H2', T) for T in T_above_K])
            valid = ~np.isnan(god_feo_fe3o4_h2)
            if np.any(valid):
                fig.add_trace(go.Scatter(
                    x=god_feo_fe3o4_h2[valid],
                    y=T_above[valid],
                    mode='lines',
                    line=dict(color=COLORS['H2_system'], width=2, dash='dash'),
                    name='System H₂/H₂O',
                    legendgroup='H2',
                    showlegend=False,
                    hovertemplate='Fe₃O₄ + H₂ = 3FeO + H₂O<br>GOD: %{x:.2f}<br>T: %{y:.0f}C<extra></extra>'
                ))

        # Fe/Fe3O4 Grenze H2 (unterhalb 570C) - NICHT senkrecht, hat Steigung
        T_below = T_range_C[T_range_C <= T_wustite_C]
        T_below_K = T_below + 273.15
        if len(T_below) > 0:
            god_fe_fe3o4_h2 = np.array([thermo.GOD('Fe3O4_Fe_H2', T) for T in T_below_K])
            valid = ~np.isnan(god_fe_fe3o4_h2)
            if np.any(valid):
                fig.add_trace(go.Scatter(
                    x=god_fe_fe3o4_h2[valid],
                    y=T_below[valid],
                    mode='lines',
                    line=dict(color=COLORS['H2_system'], width=2, dash='dash'),
                    name='System H₂/H₂O',
                    legendgroup='H2',
                    showlegend=False,
                    hovertemplate='¼Fe₃O₄ + H₂ = ¾Fe + H₂O<br>GOD: %{x:.2f}<br>T: %{y:.0f}C<extra></extra>'
                ))

    # === Boudouard-Linie (gepunktet) ===
    if show_boudouard:
        god_boud = np.array([thermo.boudouard_GOD(T) for T in T_range_K])
        fig.add_trace(go.Scatter(
            x=god_boud,
            y=T_range_C,
            mode='lines',
            line=dict(color=COLORS['boudouard'], width=1.5, dash='dot'),
            name='2CO = C + CO₂',
            hovertemplate='Boudouard: 2CO = C + CO₂<br>GOD: %{x:.2f}<br>T: %{y:.0f}C<extra></extra>'
        ))

    # === Phasenbeschriftungen ===
    if show_labels:
        # Berechne typische GOD-Werte für die Phasengrenzen zur besseren Positionierung
        T_mid_upper = (T_max_C + T_wustite_C) / 2 + 273.15  # Mitte des oberen Bereichs
        T_mid_lower = (T_wustite_C + T_min_C) / 2 + 273.15  # Mitte des unteren Bereichs

        # GOD-Werte an den Grenzen bei mittlerer Temperatur (oberer Bereich)
        GOD_fe_feo_mid = thermo.GOD('FeO_Fe_CO', T_mid_upper)
        GOD_feo_fe3o4_mid = thermo.GOD('Fe3O4_FeO_CO', T_mid_upper)

        # Fe-Bereich: links von der Fe/FeO Grenze, etwas tiefer positioniert
        fe_x = GOD_fe_feo_mid / 2  # Mitte zwischen GOD=0 und Fe/FeO-Grenze
        fe_y = T_wustite_C + (T_max_C - T_wustite_C) * 0.35  # 35% oberhalb von Wüstit-Linie
        fig.add_annotation(
            x=fe_x, y=fe_y,
            text='<b>Fe</b>',
            showarrow=False,
            font=dict(size=24, color=COLORS['Fe']),
        )

        # FeO-Bereich: zwischen Fe/FeO und FeO/Fe3O4 Grenzen
        feo_x = (GOD_fe_feo_mid + GOD_feo_fe3o4_mid) / 2
        feo_y = (T_max_C + T_wustite_C) / 2
        fig.add_annotation(
            x=feo_x, y=feo_y,
            text='<b>Fe<sub>(1-x)</sub>O</b>',
            showarrow=False,
            font=dict(size=20, color=COLORS['FeO']),
        )

        # Fe3O4-Bereich: rechts von FeO/Fe3O4 Grenze, etwas höher positioniert
        # Im unteren Bereich ist Fe3O4 der grosse Bereich rechts vom Tripelpunkt
        fe3o4_x = (GOD_triple + 0.75) / 2  # Zwischen Tripelpunkt und rechtem Rand
        fe3o4_y = T_min_C + (T_wustite_C - T_min_C) * 0.65  # 65% oberhalb von T_min
        fig.add_annotation(
            x=fe3o4_x, y=fe3o4_y,
            text='<b>Fe<sub>3</sub>O<sub>4</sub></b>',
            showarrow=False,
            font=dict(size=20, color=COLORS['Fe3O4']),
        )

        # Fe2O3-Pfeil am rechten Rand - Pfeil zeigt NACH RECHTS (Spitze bei GOD=1)
        fig.add_annotation(
            x=1.0, y=(T_max_C + T_min_C) / 2,
            text='<b>Fe<sub>2</sub>O<sub>3</sub></b>',
            showarrow=True,
            arrowhead=2,
            arrowsize=1.5,
            arrowwidth=2,
            ax=-50,  # Text links vom Pfeil
            ay=0,
            font=dict(size=16, color=COLORS['Fe2O3']),
        )

    # === Benutzerdefinierte Punkte ===
    if custom_points and len(custom_points) > 0:
        god_values = [p[0] for p in custom_points]
        temp_values = [p[1] for p in custom_points]
        labels = [str(p[2]) for p in custom_points]

        # Punkte als rote Marker
        fig.add_trace(go.Scatter(
            x=god_values,
            y=temp_values,
            mode='markers+text',
            marker=dict(
                size=12,
                color='red',
                symbol='circle',
                line=dict(width=1, color='darkred')
            ),
            text=labels,
            textposition='top right',
            textfont=dict(size=12, color='red', family='Arial Black'),
            name='Custom Points',
            hovertemplate='Point %{text}<br>GOD: %{x:.3f}<br>T: %{y:.0f}°C<extra></extra>',
            showlegend=True,
        ))

    # === Layout ===
    fig.update_layout(
        xaxis=dict(
            title='Gas oxidation degree GOD [-]',
            range=[0, 1],
            dtick=0.1,
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(200, 200, 200, 0.5)',
            zeroline=False,
        ),
        yaxis=dict(
            title='Temperature [°C]',
            range=[T_min_C, T_max_C],
            dtick=100,
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(200, 200, 200, 0.5)',
            zeroline=False,
        ),
        legend=dict(
            x=0.80,
            y=0.01,
            bgcolor='rgba(255, 255, 255, 0.9)',
            bordercolor='gray',
            borderwidth=1,
            font=dict(size=11),
        ),
        template='plotly_white',
        height=700,
        margin=dict(l=80, r=40, t=80, b=60),
        hovermode='closest',
    )

    return fig


if __name__ == '__main__':
    # Test: Diagramm erstellen und als HTML speichern
    print("Erstelle Baur-Glaessner-Diagramm...")

    thermo = BaurGlassnerThermo()

    fig = create_baur_glassner_diagram(
        thermo=thermo,
        show_CO_system=True,
        show_H2_system=True,
        show_boudouard=True,
        show_labels=True,
        show_wustite_line=True,
        T_min_C=300,
        T_max_C=1000,
    )
    fig.write_html('baur_glassner_test.html')
    print("Diagramm gespeichert als 'baur_glassner_test.html'")
