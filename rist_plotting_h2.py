"""
Plotting module for the H2 Rist Diagram.

Creates interactive Plotly diagrams for the visualization of
mass balance in hydrogen direct reduction processes.

X-Axis: O/H2 from 1 to 2 (only gas phase, no direct reduction)
- O/H2 = 1: 100% H2 (fresh reducing gas)
- O/H2 = 2: 100% H2O (fully oxidized gas)
"""

import numpy as np
import plotly.graph_objects as go
from typing import Optional, Tuple

from thermodynamics import BaurGlassnerThermo, WUSTITE_STABILITY_K
from rist_diagram_h2 import RistDiagramH2, O_FE_RATIOS, god_to_o_h2, o_h2_to_god
from mass_balance_h2 import H2MassBalanceResults, calculate_h2_operating_line_points


# Color scheme (similar to CO diagram but with H2-specific accent)
RIST_H2_COLORS = {
    'oxide_lines': '#666666',        # Gray for oxide lines
    'operating_line': '#0066CC',     # Blue for ideal operating line
    'operating_line_real': '#9933CC',  # Purple for real operating line from mass balance
    'forbidden_zone': 'rgba(100, 150, 255, 0.3)',  # Blue transparent (H2 = blue)
    'w_point': '#CC0000',            # Red for W point
    'm_point': '#009900',            # Green for M point
    'equilibrium': '#000066',        # Dark blue for equilibrium lines
    'indirect_reduction': 'rgba(100, 200, 100, 0.25)',  # Green transparent for indirect reduction
}


def create_h2_rist_diagram(
    thermo: BaurGlassnerThermo,
    T_celsius: float,
    show_forbidden_zone: bool = True,
    show_equilibrium_points: bool = True,
    mass_balance_results: Optional[H2MassBalanceResults] = None,
) -> go.Figure:
    """
    Creates the interactive H2 Rist diagram for a given temperature.

    Args:
        thermo: BaurGlassnerThermo instance
        T_celsius: Temperature in degrees Celsius
        show_forbidden_zone: Show thermodynamically forbidden zone
        show_equilibrium_points: Show W and M points
        mass_balance_results: Optional mass balance results for operating line

    Returns:
        Plotly Figure object
    """
    fig = go.Figure()

    T_kelvin = T_celsius + 273.15
    rist = RistDiagramH2(thermo)

    # === Horizontal lines for iron oxides ===
    oxide_labels = [
        ('Fe', O_FE_RATIOS['Fe'], 'Iron (Fe)'),
        ('FeO', O_FE_RATIOS['FeO'], 'Wustite (Fe₁₋ₓO)'),
        ('Fe₃O₄', O_FE_RATIOS['Fe3O4'], 'Magnetite (Fe₃O₄)'),
        ('Fe₂O₃', O_FE_RATIOS['Fe2O3'], 'Hematite (Fe₂O₃)'),
    ]

    for name, o_fe, label in oxide_labels:
        fig.add_trace(go.Scatter(
            x=[1.0, 2.0],
            y=[o_fe, o_fe],
            mode='lines',
            line=dict(color=RIST_H2_COLORS['oxide_lines'], width=1, dash='dot'),
            name=label,
            showlegend=False,
            hovertemplate=f'{label}<br>O/Fe = {o_fe:.3f}<extra></extra>',
        ))

        # Label at left edge (at x=1.0 for H2 diagram)
        fig.add_annotation(
            x=1.02,
            y=o_fe,
            text=f'<b>{name}</b>',
            showarrow=False,
            xanchor='left',
            yanchor='bottom',
            yshift=2,
            font=dict(size=12, color=RIST_H2_COLORS['oxide_lines']),
        )

    # === Vertical line at O/H2 = 1 (100% H2) ===
    fig.add_trace(go.Scatter(
        x=[1.0, 1.0],
        y=[0, O_FE_RATIOS['Fe2O3']],
        mode='lines',
        line=dict(color='rgba(100, 100, 100, 0.5)', width=1, dash='dash'),
        name='100% H₂',
        showlegend=False,
        hovertemplate='O/H₂ = 1 (100% H₂)<extra></extra>',
    ))

    # === Equilibrium curves from Baur-Glassner (H2 system) ===
    if T_kelvin >= WUSTITE_STABILITY_K:
        # Fe/FeO equilibrium (H2)
        god_fe_feo = thermo.GOD('FeO_Fe_H2', T_kelvin)
        o_h2_fe_feo = god_to_o_h2(god_fe_feo)

        # Vertical line at Fe/FeO equilibrium
        fig.add_trace(go.Scatter(
            x=[o_h2_fe_feo, o_h2_fe_feo],
            y=[O_FE_RATIOS['Fe'], O_FE_RATIOS['FeO']],
            mode='lines',
            line=dict(color=RIST_H2_COLORS['equilibrium'], width=2),
            name='Fe/FeO Equilibrium (H₂)',
            showlegend=True,
            hovertemplate=f'Fe/FeO Equilibrium (H₂)<br>O/H₂ = {o_h2_fe_feo:.3f}<br>GOD = {god_fe_feo:.3f}<extra></extra>',
        ))

        # FeO/Fe3O4 equilibrium (H2)
        god_feo_fe3o4 = thermo.GOD('Fe3O4_FeO_H2', T_kelvin)
        o_h2_feo_fe3o4 = god_to_o_h2(god_feo_fe3o4)

        fig.add_trace(go.Scatter(
            x=[o_h2_feo_fe3o4, o_h2_feo_fe3o4],
            y=[O_FE_RATIOS['FeO'], O_FE_RATIOS['Fe3O4']],
            mode='lines',
            line=dict(color=RIST_H2_COLORS['equilibrium'], width=2),
            name='FeO/Fe₃O₄ Equilibrium (H₂)',
            showlegend=True,
            hovertemplate=f'FeO/Fe₃O₄ Equilibrium (H₂)<br>O/H₂ = {o_h2_feo_fe3o4:.3f}<br>GOD = {god_feo_fe3o4:.3f}<extra></extra>',
        ))
    else:
        # Below 570C: Only Fe/Fe3O4 equilibrium
        god_fe_fe3o4 = thermo.GOD('Fe3O4_Fe_H2', T_kelvin)
        o_h2_fe_fe3o4 = god_to_o_h2(god_fe_fe3o4)

        fig.add_trace(go.Scatter(
            x=[o_h2_fe_fe3o4, o_h2_fe_fe3o4],
            y=[O_FE_RATIOS['Fe'], O_FE_RATIOS['Fe3O4']],
            mode='lines',
            line=dict(color=RIST_H2_COLORS['equilibrium'], width=2),
            name='Fe/Fe₃O₄ Equilibrium (H₂)',
            showlegend=True,
            hovertemplate=f'Fe/Fe₃O₄ Equilibrium (H₂)<br>O/H₂ = {o_h2_fe_fe3o4:.3f}<br>GOD = {god_fe_fe3o4:.3f}<extra></extra>',
        ))

    # === Thermodynamically forbidden zone (combined W and M) ===
    if show_forbidden_zone:
        vertices = rist.get_forbidden_zone_combined(T_kelvin)
        x_coords = [v[0] for v in vertices] + [vertices[0][0]]  # Close polygon
        y_coords = [v[1] for v in vertices] + [vertices[0][1]]

        fig.add_trace(go.Scatter(
            x=x_coords,
            y=y_coords,
            fill='toself',
            fillcolor=RIST_H2_COLORS['forbidden_zone'],
            line=dict(color='rgba(100, 150, 255, 0.5)', width=1),
            name='Forbidden Zone (H₂)',
            showlegend=True,
            hovertemplate='Thermodynamically forbidden zone<extra></extra>',
        ))

    # === W point (Wustite) and M point (Magnetite) ===
    if show_equilibrium_points:
        w_point = rist.get_wustite_point(T_kelvin)
        m_point = rist.get_magnetite_point(T_kelvin)

        # W point
        god_w = o_h2_to_god(w_point[0])
        god_w_str = f"{god_w:.3f}" if god_w is not None else "N/A"
        fig.add_trace(go.Scatter(
            x=[w_point[0]],
            y=[w_point[1]],
            mode='markers+text',
            marker=dict(size=14, color=RIST_H2_COLORS['w_point'], symbol='circle'),
            text=['W'],
            textposition='top center',
            textfont=dict(size=14, color=RIST_H2_COLORS['w_point'], family='Arial Black'),
            name='Wustite Point (W)',
            showlegend=True,
            hovertemplate=f'Wustite Point (W)<br>O/H₂ = {w_point[0]:.3f}<br>O/Fe = {w_point[1]:.3f}<br>GOD = {god_w_str}<extra></extra>',
        ))

        # M point
        god_m = o_h2_to_god(m_point[0])
        god_m_str = f"{god_m:.3f}" if god_m is not None else "N/A"
        fig.add_trace(go.Scatter(
            x=[m_point[0]],
            y=[m_point[1]],
            mode='markers+text',
            marker=dict(size=14, color=RIST_H2_COLORS['m_point'], symbol='circle'),
            text=['M'],
            textposition='top center',
            textfont=dict(size=14, color=RIST_H2_COLORS['m_point'], family='Arial Black'),
            name='Magnetite Point (M)',
            showlegend=True,
            hovertemplate=f'Magnetite Point (M)<br>O/H₂ = {m_point[0]:.3f}<br>O/Fe = {m_point[1]:.3f}<br>GOD = {god_m_str}<extra></extra>',
        ))

    # === Operating line from mass balance ===
    if mass_balance_results is not None:
        mb = mass_balance_results
        start_mb, end_mb = calculate_h2_operating_line_points(mb)

        # Point at O/H2 = 1 (pure H2 input)
        o_fe_at_oh1 = mb.slope_k * 1.0 + mb.intercept_d

        # === Indirect reduction zone (all reduction is indirect with H2) ===
        # Use dynamic o_fe_ore from mass balance (1.5 for Fe2O3, 1.333 for Fe3O4)
        o_fe_ore = mb.o_fe_ore
        
        if 0 <= o_fe_at_oh1 <= o_fe_ore:
            # H2 reduction is 100% indirect (gas-based)
            fig.add_trace(go.Scatter(
                x=[1.0, 1.0, end_mb[0], 1.0],
                y=[o_fe_at_oh1, o_fe_ore, end_mb[1], o_fe_at_oh1],
                fill='toself',
                fillcolor=RIST_H2_COLORS['indirect_reduction'],
                line=dict(color='rgba(100, 200, 100, 0.4)', width=1),
                name=f'Indirect Reduction (100%)',
                showlegend=True,
                hovertemplate=f'Indirect Reduction<br>100%<br>(H₂/H₂O gas)<extra></extra>',
            ))

        # Draw operating line (above the zones)
        fig.add_trace(go.Scatter(
            x=[start_mb[0], end_mb[0]],
            y=[start_mb[1], end_mb[1]],
            mode='lines+markers',
            line=dict(color=RIST_H2_COLORS['operating_line_real'], width=3),
            marker=dict(size=10, color=RIST_H2_COLORS['operating_line_real'], symbol='diamond'),
            name=f'Operating Line (H₂)',
            showlegend=True,
            hovertemplate='Operating Line (H₂)<br>O/H₂ = %{x:.3f}<br>O/Fe = %{y:.3f}<extra></extra>',
        ))

        # Mark point at O/H2 = 1
        if 0 <= o_fe_at_oh1 <= o_fe_ore:
            fig.add_trace(go.Scatter(
                x=[1.0],
                y=[o_fe_at_oh1],
                mode='markers',
                marker=dict(size=14, color=RIST_H2_COLORS['operating_line_real'], symbol='star'),
                name='H₂ Input Point (O/H₂=1)',
                showlegend=True,
                hovertemplate=f'H₂ Input Point<br>O/H₂ = 1.0<br>O/Fe = {o_fe_at_oh1:.3f}<extra></extra>',
            ))

        # Annotation for H2/Fe
        mid_x_mb = (start_mb[0] + end_mb[0]) / 2
        mid_y_mb = (start_mb[1] + end_mb[1]) / 2
        fig.add_annotation(
            x=mid_x_mb,
            y=mid_y_mb,
            text=f'H₂/Fe = {mb.slope_k:.2f}<br>H₂: {mb.h2_amount_nm3:.0f} Nm³/t',
            showarrow=True,
            arrowhead=0,
            ax=-60,
            ay=40,
            font=dict(size=10, color=RIST_H2_COLORS['operating_line_real']),
            bgcolor='white',
            bordercolor=RIST_H2_COLORS['operating_line_real'],
            borderwidth=1,
        )

    # === Layout ===
    fig.update_layout(
        title=dict(
            text=f'<b>H₂ Rist Diagram at {T_celsius:.0f}°C</b>',
            x=0.5,
            xanchor='center',
            font=dict(size=16),
        ),
        xaxis=dict(
            title='O/H₂ [-]  (Oxygen per Hydrogen)',
            range=[0.9, 2.1],
            dtick=0.25,
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(200, 200, 200, 0.5)',
            zeroline=False,
        ),
        yaxis=dict(
            title='O/Fe [-]  (Oxygen per Iron)',
            range=[-0.1, 1.6],
            dtick=0.25,
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(200, 200, 200, 0.5)',
            zeroline=True,
            zerolinewidth=2,
            zerolinecolor='rgba(0, 0, 0, 0.3)',
        ),
        legend=dict(
            x=0.75,
            y=0.1,
            bgcolor='rgba(255, 255, 255, 0.9)',
            bordercolor='gray',
            borderwidth=1,
            font=dict(size=10),
        ),
        template='plotly_white',
        height=550,
        margin=dict(l=80, r=40, t=60, b=80),
        hovermode='closest',
    )

    # Annotations for X-axis zones (only H2 and H2O, no direct reduction)
    fig.add_annotation(
        x=1.0,
        y=-0.08,
        text='100% H₂',
        showarrow=False,
        font=dict(size=10, color='gray'),
        xref='x',
        yref='y',
    )
    fig.add_annotation(
        x=2.0,
        y=-0.08,
        text='100% H₂O',
        showarrow=False,
        font=dict(size=10, color='gray'),
        xref='x',
        yref='y',
    )


    return fig


if __name__ == '__main__':
    # Test: Create diagram
    print("Creating H2 Rist Diagram...")

    thermo = BaurGlassnerThermo()

    for T_C in [700, 900, 1000]:
        fig = create_h2_rist_diagram(
            thermo=thermo,
            T_celsius=T_C,
            show_forbidden_zone=True,
            show_equilibrium_points=True,
        )
        filename = f'h2_rist_diagram_{T_C}C.html'
        fig.write_html(filename)
        print(f"Diagram saved as '{filename}'")
