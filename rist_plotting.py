"""
Plotting module for the Rist Diagram.

Creates interactive Plotly diagrams for the visualization of
mass balance in iron oxide reduction processes.

X-Axis: O/C from 0 to 2
- O/C = 0: Direct reduction (100% carbon)
- O/C = 1: 100% CO
- O/C = 2: 100% CO2
"""

import numpy as np
import plotly.graph_objects as go
from typing import Optional, Tuple

from thermodynamics import BaurGlassnerThermo, WUSTITE_STABILITY_K
from rist_diagram import RistDiagram, O_FE_RATIOS, god_to_o_c, o_c_to_god
from mass_balance import MassBalanceResults, calculate_operating_line_points


# Color scheme
RIST_COLORS = {
    'oxide_lines': '#666666',       # Gray for oxide lines
    'operating_line': '#0066CC',    # Blue for ideal operating line
    'operating_line_real': '#FF6600',  # Orange for real operating line from mass balance
    'forbidden_zone': 'rgba(255, 100, 100, 0.3)',  # Red transparent
    'w_point': '#CC0000',           # Red for W point
    'm_point': '#009900',           # Green for M point
    'equilibrium': '#000000',       # Black for equilibrium lines
    'direct_reduction': 'rgba(100, 150, 255, 0.25)',  # Blue transparent for direct reduction
    'indirect_reduction': 'rgba(100, 200, 100, 0.25)',  # Green transparent for indirect reduction
}


def create_rist_diagram(
    thermo: BaurGlassnerThermo,
    T_celsius: float,
    show_operating_line: bool = True,
    show_forbidden_zone: bool = True,
    show_equilibrium_points: bool = True,
    mass_balance_results: Optional[MassBalanceResults] = None,
) -> go.Figure:
    """
    Creates the interactive Rist diagram for a given temperature.

    Args:
        thermo: BaurGlassnerThermo instance
        T_celsius: Temperature in degrees Celsius
        show_operating_line: Show ideal operating line
        show_forbidden_zone: Show thermodynamically forbidden zone
        show_equilibrium_points: Show W and M points
        mass_balance_results: Optional mass balance results for real operating line

    Returns:
        Plotly Figure object
    """
    fig = go.Figure()

    T_kelvin = T_celsius + 273.15
    rist = RistDiagram(thermo)

    # === Horizontal lines for iron oxides ===
    oxide_labels = [
        ('Fe', O_FE_RATIOS['Fe'], 'Iron (Fe)'),
        ('FeO', O_FE_RATIOS['FeO'], 'Wustite (Fe_{1-x}O)'),
        ('Fe3O4', O_FE_RATIOS['Fe3O4'], 'Magnetite (Fe3O4)'),
        ('Fe2O3', O_FE_RATIOS['Fe2O3'], 'Hematite (Fe2O3)'),
    ]

    for name, o_fe, label in oxide_labels:
        fig.add_trace(go.Scatter(
            x=[0, 2.0],
            y=[o_fe, o_fe],
            mode='lines',
            line=dict(color=RIST_COLORS['oxide_lines'], width=1, dash='dot'),
            name=label,
            showlegend=False,
            hovertemplate=f'{label}<br>O/Fe = {o_fe:.3f}<extra></extra>',
        ))

        # Label at left edge
        fig.add_annotation(
            x=0.02,
            y=o_fe,
            text=f'<b>{name}</b>',
            showarrow=False,
            xanchor='left',
            yanchor='bottom',
            yshift=2,
            font=dict(size=12, color=RIST_COLORS['oxide_lines']),
        )

    # === Vertical line at O/C = 1 (100% CO) ===
    fig.add_trace(go.Scatter(
        x=[1.0, 1.0],
        y=[0, O_FE_RATIOS['Fe2O3']],
        mode='lines',
        line=dict(color='rgba(100, 100, 100, 0.5)', width=1, dash='dash'),
        name='100% CO',
        showlegend=False,
        hovertemplate='O/C = 1 (100% CO)<extra></extra>',
    ))

    # === Equilibrium curves from Baur-Glassner ===
    if T_kelvin >= WUSTITE_STABILITY_K:
        # Fe/FeO equilibrium
        god_fe_feo = thermo.GOD('FeO_Fe_CO', T_kelvin)
        o_c_fe_feo = god_to_o_c(god_fe_feo)

        # Vertical line at Fe/FeO equilibrium
        fig.add_trace(go.Scatter(
            x=[o_c_fe_feo, o_c_fe_feo],
            y=[O_FE_RATIOS['Fe'], O_FE_RATIOS['FeO']],
            mode='lines',
            line=dict(color=RIST_COLORS['equilibrium'], width=2),
            name='Fe/FeO Equilibrium',
            showlegend=True,
            hovertemplate=f'Fe/FeO Equilibrium<br>O/C = {o_c_fe_feo:.3f}<br>GOD = {god_fe_feo:.3f}<extra></extra>',
        ))

        # FeO/Fe3O4 equilibrium
        god_feo_fe3o4 = thermo.GOD('Fe3O4_FeO_CO', T_kelvin)
        o_c_feo_fe3o4 = god_to_o_c(god_feo_fe3o4)

        fig.add_trace(go.Scatter(
            x=[o_c_feo_fe3o4, o_c_feo_fe3o4],
            y=[O_FE_RATIOS['FeO'], O_FE_RATIOS['Fe3O4']],
            mode='lines',
            line=dict(color=RIST_COLORS['equilibrium'], width=2),
            name='FeO/Fe3O4 Equilibrium',
            showlegend=True,
            hovertemplate=f'FeO/Fe3O4 Equilibrium<br>O/C = {o_c_feo_fe3o4:.3f}<br>GOD = {god_feo_fe3o4:.3f}<extra></extra>',
        ))
    else:
        # Below 570C: Only Fe/Fe3O4 equilibrium
        god_fe_fe3o4 = thermo.GOD('Fe3O4_Fe_CO', T_kelvin)
        o_c_fe_fe3o4 = god_to_o_c(god_fe_fe3o4)

        fig.add_trace(go.Scatter(
            x=[o_c_fe_fe3o4, o_c_fe_fe3o4],
            y=[O_FE_RATIOS['Fe'], O_FE_RATIOS['Fe3O4']],
            mode='lines',
            line=dict(color=RIST_COLORS['equilibrium'], width=2),
            name='Fe/Fe3O4 Equilibrium',
            showlegend=True,
            hovertemplate=f'Fe/Fe3O4 Equilibrium<br>O/C = {o_c_fe_fe3o4:.3f}<br>GOD = {god_fe_fe3o4:.3f}<extra></extra>',
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
            fillcolor=RIST_COLORS['forbidden_zone'],
            line=dict(color='rgba(255, 100, 100, 0.5)', width=1),
            name='Forbidden Zone',
            showlegend=True,
            hovertemplate='Thermodynamically forbidden zone<extra></extra>',
        ))

    # === W point (Wustite) and M point (Magnetite) ===
    if show_equilibrium_points:
        w_point = rist.get_wustite_point(T_kelvin)
        m_point = rist.get_magnetite_point(T_kelvin)

        # W point
        god_w = o_c_to_god(w_point[0])
        god_w_str = f"{god_w:.3f}" if god_w is not None else "N/A"
        fig.add_trace(go.Scatter(
            x=[w_point[0]],
            y=[w_point[1]],
            mode='markers+text',
            marker=dict(size=14, color=RIST_COLORS['w_point'], symbol='circle'),
            text=['W'],
            textposition='top center',
            textfont=dict(size=14, color=RIST_COLORS['w_point'], family='Arial Black'),
            name='Wustite Point (W)',
            showlegend=True,
            hovertemplate=f'Wustite Point (W)<br>O/C = {w_point[0]:.3f}<br>O/Fe = {w_point[1]:.3f}<br>GOD = {god_w_str}<extra></extra>',
        ))

        # M point
        god_m = o_c_to_god(m_point[0])
        god_m_str = f"{god_m:.3f}" if god_m is not None else "N/A"
        fig.add_trace(go.Scatter(
            x=[m_point[0]],
            y=[m_point[1]],
            mode='markers+text',
            marker=dict(size=14, color=RIST_COLORS['m_point'], symbol='circle'),
            text=['M'],
            textposition='top center',
            textfont=dict(size=14, color=RIST_COLORS['m_point'], family='Arial Black'),
            name='Magnetite Point (M)',
            showlegend=True,
            hovertemplate=f'Magnetite Point (M)<br>O/C = {m_point[0]:.3f}<br>O/Fe = {m_point[1]:.3f}<br>GOD = {god_m_str}<extra></extra>',
        ))

    # === Ideal operating line ===
    if show_operating_line:
        start, end, c_per_fe = rist.get_ideal_operating_line(T_kelvin)

        fig.add_trace(go.Scatter(
            x=[start[0], end[0]],
            y=[start[1], end[1]],
            mode='lines+markers',
            line=dict(color=RIST_COLORS['operating_line'], width=3),
            marker=dict(size=8, color=RIST_COLORS['operating_line']),
            name=f'Ideal Operating Line',
            showlegend=True,
            hovertemplate='Ideal Operating Line<br>O/C = %{x:.3f}<br>O/Fe = %{y:.3f}<extra></extra>',
        ))

        # Annotation for C consumption
        mid_x = (start[0] + end[0]) / 2
        mid_y = (start[1] + end[1]) / 2
        fig.add_annotation(
            x=mid_x,
            y=mid_y,
            text=f'C/Fe = {c_per_fe:.2f}',
            showarrow=True,
            arrowhead=0,
            ax=40,
            ay=-30,
            font=dict(size=11, color=RIST_COLORS['operating_line']),
            bgcolor='white',
            bordercolor=RIST_COLORS['operating_line'],
            borderwidth=1,
        )

    # === Real operating line from mass balance ===
    if mass_balance_results is not None:
        mb = mass_balance_results
        start_mb, end_mb = calculate_operating_line_points(mb)

        # Point at O/C = 1 (boundary indirect/direct)
        o_fe_at_oc1 = mb.slope_k * 1.0 + mb.intercept_d

        # === Reduction zones based on RG point ===
        # Use dynamic o_fe_ore from mass balance (1.5 for Fe2O3, 1.333 for Fe3O4)
        o_fe_ore = mb.o_fe_ore
        
        if 0 <= o_fe_at_oc1 <= o_fe_ore:
            # Direct reduction: Zone below RG point (from Fe=0 to RG point)
            # Polygon: Start -> (1, 0) -> (1, o_fe_at_oc1) -> back along operating line
            # Simplified: Triangle from Start to O/C=1
            fig.add_trace(go.Scatter(
                x=[start_mb[0], 1.0, 1.0, start_mb[0]],
                y=[start_mb[1], 0, o_fe_at_oc1, start_mb[1]],
                fill='toself',
                fillcolor=RIST_COLORS['direct_reduction'],
                line=dict(color='rgba(100, 150, 255, 0.4)', width=1),
                name=f'Direct Reduction ({mb.direct_reduction_pct:.0f}%)',
                showlegend=True,
                hovertemplate=f'Direct Reduction<br>{mb.direct_reduction_pct:.1f}%<br>(solid carbon)<extra></extra>',
            ))

            # Indirect reduction: Zone above RG point (from RG point to o_fe_ore)
            # Polygon: (1, o_fe_at_oc1) -> (1, o_fe_ore) -> End -> back along operating line
            fig.add_trace(go.Scatter(
                x=[1.0, 1.0, end_mb[0], 1.0],
                y=[o_fe_at_oc1, o_fe_ore, end_mb[1], o_fe_at_oc1],
                fill='toself',
                fillcolor=RIST_COLORS['indirect_reduction'],
                line=dict(color='rgba(100, 200, 100, 0.4)', width=1),
                name=f'Indirect Reduction ({mb.indirect_reduction_pct:.0f}%)',
                showlegend=True,
                hovertemplate=f'Indirect Reduction<br>{mb.indirect_reduction_pct:.1f}%<br>(CO/CO2 gas)<extra></extra>',
            ))

        # Draw operating line (above the zones)
        fig.add_trace(go.Scatter(
            x=[start_mb[0], end_mb[0]],
            y=[start_mb[1], end_mb[1]],
            mode='lines+markers',
            line=dict(color=RIST_COLORS['operating_line_real'], width=3),
            marker=dict(size=10, color=RIST_COLORS['operating_line_real'], symbol='diamond'),
            name=f'Operating Line (Mass Balance)',
            showlegend=True,
            hovertemplate='Operating Line (Mass Balance)<br>O/C = %{x:.3f}<br>O/Fe = %{y:.3f}<extra></extra>',
        ))

        # Mark RG point
        if 0 <= o_fe_at_oc1 <= o_fe_ore:
            fig.add_trace(go.Scatter(
                x=[1.0],
                y=[o_fe_at_oc1],
                mode='markers',
                marker=dict(size=14, color=RIST_COLORS['operating_line_real'], symbol='star'),
                name='RG Point (O/C=1)',
                showlegend=True,
                hovertemplate=f'Reduction Degree Point<br>O/C = 1.0<br>O/Fe = {o_fe_at_oc1:.3f}<br>Indirect: {mb.indirect_reduction_pct:.1f}%<br>Direct: {mb.direct_reduction_pct:.1f}%<extra></extra>',
            ))

        # Annotation for mass balance C/Fe
        mid_x_mb = (start_mb[0] + end_mb[0]) / 2
        mid_y_mb = (start_mb[1] + end_mb[1]) / 2
        fig.add_annotation(
            x=mid_x_mb,
            y=mid_y_mb,
            text=f'C/Fe = {mb.slope_k:.2f}<br>Ind: {mb.indirect_reduction_pct:.0f}% / Dir: {mb.direct_reduction_pct:.0f}%',
            showarrow=True,
            arrowhead=0,
            ax=-60,
            ay=40,
            font=dict(size=10, color=RIST_COLORS['operating_line_real']),
            bgcolor='white',
            bordercolor=RIST_COLORS['operating_line_real'],
            borderwidth=1,
        )

    # === Layout ===
    fig.update_layout(
        title=dict(
            text=f'<b>Rist Diagram at {T_celsius:.0f}°C</b>',
            x=0.5,
            xanchor='center',
            font=dict(size=16),
        ),
        xaxis=dict(
            title='O/C [-]  (Oxygen per Carbon)',
            range=[-0.1, 2.1],
            dtick=0.25,
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(200, 200, 200, 0.5)',
            zeroline=True,
            zerolinewidth=2,
            zerolinecolor='rgba(0, 0, 0, 0.3)',
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
            x=0.13,
            y=1.0,
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

    # Annotations for X-axis zones
    fig.add_annotation(
        x=0.0,
        y=-0.08,
        text='C',
        showarrow=False,
        font=dict(size=10, color='gray'),
        xref='x',
        yref='y',
    )
    fig.add_annotation(
        x=1.0,
        y=-0.08,
        text='CO',
        showarrow=False,
        font=dict(size=10, color='gray'),
        xref='x',
        yref='y',
    )
    fig.add_annotation(
        x=2.0,
        y=-0.08,
        text='CO₂',
        showarrow=False,
        font=dict(size=10, color='gray'),
        xref='x',
        yref='y',
    )

    return fig


if __name__ == '__main__':
    # Test: Create diagram
    print("Creating Rist Diagram...")

    thermo = BaurGlassnerThermo()

    for T_C in [700, 900, 1000]:
        fig = create_rist_diagram(
            thermo=thermo,
            T_celsius=T_C,
            show_operating_line=True,
            show_forbidden_zone=True,
            show_equilibrium_points=True,
        )
        filename = f'rist_diagram_{T_C}C.html'
        fig.write_html(filename)
        print(f"Diagram saved as '{filename}'")
