"""
Baur-Glaessner Diagram - Interaktive Streamlit App

Zeigt das Phasendiagramm fuer die Systeme Fe-C-O und Fe-H2-O.
"""

import streamlit as st
import pandas as pd
from thermodynamics import BaurGlassnerThermo, SHOMATE_DATA, calc_H, calc_S, calc_G
from plotting import create_baur_glassner_diagram
from rist_plotting import create_rist_diagram
from rist_diagram import RistDiagram, O_FE_RATIOS, o_c_to_god
from mass_balance import (
    OreComposition, CokeComposition, HotBlastComposition, HotMetalComposition,
    FluxComposition, calculate_mass_balance, calculate_ore_amount_from_fe,
    calculate_operating_line_points, calculate_flux_amount
)
# H2 Rist Diagram imports
from rist_diagram_h2 import RistDiagramH2, o_h2_to_god
from rist_plotting_h2 import create_h2_rist_diagram
from mass_balance_h2 import (
    OreCompositionH2, ReductionGasComposition, DRIComposition,
    calculate_h2_mass_balance, calculate_ore_amount_from_fe_h2, calculate_h2_operating_line_points
)

# Seitenkonfiguration
st.set_page_config(
    page_title="Baur-Glaessner Diagram Generator",
    page_icon="",
    layout="wide"
)

# Titel
st.title("Baur-Glaessner Diagram Generator")
st.markdown("*Interactive phase equilibrium diagram for iron oxide reduction - by Maximilian Schaber*")
st.caption("Disclaimer: All information is provided without guarantee. No responsibility is taken for the accuracy, completeness, or timeliness of the content.")


# Thermodynamik laden
@st.cache_resource
def load_thermo():
    return BaurGlassnerThermo()


thermo = load_thermo()

# === Sidebar ===
with st.sidebar:
    st.header("Display Options")

    st.subheader("Gas Systems")
    show_CO = st.checkbox("CO/CO₂ system (solid lines)", value=True)
    show_H2 = st.checkbox("H₂/H₂O system (dashed lines)", value=True)
    show_boudouard = st.checkbox("Boudouard equilibrium (dotted)", value=True)

    st.subheader("Diagram Options")
    show_labels = st.checkbox("Show phase labels", value=True)

    st.divider()

    # === Punkt-Rechner ===
    st.subheader("Add Custom Points")

    # Gas-System Auswahl
    gas_system = st.radio("Gas System", ["CO/CO₂", "H₂/H₂O"], horizontal=True)

    # Eingabefelder
    col_gas, col_temp = st.columns(2)
    with col_gas:
        if gas_system == "CO/CO₂":
            reducing_gas_pct = st.number_input(
                "CO [%]",
                min_value=0.0,
                max_value=100.0,
                value=50.0,
                step=1.0,
                help="Percentage of CO in CO+CO₂ mixture"
            )
        else:
            reducing_gas_pct = st.number_input(
                "H₂ [%]",
                min_value=0.0,
                max_value=100.0,
                value=50.0,
                step=1.0,
                help="Percentage of H₂ in H₂+H₂O mixture"
            )

    with col_temp:
        point_temp = st.number_input(
            "Temperature [°C]",
            min_value=200,
            max_value=1200,
            value=700,
            step=10
        )

    # GOD berechnen
    god_calculated = 1.0 - (reducing_gas_pct / 100.0)

    st.info(f"**Calculated GOD:** {god_calculated:.3f}")

    # Session State für Punkte initialisieren
    if 'custom_points' not in st.session_state:
        st.session_state.custom_points = []

    # Punkt hinzufügen
    col_add, col_clear = st.columns(2)
    with col_add:
        if st.button("Add Point", use_container_width=True):
            point_number = len(st.session_state.custom_points) + 1
            gas_label = "CO" if gas_system == "CO/CO₂" else "H₂"
            st.session_state.custom_points.append({
                'number': point_number,
                'god': god_calculated,
                'temp': point_temp,
                'gas_pct': reducing_gas_pct,
                'gas_type': gas_label
            })
            st.rerun()

    with col_clear:
        if st.button("Clear All", use_container_width=True):
            st.session_state.custom_points = []
            st.rerun()

    # Tabelle mit Punkten anzeigen
    if st.session_state.custom_points:
        st.markdown("**Added Points:**")
        points_df = pd.DataFrame(st.session_state.custom_points)
        points_df.columns = ['#', 'GOD', 'T [°C]', 'Gas [%]', 'Type']
        st.dataframe(points_df, use_container_width=True, hide_index=True)

    st.divider()

    st.subheader("Temperature Range")
    T_min = st.slider("Min Temperature [C]", 200, 600, 300)
    T_max = st.slider("Max Temperature [C]", 800, 1200, 1000)

    st.divider()
    st.caption("**X-Axis:** Gas Oxidation Degree (GOD)")
    st.caption("GOD = CO2/(CO+CO2) or H2O/(H2+H2O)")

# === Punkte für Diagramm vorbereiten ===
custom_points_for_plot = None
if st.session_state.custom_points:
    custom_points_for_plot = [
        (p['god'], p['temp'], p['number'])
        for p in st.session_state.custom_points
    ]

# === Hauptdiagramm ===
fig = create_baur_glassner_diagram(
    thermo=thermo,
    show_CO_system=show_CO,
    show_H2_system=show_H2,
    show_boudouard=show_boudouard,
    show_labels=show_labels,
    T_min_C=T_min,
    T_max_C=T_max,
    custom_points=custom_points_for_plot,
)

st.plotly_chart(fig, use_container_width=True)

# === Legende ===
col1, col2 = st.columns(2)

with col1:
    st.markdown("""
    **Line Types:**
    - **Solid**: CO/CO₂ system
    - **Dashed**: H₂/H₂O system
    - **Dotted**: Boudouard equilibrium (2CO = C + CO₂)
    """)

with col2:
    st.markdown("""
    **Key Features:**
    - **570°C**: FeO stability limit (triple point)
    - **Below 570°C**: Direct Fe ↔ Fe₃O₄ transition
    - **Above 570°C**: FeO (Wustite) is stable
    """)

# === Info-Box ===
with st.expander("About the Baur-Glaessner Diagram"):
    st.markdown("""
    ### Thermodynamic Phase Diagram for Iron Oxide Reduction

    The Baur-Glaessner diagram shows the thermodynamic stability regions of iron
    and its oxides as a function of temperature and gas composition (GOD).

    **Phase Regions:**
    - **Fe**: Metallic iron (most reduced state)
    - **Fe(1-x)O (Wustite)**: Non-stoichiometric Fe(II) oxide, only stable above ~570C
    - **Fe3O4 (Magnetite)**: Mixed Fe(II)/Fe(III) oxide
    - **Fe2O3 (Hematite)**: Fe(III) oxide (far right, GOD → 1)

    **Gas Oxidation Degree (GOD):**
    - GOD = 0: Pure reducing gas (100% CO or H2)
    - GOD = 1: Pure oxidized gas (100% CO2 or H2O)
    - Phase boundaries show equilibrium compositions

    **Triple Point at 570C:**
    All three boundaries (Fe/FeO, FeO/Fe3O4, Fe/Fe3O4) meet at this temperature.
    Below 570C, FeO is thermodynamically unstable and the Fe/Fe3O4 boundary
    becomes vertical.

    **Boudouard Equilibrium:**
    The reaction 2CO = C + CO₂ defines the carbon deposition boundary.
    At low temperatures with high GOD, solid carbon can precipitate.
    """)

# === Boudouard-Gleichgewicht ===
with st.expander("Boudouard Equilibrium Data"):
    st.markdown("### 2CO = C + CO₂")

    boud_data = []
    for T_C in range(400, 1050, 50):
        T_K = T_C + 273.15
        pct_CO = thermo.boudouard_reducing_gas_fraction(T_K) * 100
        GOD = thermo.boudouard_GOD(T_K)

        boud_data.append({
            'T [C]': T_C,
            'GOD': f"{GOD:.3f}",
            '%CO': f"{pct_CO:.1f}",
            '%CO2': f"{100-pct_CO:.1f}",
        })

    df_boud = pd.DataFrame(boud_data)
    st.dataframe(df_boud, use_container_width=True, hide_index=True)

# === Datenquellen ===
with st.expander("Thermodynamic Data Sources"):
    st.markdown("""
    ### Data Sources

    The thermodynamic calculations in this application are based on **NIST-JANAF Thermochemical Tables**:

    **Primary Source:**
    - Chase, M.W., Jr. (1998). *NIST-JANAF Thermochemical Tables, Fourth Edition*.
      Journal of Physical and Chemical Reference Data, Monograph 9, 1-1951.

    **Shomate Equation Coefficients:**
    - All heat capacity, enthalpy, and entropy calculations use Shomate polynomial coefficients
      from the NIST Chemistry WebBook (https://webbook.nist.gov/chemistry/)

    **Species included:**
    - Iron phases: Fe (solid), FeO (Wustite), Fe₃O₄ (Magnetite), Fe₂O₃ (Hematite)
    - Gas species: CO, CO₂, H₂, H₂O, C (Graphite)

    **Equilibrium Constants:**
    - Phase boundary equilibrium constants are derived from the Shomate data using
      standard thermodynamic relationships (ΔG = -RT ln K)
    - Boudouard equilibrium uses empirical correlation: log₁₀(K) = -8916/T + 9.113
    """)

# =============================================================================
# RIST DIAGRAM SECTION
# =============================================================================

st.divider()
st.header("Rist Diagram - Generator")
st.markdown("*Mass balance diagram for iron oxide reduction processes - by Maximilian Schaber*")
st.caption("Disclaimer: All information is provided without guarantee. No responsibility is taken for the accuracy, completeness, or timeliness of the content.")
st.markdown("""
The **Rist Diagram** represents the mass balance of reduction processes.
The axes show the ratio O/Fe (oxygen per iron) and O/C (oxygen per carbon in the gas).
""")

# === Tabs fuer Einstellungen und Massenbilanz ===
tab_settings, tab_massenbilanz = st.tabs(["Diagram Settings", "Mass Balance Input"])

with tab_settings:
    col_rist1, col_rist2 = st.columns([1, 2])

    with col_rist1:
        rist_temp = st.number_input(
            "Temperature [°C]",
            min_value=400,
            max_value=1200,
            value=900,
            step=50,
            help="Temperature for calculating equilibrium points W and M"
        )

        show_rist_operating_line = False  # Disabled - no longer needed
        show_rist_forbidden = st.checkbox("Show forbidden zone", value=True)
        show_massenbilanz_line = st.checkbox("Show operating line from mass balance", value=True)

    with col_rist2:
        rist = RistDiagram(thermo)
        T_K = rist_temp + 273.15

        w_point = rist.get_wustite_point(T_K)
        m_point = rist.get_magnetite_point(T_K)
        _, _, c_per_fe = rist.get_ideal_operating_line(T_K)

        god_w = o_c_to_god(w_point[0])
        god_m = o_c_to_god(m_point[0])
        god_w_str = f"{god_w:.3f}" if god_w is not None else "N/A"
        god_m_str = f"{god_m:.3f}" if god_m is not None else "N/A"

        st.info(f"""
        **Equilibrium Points at {rist_temp}°C:**
        - **Wustite Point (W):** O/C = {w_point[0]:.3f}, O/Fe = {w_point[1]:.3f} (GOD = {god_w_str})
        - **Magnetite Point (M):** O/C = {m_point[0]:.3f}, O/Fe = {m_point[1]:.3f} (GOD = {god_m_str})
        """)

with tab_massenbilanz:
    st.markdown("### Input Materials")

    # === Iron Ore ===
    st.markdown("**Iron Ore**")
    col_ore1, col_ore2, col_ore3, col_ore4, col_ore5 = st.columns(5)
    with col_ore1:
        ore_fe2o3 = st.number_input("Fe₂O₃ [%]", min_value=0.0, max_value=100.0, value=94.4, step=0.1, key="ore_fe2o3",
                                     help="Hematite content")
    with col_ore2:
        ore_fe3o4 = st.number_input("Fe₃O₄ [%]", min_value=0.0, max_value=100.0, value=0.0, step=0.1, key="ore_fe3o4",
                                     help="Magnetite content")
    with col_ore3:
        ore_sio2 = st.number_input("SiO₂ [%]", min_value=0.0, max_value=100.0, value=2.5, step=0.1, key="ore_sio2")
    with col_ore4:
        ore_cao = st.number_input("CaO [%]", min_value=0.0, max_value=100.0, value=3.1, step=0.1, key="ore_cao")
    with col_ore5:
        ore_sum = ore_fe2o3 + ore_fe3o4 + ore_sio2 + ore_cao
        st.metric("Sum", f"{ore_sum:.1f}%")

    # === Coke ===
    st.markdown("**Coke**")
    col_coke1, col_coke2, col_coke3, col_coke4 = st.columns(4)
    with col_coke1:
        coke_c = st.number_input("C [%]", min_value=0.0, max_value=100.0, value=87.0, step=0.1, key="coke_c")
    with col_coke2:
        coke_sio2 = st.number_input("SiO2 [%]", min_value=0.0, max_value=100.0, value=13.0, step=0.1, key="coke_sio2")
    with col_coke3:
        coke_amount = st.number_input("Amount [kg/tHM]", min_value=100.0, max_value=1000.0, value=500.0, step=10.0, key="coke_amount")
    with col_coke4:
        coke_sum = coke_c + coke_sio2
        st.metric("Sum", f"{coke_sum:.1f}%")

    # === Hot Blast ===
    st.markdown("**Hot Blast**")
    col_wind1, col_wind2, col_wind3, col_wind4 = st.columns(4)
    with col_wind1:
        wind_o2 = st.number_input("O2 [Vol-%]", min_value=0.0, max_value=100.0, value=21.0, step=0.5, key="wind_o2")
    with col_wind2:
        wind_n2 = st.number_input("N2 [Vol-%]", min_value=0.0, max_value=100.0, value=79.0, step=0.5, key="wind_n2")
    with col_wind3:
        wind_amount = st.number_input("Amount [Nm³/tHM]", min_value=500.0, max_value=2000.0, value=1200.0, step=50.0, key="wind_amount")
    with col_wind4:
        wind_sum = wind_o2 + wind_n2
        st.metric("Sum", f"{wind_sum:.1f}%")

    # === Flux ===
    st.markdown("**Flux (Limestone)**")
    col_flux1, col_flux2, col_flux3, col_flux4 = st.columns(4)
    with col_flux1:
        flux_caco3 = st.number_input("CaCO3 [%]", min_value=0.0, max_value=100.0, value=100.0, step=1.0, key="flux_caco3")
    with col_flux2:
        flux_cao = st.number_input("CaO [%]", min_value=0.0, max_value=100.0, value=0.0, step=1.0, key="flux_cao")
    with col_flux3:
        slag_basicity = st.number_input("Basicity B2", min_value=0.5, max_value=2.0, value=1.0, step=0.1, key="slag_basicity",
                                        help="Target basicity B2 = CaO/SiO2 in the slag")
    with col_flux4:
        flux_sum = flux_caco3 + flux_cao
        st.metric("Sum", f"{flux_sum:.1f}%")

    st.markdown("### Products")

    # === Hot Metal ===
    st.markdown("**Hot Metal**")
    col_hm1, col_hm2, col_hm3, col_hm4 = st.columns(4)
    with col_hm1:
        hm_fe = st.number_input("Fe [%]", min_value=90.0, max_value=100.0, value=95.7, step=0.1, key="hm_fe")
    with col_hm2:
        hm_si = st.number_input("Si [%]", min_value=0.0, max_value=5.0, value=0.3, step=0.1, key="hm_si")
    with col_hm3:
        hm_c = st.number_input("C [%]", min_value=0.0, max_value=10.0, value=4.0, step=0.1, key="hm_c")
    with col_hm4:
        hm_sum = hm_fe + hm_si + hm_c
        st.metric("Sum", f"{hm_sum:.1f}%")

# === Massenbilanz berechnen ===
ore = OreComposition(Fe2O3=ore_fe2o3, Fe3O4=ore_fe3o4, SiO2=ore_sio2, CaO=ore_cao)
coke = CokeComposition(C=coke_c, SiO2=coke_sio2)
hot_blast = HotBlastComposition(O2=wind_o2, N2=wind_n2)
hot_metal = HotMetalComposition(Fe=hm_fe, Si=hm_si, C=hm_c)
flux = FluxComposition(CaCO3=flux_caco3, CaO=flux_cao)

# Calculate ore amount from Fe balance
ore_amount = calculate_ore_amount_from_fe(ore, hot_metal)

# Calculate mass balance (flux amount is calculated automatically from basicity)
mass_balance_results = None
try:
    mass_balance_results = calculate_mass_balance(
        ore=ore,
        ore_amount=ore_amount,
        coke=coke,
        coke_amount=coke_amount,
        hot_blast=hot_blast,
        hot_blast_amount=wind_amount,
        hot_metal=hot_metal,
        slag_basicity_B2=slag_basicity,
        flux=flux,
        flux_amount=None,  # Automatisch berechnen
    )
except Exception as e:
    st.error(f"Error in mass balance calculation: {e}")

# === Display results ===
if mass_balance_results is not None:
    st.markdown("### Calculated Values")
    col_res1, col_res2, col_res3, col_res4 = st.columns(4)

    with col_res1:
        st.metric("Ore Amount", f"{ore_amount:.1f} kg/tHM")
        st.metric("Flux Amount", f"{mass_balance_results.flux_amount:.1f} kg/tHM")

    with col_res2:
        st.metric("Slag Amount", f"{mass_balance_results.slag_amount:.1f} kg/tHM")
        st.metric("Basicity B2", f"{mass_balance_results.slag_basicity_B2:.2f}")

    with col_res3:
        st.metric("C/Fe (Slope k)", f"{mass_balance_results.slope_k:.3f}")
        st.metric("Intercept d", f"{mass_balance_results.intercept_d:.3f}")

    with col_res4:
        st.metric("Indirect Reduction", f"{mass_balance_results.indirect_reduction_pct:.1f}%")
        st.metric("Direct Reduction", f"{mass_balance_results.direct_reduction_pct:.1f}%")

    # Slag details and top gas
    st.info(f"""
    **Slag:** CaO = {mass_balance_results.slag_cao:.1f} kg/tHM, SiO₂ = {mass_balance_results.slag_sio2:.1f} kg/tHM

    **Top Gas (Gichtgas):** {mass_balance_results.gichtgas_volume_nm3:.0f} Nm³/tHM → **{mass_balance_results.gichtgas_co_pct:.1f}% CO, {mass_balance_results.gichtgas_co2_pct:.1f}% CO₂, {mass_balance_results.gichtgas_n2_pct:.1f}% N₂**
    
    **Top Gas GOD:** O/C = {mass_balance_results.o_c_gichtgas:.3f}, GOD = {mass_balance_results.god_gichtgas:.1f}%

    **Line Equation:** O/Fe = {mass_balance_results.slope_k:.3f} × (O/C) + ({mass_balance_results.intercept_d:.3f})
    """)

# === Rist-Diagramm erstellen und anzeigen ===
fig_rist = create_rist_diagram(
    thermo=thermo,
    T_celsius=rist_temp,
    show_operating_line=show_rist_operating_line,
    show_forbidden_zone=show_rist_forbidden,
    show_equilibrium_points=True,
    mass_balance_results=mass_balance_results if show_massenbilanz_line else None,
)

st.plotly_chart(fig_rist, use_container_width=True)

# Explanation of the Rist Diagram
with st.expander("About the Rist Diagram"):
    st.markdown("""
    ### The Rist Diagram

    The Rist Diagram was developed in the 1960s and is an important
    tool for analyzing reduction processes in blast furnaces and direct reduction.

    **Axes:**
    - **Y-Axis (O/Fe):** Oxygen per mole of iron
      - Fe (Iron): O/Fe = 0
      - FeO (Wustite): O/Fe = 1.05 (non-stoichiometric)
      - Fe₃O₄ (Magnetite): O/Fe = 1.33
      - Fe₂O₃ (Hematite): O/Fe = 1.5

    - **X-Axis (O/C):** Oxygen per mole of carbon
      - O/C = 0: Direct reduction (carbon only)
      - O/C = 1: 100% CO
      - O/C = 2: 100% CO₂

    **Operating Line (orange):**
    - The orange line shows the operating line based on mass balance.
    - The slope k = C/Fe represents carbon consumption per mole of iron.
    - The intercept d = -O_external/Fe accounts for external oxygen (blast + Si reduction).
    - The star marker at O/C = 1 shows the split between direct/indirect reduction.

    **Key Points:**
    - **W (Wustite Point):** Fe/FeO equilibrium from the Baur-Glassner diagram.
    - **M (Magnetite Point):** FeO/Fe₃O₄ equilibrium.

    **Reduction Fractions:**
    - **Indirect Reduction:** Reduction by CO gas (O/C ≥ 1)
    - **Direct Reduction:** Reduction by solid carbon (O/C < 1)
    """)

# Table with values at different temperatures
with st.expander("Equilibrium Values at Different Temperatures"):
    rist_data = []
    for T_C in range(600, 1100, 50):
        T_K = T_C + 273.15
        w = rist.get_wustite_point(T_K)
        m = rist.get_magnetite_point(T_K)
        _, _, c_fe = rist.get_ideal_operating_line(T_K)

        god_w = o_c_to_god(w[0])
        god_m = o_c_to_god(m[0])
        rist_data.append({
            'T [C]': T_C,
            'W: O/C': f"{w[0]:.3f}",
            'W: GOD': f"{god_w:.3f}" if god_w is not None else "N/A",
            'W: %CO': f"{(1-god_w)*100:.1f}" if god_w is not None else "N/A",
            'M: O/C': f"{m[0]:.3f}",
            'M: GOD': f"{god_m:.3f}" if god_m is not None else "N/A",
            'C/Fe': f"{c_fe:.3f}",
        })

    df_rist = pd.DataFrame(rist_data)
    st.dataframe(df_rist, use_container_width=True, hide_index=True)

# =============================================================================
# ENDE RIST-DIAGRAMM SEKTION
# =============================================================================

# =============================================================================
# H2 RIST DIAGRAM SECTION
# =============================================================================

st.divider()
st.header("H₂ Rist Diagram - Generator")
st.markdown("*Mass balance diagram for hydrogen direct reduction (H₂-DR) - by Maximilian Schaber*")
st.caption("Disclaimer: All information is provided without guarantee. No responsibility is taken for the accuracy, completeness, or timeliness of the content.")
st.markdown("""
The **H₂ Rist Diagram** represents the mass balance of hydrogen-based reduction processes.
The axes show the ratio O/Fe (oxygen per iron) and O/H₂ (oxygen per hydrogen).
- **X-Axis (O/H₂):** O/H₂ = 1 means 100% H₂, O/H₂ = 2 means 100% H₂O
- **Y-Axis (O/Fe):** Same as the CO Rist diagram (Fe₂O₃ = 1.5, Fe = 0)
""")

# === Tabs for H2 Settings and Mass Balance ===
tab_h2_settings, tab_h2_massenbilanz = st.tabs(["Diagram Settings", "Mass Balance Input"])

with tab_h2_settings:
    col_h2_1, col_h2_2 = st.columns([1, 2])

    with col_h2_1:
        h2_rist_temp = st.number_input(
            "Temperature [°C]",
            min_value=400,
            max_value=1200,
            value=900,
            step=50,
            key="h2_rist_temp",
            help="Temperature for calculating equilibrium points W and M (H₂ system)"
        )

        show_h2_forbidden = st.checkbox("Show forbidden zone", value=True, key="h2_forbidden")
        show_h2_massenbilanz_line = st.checkbox("Show operating line from mass balance", value=True, key="h2_massenbilanz_line")

    with col_h2_2:
        rist_h2 = RistDiagramH2(thermo)
        T_K_h2 = h2_rist_temp + 273.15

        w_point_h2 = rist_h2.get_wustite_point(T_K_h2)
        m_point_h2 = rist_h2.get_magnetite_point(T_K_h2)

        god_w_h2 = o_h2_to_god(w_point_h2[0])
        god_m_h2 = o_h2_to_god(m_point_h2[0])
        god_w_h2_str = f"{god_w_h2:.3f}" if god_w_h2 is not None else "N/A"
        god_m_h2_str = f"{god_m_h2:.3f}" if god_m_h2 is not None else "N/A"

        st.info(f"""
        **Equilibrium Points (H₂ system) at {h2_rist_temp}°C:**
        - **Wustite Point (W):** O/H₂ = {w_point_h2[0]:.3f}, O/Fe = {w_point_h2[1]:.3f} (GOD = {god_w_h2_str})
        - **Magnetite Point (M):** O/H₂ = {m_point_h2[0]:.3f}, O/Fe = {m_point_h2[1]:.3f} (GOD = {god_m_h2_str})
        """)

with tab_h2_massenbilanz:
    st.markdown("### Input Materials")

    # === Iron Ore (same structure as BF) ===
    st.markdown("**Iron Ore**")
    col_h2_ore1, col_h2_ore2, col_h2_ore3, col_h2_ore4, col_h2_ore5 = st.columns(5)
    with col_h2_ore1:
        h2_ore_fe2o3 = st.number_input("Fe₂O₃ [%]", min_value=0.0, max_value=100.0, value=94.4, step=0.1, key="h2_ore_fe2o3",
                                        help="Hematite content")
    with col_h2_ore2:
        h2_ore_fe3o4 = st.number_input("Fe₃O₄ [%]", min_value=0.0, max_value=100.0, value=0.0, step=0.1, key="h2_ore_fe3o4",
                                        help="Magnetite content")
    with col_h2_ore3:
        h2_ore_sio2 = st.number_input("SiO₂ [%]", min_value=0.0, max_value=100.0, value=2.5, step=0.1, key="h2_ore_sio2")
    with col_h2_ore4:
        h2_ore_cao = st.number_input("CaO [%]", min_value=0.0, max_value=100.0, value=3.1, step=0.1, key="h2_ore_cao")
    with col_h2_ore5:
        h2_ore_sum = h2_ore_fe2o3 + h2_ore_fe3o4 + h2_ore_sio2 + h2_ore_cao
        st.metric("Sum", f"{h2_ore_sum:.1f}%")

    # === Reduction Gas ===
    st.markdown("**Reduction Gas (H₂/H₂O)**")
    col_h2_gas1, col_h2_gas2, col_h2_gas3, col_h2_gas4 = st.columns(4)
    with col_h2_gas1:
        h2_gas_h2 = st.number_input("H₂ [Vol-%]", min_value=0.0, max_value=100.0, value=95.0, step=1.0, key="h2_gas_h2")
    with col_h2_gas2:
        h2_gas_h2o = st.number_input("H₂O [Vol-%]", min_value=0.0, max_value=100.0, value=5.0, step=1.0, key="h2_gas_h2o")
    with col_h2_gas3:
        h2_gas_amount = st.number_input("Amount [Nm³/tDRI]", min_value=100.0, max_value=3000.0, value=800.0, step=50.0, key="h2_gas_amount",
                                        help="Total reduction gas input per tonne DRI")
    with col_h2_gas4:
        h2_gas_god = h2_gas_h2o / (h2_gas_h2 + h2_gas_h2o) * 100 if (h2_gas_h2 + h2_gas_h2o) > 0 else 0
        st.metric("Input GOD", f"{h2_gas_god:.1f}%")

    st.markdown("### Product")

    # === DRI/HBI Product ===
    st.markdown("**DRI / HBI (Direct Reduced Iron)**")
    col_h2_dri1, col_h2_dri2, col_h2_dri3, col_h2_dri4 = st.columns(4)
    with col_h2_dri1:
        h2_dri_fe = st.number_input("Fe metallic [%]", min_value=0.0, max_value=100.0, value=92.0, step=0.5, key="h2_dri_fe")
    with col_h2_dri2:
        h2_dri_feo = st.number_input("FeO [%]", min_value=0.0, max_value=20.0, value=3.0, step=0.5, key="h2_dri_feo")
    with col_h2_dri3:
        h2_dri_gangue = st.number_input("Gangue [%]", min_value=0.0, max_value=20.0, value=5.0, step=0.5, key="h2_dri_gangue")
    with col_h2_dri4:
        h2_dri_sum = h2_dri_fe + h2_dri_feo + h2_dri_gangue
        st.metric("Sum", f"{h2_dri_sum:.1f}%")

    # Calculate metallization
    fe_in_feo = h2_dri_feo * 56.0 / 72.0
    fe_total = h2_dri_fe + fe_in_feo
    metallization = (h2_dri_fe / fe_total * 100) if fe_total > 0 else 0
    st.caption(f"**Metallization:** {metallization:.1f}%")

# === H2 Mass Balance Calculation ===
h2_ore = OreCompositionH2(Fe2O3=h2_ore_fe2o3, Fe3O4=h2_ore_fe3o4, SiO2=h2_ore_sio2, CaO=h2_ore_cao)
h2_input_gas = ReductionGasComposition(H2=h2_gas_h2, H2O=h2_gas_h2o)
h2_dri = DRIComposition(Fe_met=h2_dri_fe, FeO=h2_dri_feo, C=0.0, Gangue=h2_dri_gangue)

# Calculate ore amount from Fe balance
h2_ore_amount = calculate_ore_amount_from_fe_h2(h2_ore, h2_dri)

# Calculate H2 mass balance
h2_mass_balance_results = None
try:
    h2_mass_balance_results = calculate_h2_mass_balance(
        ore=h2_ore,
        ore_amount=h2_ore_amount,
        input_gas=h2_input_gas,
        gas_amount_nm3=h2_gas_amount,
        dri=h2_dri,
    )
except Exception as e:
    st.error(f"Error in H₂ mass balance calculation: {e}")

# === Display H2 results ===
if h2_mass_balance_results is not None:
    st.markdown("### Calculated Values")
    col_h2_res1, col_h2_res2, col_h2_res3, col_h2_res4 = st.columns(4)

    with col_h2_res1:
        st.metric("Ore Amount", f"{h2_ore_amount:.1f} kg/tDRI")
        st.metric("H₂ Consumption", f"{h2_mass_balance_results.h2_amount_nm3:.0f} Nm³/tDRI")

    with col_h2_res2:
        st.metric("H₂/Fe (Slope k)", f"{h2_mass_balance_results.slope_k:.3f}")
        st.metric("Gas Utilization", f"{h2_mass_balance_results.gas_utilization:.1f}%")

    with col_h2_res3:
        st.metric("Output GOD", f"{h2_mass_balance_results.god_output*100:.1f}%")
        st.metric("Output O/H₂", f"{h2_mass_balance_results.o_h2_output:.3f}")

    with col_h2_res4:
        st.metric("Indirect Reduction", f"{h2_mass_balance_results.indirect_reduction_pct:.0f}%")
        st.metric("Metallization", f"{metallization:.1f}%")

    # Warning if not enough gas
    if h2_mass_balance_results.gas_utilization > 100:
        st.warning(f"⚠️ Not enough reduction gas! Gas utilization is {h2_mass_balance_results.gas_utilization:.1f}%. Increase gas amount or reduce product metallization.")

    # Calculate output gas composition
    h2_output_pct = (1 - h2_mass_balance_results.god_output) * 100
    h2o_output_pct = h2_mass_balance_results.god_output * 100

    st.info(f"""
    **Gas Input:** {h2_mass_balance_results.gas_input_nm3:.0f} Nm³/tDRI ({h2_gas_h2:.0f}% H₂, {h2_gas_h2o:.0f}% H₂O)
    
    **Gas Output:** GOD = {h2_mass_balance_results.god_output*100:.1f}% → **{h2_output_pct:.1f}% H₂, {h2o_output_pct:.1f}% H₂O**
    
    **Line Equation:** O/Fe = {h2_mass_balance_results.slope_k:.3f} × (O/H₂) + ({h2_mass_balance_results.intercept_d:.3f})
    """)

# === H2 Rist Diagram Plot ===
fig_h2_rist = create_h2_rist_diagram(
    thermo=thermo,
    T_celsius=h2_rist_temp,
    show_forbidden_zone=show_h2_forbidden,
    show_equilibrium_points=True,
    mass_balance_results=h2_mass_balance_results if show_h2_massenbilanz_line else None,
)

st.plotly_chart(fig_h2_rist, use_container_width=True)

# Explanation of the H2 Rist Diagram
with st.expander("About the H₂ Rist Diagram"):
    st.markdown("""
    ### The H₂ Rist Diagram

    The H₂ Rist Diagram is an adaptation of the traditional Rist Diagram for
    **hydrogen-based direct reduction** processes (H₂-DR).

    **Axes:**
    - **Y-Axis (O/Fe):** Oxygen per mole of iron (same as CO Rist)
      - Fe (Iron): O/Fe = 0
      - FeO (Wustite): O/Fe = 1.05
      - Fe₃O₄ (Magnetite): O/Fe = 1.33
      - Fe₂O₃ (Hematite): O/Fe = 1.5

    - **X-Axis (O/H₂):** Oxygen per mole of hydrogen
      - O/H₂ = 1: 100% H₂ (fresh reducing gas)
      - O/H₂ = 2: 100% H₂O (fully oxidized gas)

    **Key Differences from CO Rist:**
    - **No direct reduction:** H₂ reduction is 100% indirect (gas-based)
    - **No carbon input:** Pure H₂ reduction produces no CO₂ emissions
    - **Water as byproduct:** Instead of CO₂, water vapor (H₂O) is produced

    **Operating Line:**
    - Slope k = H₂/Fe represents hydrogen consumption per mole of iron
    - Intercept d = 0 (no external oxygen source in pure H₂ reduction)

    **Advantages of H₂ Reduction:**
    - CO₂-free iron production (green steel)
    - Lower energy consumption at high temperatures
    - Faster reduction kinetics above ~800°C
    """)

# H2 Equilibrium values table
with st.expander("H₂ Equilibrium Values at Different Temperatures"):
    h2_rist_data = []
    for T_C in range(600, 1100, 50):
        T_K = T_C + 273.15
        w_h2 = rist_h2.get_wustite_point(T_K)
        m_h2 = rist_h2.get_magnetite_point(T_K)

        god_w_h2 = o_h2_to_god(w_h2[0])
        god_m_h2 = o_h2_to_god(m_h2[0])
        h2_rist_data.append({
            'T [°C]': T_C,
            'W: O/H₂': f"{w_h2[0]:.3f}",
            'W: GOD': f"{god_w_h2:.3f}" if god_w_h2 is not None else "N/A",
            'W: %H₂': f"{(1-god_w_h2)*100:.1f}" if god_w_h2 is not None else "N/A",
            'M: O/H₂': f"{m_h2[0]:.3f}",
            'M: GOD': f"{god_m_h2:.3f}" if god_m_h2 is not None else "N/A",
        })

    df_h2_rist = pd.DataFrame(h2_rist_data)
    st.dataframe(df_h2_rist, use_container_width=True, hide_index=True)

# =============================================================================
# END H2 RIST DIAGRAM SECTION
# =============================================================================
