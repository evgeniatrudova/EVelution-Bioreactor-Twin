# ==============================================================================
# FED-BATCH BIOREACTOR OPTIMISATION FOR EV YIELD PRODUCTION
# Based on the Multi-Machinery Model of Biogenesis
#
# Author / Lead Software Engineer: Evgenia Trudova
# Academic Origin: Swedish University of Agricultural Sciences (SLU), 2026
#
# COPYRIGHT NOTICE & TERMS OF USE
# Copyright (c) 2026 Evgenia Trudova. All Rights Reserved.
# UNAUTHORIZED COMMERCIAL USE STRICTLY PROHIBITED.
# ==============================================================================
import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

# --- 1. THEME & STYLING ---
st.set_page_config(page_title="EVelution Bioreactor Twin", layout="wide")

st.markdown("""
<style>
    /* Blue-centric theme enforcement */
    :root { --primary-color: #779ECB; }
    .stSlider [data-baseweb="slider"] > div > div > div > div { background-color: #779ECB !important; }
    #MainMenu {visibility: hidden;} footer {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# --- 2. GLOBAL CONSTANTS ---
fixed_height = 400
fixed_margin = dict(t=40, b=10, l=20, r=20)
C_GREEN, C_BLUE, C_PURPLE, C_STAR = "#77DD77", "#779ECB", "#B39EB5", "#E39777"

# --- 3. BIOPHYSICAL ENGINE ---
class BiogenesisEngine:
    @staticmethod
    def calc_flux(o2, temp, ph, s_o2, s_temp, s_ph, base_rate):
        T0, R_gas = 310.15, 8.314
        Tk = temp + 273.15
        K, n, x0 = 10.0, 2.0, 21.0
        lambda_hyp = (K**n + x0**n) / (K**n + o2**n)
        Ea = 50000
        arrhenius = np.exp(-(Ea/R_gas) * (1/Tk - 1/T0))
        T_diff = temp - 37.0
        A_stress = (max(0, T_diff)**2) / (1.0**2 + max(0, T_diff)**2)
        thermal_flux = arrhenius * (1 + A_stress)
        pH0 = 7.4
        dG = 2.303 * R_gas * T0 * (ph - pH0)
        gibbs = np.exp(-dG / (R_gas * T0))
        pH_mod = 1 + (0.5 * (10**(-ph) / 10**(-pH0))**2) / (0.1**2 + (10**(-ph) / 10**(-pH0))**2)
        return max(0, base_rate * (lambda_hyp**s_o2) * (thermal_flux**s_temp) * ((gibbs * pH_mod)**s_ph))

class FedBatchBioreactorModel:
    def __init__(self): self.base_rate = 1.2e9
    def run_simulation(self, init_o2, target_temp, target_ph, mixing_homogeneity, duration_hours, s_o2, s_temp, s_ph):
        viability = 100.0
        history = {"Hour": [], "Therapeutic EVs": [], "Stress-Altered EVs": [], "Apoptotic Impurities": [], "Cell Viability (%)": []}
        for hour in range(1, duration_hours + 1):
            rate = BiogenesisEngine.calc_flux(init_o2, target_temp, target_ph, s_o2, s_temp, s_ph, self.base_rate)
            shear = ((mixing_homogeneity - 85.0)**1.5 * 0.1) if mixing_homogeneity > 85 else 0
            viability = max(0, viability - (0.5 + shear))
            history["Hour"].append(hour); history["Therapeutic EVs"].append(rate * (viability / 100))
            history["Stress-Altered EVs"].append(rate * ((100 - viability) / 100) * 0.3)
            history["Apoptotic Impurities"].append(self.base_rate * ((100 - viability) / 100) * 5)
            history["Cell Viability (%)"].append(viability)
        return pd.DataFrame(history)

# --- 4. APP UI ---
st.title("EVelution Bioreactor Optimisation")
st.caption("Multi-Machinery Model (MMModel) | Default Cell Line: MSC | Author: Evgenia Trudova")

# Metric Definition Expander
 with st.expander("Formula Library", expanded=False):
    # Tabbed view keeps the UI clean and avoids 'wall of text'
    tab_bio, tab_math = st.tabs(["Biology (Context)", "Model (Manual Calc)"])

    with tab_bio:
        st.markdown("""
        * **Yield Performance**: Total therapeutic EV quantity projected post-DSP. Calculated by integrating hourly production flux ($\Phi$) over the culture duration, scaled by volume.
        * **Harvest Concentration**: Instantaneous density of EVs in the bioreactor at the moment of harvest ($\Phi_{final}$). Dictates the upstream input load for downstream purification.
        * **Downstream Purity**: Represents the success rate of the purification workflow (TFF/Chromatography). Modeled as a recovery coefficient ($\eta_{purity}$).
        * **Cargo Consistency**: Quality index denoting the percentage of EVs that contain the active therapeutic payload ($\phi_{consistency}$).
        """)

    with tab_math:
        st.markdown("### Student Calculation Guide")
        st.latex(r"1. \text{ Yield Performance: } Y_{perf} = \left( \sum_{t=1}^{T} \Phi(t) \cdot V_{react} \right) \cdot \eta_{purity} \cdot \phi_{consistency}")
        st.latex(r"2. \text{ Harvest Conc: } \Phi_{final} = \Phi(t_{harvest})")
        st.latex(r"3. \text{ Purity: } \eta_{purity} = \frac{EV_{purified}}{EV_{crude}}")
        st.latex(r"4. \text{ Consistency: } \phi_{consistency} = \frac{EV_{loaded}}{EV_{total}}")

# Sidebar
with st.sidebar:
    st.header("Calibration")
    s_o2 = st.slider("Hypoxia Sensitivity", 0.0, 2.0, 1.2)
    s_temp = st.slider("Thermal Sensitivity", 0.0, 2.0, 0.8)
    s_ph = st.slider("pH Sensitivity", 0.0, 2.0, 0.9)
    st.header("Experimental")
    vol = st.number_input("Volume (L)", value=50.0)
    o2, temp, ph = st.slider("Oxygen (%)", 0, 21, 21), st.slider("Temp (°C)", 30, 45, 37), st.slider("pH", 6.0, 8.0, 7.4)
    mix, dur = st.slider("Mixing (%)", 50, 100, 85), st.slider("Duration (h)", 12, 72, 48)
    target = st.number_input("Target Yield", value=1e15, format="%.1e")
    st.divider()
    if st.button("Export to PDF"): st.warning("Export functionality requires FPDF integration.")
    st.file_uploader("Upload History File", type="csv")

# Simulation
model = FedBatchBioreactorModel()
df = model.run_simulation(o2, temp, ph, mix, dur, s_o2, s_temp, s_ph)
true_val = df["Therapeutic EVs"].sum() * vol * 1000 * 0.78 * 0.62

# Metrics
m1, m2, m3, m4 = st.columns(4)
m1.metric("Yield Performance", f"{true_val:.2e}")
m2.metric("Harvest Conc", f"{df['Therapeutic EVs'].iloc[-1]:.2e} ev/mL")
m3.metric("Downstream Purity", "78.0%")
m4.metric("Cargo Consistency", "62.0%")

# --- 5. ANALYTICS GRID ---
st.divider()
st.subheader("Analytics")
r1c1, r1c2 = st.columns(2)
r2c1, r2c2 = st.columns(2)

with r1c1:
    st.markdown("### Process Accumulation")
    fig = px.line(df, x="Hour", y=["Therapeutic EVs", "Stress-Altered EVs", "Apoptotic Impurities"], log_y=True,
                  color_discrete_map={"Therapeutic EVs": C_GREEN, "Stress-Altered EVs": C_PURPLE, "Apoptotic Impurities": C_BLUE})
    fig.update_layout(height=fixed_height, margin=fixed_margin)
    st.plotly_chart(fig, use_container_width=True)
    with st.expander("Explore Logic"):
        t1, t2 = st.tabs(["Biology", "Model"])
        t1.markdown("This graph tracks therapeutic accumulation versus impurity buildup. The harvest window is constrained by the impurity crossover point. Harvesting before this ensures optimal purity levels.")
        t2.latex(r"\Phi_{total} = \Phi_{Therapeutic} + \Phi_{Stress} + \Phi_{Apoptotic}")

with r1c2:
    st.markdown("### Cellular Viability")
    crit = df[df["Cell Viability (%)"] < 50.0]
    fig = px.line(df, x="Hour", y="Cell Viability (%)", color_discrete_sequence=[C_BLUE])
    if not crit.empty: fig.add_trace(go.Scatter(x=[crit.iloc[0]["Hour"]], y=[crit.iloc[0]["Cell Viability (%)"]], mode='markers', marker=dict(size=12, color=C_STAR, symbol='star')))
    fig.update_layout(height=fixed_height, margin=fixed_margin, showlegend=False)
    st.plotly_chart(fig, use_container_width=True)
    with st.expander("Explore Logic"):
        t1, t2 = st.tabs(["Biology", "Model"])
        t1.markdown("Viability decline signifies the transition from growth to apoptotic phase. The 'Death Cliff' marker identifies when cell repair mechanisms collapse under toxic stress.")
        t2.latex(r"\frac{dV}{dt} = -(\kappa_{tox} + \tau_{shear})")

with r2c1:
    st.markdown("### Yield-to-Value Bridge")
    fig = go.Figure(go.Funnel(y=["Raw Yield", "Intact (Purity)", "Functional"], x=[true_val/(0.78*0.62), true_val/0.62, true_val], textinfo="value+percent previous", marker={"color": [C_BLUE, C_PURPLE, C_GREEN]}))
    fig.update_layout(height=fixed_height, margin=fixed_margin)
    st.plotly_chart(fig, use_container_width=True)
    with st.expander("Explore Logic"):
        t1, t2 = st.tabs(["Biology", "Model"])
        t1.markdown("This funnel visualizes the mass balance across downstream stages. Narrowing segments highlight cumulative yield loss during purification. It serves as a diagnostic for loading/recovery efficiency.")
        t2.latex(r"V_{final} = Yield_{raw} \cdot \eta_{purity} \cdot \phi_{consistency}")

with r2c2:
    st.markdown("### Yield Sensitivity Analysis")
    sens = [model.run_simulation(o2, temp, ph, mix, d, s_o2, s_temp, s_ph)["Therapeutic EVs"].sum() * vol * 1000 * 0.78 * 0.62 for d in range(12, 96, 6)]
    fig = px.line(x=range(12, 96, 6), y=sens, color_discrete_sequence=[C_GREEN])
    idx = np.argmax(sens)
    fig.add_trace(go.Scatter(x=[list(range(12, 96, 6))[idx]], y=[sens[idx]], mode='markers', marker=dict(size=12, color=C_STAR, symbol='star')))
    fig.update_layout(height=fixed_height, margin=fixed_margin, showlegend=False)
    st.plotly_chart(fig, use_container_width=True)
    with st.expander("Explore Logic"):
        t1, t2 = st.tabs(["Biology", "Model"])
        t1.markdown("Determines the 'sweet spot' for harvest duration. The inflection point occurs where incremental EV gain is offset by culture necrosis and byproduct toxicity.")
        t2.latex(r"\frac{d}{dt}Yield(t) = 0 \quad at \quad t_{optimal}")
