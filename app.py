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

# --- 1. THEME OVERRIDE ---
st.markdown("""
<style>
    :root {
        --primary-color: #77DD77;
        --background-color: #0E1117;
        --secondary-background-color: #1E1E2E;
        --text-color: #E0E0E0;
    }
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# --- 2. LAYOUT & COLOR CONFIGURATION ---
fixed_height = 400
fixed_margin = dict(t=40, b=10, l=20, r=20)
C_GREEN = "#77DD77"
C_BLUE = "#779ECB"
C_PURPLE = "#B39EB5"
C_STAR = "#E39777"

# --- CORE BIOPHYSICAL ENGINE ---
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
        flux = base_rate * (lambda_hyp**s_o2) * (thermal_flux**s_temp) * ((gibbs * pH_mod)**s_ph)
        return max(0, flux)

# --- BIOREACTOR SIMULATION MODEL ---
class FedBatchBioreactorModel:
    def __init__(self):
        # Base rate calibrated for MSC metabolic flux
        self.base_rate = 1.2e9

    def run_simulation(self, init_o2, target_temp, target_ph, mixing_homogeneity, duration_hours, s_o2, s_temp, s_ph):
        viability = 100.0
        history = {"Hour": [], "Therapeutic EVs": [], "Stress-Altered EVs": [], "Apoptotic Impurities": [], "Cell Viability (%)": []}
        for hour in range(1, duration_hours + 1):
            rate = BiogenesisEngine.calc_flux(init_o2, target_temp, target_ph, s_o2, s_temp, s_ph, self.base_rate)
            shear = ((mixing_homogeneity - 85.0)**1.5 * 0.1) if mixing_homogeneity > 85 else 0
            viability = max(0, viability - (0.5 + shear))
            thera = rate * (viability / 100)
            stress = rate * ((100 - viability) / 100) * 0.3
            apop = self.base_rate * ((100 - viability) / 100) * 5
            history["Hour"].append(hour)
            history["Therapeutic EVs"].append(thera)
            history["Stress-Altered EVs"].append(stress)
            history["Apoptotic Impurities"].append(apop)
            history["Cell Viability (%)"].append(viability)
        return pd.DataFrame(history)

# --- STREAMLIT FRONTEND UI ---
st.set_page_config(page_title="EVelution Digital Twin", layout="wide")
st.title("Bioreactor Optimisation")
st.caption("Multi-Machinery Model (MMModel) Biophysical Twin | **Cell Line: MSC (Mesenchymal Stem Cell)** | Author: Evgenia Trudova")

# --- EXPANDERS ---
with st.expander("🔬 MSC Model Foundation & Biophysical Formulas"):
    st.markdown("The default simulation parameters are calibrated for **Mesenchymal Stem Cells (MSCs)**. Sensitivity coefficients reflect MSC metabolic flux and membrane vulnerability to sheer stress.")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Hypoxia (Hill Sigmoid)**"); st.latex(r"\lambda_{hyp} = \frac{K^n + x_0^n}{K^n + x^n}")
        st.markdown("**Thermal (Arrhenius)**"); st.latex(r"R = r^* \cdot \exp\left[-\frac{E_a}{R}\left(\frac{1}{T_k} - \frac{1}{T_0}\right)\right] \cdot (1 + A_{stress})")
    with c2:
        st.markdown("**pH (Electrochemical)**"); st.latex(r"R_v = r^* \cdot e^{-\frac{\Delta G}{RT}} \cdot [1 + \frac{A_0 (H_e/H_0)^n}{K_{pH}^n + (H_e/H_0)^n}]")
        st.markdown("**Functional Value Index**"); st.latex(r"V_{final} = \text{Yield} \times \text{Purity} \times \text{Consistency}")

with st.expander("📊 Dashboard Metrics Breakdown"):
    tab_bio, tab_mod = st.tabs(["Biology", "Model"])
    tab_bio.markdown("""
    * **Yield Performance:** Total projected therapeutic EV quantity post-DSP. Calculated by summing the hourly production flux ($\Phi$) over duration, scaled by volume, and adjusted for recovery/loading efficiency.
    * **Harvest Conc:** Concentration of EVs in the bioreactor at harvest. Calculated as the instantaneous flux ($\Phi$) at the final simulation hour.
    * **Downstream Purity:** Efficiency index for purification workflow. Calculated as a fixed recovery coefficient ($\eta_{purity}$).
    * **Cargo Consistency:** Quality index representing therapeutic payload loading. Calculated as a fixed loading coefficient ($\phi_{consistency}$).
    """)
    tab_mod.latex(r"Y_{perf} = \left( \sum_{t=1}^{t_{dur}} \Phi_{thera}(t) \cdot V_{react} \right) \cdot \eta_{purity} \cdot \phi_{consistency}")

# --- SIDEBAR ---
with st.sidebar:
    st.header("1. MSC Sensitivity Calibration")
    s_o2 = st.slider("Hypoxia Sensitivity", 0.0, 2.0, 1.2)
    s_temp = st.slider("Thermal Sensitivity", 0.0, 2.0, 0.8)
    s_ph = st.slider("pH Sensitivity", 0.0, 2.0, 0.9)
    st.header("2. Experimental Parameters")
    vol = st.number_input("Volume (L)", value=50.0)
    o2 = st.slider("Oxygen (%)", 0.0, 21.0, 21.0)
    temp = st.slider("Temp (°C)", 30.0, 45.0, 37.0)
    ph = st.slider("pH", 6.0, 8.0, 7.4)
    mix = st.slider("Mixing (%)", 50.0, 100.0, 85.0)
    dur = st.slider("Duration (h)", 12, 72, 48)
    st.header("3. Target Yield")
    target_yield = st.number_input("Desired Target Yield", value=1e15, format="%.1e")

# --- SIMULATION ---
model = FedBatchBioreactorModel()
df = model.run_simulation(o2, temp, ph, mix, dur, s_o2, s_temp, s_ph)
total_prod = df["Therapeutic EVs"].sum()
true_val = total_prod * vol * 1000 * 0.78 * 0.62
completion = (true_val / target_yield) * 100

# --- METRIC DISPLAY ---
m1, m2, m3, m4 = st.columns(4)
m1.metric("Yield Performance", f"{true_val:.2e}", f"{true_val - target_yield:+.1e}")
m2.metric("Harvest Conc", f"{df['Therapeutic EVs'].iloc[-1]:.2e} ev/mL")
m3.metric("Downstream Purity", "78.0%")
m4.metric("Cargo Consistency", "62.0%")
st.progress(min(completion / 100, 1.0))
st.markdown(f"**Goal Progress: {completion:.1f}%**")

# --- ANALYTICS GRID ---
st.divider()
st.subheader("Analytics")
row1_col1, row1_col2 = st.columns(2)
row2_col1, row2_col2 = st.columns(2)

with row1_col1:
    st.markdown("### Process Accumulation")
    fig_acc = px.line(df, x="Hour", y=["Therapeutic EVs", "Stress-Altered EVs", "Apoptotic Impurities"], log_y=True,
                      color_discrete_map={"Therapeutic EVs": C_GREEN, "Stress-Altered EVs": C_PURPLE, "Apoptotic Impurities": C_BLUE})
    fig_acc.update_layout(height=fixed_height, margin=fixed_margin)
    st.plotly_chart(fig_acc, use_container_width=True)
    with st.expander("Explore Logic"):
        tab_bio, tab_mod = st.tabs(["Biology", "Model"])
        tab_bio.markdown("Tracks the crossover point between therapeutic EVs and necrotic debris. Harvesting past this point drastically increases purification bottlenecks.")
        tab_mod.latex(r"\Phi_{total} = \Phi_{Therapeutic} + \Phi_{Stress} + \Phi_{Apoptotic}")

with row1_col2:
    st.markdown("### Cellular Viability")
    crit = df[df["Cell Viability (%)"] < 50.0]
    fig_via = px.line(df, x="Hour", y="Cell Viability (%)", color_discrete_sequence=[C_BLUE])
    if not crit.empty:
        fig_via.add_trace(go.Scatter(x=[crit.iloc[0]["Hour"]], y=[crit.iloc[0]["Cell Viability (%)"]], mode='markers', marker=dict(size=12, color=C_STAR, symbol='star')))
    fig_via.update_layout(height=fixed_height, margin=fixed_margin, showlegend=False)
    st.plotly_chart(fig_via, use_container_width=True)
    with st.expander("Explore Logic"):
        tab_bio, tab_mod = st.tabs(["Biology", "Model"])
        tab_bio.markdown("The 'Death Cliff': Viability collapses once cumulative stress parameters surpass homeostatic repair capacity, signaling an apoptotic harvest profile.")
        tab_mod.latex(r"\frac{dV}{dt} = -(\kappa_{tox} + \tau_{shear})")

with row2_col1:
    st.markdown("### Yield-to-Value Bridge")
    raw = true_val / (0.78 * 0.62)
    fig_funnel = go.Figure(go.Funnel(y=["Raw Yield", "Intact (Purity)", "Functional"], x=[raw, true_val/0.62, true_val], textinfo="value+percent previous", marker={"color": [C_BLUE, C_PURPLE, C_GREEN]}))
    fig_funnel.update_layout(height=fixed_height, margin=fixed_margin)
    st.plotly_chart(fig_funnel, use_container_width=True)
    with st.expander("Explore Logic"):
        tab_bio, tab_mod = st.tabs(["Biology", "Model"])
        tab_bio.markdown("Visualizes mass balance cascade. Horizontal width represents yield per stage; narrowing intervals illustrate downstream processing losses.")
        tab_mod.latex(r"V_{final} = Yield_{raw} \cdot \eta_{purity} \cdot \phi_{consistency}")

with row2_col2:
    st.markdown("### Yield Sensitivity Analysis")
    dur_rng = range(12, 96, 6)
    sens_data = [model.run_simulation(o2, temp, ph, mix, d, s_o2, s_temp, s_ph)["Therapeutic EVs"].sum() * vol * 1000 * 0.78 * 0.62 for d in dur_rng]
    fig_sens = px.line(x=list(dur_rng), y=sens_data, labels={'x': 'Duration (h)', 'y': 'Total Yield'}, color_discrete_sequence=[C_GREEN])
    idx = np.argmax(sens_data)
    fig_sens.add_trace(go.Scatter(x=[list(dur_rng)[idx]], y=[sens_data[idx]], mode='markers', marker=dict(size=12, color=C_STAR, symbol='star')))
    fig_sens.update_layout(height=fixed_height, margin=fixed_margin, showlegend=False)
    st.plotly_chart(fig_sens, use_container_width=True)
    with st.expander("Explore Logic"):
        tab_bio, tab_mod = st.tabs(["Biology", "Model"])
        tab_bio.markdown("Identifies the metabolic harvest inflection point. Optimal duration is where incremental production is balanced by cell death.")
        tab_mod.latex(r"\frac{d}{dt}Yield(t) = 0 \text{ at } t_{optimal}")
