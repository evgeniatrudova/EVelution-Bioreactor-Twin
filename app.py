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
import plotly.express as px
import plotly.graph_objects as go

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

# --- STREAMLIT FRONTEND ---
st.set_page_config(page_title="EVelution Digital Twin", layout="wide")
st.title("Bioreactor Optimisation")
st.caption("Multi-Machinery Model Bioreactor Twin")

# 1. Global Expander: MMModel
with st.expander("MMModel Formulas"):
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Hypoxia (Hill Sigmoid)**")
        st.latex(r"\lambda_{hyp} = \frac{K^n + x_0^n}{K^n + x^n}")
        st.markdown("**Thermal (Arrhenius)**")
        st.latex(r"R = r^* \cdot \exp\left[-\frac{E_a}{R}\left(\frac{1}{T_k} - \frac{1}{T_0}\right)\right] \cdot (1 + A_{stress})")
    with c2:
        st.markdown("**pH (Electrochemical)**")
        st.latex(r"R_v = r^* \cdot e^{-\frac{\Delta G}{RT}} \cdot [1 + \frac{A_0 (H_e/H_0)^n}{K_{pH}^n + (H_e/H_0)^n}]")

with st.expander("Dashboard Metrics"):
    tab_bio, tab_mod = st.tabs(["Biology", "Model"])
    
    with tab_bio:
        st.markdown("""
        * **Yield Performance:** The total quantity of therapeutic EVs projected at the end of the full process. It is calculated by summing the cumulative production flux over the bioreactor run, scaled by vessel volume, and adjusted for recovery and loading efficiencies.
        * **Harvest Conc:** The concentration of EVs present in the bioreactor at the moment of harvest. It is calculated as the instantaneous production flux observed at the final simulation hour, prior to any purification.
        * **Downstream Purity:** An efficiency index for the purification workflow (e.g., TFF or chromatography). It is calculated as a fixed recovery coefficient representing the fraction of vesicles retained during processing.
        * **Cargo Consistency:** A quality index denoting the percentage of EVs that contain the target therapeutic payload. It is calculated as a fixed coefficient applied to the final yield to account for empty or mis-loaded vesicles.
        """)
    
    with tab_mod:
        st.markdown("### Process Formulas")
        st.latex(r"Y_{perf} = \left( \sum_{t=1}^{t_{dur}} \Phi_{thera}(t) \cdot V_{react} \right) \cdot \eta_{purity} \cdot \phi_{consistency}")
        st.markdown("---")
        st.latex(r"C_{harvest} = \Phi_{thera}(t_{final})")
        st.markdown("---")
        st.latex(r"\text{Where } \eta_{purity} \text{ is the recovery efficiency and } \phi_{consistency} \text{ is the loading factor.}")
        
st.divider()
# --- SIDEBAR ---
with st.sidebar:
    st.header("Cell Line Sensitivity")
    s_o2, s_temp, s_ph = st.slider("Hypoxia", 0.0, 2.0, 1.2), st.slider("Temperature", 0.0, 2.0, 0.8), st.slider("pH", 0.0, 2.0, 0.9)
    st.header("Experimental")
    vol, o2, temp, ph, mix, dur = st.number_input("Volume (L)", 50.0), st.slider("Oxygen (%)", 0.0, 21.0, 21.0), st.slider("Temp (°C)", 30.0, 45.0, 37.0), st.slider("pH", 6.0, 8.0, 7.4), st.slider("Mixing (%)", 50.0, 100.0, 85.0), st.slider("Duration (h)", 12, 72, 48)
    st.header("Yield")
    target_yield = st.number_input("Desired Yield", value=1e15, format="%.1e")

# --- SIMULATION ---
model = FedBatchBioreactorModel()
df = model.run_simulation(o2, temp, ph, mix, dur, s_o2, s_temp, s_ph)
true_val = df["Therapeutic EVs"].sum() * vol * 1000 * 0.78 * 0.62
completion = (true_val / target_yield) * 100

# --- METRIC DISPLAY ---
m1, m2, m3, m4 = st.columns(4)
m1.metric("Yield Performance", f"{true_val:.2e}", f"{true_val - target_yield:+.1e}")
m2.metric("Harvest Conc", f"{df['Therapeutic EVs'].iloc[-1]:.2e} ev/mL")
m3.metric("Downstream Purity", "78.0%")
m4.metric("Cargo Consistency", "62.0%")
st.progress(min(completion / 100, 1.0))
st.write(f"**Goal Progress: {completion:.1f}%**")

# --- DASHBOARD VISUALIZATION GRID ---
st.divider()
st.subheader("Analytics")

# Create two rows of two columns each
row1_col1, row1_col2 = st.columns(2)
row2_col1, row2_col2 = st.columns(2)

# 1. Process Accumulation (Row 1, Left)
with row1_col1:
    st.markdown("### Process ")
    fig_acc = px.line(df, x="Hour", y=["Therapeutic EVs", "Stress-Altered EVs", "Apoptotic Impurities"], log_y=True)
    st.plotly_chart(fig_acc, use_container_width=True)
    with st.expander("Explore Logic"):
        tab_bio, tab_mod = st.tabs(["Biology", "Model"])
        tab_bio.markdown("Tracks therapeutic product and byproduct accumulation over time to identify the optimal harvest point before impurity dominance.")
        tab_mod.latex(r"\Phi_{total} = \Phi_{Therapeutic} + \Phi_{Stress} + \Phi_{Apoptotic}")

# 2. Cellular Viability (Row 1, Right)
with row1_col2:
    st.markdown("### Cellular Viability")
    fig_via = px.line(df, x="Hour", y="Cell Viability (%)")
    st.plotly_chart(fig_via, use_container_width=True)
    with st.expander("Explore Logic"):
        tab_bio, tab_mod = st.tabs(["Biology", "Model"])
        tab_bio.markdown("Viability declines due to cumulative metabolic stress and mechanical impeller shear, identifying the culture 'Death Cliff'.")
        tab_mod.latex(r"\frac{dV}{dt} = -(\kappa_{tox} + \tau_{shear})")


# Define common layout parameters
fixed_height = 400
fixed_margin = dict(t=30, b=0, l=10, r=10)

# 3. Yield-to-Value Bridge (Row 2, Left)
with row2_col1:
    st.markdown("### Yield-to-Value Bridge")
    raw_yield = true_val / (0.78 * 0.62)
    purity_yield = true_val / 0.62
    
    fig_funnel = go.Figure(go.Funnel(
        y=["Raw Target Yield", "Intact EVs (Purity)", "Functional Value"],
        x=[raw_yield, purity_yield, true_val],
        textinfo="value+percent previous",
        marker={"color": ["#636EFA", "#EF553B", "#00CC96"]}
    ))
    
    # Enforce standard dimensions
    fig_funnel.update_layout(height=fixed_height, margin=fixed_margin)
    st.plotly_chart(fig_funnel, use_container_width=True)
    
    with st.expander("Explore Logic"):
        tab_bio, tab_mod = st.tabs(["Biology", "Model"])
        tab_bio.markdown("Visualizes the loss cascade from crude bioreactor harvest to the final functional therapeutic product.")
        tab_mod.latex(r"V_{final} = Yield_{raw} \cdot \eta_{purity} \cdot \phi_{consistency}")

# 4. Yield Sensitivity Analysis (Row 2, Right)
with row2_col2:
    st.markdown("### Yield Sensitivity")
    dur_range = range(12, 96, 6)
    sens_data = [
        model.run_simulation(o2, temp, ph, mix, d, s_o2, s_temp, s_ph)["Therapeutic EVs"].sum() * vol * 1000 * 0.78 * 0.62 
        for d in dur_range
    ]
    
    fig_sens = px.line(x=list(dur_range), y=sens_data, labels={'x': 'Duration (h)', 'y': 'Total Yield'})
    
    # Enforce standard dimensions to match the funnel
    fig_sens.update_layout(height=fixed_height, margin=fixed_margin)
    st.plotly_chart(fig_sens, use_container_width=True)
    
    with st.expander("Explore Logic"):
        tab_bio, tab_mod = st.tabs(["Biology", "Model"])
        tab_bio.markdown("Maps the batch duration 'sweet spot' where yield is maximized before cell necrosis dominates.")
        tab_mod.latex(r"Yield_{target} = \int_{0}^{t} \Phi(t) dt")
