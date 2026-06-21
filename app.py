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
#
# ==============================================================================
# ==============================================================================

import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px

# --- GENERALIZED BIOPHYSICAL ENGINE ---
class BiogenesisEngine:
    """Universal MMModel engine using sensitivity coefficients."""
    @staticmethod
    def calc_flux(o2, temp, ph, s_o2, s_temp, s_ph, base_rate):
        # Universal normalization: 0 = Homeostasis, 1 = Lethal Stress
        # Hypoxia: Hill-type inhibition
        lambda_hyp = (10.0**2 + 21.0**2) / (10.0**2 + o2**2)
        
        # Thermal: Arrhenius scaling
        Tk = temp + 273.15
        T0 = 310.15
        arrhenius = np.exp(-(50000/8.314) * (1/Tk - 1/T0))
        
        # pH: Electrochemical potential
        pH0 = 7.4
        pH_mod = 1 + (0.5 * (10**(-ph) / 10**(-pH0))**2) / (0.1**2 + (10**(-ph) / 10**(-pH0))**2)
        
        # Weighted Sensitivity Integration
        flux = base_rate * (lambda_hyp**s_o2) * (arrhenius**s_temp) * (pH_mod**s_ph)
        return max(0, flux)

# --- STREAMLIT FRONTEND ---
st.set_page_config(page_title="EVelution Engine", layout="wide")
st.title("Bioreactor Optimisation")
st.caption("Universal MMModel Biophysical Twin | Author: Evgenia Trudova")

with st.expander("Model Foundation & Biophysical Formulas"):
    st.markdown("This engine utilizes phenotypic sensitivity coefficients to model EV biogenesis across any cell line.")
    col_math1, col_math2 = st.columns(2)
    with col_math1:
        st.markdown("**Hypoxia (Hill Sigmoid)**")
        st.latex(r"\lambda_{hyp} = \frac{K^n + x_0^n}{K^n + x^n}")
        st.markdown("**Thermal (Arrhenius)**")
        st.latex(r"R = r^* \cdot \exp\left[-\frac{E_a}{R}\left(\frac{1}{T_k} - \frac{1}{T_0}\right)\right]")
    with col_math2:
        st.markdown("**pH (Electrochemical)**")
        st.latex(r"R_v = r^* \cdot [1 + \frac{A_0 (H_e/H_0)^n}{K_{pH}^n + (H_e/H_0)^n}]")
        st.markdown("**Functional Value Index**")
        st.latex(r"V_{final} = \text{Yield} \times \text{Purity} \times \text{Consistency}")

st.divider()

# --- INPUTS ---
col_sidebar1, col_sidebar2 = st.columns([1, 2])
with col_sidebar1:
    st.header("Calibration")
    s_o2 = st.slider("Hypoxia Sensitivity", 0.0, 2.0, 1.0, help="Higher = faster collapse under low O2")
    s_temp = st.slider("Thermal Sensitivity", 0.0, 2.0, 1.0)
    s_ph = st.slider("pH Sensitivity", 0.0, 2.0, 1.0)
    st.divider()
    st.file_uploader("Upload Empirical Data (CSV)")
    st.button("Export PDF Report")

with col_sidebar2:
    st.header("Lab Parameters")
    cols = st.columns(3)
    vol = cols[0].number_input("Volume (L)", value=50.0)
    o2 = cols[1].slider("Oxygen (%)", 0.0, 21.0, 15.0)
    temp = cols[2].slider("Temp (°C)", 30.0, 45.0, 38.0)
    
    cols2 = st.columns(3)
    ph = cols2[0].slider("pH", 6.0, 8.0, 7.0)
    mix = cols2[1].slider("Mixing (%)", 50.0, 100.0, 85.0)
    dur = cols2[2].slider("Duration (h)", 12, 72, 48)

# --- SIMULATION ---
engine = BiogenesisEngine()
data = [engine.calc_flux(o2, temp, ph, s_o2, s_temp, s_ph, 1.2e9) for h in range(dur)]
df = pd.DataFrame({"Hour": range(dur), "EV Flux": data})

# --- DASHBOARD ---
metric_col1, metric_col2, metric_col3 = st.columns(3)
metric_col1.metric("Predicted Yield", f"{df['EV Flux'].sum()*vol:.2e}")
metric_col2.metric("Harvest Conc", f"{df['EV Flux'].iloc[-1]:.2e} ev/mL")
metric_col3.metric("Feasibility", "Optimal" if mix < 90 else "Shear Risk")

st.divider()

# Visuals
st.markdown("### Process Accumulation")
fig = px.line(df, x="Hour", y="EV Flux", log_y=True)
st.plotly_chart(fig, use_container_width=True)

with st.expander("Explore Logic"):
    tab_bio, tab_mod = st.tabs(["Biology", "Model"])
    tab_bio.markdown("Dynamics based on cell membrane stability and metabolic stress response.")
    tab_mod.markdown(r"Calculated using non-linear sensitivity coefficients $\sigma_{o2}, \sigma_{temp}, \sigma_{ph}$.")
