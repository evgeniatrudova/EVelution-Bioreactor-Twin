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

# --- CORE BIOPHYSICAL ENGINE (Thesis Formulas) ---
class BiogenesisEngine:
    """Implements the MMModel differential equations."""
    
    @staticmethod
    def calc_flux(o2, temp, ph, s_o2, s_temp, s_ph, base_rate):
        # Constants
        T0, R_gas = 310.15, 8.314
        Tk = temp + 273.15
        
        # 1. Hypoxia (Hill-type Sigmoid)
        K, n, x0 = 10.0, 2.0, 21.0
        lambda_hyp = (K**n + x0**n) / (K**n + o2**n)
        
        # 2. Thermal (Arrhenius-inspired)
        Ea = 50000
        arrhenius = np.exp(-(Ea/R_gas) * (1/Tk - 1/T0))
        T_diff = temp - 37.0
        A_stress = (max(0, T_diff)**2) / (1.0**2 + max(0, T_diff)**2)
        thermal_flux = arrhenius * (1 + A_stress)
        
        # 3. pH (Gibbs electrochemical shift)
        pH0 = 7.4
        dG = 2.303 * R_gas * T0 * (ph - pH0)
        gibbs = np.exp(-dG / (R_gas * T0))
        pH_mod = 1 + (0.5 * (10**(-ph) / 10**(-pH0))**2) / (0.1**2 + (10**(-ph) / 10**(-pH0))**2)
        
        # Integration
        flux = base_rate * (lambda_hyp**s_o2) * (thermal_flux**s_temp) * ((gibbs * pH_mod)**s_ph)
        return max(0, flux)

# --- BIOREACTOR SIMULATION MODEL ---
class FedBatchBioreactorModel:
    def __init__(self):
        self.base_rate = 1.2e9

    def run_simulation(self, init_o2, target_temp, target_ph, mixing_homogeneity, duration_hours, s_o2, s_temp, s_ph):
        viability = 100.0
        history = {"Hour": [], "Therapeutic EVs": [], "Apoptotic Impurities": [], "Cell Viability (%)": []}
        
        for hour in range(1, duration_hours + 1):
            rate = BiogenesisEngine.calc_flux(init_o2, target_temp, target_ph, s_o2, s_temp, s_ph, self.base_rate)
            shear = ((mixing_homogeneity - 85.0)**1.5 * 0.1) if mixing_homogeneity > 85 else 0
            viability = max(0, viability - (0.5 + shear))
            
            history["Hour"].append(hour)
            history["Therapeutic EVs"].append(rate * (viability/100))
            history["Apoptotic Impurities"].append(self.base_rate * (1 - viability/100) * 10)
            history["Cell Viability (%)"].append(viability)
            
        return pd.DataFrame(history)

# --- STREAMLIT FRONTEND UI ---
st.set_page_config(page_title="EVelution Engine", layout="wide")

st.title("Bioreactor Optimisation")
st.caption("Default Configuration: Mesenchymal Stromal Cells | MMModel ")

# MMModel Explanation
with st.expander("Model Foundation & Biophysical Formulas"):
    st.markdown("This engine utilizes phenotypic sensitivity coefficients to model EV biogenesis across any cell line.")
    col_math1, col_math2 = st.columns(2)
    with col_math1:
        st.markdown("**Hypoxia (Hill Sigmoid)**")
        st.latex(r"\lambda_{hyp} = \frac{K^n + x_0^n}{K^n + x^n}")
        st.markdown("**Thermal (Arrhenius)**")
        st.latex(r"R = r^* \cdot \exp\left[-\frac{E_a}{R}\left(\frac{1}{T_k} - \frac{1}{T_0}\right)\right] \cdot (1 + A_{stress})")
    with col_math2:
        st.markdown("**pH (Electrochemical)**")
        st.latex(r"R_v = r^* \cdot e^{-\frac{\Delta G}{RT}} \cdot [1 + \frac{A_0 (H_e/H_0)^n}{K_{pH}^n + (H_e/H_0)^n}]")
        st.markdown("**Functional Value Index**")
        st.latex(r"V_{final} = \text{Yield} \times \text{Purity} \times \text{Consistency}")

st.divider()

# --- SIDEBAR: ORGANIZED HIERARCHY WITH MSC DEFAULTS ---
with st.sidebar:
    st.header("Cell Line")
    s_o2 = st.slider("Hypoxia Sensitivity", 0.0, 2.0, 1.2, help="MSC baseline sensitivity to O2 tension")
    s_temp = st.slider("Thermal Sensitivity", 0.0, 2.0, 0.8)
    s_ph = st.slider("pH Sensitivity", 0.0, 2.0, 0.9)
    
    st.divider()
    
    st.header("Experimental")
    vol = st.number_input("Volume (L)", value=50.0)
    o2 = st.slider("Oxygen (%)", 0.0, 21.0, 21.0) # Standard Normoxia
    temp = st.slider("Temp (°C)", 30.0, 45.0, 37.0) # MSC Optimal
    ph = st.slider("pH", 6.0, 8.0, 7.4) # MSC Optimal
    mix = st.slider("Mixing (%)", 50.0, 100.0, 85.0)
    dur = st.slider("Duration (h)", 12, 72, 48)
    
    st.divider()
    
    st.header("File Managment")
    target_clinical_yield = st.number_input("Desired Yield (EVs)", value=1e15, format="%.1e")
    uploaded_file = st.file_uploader("Upload CSV", type=["csv"])
    if st.button("Export PDF Report"):
        st.info("Generating PDF report...")

# --- SIMULATION & DASHBOARD ---
model = FedBatchBioreactorModel()
df = model.run_simulation(o2, temp, ph, mix, dur, s_o2, s_temp, s_ph)

col1, col2, col3 = st.columns(3)
col1.metric("Predicted Yield", f"{df['Therapeutic EVs'].sum()*vol:.2e}")
col2.metric("Harvest Concentration", f"{df['Therapeutic EVs'].iloc[-1]:.2e} ev/mL")
col3.metric("Viability", f"{df['Cell Viability (%)'].iloc[-1]:.1f}%")

st.divider()

col_left, col_right = st.columns(2)
with col_left:
    st.markdown("### Process Accumulation")
    fig = px.line(df, x="Hour", y=["Therapeutic EVs", "Apoptotic Impurities"], log_y=True)
    st.plotly_chart(fig, use_container_width=True)
with col_right:
    st.markdown("### Cellular Viability")
    fig2 = px.line(df, x="Hour", y="Cell Viability (%)")
    st.plotly_chart(fig2, use_container_width=True)
