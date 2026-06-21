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

# --- CORE BIOPHYSICAL ENGINE ---
class BiogenesisEngine:
    """Implements the MMModel differential equations."""
    
    @staticmethod
    def calc_hypoxia_lambda(o2):
        K, n, x0 = 10.0, 2.0, 21.0
        return (K**n + x0**n) / (K**n + o2**n)

    @staticmethod
    def calc_thermal_flux(temp):
        T0 = 310.15 # 37C in Kelvin
        Tk = temp + 273.15
        Ea, R_gas = 50000, 8.314
        arrhenius = np.exp(-(Ea/R_gas) * (1/Tk - 1/T0))
        # Stress function (A_stress)
        T_diff = temp - 37.0
        A_stress = (max(0, T_diff)**2) / (1.0**2 + max(0, T_diff)**2)
        return arrhenius * (1 + A_stress)

    @staticmethod
    def calc_ph_flux(ph):
        pH0 = 7.4
        R_gas, T0 = 8.314, 310.15
        dG = 2.303 * R_gas * T0 * (ph - pH0)
        gibbs = np.exp(-dG / (R_gas * T0))
        # Electrochemical modulator
        pH_mod = 1 + (0.5 * (10**(-ph) / 10**(-pH0))**2) / (0.1**2 + (10**(-ph) / 10**(-pH0))**2)
        return gibbs * pH_mod

    @staticmethod
    def get_ev_rate(o2, temp, ph, base_rate, s_o2, s_temp, s_ph):
        # Thesis-derived integrative rate formula
        l_hyp = BiogenesisEngine.calc_hypoxia_lambda(o2)
        t_flux = BiogenesisEngine.calc_thermal_flux(temp)
        p_flux = BiogenesisEngine.calc_ph_flux(ph)
        
        # Weighted sensitivity integration
        rate = base_rate * (l_hyp**s_o2) * (t_flux**s_temp) * (p_flux**s_ph)
        return max(0, rate)

# --- BIOREACTOR SIMULATION ---
class FedBatchBioreactorModel:
    def __init__(self):
        self.base_rate = 1.2e9

    def run_simulation(self, init_o2, target_temp, target_ph, mixing_homogeneity, duration_hours, s_o2, s_temp, s_ph):
        history = {"Hour": [], "Therapeutic EVs": [], "Apoptotic Impurities": [], "Cell Viability (%)": []}
        viability = 100.0
        
        for hour in range(1, duration_hours + 1):
            # Calculate flux using Thesis formulas
            rate = BiogenesisEngine.get_ev_rate(init_o2, target_temp, target_ph, self.base_rate, s_o2, s_temp, s_ph)
            
            # Damage logic (simplified for demonstration)
            shear = ((mixing_homogeneity - 85.0)**1.5 * 0.1) if mixing_homogeneity > 85 else 0
            viability = max(0, viability - (0.5 + shear))
            
            history["Hour"].append(hour)
            history["Therapeutic EVs"].append(rate * (viability/100))
            history["Apoptotic Impurities"].append(self.base_rate * (1 - viability/100) * 10)
            history["Cell Viability (%)"].append(viability)
            
        return pd.DataFrame(history)

# --- STREAMLIT FRONTEND ---
st.set_page_config(page_title="EVelution Engine", layout="wide")
st.title("Bioreactor Optimisation")
st.caption("Multi-Machinery Model (MMModel) Biophysical Twin | Author: Evgenia Trudova")

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

# --- SIDEBAR ---
with st.sidebar:
    st.header("Calibration")
    s_o2 = st.slider("Hypoxia Sensitivity", 0.0, 2.0, 1.0)
    s_temp = st.slider("Thermal Sensitivity", 0.0, 2.0, 1.0)
    s_ph = st.slider("pH Sensitivity", 0.0, 2.0, 1.0)
    st.divider()
    st.header("Lab Parameters")
    vol = st.number_input("Volume (L)", value=50.0)
    o2 = st.slider("Oxygen (%)", 0.0, 21.0, 15.0)
    temp = st.slider("Temp (°C)", 30.0, 45.0, 38.0)
    ph = st.slider("pH", 6.0, 8.0, 7.0)
    mix = st.slider("Mixing (%)", 50.0, 100.0, 85.0)
    dur = st.slider("Duration (h)", 12, 72, 48)

# --- RUN MODEL ---
model = FedBatchBioreactorModel()
df = model.run_simulation(o2, temp, ph, mix, dur, s_o2, s_temp, s_ph)

# --- DASHBOARD ---
col1, col2, col3 = st.columns(3)
col1.metric("Predicted Yield", f"{df['Therapeutic EVs'].sum()*vol:.2e}")
col2.metric("Harvest Conc", f"{df['Therapeutic EVs'].iloc[-1]:.2e} ev/mL")
col3.metric("Viability", f"{df['Cell Viability (%)'].iloc[-1]:.1f}%")

st.divider()

col_left, col_right = st.columns(2)
with col_left:
    st.markdown("### Accumulation of Particles")
    fig = px.line(df, x="Hour", y=["Therapeutic EVs", "Apoptotic Impurities"], log_y=True)
    st.plotly_chart(fig, use_container_width=True)
with col_right:
    st.markdown("### Cellular Viability")
    fig2 = px.line(df, x="Hour", y="Cell Viability (%)")
    st.plotly_chart(fig2, use_container_width=True)
