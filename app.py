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
    @staticmethod
    def calc_flux(o2, temp, ph, s_o2, s_temp, s_ph, base_rate):
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
        history = {"Hour": [], "Therapeutic EVs": [], "Stress-Altered EVs": [], "Apoptotic Impurities": [], "Cell Viability (%)": []}
        
        for hour in range(1, duration_hours + 1):
            # Calculate flux using Thesis-derived Engine
            rate = BiogenesisEngine.calc_flux(init_o2, target_temp, target_ph, s_o2, s_temp, s_ph, self.base_rate)
            
            # Damage calculation
            shear = ((mixing_homogeneity - 85.0)**1.5 * 0.1) if mixing_homogeneity > 85 else 0
            viability = max(0, viability - (0.5 + shear))
            
            # Partitioning fluxes
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
st.caption("Multi-Machinery Model (MMModel) Biophysical Twin | Author: Evgenia Trudova")

# The MMModel Dropdown
with st.expander("Model Foundation & Biophysical Formulas"):
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

# --- SIDEBAR: REORGANIZED HIERARCHY ---
with st.sidebar:
    st.header("Cell Line Sensitivity")
    s_o2 = st.slider("Hypoxia", 0.0, 2.0, 1.2)
    s_temp = st.slider("Temperature", 0.0, 2.0, 0.8)
    s_ph = st.slider("pH", 0.0, 2.0, 0.9)
    
    st.divider()
    
    st.header("Experimental")
    vol = st.number_input("Volume (L)", value=50.0)
    o2 = st.slider("Oxygen (%)", 0.0, 21.0, 21.0)
    temp = st.slider("Temp (°C)", 30.0, 45.0, 37.0)
    ph = st.slider("pH", 6.0, 8.0, 7.4)
    mix = st.slider("Mixing (%)", 50.0, 100.0, 85.0)
    dur = st.slider("Duration (h)", 12, 72, 48)
    
    st.divider()
    
    st.header("Desired Yield")
    target_clinical_yield = st.number_input("Desired Target Yield (EVs)", value=1e15, step=1e14, format="%.1e")
    
    st.divider()
    
    st.header("Data Management")
    st.file_uploader("Upload Run Data (CSV)", type=["csv"])
    if st.button("Export to PDF Report"):
        st.info("Generating report...")

# --- SIMULATION & CALCULATIONS ---
model = FedBatchBioreactorModel()
df = model.run_simulation(o2, temp, ph, mix, dur, s_o2, s_temp, s_ph)

# Final Metrics
final_thera = df["Therapeutic EVs"].iloc[-1] * vol * 1000
true_val = final_thera * 0.78 * 0.62
completion_pct = (true_val / target_clinical_yield) * 100
goal_delta = true_val - target_clinical_yield

# --- DASHBOARD RENDER ---
col1, col2, col3, col4 = st.columns(4)
col1.metric("Yield Performance", f"{true_val:.2e}", f"{goal_delta:+.1e} vs Goal")
col2.metric("Harvest Conc", f"{df['Therapeutic EVs'].iloc[-1]:.2e} ev/mL")
col3.metric("Downstream Purity", "78.0%")
col4.metric("Cargo Consistency", "62.0%")

# Progress Bar & Status Banner
st.markdown(f"**Goal Progress: {completion_pct:.1f}%**")
st.progress(min(completion_pct / 100, 1.0))

if completion_pct >= 100:
    st.success(f"Success: Configuration meets clinical yield target of {target_clinical_yield:.1e} EVs.")
else:
    st.warning(f"Configuration shortfall: {abs(goal_delta):.1e} EVs below target. Adjust temperature or duration to optimize yield.")

st.divider()

# --- GRAPHS ---
col_left, col_right = st.columns(2)
with col_left:
    st.markdown("### Process Accumulation")
    fig = px.line(df, x="Hour", y=["Therapeutic EVs", "Stress-Altered EVs", "Apoptotic Impurities"], log_y=True)
    st.plotly_chart(fig, use_container_width=True)
    
    with st.expander("Explore Logic"):
        tab_bio, tab_mod = st.tabs(["Biology", "Model"])
        tab_bio.markdown("Therapeutic yield decays as metabolic stress increases. Apoptotic debris eventually dominates the harvest profile.")
        tab_mod.markdown(r"$\Phi_{total} = \Phi_{Therapeutic} + \Phi_{Stress} + \Phi_{Apoptotic}$")

with col_right:
    st.markdown("### Cellular Viability")
    fig2 = px.line(df, x="Hour", y="Cell Viability (%)")
    st.plotly_chart(fig2, use_container_width=True)
    
    with st.expander("Explore Membrane Degradation"):
        tab_bio, tab_mod = st.tabs(["Biology", "Model"])
        tab_bio.markdown("Cumulative necrosis results from toxic byproduct accumulation and physical impeller stress. The threshold represents the 'Death Cliff'.")
        tab_mod.markdown(r"$\frac{dV}{dt} = -(\kappa_{tox} + \tau_{shear})$")
