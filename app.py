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
# CORE MATHEMATICAL FORMULAS (The Multi-Machinery Model)
# ------------------------------------------------------------------------------
# 1. ACUTE STRESS BURSTS (Vesicle release triggered by environmental panic)
#    - Hypoxia Burst = ((Opt_O2 - Hypoxia_Thresh) - Current_O2)^2 * Hypoxia_Mult
#    - Thermal Burst = (|Opt_Temp - Current_Temp|^1.5) * Thermal_Mult
#    - pH Burst      = ((Opt_pH - Current_pH)^3) * pH_Mult
#
# 2. HOURLY CELLULAR DAMAGE (Cell death from toxic environments)
#    - Anoxia Damage = ((2.0 - Current_O2) * O2_Death_Rate) / Robustness
#    - Heat Damage   = ((Current_Temp - 41.0) * Temp_Death_Rate) / Robustness
#    - Acid Damage   = ((6.4 - Current_pH) * pH_Death_Rate) / Robustness
#
# 3. FUNCTIONAL VALUE INDEX (The Yield-to-Value Bridge)
#    - Total_Yield   = Therapeutic_EVs + Stress_EVs + Apoptotic_Debris
#    - Consistency   = Therapeutic_EVs / (Therapeutic_EVs + Stress_EVs)
#    - Purity        = (Therapeutic_EVs + Stress_EVs) / Total_Yield
#    - TRUE_VALUE    = Total_Yield * Consistency * Purity
# ==============================================================================
# ==============================================================================
# FED-BATCH BIOREACTOR OPTIMISATION FOR EV YIELD PRODUCTION
# ==============================================================================

import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# --- CORE BIOLOGICAL MODEL ---
class CellLineKinetics:
    def __init__(self, name, base_ev_rate, optimal_env, kinetic_params, robustness):
        self.name = name
        self.base_ev_rate = base_ev_rate
        self.optimal = optimal_env       
        self.kinetics = kinetic_params   
        self.robustness = robustness     

class FedBatchBioreactorModel:
    def __init__(self, cell_line):
        self.cell = cell_line

    def _calc_gradient_baseline(self, o2, temp, ph):
        dev_o2 = abs(self.cell.optimal['o2'] - o2)
        dev_temp = abs(self.cell.optimal['temp'] - temp)
        dev_ph = abs(self.cell.optimal['ph'] - ph)
        penalty = (dev_o2 * self.cell.kinetics['o2_grad_penalty']) + \
                  (dev_temp * self.cell.kinetics['temp_grad_penalty']) + \
                  (dev_ph * self.cell.kinetics['ph_grad_penalty'])
        adjusted_base = self.cell.base_ev_rate * max(0.1, (1.0 - penalty))
        return max(0, np.random.normal(loc=adjusted_base, scale=(adjusted_base * 0.05)))

    def _calc_acute_stress(self, o2, temp, ph):
        burst = 0.0
        if o2 < (self.cell.optimal['o2'] - self.cell.kinetics['hypoxia_thresh']):
            burst += ((self.cell.optimal['o2'] - self.cell.kinetics['hypoxia_thresh']) - o2)**2 * self.cell.kinetics['hypoxia_mult']
        temp_diff = abs(self.cell.optimal['temp'] - temp)
        if temp_diff > self.cell.kinetics['thermal_thresh']:
            burst += (temp_diff**1.5) * self.cell.kinetics['thermal_mult']
        if ph < (self.cell.optimal['ph'] - self.cell.kinetics['ph_thresh']):
            burst += ((self.cell.optimal['ph'] - ph)**3) * self.cell.kinetics['ph_mult']
        return burst

    def _calc_hourly_damage(self, o2, temp, ph):
        dmg = 0.0
        if o2 < 2.0: dmg += ((2.0 - o2) * self.cell.kinetics['o2_death_rate']) / self.cell.robustness
        if temp > 41.0: dmg += ((temp - 41.0) * self.cell.kinetics['temp_death_rate']) / self.cell.robustness
        if ph < 6.4: dmg += ((6.4 - ph) * self.cell.kinetics['ph_death_rate']) / self.cell.robustness
        return dmg

    def run_simulation(self, init_o2, target_temp, target_ph, mixing_homogeneity, duration_hours):
        main_weight = mixing_homogeneity / 100.0
        dead_weight = 1.0 - main_weight
        current_o2 = init_o2
        viability = {"main": 100.0, "dead": 100.0}
        
        yields = {"therapeutic": 0.0, "stress": 0.0, "apoptotic": 0.0}
        history = {"Hour": [], "Therapeutic EVs": [], "Stress-Altered EVs": [], "Apoptotic Impurities": [], "Cell Viability (%)": []}
        
        for hour in range(1, duration_hours + 1):
            current_o2 = max(0.5, current_o2 - 0.1) 
            dead_o2 = max(0.0, current_o2 * 0.1) 
            dead_ph = max(6.0, target_ph - (hour * 0.02)) 
            
            if viability["main"] > 0:
                yields["therapeutic"] += self._calc_gradient_baseline(current_o2, target_temp, target_ph) * main_weight * (viability["main"]/100)
                yields["stress"] += self._calc_acute_stress(current_o2, target_temp, target_ph) * main_weight * (viability["main"]/100)
                dmg = self._calc_hourly_damage(current_o2, target_temp, target_ph)
                viability["main"] = max(0.0, viability["main"] - dmg)
                yields["apoptotic"] += (dmg * self.cell.base_ev_rate * self.cell.kinetics['apoptotic_mult']) * main_weight
                
            if viability["dead"] > 0:
                yields["therapeutic"] += self._calc_gradient_baseline(dead_o2, target_temp, dead_ph) * dead_weight * (viability["dead"]/100)
                yields["stress"] += self._calc_acute_stress(dead_o2, target_temp, dead_ph) * dead_weight * (viability["dead"]/100)
                dmg = self._calc_hourly_damage(dead_o2, target_temp, dead_ph)
                viability["dead"] = max(0.0, viability["dead"] - dmg)
                yields["apoptotic"] += (dmg * self.cell.base_ev_rate * self.cell.kinetics['apoptotic_mult']) * dead_weight

            history["Hour"].append(hour)
            history["Therapeutic EVs"].append(yields["therapeutic"])
            history["Stress-Altered EVs"].append(yields["stress"])
            history["Apoptotic Impurities"].append(yields["apoptotic"])
            history["Cell Viability (%)"].append((viability["main"] * main_weight) + (viability["dead"] * dead_weight))

        return history

# --- STREAMLIT FRONTEND UI ---
st.set_page_config(page_title="EV Bioreactor Digital Twin", layout="wide")

# Minimal Header
st.title("Bioreactor Optimisation")
st.caption("Based on the Multi-Machinery Model of Biogenesis | Author: Evgenia Trudova")
st.divider()

# Kinetics Initialization
msc_kinetics = {
    'o2_grad_penalty': 0.02, 'temp_grad_penalty': 0.05, 'ph_grad_penalty': 0.1,
    'hypoxia_thresh': 10.0, 'hypoxia_mult': 5e8, 'thermal_thresh': 1.0, 'thermal_mult': 8e8,
    'ph_thresh': 0.2, 'ph_mult': 2e9, 'o2_death_rate': 2.0, 'temp_death_rate': 2.5, 'ph_death_rate': 4.0,
    'apoptotic_mult': 1.5e9
}
msc_cell = CellLineKinetics("Standard EV Cell Line", 1.2e9, {'o2': 21.0, 'temp': 37.0, 'ph': 7.4}, msc_kinetics, 1.0)
bioreactor = FedBatchBioreactorModel(msc_cell)

# Sidebar: UX Optimized for Goal Seeking
st.sidebar.header("1. Clinical Goal")
target_clinical_yield = st.sidebar.number_input("Desired Target Yield (EVs)", min_value=1e10, value=1e15, step=1e13, format="%.1e")

st.sidebar.header("2. Wet Lab Input")
reactor_vol = st.sidebar.number_input("Reactor Volume (L)", min_value=0.1, max_value=2000.0, value=50.0, step=0.5)
init_o2 = st.sidebar.slider("Oxygen (%)", 0.0, 21.0, 15.0, 0.5)
target_temp = st.sidebar.slider("Temperature (°C)", 30.0, 45.0, 38.0, 0.5)
target_ph = st.sidebar.slider("pH", 6.0, 8.0, 7.0, 0.1)
mixing_homo = st.sidebar.slider("Mixing Homogeneity (%)", 50.0, 100.0, 85.0, 1.0)
duration = st.sidebar.slider("Batch Duration (Hours)", 12, 72, 48, 2)

# Run Simulation
history = bioreactor.run_simulation(init_o2, target_temp, target_ph, mixing_homo, duration)
df = pd.DataFrame(history)

# Final Math & Conversions
reactor_vol_ml = reactor_vol * 1000
final_thera_ml = df["Therapeutic EVs"].iloc[-1]
final_stress_ml = df["Stress-Altered EVs"].iloc[-1]
final_apop_ml = df["Apoptotic Impurities"].iloc[-1]

total_evs_ml = final_thera_ml + final_stress_ml
total_particles_ml = total_evs_ml + final_apop_ml 

target_yield_predicted = total_particles_ml * reactor_vol_ml
functional_yield_total = total_evs_ml * reactor_vol_ml

cargo_consistency = (final_thera_ml / total_evs_ml * 100) if total_evs_ml > 0 else 0
downstream_purity = (total_evs_ml / total_particles_ml * 100) if total_particles_ml > 0 else 0

true_functional_value_total = target_yield_predicted * (cargo_consistency / 100) * (downstream_purity / 100)

# Goal Seeking Metric
goal_achieved_pct = (true_functional_value_total / target_clinical_yield) * 100

# --- DASHBOARD RENDER ---
# TOP METRICS (Responsive row)
col1, col2, col3, col4 = st.columns(4)
col1.metric("Predicted Functional Value", f"{true_functional_value_total:.2e}", f"{goal_achieved_pct:.1f}% of Target Goal")
col2.metric("Harvest Concentration", f"{total_particles_ml:.2e} ev/mL", f"In {reactor_vol}L Tank")
col3.metric("Downstream Purity", f"{downstream_purity:.1f}%")
col4.metric("Cargo Consistency", f"{cargo_consistency:.1f}%")

st.divider()

# --- GRAPH 1: YIELD TO VALUE FUNNEL (Plotly) ---
st.markdown("### The Yield-to-Value Bridge")
funnel_data = dict(
    stage=["1. Raw Target Yield", "2. Intact EVs (Purity)", "3. True Functional Value"],
    count=[target_yield_predicted, functional_yield_total, true_functional_value_total]
)
fig_funnel = px.funnel(funnel_data, x='count', y='stage', color_discrete_sequence=['#4C78A8'])
fig_funnel.update_layout(margin=dict(l=20, r=20, t=20, b=20), height=300)
st.plotly_chart(fig_funnel, use_container_width=True)

with st.expander("Formula & Explanation"):
    st.markdown(r"$$ V_{functional} = (C_{harvest} \times V_{reactor}) \times \left(\frac{Purity}{100}\right) \times \left(\frac{Consistency}{100}\right) $$")
    st.markdown("Deconstructs raw NTA scatter count by penalizing for apoptotic debris and stress-altered vesicles.")

st.divider()

col_left, col_right = st.columns(2)

# --- GRAPH 2: ACCUMULATION (Plotly) ---
with col_left:
    st.markdown("### Accumulation of Particles")
    
    # Restructure dataframe for Plotly
    df_melted = df.melt(id_vars=['Hour'], value_vars=['Therapeutic EVs', 'Stress-Altered EVs', 'Apoptotic Impurities'], 
                        var_name='Particle Type', value_name='Count (per mL)')
    
    # Colorblind safe palette
    fig_acc = px.line(df_melted, x='Hour', y='Count (per mL)', color='Particle Type',
                      color_discrete_map={"Therapeutic EVs": "#4C78A8", "Stress-Altered EVs": "#F58518", "Apoptotic Impurities": "#E45756"})
    
    fig_acc.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1), margin=dict(l=0, r=0, t=0, b=0))
    # Enforce scientific notation on Y axis for readability
    fig_acc.update_yaxes(type="log", tickformat=".1e")
    fig_acc.update_yaxes(tickformat=".1e")
    
    st.plotly_chart(fig_acc, use_container_width=True)
    
    with st.expander("Formula & Explanation"):
        st.markdown(r"$$ \frac{d(EV)}{dt} = k_{base}(1 - P_{grad}) + \sum_{s=1}^{n} k_{s}(|E_{opt} - E_{current}|)^m $$")
        st.markdown("Tracks biogenesis baseline against acute stress bursts triggered by environmental toxicity.")

# --- GRAPH 3: VIABILITY (Plotly) ---
with col_right:
    st.markdown("### Cellular Viability Curve")
    fig_viab = px.line(df, x='Hour', y='Cell Viability (%)', color_discrete_sequence=['#72B7B2'])
    fig_viab.update_layout(margin=dict(l=0, r=0, t=0, b=0), yaxis_range=[0, 100])
    st.plotly_chart(fig_viab, use_container_width=True)
    
    with st.expander("Formula & Explanation"):
        st.markdown(r"$$ Viability_{t} = Viability_{t-1} - \left( \frac{\mu_{o2} + \mu_{temp} + \mu_{pH}}{\rho_{cell}} \right) $$")
        st.markdown("Models active cell death based on cumulative environmental toxicity and unmixed hydrodynamic dead zones.")
