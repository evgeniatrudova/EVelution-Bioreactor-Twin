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
            
            # Main Zone
            if viability["main"] > 0:
                yields["therapeutic"] += self._calc_gradient_baseline(current_o2, target_temp, target_ph) * main_weight * (viability["main"]/100)
                yields["stress"] += self._calc_acute_stress(current_o2, target_temp, target_ph) * main_weight * (viability["main"]/100)
                dmg = self._calc_hourly_damage(current_o2, target_temp, target_ph)
                viability["main"] = max(0.0, viability["main"] - dmg)
                yields["apoptotic"] += (dmg * self.cell.base_ev_rate * self.cell.kinetics['apoptotic_mult']) * main_weight
                
            # Dead Zone
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

st.title("Fed-Batch Bioreactor Optimisation for EV Yield")
st.subheader("Based on the Multi-Machinery Model of Biogenesis | Author: Evgenia Trudova")
st.caption("© 2026 Evgenia Trudova. All Rights Reserved. For academic demonstration only.")

# Initialize Cell Line
msc_kinetics = {
    'o2_grad_penalty': 0.02, 'temp_grad_penalty': 0.05, 'ph_grad_penalty': 0.1,
    'hypoxia_thresh': 10.0, 'hypoxia_mult': 10.0, 'thermal_thresh': 1.0, 'thermal_mult': 15.0,
    'ph_thresh': 0.2, 'ph_mult': 50.0, 'o2_death_rate': 2.0, 'temp_death_rate': 2.5, 'ph_death_rate': 4.0,
    'apoptotic_mult': 3.5
}
msc_cell = CellLineKinetics("Standard EV Cell Line", 200, {'o2': 21.0, 'temp': 37.0, 'ph': 7.4}, msc_kinetics, 1.0)
bioreactor = FedBatchBioreactorModel(msc_cell)

# Sidebar Controls
st.sidebar.header("Bioreactor Controls")
init_o2 = st.sidebar.slider("Initial Oxygen (%)", 0.0, 21.0, 15.0, 0.5)
target_temp = st.sidebar.slider("Target Temperature (°C)", 30.0, 45.0, 38.0, 0.5)
target_ph = st.sidebar.slider("Target pH", 6.0, 8.0, 7.0, 0.1)
mixing_homo = st.sidebar.slider("Mixing Homogeneity (%)", 50.0, 100.0, 85.0, 1.0)
duration = st.sidebar.slider("Batch Duration (Hours)", 12, 72, 48, 2)

# Run Simulation
history = bioreactor.run_simulation(init_o2, target_temp, target_ph, mixing_homo, duration)
df = pd.DataFrame(history).set_index("Hour")

# Calculate KPIs
final_thera = df["Therapeutic EVs"].iloc[-1]
final_stress = df["Stress-Altered EVs"].iloc[-1]
final_apop = df["Apoptotic Impurities"].iloc[-1]
final_viab = df["Cell Viability (%)"].iloc[-1]

total_evs = final_thera + final_stress
cargo_consistency = (final_thera / total_evs * 100) if total_evs > 0 else 0
downstream_purity = (total_evs / (total_evs + final_apop) * 100) if (total_evs + final_apop) > 0 else 0

# Display KPIs
col1, col2, col3, col4 = st.columns(4)
col1.metric("Final Cell Viability", f"{final_viab:.1f}%")
col2.metric("Total EV Yield", f"{total_evs:,.0f}")
col3.metric("Cargo Consistency", f"{cargo_consistency:.1f}%")
col4.metric("Downstream Purity", f"{downstream_purity:.1f}%")

# Display Charts
st.markdown("### Accumulation of Bioreactor Particles Over Time")
st.line_chart(df[["Therapeutic EVs", "Stress-Altered EVs", "Apoptotic Impurities"]])

st.markdown("### Cellular Viability Curve")
st.line_chart(df[["Cell Viability (%)"]], color="#ff0000")
