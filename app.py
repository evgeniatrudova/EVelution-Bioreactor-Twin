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

import streamlit as st
import numpy as np
import pandas as pd

# --- CORE BIOLOGICAL MODEL ---
class CellLineKinetics:
    def __init__(self, name, base_ev_rate, optimal_env, kinetic_params, robustness):
        self.name = name
        self.base_ev_rate = base_ev_rate  # Baseline EVs per mL per hour
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
        
        # Tracking yields per mL
        yields = {"therapeutic": 0.0, "stress": 0.0, "apoptotic": 0.0}
        history = {"Hour": [], "Therapeutic EVs (per mL)": [], "Stress-Altered EVs (per mL)": [], "Apoptotic Impurities (per mL)": [], "Cell Viability (%)": []}
        
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
            history["Therapeutic EVs (per mL)"].append(yields["therapeutic"])
            history["Stress-Altered EVs (per mL)"].append(yields["stress"])
            history["Apoptotic Impurities (per mL)"].append(yields["apoptotic"])
            history["Cell Viability (%)"].append((viability["main"] * main_weight) + (viability["dead"] * dead_weight))

        return history

# --- STREAMLIT FRONTEND UI ---
st.set_page_config(page_title="EV Bioreactor Digital Twin", layout="wide")

st.title("Fed-Batch Bioreactor Optimisation")
st.markdown("**Based on the Multi-Machinery Model of Biogenesis | Author: Evgenia Trudova**")

# Initialize Cell Line (Scaled to realistic realistic billions of EVs/mL)
msc_kinetics = {
    'o2_grad_penalty': 0.02, 'temp_grad_penalty': 0.05, 'ph_grad_penalty': 0.1,
    'hypoxia_thresh': 10.0, 'hypoxia_mult': 5e8, 'thermal_thresh': 1.0, 'thermal_mult': 8e8,
    'ph_thresh': 0.2, 'ph_mult': 2e9, 'o2_death_rate': 2.0, 'temp_death_rate': 2.5, 'ph_death_rate': 4.0,
    'apoptotic_mult': 1.5e9
}
# Base rate set to 1.2 billion EVs produced per mL per hour under ideal conditions
msc_cell = CellLineKinetics("Standard EV Cell Line", 1.2e9, {'o2': 21.0, 'temp': 37.0, 'ph': 7.4}, msc_kinetics, 1.0)
bioreactor = FedBatchBioreactorModel(msc_cell)

# Sidebar Controls (The Wet Lab Interface)
st.sidebar.header("Wet Lab Input")
init_o2 = st.sidebar.slider("Oxygen (%)", 0.0, 21.0, 15.0, 0.5)
target_temp = st.sidebar.slider("Temperature (°C)", 30.0, 45.0, 38.0, 0.5)
target_ph = st.sidebar.slider("pH", 6.0, 8.0, 7.0, 0.1)
mixing_homo = st.sidebar.slider("Mixing Homogeneity (%)", 50.0, 100.0, 85.0, 1.0)
duration = st.sidebar.slider("Batch Duration (Hours)", 12, 72, 48, 2)
reactor_vol = st.sidebar.number_input("Reactor Volume (L)", min_value=0.1, max_value=2000.0, value=50.0, step=0.5)

# Run Simulation
history = bioreactor.run_simulation(init_o2, target_temp, target_ph, mixing_homo, duration)
df = pd.DataFrame(history).set_index("Hour")

# Extract final harvest numbers (per mL)
final_thera_ml = df["Therapeutic EVs (per mL)"].iloc[-1]
final_stress_ml = df["Stress-Altered EVs (per mL)"].iloc[-1]
final_apop_ml = df["Apoptotic Impurities (per mL)"].iloc[-1]
final_viab = df["Cell Viability (%)"].iloc[-1]

# --- THE MATH & CONVERSIONS ---
reactor_vol_ml = reactor_vol * 1000

# Per mL Concentrations
total_evs_ml = final_thera_ml + final_stress_ml
total_particles_ml = total_evs_ml + final_apop_ml # What the NTA actually counts per mL

# Total Tank Yields (Absolute Numbers)
target_yield_total = total_particles_ml * reactor_vol_ml
functional_yield_total = total_evs_ml * reactor_vol_ml

# Quality Coefficients
cargo_consistency = (final_thera_ml / total_evs_ml * 100) if total_evs_ml > 0 else 0
downstream_purity = (total_evs_ml / total_particles_ml * 100) if total_particles_ml > 0 else 0

# True Functional Value (Scale adjusted for readability on graph)
true_functional_value_total = target_yield_total * (cargo_consistency / 100) * (downstream_purity / 100)

# --- DASHBOARD RENDER ---
# TOP METRICS
col1, col2, col3, col4 = st.columns(4)
col1.metric("Harvest Concentration", f"{total_particles_ml:.2e} ev/mL", "Total particles measured by NTA")
col2.metric("Target Yield (Total EVs)", f"{target_yield_total:.2e}", f"Across {reactor_vol}L Volume")
col3.metric("Downstream Purity", f"{downstream_purity:.1f}%", "Ratio of intact vesicles to debris")
col4.metric("Cargo Consistency", f"{cargo_consistency:.1f}%", "Ratio of therapeutic to stress EVs")

st.divider()

# --- GRAPH 1: YIELD TO VALUE BRIDGE ---
st.markdown("### The Yield-to-Value Bridge")
# Create a dataframe specifically for the bar chart comparison
bridge_data = pd.DataFrame({
    "Metric": ["1. Raw Target Yield (Total Particles)", "2. Total Intact EVs (Purity Applied)", "3. True Functional Value (Consistency Applied)"],
    "Total Count": [target_yield_total, functional_yield_total, true_functional_value_total]
}).set_index("Metric")

st.bar_chart(bridge_data, color="#4CAF50")

with st.expander("Mathematical Formula & Academic Explanation"):
    st.markdown(r"$$ V_{functional} = (C_{harvest} \times V_{reactor}) \times \left(\frac{Purity}{100}\right) \times \left(\frac{Consistency}{100}\right) $$")
    st.markdown("""
    Standard bioprocessing often conflates total nanoparticle count with successful biogenesis. The Yield-to-Value Bridge mathematically deconstructs the raw NTA scatter count ($C_{harvest}$). 
    
    It calculates the absolute particle mass by multiplying concentration by the total reactor volume ($V_{reactor}$). It then applies a fractional purity coefficient to subtract apoptotic bodies and genomic debris caused by hydrodynamic shear stress and necrosis. Finally, it applies a consistency penalty, filtering out structurally intact but functionally altered "stress vesicles" generated by extreme metabolic shifts, revealing the true therapeutic lot size.
    """)

st.divider()

col_left, col_right = st.columns(2)

# --- GRAPH 2: ACCUMULATION OF PARTICLES ---
with col_left:
    st.markdown("### Accumulation of Particles Over Time")
    st.line_chart(df[["Therapeutic EVs (per mL)", "Stress-Altered EVs (per mL)", "Apoptotic Impurities (per mL)"]])
    
    with st.expander("Mathematical Formula & Academic Explanation"):
        st.markdown(r"$$ \frac{d(EV)}{dt} = k_{base}(1 - P_{grad}) + \sum_{s=1}^{n} k_{s}(|E_{opt} - E_{current}|)^m $$")
        st.markdown("""
        
        Vesicle accumulation is modeled as a set of non-linear kinetic equations. The baseline biogenesis rate ($k_{base}$) occurs during cellular homeostasis but is suppressed by a physiological gradient penalty ($P_{grad}$) when conditions drift. 
        
        Acute environmental stressors—such as dropping below the hypoxia threshold or exceeding physiological temperature—trigger exponential stress bursts (the summation term). The exponent ($m$) dictates the severity of the cellular response, simulating how cells rapidly shed lipid bilayers to eject toxic proteins or communicate distress to adjacent cells before impending apoptosis.
        """)

# --- GRAPH 3: CELLULAR VIABILITY ---
with col_right:
    st.markdown("### Cellular Viability Curve")
    st.line_chart(df[["Cell Viability (%)"]], color="#ff4b4b")
    
    with st.expander("Mathematical Formula & Academic Explanation"):
        st.markdown(r"$$ Viability_{t} = Viability_{t-1} - \left( \frac{\mu_{o2} + \mu_{temp} + \mu_{pH}}{\rho_{cell}} \right) $$")
        st.markdown("""
        
        Viability is not modeled as a static timeline, but as an actively depleting resource dependent on cumulative environmental toxicity ($\mu$). 
        
        The model uniquely accounts for the physical reality of large-scale bioreactors by integrating the Mixing Homogeneity factor. Stagnant hydrodynamic "dead zones" experience rapid localized oxygen depletion and severe lactate-induced pH drops. The total fractional cell death is determined by the cellular membrane robustness ($\rho_{cell}$). When viability crashes, the corresponding rate of Apoptotic Impurities (shown in the adjacent graph) spikes exponentially.
        """)
