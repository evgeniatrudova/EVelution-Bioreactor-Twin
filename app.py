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
    """Stores the specific biological traits and tolerances of the EV-producing cell line."""
    def __init__(self, name, base_ev_rate, optimal_env, kinetic_params, robustness):
        self.name = name
        self.base_ev_rate = base_ev_rate  # Baseline EVs produced per hour under perfect conditions
        self.optimal = optimal_env        # Dictionary of perfect homeostatic conditions (O2, Temp, pH)
        self.kinetics = kinetic_params    # The multipliers for how severely the cell reacts to stress
        self.robustness = robustness      # Resistance to lysis (higher = tougher cell membrane)

class FedBatchBioreactorModel:
    """Simulates the physical bioreactor environment and its impact on the cells over time."""
    def __init__(self, cell_line):
        self.cell = cell_line

    def _calc_gradient_baseline(self, o2, temp, ph):
        """Calculates normal, healthy 'Therapeutic EV' production, penalized by mild sub-optimal conditions."""
        dev_o2 = abs(self.cell.optimal['o2'] - o2)
        dev_temp = abs(self.cell.optimal['temp'] - temp)
        dev_ph = abs(self.cell.optimal['ph'] - ph)
        
        # Calculate how far we are from homeostatic perfection
        penalty = (dev_o2 * self.cell.kinetics['o2_grad_penalty']) + \
                  (dev_temp * self.cell.kinetics['temp_grad_penalty']) + \
                  (dev_ph * self.cell.kinetics['ph_grad_penalty'])
                  
        # Cells reduce baseline production when mildly stressed, down to a minimum of 10%
        adjusted_base = self.cell.base_ev_rate * max(0.1, (1.0 - penalty))
        
        # Add slight biological noise (Gaussian distribution) for realism
        return max(0, np.random.normal(loc=adjusted_base, scale=(adjusted_base * 0.05)))

    def _calc_acute_stress(self, o2, temp, ph):
        """Calculates 'Stress-Altered EVs' triggered by severe, acute environmental panic."""
        burst = 0.0
        
        # Hypoxic Trigger (Exponential increase as O2 drops below threshold)
        if o2 < (self.cell.optimal['o2'] - self.cell.kinetics['hypoxia_thresh']):
            burst += ((self.cell.optimal['o2'] - self.cell.kinetics['hypoxia_thresh']) - o2)**2 * self.cell.kinetics['hypoxia_mult']
            
        # Thermal Shock Trigger (Non-linear increase based on temp deviation)
        temp_diff = abs(self.cell.optimal['temp'] - temp)
        if temp_diff > self.cell.kinetics['thermal_thresh']:
            burst += (temp_diff**1.5) * self.cell.kinetics['thermal_mult']
            
        # Acidic Shock Trigger (Cubic increase as pH drops dangerously low)
        if ph < (self.cell.optimal['ph'] - self.cell.kinetics['ph_thresh']):
            burst += ((self.cell.optimal['ph'] - ph)**3) * self.cell.kinetics['ph_mult']
            
        return burst

    def _calc_hourly_damage(self, o2, temp, ph):
        """Calculates percentage of cell death based on extreme toxicity or dead zones."""
        dmg = 0.0
        if o2 < 2.0: dmg += ((2.0 - o2) * self.cell.kinetics['o2_death_rate']) / self.cell.robustness
        if temp > 41.0: dmg += ((temp - 41.0) * self.cell.kinetics['temp_death_rate']) / self.cell.robustness
        if ph < 6.4: dmg += ((6.4 - ph) * self.cell.kinetics['ph_death_rate']) / self.cell.robustness
        return dmg

    def run_simulation(self, init_o2, target_temp, target_ph, mixing_homogeneity, duration_hours):
        """Runs the step-by-step physical simulation over the specified batch duration."""
        
        # Split the reactor volume into a perfectly mixed zone and a stagnant "dead" zone
        main_weight = mixing_homogeneity / 100.0
        dead_weight = 1.0 - main_weight
        
        current_o2 = init_o2
        viability = {"main": 100.0, "dead": 100.0}
        yields = {"therapeutic": 0.0, "stress": 0.0, "apoptotic": 0.0}
        
        # Track data for Streamlit charts
        history = {"Hour": [], "Therapeutic EVs": [], "Stress-Altered EVs": [], "Apoptotic Impurities": [], "Cell Viability (%)": []}
        
        for hour in range(1, duration_hours + 1):
            # Simulate natural physical decay in the tank
            current_o2 = max(0.5, current_o2 - 0.1) # Natural O2 consumption
            dead_o2 = max(0.0, current_o2 * 0.1)    # Dead zones lack oxygen transfer
            dead_ph = max(6.0, target_ph - (hour * 0.02)) # Lactate buildup in stagnant zones
            
            # --- MAIN MIXED ZONE CALCULATIONS ---
            if viability["main"] > 0:
                # Add EVs based on remaining living cells in the main volume
                yields["therapeutic"] += self._calc_gradient_baseline(current_o2, target_temp, target_ph) * main_weight * (viability["main"]/100)
                yields["stress"] += self._calc_acute_stress(current_o2, target_temp, target_ph) * main_weight * (viability["main"]/100)
                
                # Apply damage, subtract living cells, and add the resulting debris to Apoptotic yield
                dmg = self._calc_hourly_damage(current_o2, target_temp, target_ph)
                viability["main"] = max(0.0, viability["main"] - dmg)
                yields["apoptotic"] += (dmg * self.cell.base_ev_rate * self.cell.kinetics['apoptotic_mult']) * main_weight
                
            # --- DEAD ZONE CALCULATIONS ---
            if viability["dead"] > 0:
                yields["therapeutic"] += self._calc_gradient_baseline(dead_o2, target_temp, dead_ph) * dead_weight * (viability["dead"]/100)
                yields["stress"] += self._calc_acute_stress(dead_o2, target_temp, dead_ph) * dead_weight * (viability["dead"]/100)
                
                dmg = self._calc_hourly_damage(dead_o2, target_temp, dead_ph)
                viability["dead"] = max(0.0, viability["dead"] - dmg)
                yields["apoptotic"] += (dmg * self.cell.base_ev_rate * self.cell.kinetics['apoptotic_mult']) * dead_weight

            # Log hourly snapshot for charts
            history["Hour"].append(hour)
            history["Therapeutic EVs"].append(yields["therapeutic"])
            history["Stress-Altered EVs"].append(yields["stress"])
            history["Apoptotic Impurities"].append(yields["apoptotic"])
            history["Cell Viability (%)"].append((viability["main"] * main_weight) + (viability["dead"] * dead_weight))

        return history

# --- STREAMLIT FRONTEND UI ---
st.set_page_config(page_title="EV Bioreactor Digital Twin", layout="wide")

st.title("Fed-Batch Bioreactor Optimisation: Yield-to-Value Bridge")
st.markdown("**Based on the Multi-Machinery Model of Biogenesis | Author: Evgenia Trudova**")

# Define the baseline metabolic traits for our theoretical cell line
msc_kinetics = {
    'o2_grad_penalty': 0.02, 'temp_grad_penalty': 0.05, 'ph_grad_penalty': 0.1,
    'hypoxia_thresh': 10.0, 'hypoxia_mult': 10.0, 'thermal_thresh': 1.0, 'thermal_mult': 15.0,
    'ph_thresh': 0.2, 'ph_mult': 50.0, 'o2_death_rate': 2.0, 'temp_death_rate': 2.5, 'ph_death_rate': 4.0,
    'apoptotic_mult': 3.5
}
msc_cell = CellLineKinetics("Standard EV Cell Line", 200, {'o2': 21.0, 'temp': 37.0, 'ph': 7.4}, msc_kinetics, 1.0)
bioreactor = FedBatchBioreactorModel(msc_cell)

# Sidebar Controls (The Wet Lab Interface)
st.sidebar.header("Wet Lab Constraints (Inputs)")
init_o2 = st.sidebar.slider("Initial Oxygen (%)", 0.0, 21.0, 15.0, 0.5)
target_temp = st.sidebar.slider("Target Temperature (°C)", 30.0, 45.0, 38.0, 0.5)
target_ph = st.sidebar.slider("Target pH", 6.0, 8.0, 7.0, 0.1)
mixing_homo = st.sidebar.slider("Mixing Homogeneity (%)", 50.0, 100.0, 85.0, 1.0)
duration = st.sidebar.slider("Batch Duration (Hours)", 12, 72, 48, 2)

# Run Simulation and process data
history = bioreactor.run_simulation(init_o2, target_temp, target_ph, mixing_homo, duration)
df = pd.DataFrame(history).set_index("Hour")

# Extract final harvest numbers
final_thera = df["Therapeutic EVs"].iloc[-1]
final_stress = df["Stress-Altered EVs"].iloc[-1]
final_apop = df["Apoptotic Impurities"].iloc[-1]
final_viab = df["Cell Viability (%)"].iloc[-1]

# --- THE YIELD-TO-VALUE MATH ---
# Calculate total vesicles (Therapeutic + Stress induced)
total_evs = final_thera + final_stress
# Add cellular debris to get the raw number the NTA machine will read
total_particles = total_evs + final_apop 

# Calculate quality coefficients
cargo_consistency = (final_thera / total_evs * 100) if total_evs > 0 else 0
downstream_purity = (total_evs / total_particles * 100) if total_particles > 0 else 0

# The ultimate formula: Total Yield filtered by Purity and Quality constraints
functional_value = total_particles * (cargo_consistency / 100) * (downstream_purity / 100)

# --- DASHBOARD RENDER ---
st.markdown("### 1. The Yield-to-Value Bridge")
st.info("Traditional bioprocessing optimizes for **Total Particle Yield** (what an NTA counts). This model translates that raw count into true **Functional Value** by mathematically penalizing the batch for empty cargo and apoptotic cellular debris.")

# Render KPI Columns
col1, col2, col3, col4 = st.columns(4)
col1.metric("Raw Particle Yield (NTA)", f"{total_particles:,.0f}", "Total entities counted")
col2.metric("Downstream Purity", f"{downstream_purity:.1f}%", "Free of apoptotic debris")
col3.metric("Cargo Consistency", f"{cargo_consistency:.1f}%", "Therapeutic alignment")
col4.metric("TRUE FUNCTIONAL VALUE", f"{functional_value:,.0f}", "Optimized target metric")

st.divider()

# Render Charts side-by-side
col_left, col_right = st.columns(2)

with col_left:
    st.markdown("### 2. Accumulation of Particles Over Time")
    st.line_chart(df[["Therapeutic EVs", "Stress-Altered EVs", "Apoptotic Impurities"]])
    with st.expander("How to read this chart"):
        st.write("""
        **What you are seeing:** This chart breaks your raw NTA yield into its actual biological components. 
        - **Therapeutic EVs:** Steady baseline biogenesis under optimal conditions.
        - **Stress-Altered EVs:** Rapid spikes when the bioreactor environment triggers cellular panic (e.g., hypoxia).
        - **Apoptotic Impurities:** Accumulation of cellular debris as dead zones or toxic parameters cause cell lysis.
        **Optimization Goal:** Maximize the gap between Therapeutic EVs and Apoptotic Impurities.
        """)

with col_right:
    st.markdown("### 3. Cellular Viability Curve")
    st.line_chart(df[["Cell Viability (%)"]], color="#ff0000")
    with st.expander("How to read this chart"):
        st.write("""
        **What you are seeing:** The holistic health of the cell culture, factoring in both the main mixed zone and the unmixed 'dead zones'.
        **Optimization Goal:** While dropping viability often causes an initial spike in stress-altered vesicles, a complete crash (below 60%) will flood the harvest with genomic DNA and fragmented lipid bilayers, permanently destroying downstream purity.
        """)
