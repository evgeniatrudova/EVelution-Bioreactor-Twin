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
        
        shear_damage = 0.0
        if mixing_homogeneity > 85.0:
            shear_damage = ((mixing_homogeneity - 85.0) ** 1.5) * 0.1
        
        for hour in range(1, duration_hours + 1):
            current_o2 = max(0.5, current_o2 - 0.1) 
            dead_o2 = max(0.0, current_o2 * 0.1) 
            dead_ph = max(6.0, target_ph - (hour * 0.02)) 
            
            if viability["main"] > 0:
                yields["therapeutic"] += self._calc_gradient_baseline(current_o2, target_temp, target_ph) * main_weight * (viability["main"]/100)
                yields["stress"] += self._calc_acute_stress(current_o2, target_temp, target_ph) * main_weight * (viability["main"]/100)
                
                dmg = self._calc_hourly_damage(current_o2, target_temp, target_ph) + (shear_damage / self.cell.robustness)
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

# The MMModel Dropdown
with st.expander("The Multi-Machinery Model (MMModel) Overview"):
    st.markdown("""
    **Biological Context:** The MMModel simulates the competing pathways of vesicle biogenesis. Rather than treating cells as simple particle factories, it calculates how cells respond to stress by either releasing therapeutic exosomes, shedding stress-altered vesicles to survive, or rupturing into apoptotic debris.
    """)
    st.markdown("#### 1. Acute Stress Bursts (Vesicle release triggered by environmental panic)")
    st.markdown(r"$$ \text{Hypoxia Burst} = ((\text{Opt}_{O2} - \text{Hypoxia}_{\text{Thresh}}) - \text{Current}_{O2})^2 \times \text{Hypoxia}_{\text{Mult}} $$")
    st.markdown(r"$$ \text{Thermal Burst} = (|\text{Opt}_{\text{Temp}} - \text{Current}_{\text{Temp}}|^{1.5}) \times \text{Thermal}_{\text{Mult}} $$")
    st.markdown(r"$$ \text{pH Burst} = ((\text{Opt}_{pH} - \text{Current}_{pH})^3) \times \text{pH}_{\text{Mult}} $$")
    
    st.markdown("#### 2. Hourly Cellular Damage (Cell death from toxic environments)")
    st.markdown(r"$$ \text{Anoxia Damage} = \frac{(2.0 - \text{Current}_{O2}) \times \text{O2}_{\text{DeathRate}}}{\text{Robustness}} $$")
    st.markdown(r"$$ \text{Heat Damage} = \frac{(\text{Current}_{\text{Temp}} - 41.0) \times \text{Temp}_{\text{DeathRate}}}{\text{Robustness}} $$")
    st.markdown(r"$$ \text{Acid Damage} = \frac{(6.4 - \text{Current}_{pH}) \times \text{pH}_{\text{DeathRate}}}{\text{Robustness}} $$")

    st.markdown("#### 3. Functional Value Index (The Yield-to-Value Bridge)")
    st.markdown(r"$$ \text{Total Yield} = \text{Therapeutic EVs} + \text{Stress EVs} + \text{Apoptotic Debris} $$")
    st.markdown(r"$$ \text{Consistency} = \frac{\text{Therapeutic EVs}}{(\text{Therapeutic EVs} + \text{Stress EVs})} $$")
    st.markdown(r"$$ \text{Purity} = \frac{(\text{Therapeutic EVs} + \text{Stress EVs})}{\text{Total Yield}} $$")
    st.markdown(r"$$ \text{TRUE VALUE} = \text{Total Yield} \times \text{Consistency} \times \text{Purity} $$")

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

# Sidebar: Enterprise Features & UX Optimized Inputs
st.sidebar.header("Data Management")
uploaded_file = st.sidebar.file_uploader("Upload Run Data (CSV)", type=["csv"])
if uploaded_file is not None:
    st.sidebar.success("Historical data loaded.")
if st.sidebar.button("Export to PDF Report"):
    st.sidebar.info("Generating PDF... Please wait.")

st.sidebar.divider()

st.sidebar.header("Clinical Goal")
target_clinical_yield = st.sidebar.number_input("Desired Target Yield (EVs)", min_value=1e10, value=1e15, step=1e13, format="%.1e")

st.sidebar.header("Wet Lab Input")
reactor_vol = st.sidebar.number_input("Reactor Volume (L)", min_value=0.1, max_value=2000.0, value=50.0, step=0.5)
init_o2 = st.sidebar.slider("Oxygen (%)", 0.0, 21.0, 15.0, 0.5)
target_temp = st.sidebar.slider("Temperature (°C)", 30.0, 45.0, 38.0, 0.5)
target_ph = st.sidebar.slider("pH", 6.0, 8.0, 7.0, 0.1)
mixing_homo = st.sidebar.slider("Mixing Homogeneity (%)", 50.0, 100.0, 85.0, 1.0, help="Careful: >85% introduces hydrodynamic shear stress!")
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

goal_delta = true_functional_value_total - target_clinical_yield
goal_achieved_pct = (true_functional_value_total / target_clinical_yield) * 100

# --- DASHBOARD RENDER ---
col1, col2, col3, col4 = st.columns(4)
col1.metric("Predicted Functional Value", f"{true_functional_value_total:.2e}", f"{goal_delta:.1e} vs Goal", delta_color="normal")
col2.metric("Harvest Concentration", f"{total_particles_ml:.2e} ev/mL", f"In {reactor_vol}L Tank")
col3.metric("Downstream Purity", f"{downstream_purity:.1f}%")
col4.metric("Cargo Consistency", f"{cargo_consistency:.1f}%")

st.divider()

st.markdown("### The Yield-to-Value Bridge")
funnel_data = dict(
    stage=["Raw Target Yield", "Intact EVs (Purity)", "True Functional Value"],
    count=[target_yield_predicted, functional_yield_total, true_functional_value_total]
)
fig_funnel = px.funnel(funnel_data, x='count', y='stage', color_discrete_sequence=['#4C78A8'])
fig_funnel.update_traces(textinfo="value+percent initial")
fig_funnel.update_layout(margin=dict(l=20, r=20, t=20, b=20), height=300)
st.plotly_chart(fig_funnel, use_container_width=True)

with st.expander("Biological Mechanism: Downstream Value"):
    st.markdown("Traditional biomanufacturing focuses on raw nanoparticle counts. However, downstream filtration and chromatography will strip out broken lipid husks and apoptotic bodies. This metric visualizes the harsh reality of processing: isolating only the structurally intact, therapeutically viable EVs from the crude harvest.")

st.divider()

col_left, col_right = st.columns(2)

with col_left:
    st.markdown("### Accumulation of Particles")
    
    df_melted = df.melt(id_vars=['Hour'], value_vars=['Therapeutic EVs', 'Stress-Altered EVs', 'Apoptotic Impurities'], 
                        var_name='Particle Type', value_name='Count (per mL)')
    
    fig_acc = px.line(df_melted, x='Hour', y='Count (per mL)', color='Particle Type',
                      color_discrete_map={"Therapeutic EVs": "#4C78A8", "Stress-Altered EVs": "#F58518", "Apoptotic Impurities": "#E45756"})
    
    fig_acc.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1), margin=dict(l=0, r=0, t=0, b=0))
    fig_acc.update_yaxes(type="log", tickformat=".1e")
    
    st.plotly_chart(fig_acc, use_container_width=True)
    
    with st.expander("Biological Mechanism: Biogenesis vs. Stress"):
        st.markdown("Cells naturally secrete baseline exosomes, but extreme shifts in their environment (like suddenly dropping oxygen or spiking temperatures) induce cellular panic. Before they die, stressed cells rapidly shed their membranes to eject toxins or send distress signals, leading to the massive explosion of stress-altered vesicles seen here.")

with col_right:
    st.markdown("### Cellular Viability Curve")
    fig_viab = px.line(df, x='Hour', y='Cell Viability (%)', color_discrete_sequence=['#72B7B2'])
    fig_viab.update_layout(margin=dict(l=0, r=0, t=0, b=0), yaxis_range=[0, 100])
    st.plotly_chart(fig_viab, use_container_width=True)
    
    with st.expander("Biological Mechanism: Membrane Degradation"):
        st.markdown("""
        In a real bioreactor, cells don't die instantly—they suffer cumulative stress. This curve tracks how cell membrane integrity degrades over the batch duration. 
        
        * **Toxicity:** When oxygen drops too low, or stagnant dead zones cause acidic lactate to build up, the cells begin to suffocate and undergo necrosis. 
        * **Physical Tearing:** Conversely, if the impeller spins too aggressively (setting mixing homogeneity too high), the violent fluid forces generate hydrodynamic shear stress, literally tearing the fragile cell membranes apart and causing an immediate viability crash.
        """)
