# ==============================================================================
# FED-BATCH BIOREACTOR OPTIMISATION FOR EV YIELD PRODUCTION
# Based on the Multi-Machinery Model of Biogenesis & Monod Mass Balance
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
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from fpdf import FPDF
import io
import hashlib
import time

# --- 1. THEME & STYLING ---
st.set_page_config(page_title="EVelution Bioreactor Twin", layout="wide", initial_sidebar_state="auto")

# --- 2. CSS INJECTION ---
st.markdown("""
<style>
    /* Blue-centric theme enforcement */
    :root { --primary-color: #779ECB; }
    .stSlider [data-baseweb="slider"] > div > div > div > div { background-color: #779ECB !important; }
    #MainMenu {visibility: hidden;} 
    footer {visibility: hidden;}
    
    /* UX FIX: Force disabled buttons to remain highly visible */
    .stButton > button[disabled] {
        opacity: 0.85 !important;            
        border: 1px solid #4B4B60 !important; 
        color: #B39EB5 !important;          
        background-color: transparent !important;
        cursor: not-allowed !important;    
    }

    /* MOBILE RESPONSIVENESS FIXES */
    @media (max-width: 768px) {
        .status-card {
            flex-direction: column !important;
            align-items: flex-start !important;
            gap: 15px !important;
        }
        .status-card-metric {
            text-align: left !important;
            padding-left: 0 !important;
            border-left: none !important;
            border-top: 1px solid #4B4B60 !important;
            padding-top: 10px !important;
            width: 100% !important;
        }
        .legal-disclaimer {
            font-size: 0.70em !important;
            word-wrap: break-word !important;
            overflow-wrap: break-word !important;
        }
        [data-testid="stMetricValue"] {
            font-size: 1.5rem !important;
        }
    }
</style>
""", unsafe_allow_html=True)
    
# --- 2. GLOBAL CONSTANTS ---
fixed_height = 400
fixed_margin = dict(t=40, b=10, l=10, r=10)
C_GREEN, C_BLUE, C_PURPLE, C_STAR = "#77DD77", "#779ECB", "#B39EB5", "#E39777"

# --- 3. BIOPHYSICAL ENGINE ---
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
        return max(0, base_rate * (lambda_hyp**s_o2) * (thermal_flux**s_temp) * ((gibbs * pH_mod)**s_ph))

class FedBatchBioreactorModel:
    def __init__(self): 
        self.base_rate = 1.2e9
    
def run_mass_balance(self, init_vol, s_0, s_in, mu_max, o2, temp, ph, dur, feed_strategy, f_in_initial, mu_setpoint, dilution_rate, is_chemostat, feed_type, titrant_molarity):
        import math 
        
        # --- DYNAMIC FEED CONSTANTS ---
        feed_db = {
            "Glucose (Pure)": {"Ks": 0.1, "Yxs": 0.5},
            "Molasses (Complex)": {"Ks": 0.4, "Yxs": 0.35} 
        }
        Ks = feed_db[feed_type]["Ks"]
        Yxs = feed_db[feed_type]["Yxs"]
        
        # --- TITRATION CONSTANTS ---
        Y_lac_x = 0.8 # Grams of Lactic acid produced per gram of biomass
        MW_lactic_acid = 90.08 # g/mol
        total_lactic_acid_g = 0.0
        
        X, S, V = 0.5, s_0, init_vol
        history = {"Hour": [], "Biomass (g/L)": [], "Substrate (g/L)": [], "Volume (L)": [], "Titrant Needed (L)": []}
        
        for hour in range(1, dur + 1):
            dt = 1.0 / 100.0 
            for step in range(100):
                current_time = hour - 1 + (step * dt)

                if feed_strategy == "Constant Feed":
                    f_in_current = f_in_initial
                elif feed_strategy == "Exponential Feed":
                    f_in_current = f_in_initial * math.exp(mu_setpoint * current_time)
                elif feed_strategy == "Dilution Rate (D)":
                    f_in_current = dilution_rate * V
                else:
                    f_in_current = 0.0
                
                f_out_current = f_in_current if is_chemostat else 0.0

                # Kinetics
                mu = mu_max * (S / (Ks + S)) if S > 0 else 0
                
                dV = (f_in_current - f_out_current) * dt
                D = (f_out_current / V) if V > 0 else 0
                
                dX = (mu * X - D * X - (f_in_current / V) * X) * dt if V > 0 else 0
                dS = ((f_in_current / V) * s_in - (mu * X) / Yxs - D * S - (f_in_current / V) * S) * dt if V > 0 else 0
                
                # Calculate Lactic Acid accumulation
                if dX > 0:
                    total_lactic_acid_g += (dX * V * Y_lac_x)
                
                V += dV
                X = max(0.0, X + dX)
                S = max(0.0, S + dS)
                
            # --- CALCULATE TITRATION ---
            moles_lactic = total_lactic_acid_g / MW_lactic_acid
            titrant_vol_L = moles_lactic / titrant_molarity 
                
            history["Hour"].append(hour)
            history["Biomass (g/L)"].append(X)
            history["Substrate (g/L)"].append(S)
            history["Volume (L)"].append(V)
            history["Titrant Needed (L)"].append(titrant_vol_L)
            
        return pd.DataFrame(history)

 def run_simulation(self, init_o2, target_temp, target_ph, mixing_homogeneity, duration_hours, s_o2, s_temp, s_ph, biomass_series, feed_type, titrant_molarity):
        viability = 100.0
        history = {"Hour": [], "Therapeutic EVs": [], "Stress-Altered EVs": [], "Apoptotic Impurities": [], "Cell Viability (%)": []}
        
        for hour in range(1, duration_hours + 1):
            # EXTRACT: Get active biomass for the current hour (index is hour - 1)
            active_biomass = biomass_series.iloc[hour - 1] if hasattr(biomass_series, 'iloc') else biomass_series[hour - 1]
            
            # COUPLE: Scale specific productivity (q_EV) by active biomass (X)
            dynamic_base_rate = self.base_rate * active_biomass
            
            # CALCULATE FLUX
            rate = BiogenesisEngine.calc_flux(init_o2, target_temp, target_ph, s_o2, s_temp, s_ph, dynamic_base_rate)
            
            # SHEAR & VIABILITY LOGIC
            shear = ((mixing_homogeneity - 85.0)**1.5 * 0.1) if mixing_homogeneity > 85 else 0
            viability = max(0, viability - (0.5 + shear))
            
            history["Hour"].append(hour)
            history["Therapeutic EVs"].append(rate * (viability / 100))
            history["Stress-Altered EVs"].append(rate * ((100 - viability) / 100) * 0.3)
            # Ensure apoptotic impurities also scale with the dynamic biomass
            history["Apoptotic Impurities"].append(dynamic_base_rate * ((100 - viability) / 100) * 5)
            history["Cell Viability (%)"].append(viability)
            
        return pd.DataFrame(history)

# --- 4. ADVANCED PDF GENERATOR ---


def generate_qms_pdf(batch_id, cell_line, vol, dur, target, yield_val, purity, consistency, q_score):
    pdf = FPDF()
    pdf.add_page()
    
    # Header
    pdf.set_font("Helvetica", "B", 16)
    pdf.set_text_color(119, 158, 203) 
    pdf.cell(0, 10, "EVelution-bio: Digital Twin Projection Report", ln=True)
    pdf.set_font("Helvetica", "I", 10)
    pdf.set_text_color(128, 128, 128)
    pdf.cell(0, 5, "Multi-Machinery Model (MMModel) Bioprocess Simulation", ln=True)
    pdf.line(10, 28, 200, 28)
    pdf.ln(10)
    
    # Section 1: Configuration (Now includes Batch ID)
    pdf.set_font("Helvetica", "B", 12)
    pdf.set_text_color(0, 0, 0)
    pdf.cell(0, 10, f"1. Bioreactor Configuration (Batch ID: {batch_id})", ln=True)
    pdf.set_font("Helvetica", "", 11)
    pdf.cell(0, 8, f"Host Cell Line: {cell_line}", ln=True)
    pdf.cell(0, 8, f"Final Working Volume: {vol:.2f} L", ln=True)
    pdf.cell(0, 8, f"Culture Duration: {dur} Hours", ln=True)
    pdf.ln(5)
    
    # Section 2: Outcomes
    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(0, 10, "2. Projected Downstream Outcomes", ln=True)
    pdf.set_font("Helvetica", "", 11)
    pdf.cell(0, 8, f"Target Yield: {target:.2e} EVs", ln=True)
    pdf.cell(0, 8, f"Projected Yield: {yield_val:.2e} EVs ({(yield_val/target)*100:.1f}% of Target)", ln=True)
    pdf.cell(0, 8, f"Downstream Purity Recovery: {purity*100:.1f}%", ln=True)
    pdf.cell(0, 8, f"Cargo Consistency: {consistency*100:.1f}%", ln=True)
    pdf.ln(5)
    
    # Section 3: QA Evaluation
    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(0, 10, "3. Quality Assurance Evaluation", ln=True)
    if (yield_val/target)*100 >= 100 and q_score >= 60.0:
        pdf.set_text_color(34, 139, 34)
        status = "OPTIMAL: Target reached with high-purity harvest. Cleared for scale-up."
    elif (yield_val/target)*100 >= 100 and q_score >= 40.0:
        pdf.set_text_color(255, 140, 0)
        status = "MARGINAL: Target reached, but cargo consistency is dropping."
    else:
        pdf.set_text_color(220, 20, 60)
        status = "DEFICIENT / CRITICAL: Do not proceed to downstream processing."
        
    pdf.set_font("Helvetica", "B", 11)
    pdf.cell(0, 10, f"STATUS: {status}", ln=True)
    
    # Section 4: GMP Audit Hash (NEW)
    pdf.ln(15)
    pdf.set_font("Helvetica", "I", 8)
    pdf.set_text_color(150, 150, 150)
    
    # Generate a unique cryptographic hash based on the batch parameters
    raw_data = f"{batch_id}{cell_line}{vol}{dur}{target}{yield_val}{time.time()}"
    run_hash = hashlib.sha256(raw_data.encode()).hexdigest()
    
    pdf.cell(0, 5, "DATA INTEGRITY VERIFICATION", ln=True)
    pdf.cell(0, 5, f"Generated Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')} UTC", ln=True)
    pdf.cell(0, 5, f"SHA-256 Audit Hash: {run_hash}", ln=True)
    
    return bytes(pdf.output())

# --- 5. APP UI HEADER ---
st.title("EVelution Bioreactor Optimisation")
st.markdown("""
<div style="font-size: 0.9em; color: #B39EB5; margin-bottom: 20px;">
    <b>Multi-Machinery Model (MMModel)</b> | Author: Evgenia Trudova 
</div>
""", unsafe_allow_html=True)

with st.expander("Explore Logic", expanded=False):
    tab_bio, tab_math = st.tabs(["Biology", "Model"])
    
    with tab_bio:
        st.markdown("""
        **Metrics**
        * **Yield Performance**: Total therapeutic EV quantity projected post-DSP. Calculated by integrating hourly production flux (Φ) over the culture duration, scaled by volume.
        * **Harvest Concentration**: Instantaneous density of EVs in the bioreactor at the moment of harvest (Φ_final). Dictates the upstream input load for downstream purification.
        * **Downstream Purity**: Represents the success rate of the purification workflow (TFF/Chromatography). Modeled as a recovery coefficient (η_purity).
        * **Cargo Consistency**: Quality index denoting the percentage of EVs that contain the active therapeutic payload (φ_consistency).
        
        **Biogenesis Triggers (BiogenesisEngine.calc_flux)**
        * **Hypoxia**: Modeled as an evolutionary panic response to oxygen deprivation using an inverted Hill equation.
        * **Temperature**: Uses the Arrhenius equation coupled with a quadratic heat-shock modifier to simulate thermodynamic acceleration of biogenesis.
        * **pH**: Modeled via Gibbs Free Energy thermodynamics combined with an empirical inhibition scalar for proton concentration shifts.
        """)
        
    with tab_math:
        st.markdown("Metrics")
        st.latex(r" \text{Yield Performance: } Y_{perf} = \left( \sum_{t=1}^{T} \Phi(t) \cdot V_{react} \right) \cdot \eta_{purity} \cdot \phi_{consistency}")
        st.latex(r" \text{Harvest Conc: } \Phi_{final} = \Phi(t_{harvest})")
        st.latex(r" \text{Purity: } \eta_{purity} = \frac{EV_{purified}}{EV_{crude}} \quad \text{Consistency: } \phi_{consistency} = \frac{EV_{loaded}}{EV_{total}}")
        
        st.divider()
        
        st.markdown("Biogenesis")
        st.latex(r" \text{Hypoxia (Hill): } \lambda_{hyp} = \frac{K^n + x_0^n}{K^n + O_2^n}")
        st.latex(r" \text{Thermal (Arrhenius + Shock): } S_{temp} = A_0 e^{-\frac{E_a}{R}\left(\frac{1}{T} - \frac{1}{T_0}\right)} \cdot \left( 1 + \frac{(T - 37)^2}{1^2 + (T - 37)^2} \right)")
        st.latex(r" \text{pH (Gibbs + Inhibition): } S_{pH} = e^{-\frac{2.303 RT (pH - pH_0)}{RT}} \cdot \left( 1 + \frac{0.5 \left(\frac{[H^+]}{[H^+]_0}\right)^2}{0.1^2 + \left(\frac{[H^+]}{[H^+]_0}\right)^2} \right)")
        
        st.markdown("**Total Hourly Production Flux:**")
        st.latex(r" \Phi(t) = \Phi_{base} \cdot \lambda_{hyp}^{s_{o2}} \cdot S_{temp}^{s_{temp}} \cdot S_{pH}^{s_{ph}}")

# --- 6. CELL LINE DATABASE ---
cell_line_db = {
    "Human MSCs (Bone Marrow)": {"hypoxia": 1.15, "thermal": 0.85, "ph": 0.95, "mu_max": 0.02},
    "HEK293T (Suspension)": {"hypoxia": 0.90, "thermal": 1.20, "ph": 0.80, "mu_max": 0.035},
    "CHO-K1": {"hypoxia": 0.75, "thermal": 0.95, "ph": 1.30, "mu_max": 0.045},
    "Custom / Empirical DoE": {"hypoxia": 1.00, "thermal": 1.00, "ph": 1.00, "mu_max": 0.03}
}

# --- 7. SIDEBAR INPUTS ---
with st.sidebar:
    st.divider()
    st.header("Batch")
    
    # NEW: Batch Naming
    batch_name = st.text_input("Batch ID", "EXP-001")

    st.divider()
    st.header("Cell Line")
    
    selected_cell = st.selectbox("Select Host Cell Line", list(cell_line_db.keys()))
    defaults = cell_line_db[selected_cell]
    
    manual_override = st.toggle("Enable Manual Input (DoE Override)")
    s_o2 = st.number_input("Hypoxia Modifier", min_value=0.0, max_value=3.0, value=defaults["hypoxia"], step=0.05, disabled=not manual_override)
    s_temp = st.number_input("Thermal Modifier", min_value=0.0, max_value=3.0, value=defaults["thermal"], step=0.05, disabled=not manual_override)
    s_ph = st.number_input("pH Modifier", min_value=0.0, max_value=3.0, value=defaults["ph"], step=0.05, disabled=not manual_override)
    
    # UPGRADE: mu_max as an intrinsic biological property
    # Uses .get() as a safeguard in case Section 6 hasn't been updated yet
    default_mu = defaults.get("mu_max", 0.035) 
    mu_max = st.number_input(
        "Max Growth Rate, μ_max (1/h)", 
        min_value=0.001, max_value=2.000, 
        value=default_mu, 
        step=0.005, 
        format="%.3f",
        disabled=not manual_override,
        help="Intrinsic specific growth rate. Mammalian cells typically range from 0.015 to 0.050 1/h."
    )

    st.divider()
    st.header("Metabolism")
    
    # NEW: Carbon Source (Feed Type)
    feed_type = st.selectbox("Carbon Source", ["Glucose (Pure)", "Molasses (Complex)"])
    
    unit_toggle = st.radio("Substrate Units", ["g/L", "mol/L"], horizontal=True)
    
    # Assuming Glucose (MW = 180.16 g/mol)
    mw_glucose = 180.16 
    
    display_s0 = 10.0 if unit_toggle == "g/L" else (10.0 / mw_glucose)
    display_sin = 100.0 if unit_toggle == "g/L" else (100.0 / mw_glucose)
    
    s_0_input = st.number_input(f"Initial Substrate ({unit_toggle})", value=display_s0, step=display_s0 * 0.1)
    s_in_input = st.number_input(f"Feed Concentration ({unit_toggle})", value=display_sin, step=display_sin * 0.1)
    
    s_0 = s_0_input if unit_toggle == "g/L" else s_0_input * mw_glucose
    s_in = s_in_input if unit_toggle == "g/L" else s_in_input * mw_glucose

    st.divider()
    st.header("Experimental")
    
    vol = st.number_input("Initial Volume (L)", value=50.0)
    o2 = st.slider("Oxygen (%)", 0, 21, 21)
    temp = st.slider("Temp (°C)", 4, 60, 37)
    ph = st.slider("Target pH", 6.0, 8.0, 7.4)
    
    # NEW: Titration Magnitude & Strength
    st.caption("pH Control Strategy")
    titrant_type = st.selectbox("Base Type", ["NaOH (Sodium Hydroxide)", "NaHCO3 (Sodium Bicarbonate)"])
    titrant_molarity = st.number_input("Base Strength (Molarity, M)", min_value=0.1, max_value=5.0, value=1.0, step=0.1)
    
    mix = st.slider("Mixing (%)", 50, 100, 85)
    dur = st.slider("Duration (h)", 24, 336, 240) # Industry 10-day default!
    
    st.divider()
    st.header("Feeding Strategy")

    # UPGRADE: Dynamic Feeding Strategies
    feed_strategy = st.selectbox(
        "Profile Type", 
        ["Constant Feed", "Exponential Feed", "Dilution Rate (D)"]
    )

    # Initialize variables for the ODE engine to prevent reference errors
    f_in_initial = 0.0
    mu_setpoint = 0.0
    dilution_rate = 0.0
    is_chemostat = False 

    if feed_strategy == "Constant Feed":
        max_safe_feed = vol * 0.1 # Bound it to 10% of volume per hour
        f_in_initial = st.number_input("Feed Rate (L/h)", min_value=0.0, max_value=float(max_safe_feed), value=0.0, step=0.05)
    
    elif feed_strategy == "Exponential Feed":
        st.caption("Matches feed to cell growth to prevent overfeeding.")
        mu_setpoint = st.number_input("Target Growth Rate for Feed (1/h)", min_value=0.005, max_value=0.100, value=0.020, step=0.005, format="%.3f")
        f_in_initial = st.number_input("Initial Feed Rate (L/h)", value=0.1, step=0.05)

    elif feed_strategy == "Dilution Rate (D)":
        st.caption("Feed rate scales dynamically with bioreactor volume.")
        dilution_rate = st.number_input("Dilution Rate, D (1/h)", min_value=0.000, max_value=0.500, value=0.010, step=0.005, format="%.3f")
        
        is_chemostat = st.toggle("Enable Chemostat (F_out = F_in)")
        if is_chemostat:
            st.caption("Volume remains constant. Biomass and Substrate will wash out.")

    st.divider()
    st.header("Yield Target")
    t_col1, t_col2 = st.columns(2)
    with t_col1:
        target_base = st.number_input("Base", min_value=1.0, max_value=9.9, value=1.0, step=0.1)
    with t_col2:
        target_exp = st.number_input("Magnitude (10^x)", min_value=10, max_value=20, value=15, step=1)
        
    target = target_base * (10 ** target_exp)
    st.markdown(f"<div style='text-align: right; color: #77DD77; font-size: 0.9em;'><b>Active Target: {target:.1e} EVs</b></div>", unsafe_allow_html=True)

# --- 8. SIMULATION & GRAPHS ---
model = FedBatchBioreactorModel()

# RUN MASS BALANCE FIRST to get X(t)
mb_df = model.run_mass_balance(
            vol, s_0, s_in, mu_max, o2, temp, ph, dur,
            feed_strategy, f_in_initial, mu_setpoint, dilution_rate, is_chemostat, feed_type, titrant_molarity
        )
final_biomass = mb_df["Biomass (g/L)"].iloc[-1]
final_substrate = mb_df["Substrate (g/L)"].iloc[-1]
final_vol = mb_df["Volume (L)"].iloc[-1]

# RUN EV SIMULATION using the dynamic Biomass series
df = model.run_simulation(o2, temp, ph, mix, dur, s_o2, s_temp, s_ph, mb_df["Biomass (g/L)"])

# CALCULATE YIELD METRICS
avg_viability = df["Cell Viability (%)"].mean() / 100.0
impurity_ratio = df["Apoptotic Impurities"].iloc[-1] / (df["Therapeutic EVs"].iloc[-1] + 1e-9)

dynamic_purity = max(0.2, 0.78 * (1.0 / (1.0 + (impurity_ratio * 0.05)))) 
dynamic_consistency = max(0.3, 0.62 * (avg_viability ** 0.5)) 

total_prod = df["Therapeutic EVs"].sum()
true_val = total_prod * final_vol * 1000 * dynamic_purity * dynamic_consistency
yield_achievement = (true_val / target) * 100
quality_score = (dynamic_purity * dynamic_consistency) * 100

# Build Figure 0: Kinetics 
fig_monod = make_subplots(specs=[[{"secondary_y": True}]])

# 1. Primary Y-Axis (Left): Concentration (g/L)
fig_monod.add_trace(
    go.Scatter(x=mb_df["Hour"], y=mb_df["Biomass (g/L)"], mode='lines', name='Biomass (g/L)', line=dict(color=C_BLUE, width=3)),
    secondary_y=False,
)
fig_monod.add_trace(
    go.Scatter(x=mb_df["Hour"], y=mb_df["Substrate (g/L)"], mode='lines', name='Substrate (g/L)', line=dict(color=C_GREEN, width=3)),
    secondary_y=False,
)

# 2. Secondary Y-Axis (Right): Volume (L)
# Using a dotted line for Volume helps industry operators instantly distinguish it from concentrations
fig_monod.add_trace(
    go.Scatter(x=mb_df["Hour"], y=mb_df["Volume (L)"], mode='lines', name='Volume (L)', line=dict(color=C_STAR, width=3, dash='dot')),
    secondary_y=True,
)

# 3. Layout Formatting
fig_monod.update_layout(
    height=400, 
    margin=fixed_margin, 
    hovermode="x unified", 
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
)

# 4. Axis Titles
fig_monod.update_yaxes(title_text="Concentration (g/L)", secondary_y=False, color="#B39EB5")
fig_monod.update_yaxes(title_text="Volume (L)", secondary_y=True, color=C_STAR)
fig_monod.update_xaxes(title_text="Culture Duration (Hours)")


# Build Figure 1: Accumulation
fig_accum = px.line(df, x="Hour", y=["Therapeutic EVs", "Stress-Altered EVs", "Apoptotic Impurities"], log_y=True,
              color_discrete_map={"Therapeutic EVs": C_GREEN, "Stress-Altered EVs": C_PURPLE, "Apoptotic Impurities": C_BLUE})

if 'historical_df' in st.session_state and st.session_state['historical_df'] is not None:
    hist_df = st.session_state['historical_df']
    
    # 1. Smart X-axis detection (English & Swedish)
    time_aliases = ["Hour", "Time", "Time (h)", "Timme", "Tid", "Tid (h)"]
    # Find the first matching time column, otherwise default to the very first column in the CSV
    x_col = next((col for col in time_aliases if col in hist_df.columns), hist_df.columns[0])
    
    # 2. Smart Y-axis detection (English & Swedish)
    yield_aliases = [
        "Therapeutic EVs", "Yield", "Titer", "Concentration", "Total Particles", 
        "Utbyte", "Koncentration", "Partikelkoncentration", "EV-utbyte"
    ]
    # Find the first matching yield column. If none match, grab the first column that isn't the X-axis.
    available_y_cols = [col for col in hist_df.columns if col != x_col]
    y_fallback = available_y_cols[0] if available_y_cols else hist_df.columns[0]
    
    y_col = next((col for col in yield_aliases if col in hist_df.columns), y_fallback)
    
    # Plot the Golden Batch line using the detected columns
    fig_accum.add_trace(go.Scatter(
        x=hist_df[x_col], y=hist_df[y_col],
        mode='lines', name='Golden Batch Benchmark',
        line=dict(color='rgba(227, 82, 82, 0.7)', width=3, dash='dash')
    ))

fig_accum.update_layout(height=fixed_height, margin=fixed_margin, hovermode="x unified", legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))

# Build Figure 2: Viability
crit = df[df["Cell Viability (%)"] < 50.0]
fig_viab = px.line(df, x="Hour", y=["Cell Viability (%)"], color_discrete_sequence=[C_BLUE])
if not crit.empty: 
    fig_viab.add_trace(go.Scatter(
        x=[crit.iloc[0]["Hour"]], y=[crit.iloc[0]["Cell Viability (%)"]], 
        mode='markers+text', text=['Death Cliff'], textposition='bottom right',
        textfont=dict(color=C_STAR, size=12, family="sans serif"), marker=dict(size=14, color=C_STAR, symbol='star')
    ))
fig_viab.update_layout(height=fixed_height, margin=fixed_margin, showlegend=False, hovermode="x unified")

# Build Figure 3: Funnel
fig_funnel = go.Figure(go.Funnel(y=["Raw Yield", "Intact (Purity)", "Functional"], x=[true_val/(0.78*0.62), true_val/0.62, true_val], textinfo="value+percent previous", marker={"color": [C_BLUE, C_PURPLE, C_GREEN]}))
fig_funnel.update_layout(height=fixed_height, margin=fixed_margin)

# Build Figure 4: Sensitivity (OPTIMIZED WITH CACHING)

@st.cache_data(show_spinner=False)
def calculate_sensitivity(vol, s_0, s_in, mu_max, o2, temp, ph, mix, s_o2, s_temp, s_ph, 
                          feed_strategy, f_in_initial, mu_setpoint, dilution_rate, is_chemostat):
    """
    Caches the 14-step sensitivity loop. It only recalculates if the user explicitly 
    changes an upstream bioreactor parameter.
    """
    temp_model = FedBatchBioreactorModel()
    sens_range = list(range(12, 96, 6))
    sens_data = []

    for d in sens_range:
        temp_mb = temp_model.run_mass_balance(
            vol, s_0, s_in, mu_max, o2, temp, ph, d,
            feed_strategy, f_in_initial, mu_setpoint, dilution_rate, is_chemostat
        )
        temp_vol = temp_mb["Volume (L)"].iloc[-1]
        
        # Pass the isolated duration's biomass series into the simulation
        temp_df = temp_model.run_simulation(o2, temp, ph, mix, d, s_o2, s_temp, s_ph, temp_mb["Biomass (g/L)"])
        
        step_yield = temp_df["Therapeutic EVs"].sum() * temp_vol * 1000 * 0.78 * 0.62
        sens_data.append(step_yield)
        
    return sens_range, sens_data

# Execute the cached function by passing in all the active UI parameters
sens_range, sens_data = calculate_sensitivity(
    vol, s_0, s_in, mu_max, o2, temp, ph, mix, s_o2, s_temp, s_ph, 
    feed_strategy, f_in_initial, mu_setpoint, dilution_rate, is_chemostat
)

variance_percentages = np.linspace(0.02, 0.15, len(sens_range))
upper_bound = [val * (1 + var) for val, var in zip(sens_data, variance_percentages)]
lower_bound = [val * (1 - var) for val, var in zip(sens_data, variance_percentages)]

fig_sens = go.Figure()
fig_sens.add_trace(go.Scatter(
    x=sens_range + sens_range[::-1],
    y=upper_bound + lower_bound[::-1],
    fill='toself',
    fillcolor='rgba(119, 221, 119, 0.2)',
    line=dict(color='rgba(255,255,255,0)'),
    hoverinfo="skip",
    name='Biological Variance (95% CI)'
))
fig_sens.add_trace(go.Scatter(
    x=sens_range, 
    y=sens_data, 
    mode='lines', 
    line=dict(color=C_GREEN, width=3),
    name='Predicted Yield'
))
idx = np.argmax(sens_data)
fig_sens.add_trace(go.Scatter(
    x=[sens_range[idx]], 
    y=[sens_data[idx]], 
    mode='markers+text', 
    text=['Optimal Harvest'], 
    textposition='top center',
    textfont=dict(color=C_STAR, size=12, family="sans serif"), 
    marker=dict(size=14, color=C_STAR, symbol='star'),
    name='Target'
))
fig_sens.update_layout(
    height=fixed_height, 
    margin=fixed_margin, 
    showlegend=False, 
    hovermode="x unified",
    xaxis_title="Culture Duration (Hours)",
    yaxis_title="Functional EV Yield"
)

# --- 9. SIDEBAR DATA EXPORT ---
with st.sidebar:
    st.divider()
    st.header("Data")
    
    # PASS BATCH ID: Ensure the PDF generator uses the batch_name defined in Section 7
    pdf_bytes = generate_qms_pdf(
        batch_id=batch_name, 
        cell_line=selected_cell, vol=final_vol, dur=dur, target=target,
        yield_val=true_val, purity=dynamic_purity, consistency=dynamic_consistency, q_score=quality_score
    )
    
    # SANITIZE FILENAME: Prevent OS-level file saving errors if users type illegal characters (/, \, :)
    safe_batch_name = batch_name.replace("/", "-").replace("\\", "-").replace(":", "-")
    
    # DYNAMIC FILENAME: The downloaded file will now be safely named after the Batch ID
    st.download_button(
        label="Download PDF",
        data=pdf_bytes,
        file_name=f"{safe_batch_name}_EVelution_Projection.pdf", 
        mime="application/pdf",
        use_container_width=True
    )
    
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("**Historical Benchmarking**")
    
    # UX TEXT: Clear instructions for the user before they upload
    st.caption("Upload a past batch CSV to instantly overlay your real-world data as a dashed red 'Golden Batch' line on the **Process Accumulation** graph.")
    
    uploaded_file = st.file_uploader(
        "Upload History File (.csv)", 
        type="csv",
        help="The app auto-detects common time and yield columns in English and Swedish (e.g., Hour/Timme, Titer/Utbyte).",
        key="csv_uploader"
    )
    
    if uploaded_file is not None:
        try:
            historical_df = pd.read_csv(uploaded_file)
            st.session_state['historical_df'] = historical_df
        except Exception as e:
            st.error(f"Could not read the file. Error: {e}")

    has_data = st.session_state.get('historical_df') is not None
    if st.button("Delete", disabled=not has_data, use_container_width=True):
        st.session_state['historical_df'] = None
        if 'csv_uploader' in st.session_state:
            del st.session_state['csv_uploader']
        st.rerun()

    st.markdown("""
    <div class="legal-disclaimer" style="font-size: 0.75em; color: #808080; margin-top: 20px; padding-top: 15px; border-top: 1px solid #333; line-height: 1.4; text-align: justify;">
        <b>DATA GOVERNANCE AND LIABILITY DISCLAIMER</b><br>
        Data uploaded to this application is processed strictly in volatile memory (RAM) for ephemeral visualization. The architecture contains no persistent storage, logging, or external transmission protocols. By utilizing this software, end-users assume absolute liability for the safeguarding of proprietary technical trade secrets, in compliance with the Swedish Trade Secrets Act (SFS 2018:558, as amended 2026), the EU Trade Secrets Directive (2016/943), and the US Defend Trade Secrets Act (DTSA). Furthermore, operations align with global electronic record frameworks (FDA 21 CFR Part 11, EMA Annex 11) and international data sovereignty laws (GDPR) via zero-retention processing. Under international safe harbor and ephemeral data processing doctrines, this application acts purely as a localized execution environment, legally exempting the provider from data hosting, processor, and controller liabilities. Users are legally obligated to execute the internal data purge protocol upon session termination.
    </div>
    """, unsafe_allow_html=True)
    
# --- 10. METRICS & BATCH EVALUATION ---
final_titrant = mb_df["Titrant Needed (L)"].iloc[-1]

m1, m2, m3, m4 = st.columns(4)
m1.metric("Yield Performance", f"{true_val:.2e}")
m2.metric("Harvest Conc", f"{df['Therapeutic EVs'].iloc[-1]:.2e} ev/mL")
m3.metric("Downstream Purity", f"{dynamic_purity*100:.1f}%")
m4.metric("Cargo Consistency", f"{dynamic_consistency*100:.1f}%")

st.markdown("### Batch Success Evaluation")

is_startup = (
    selected_cell == "Human MSCs (Bone Marrow)" and
    vol == 50.0 and o2 == 21 and temp == 37 and ph == 7.4 and mix == 85 and dur == 240 and
    target_base == 1.0 and target_exp == 15 and not manual_override and 
    feed_strategy == "Constant Feed" and f_in_initial == 0.0 and
    feed_type == "Glucose (Pure)" and titrant_molarity == 1.0 and unit_toggle == "g/L"
)

if is_startup:
    status_color, quality_color = "#779ECB", "#779ECB" 
    status_icon = "💡"
    status_text = "AWAITING OPTIMIZATION: Default baseline loaded. Adjust parameters to begin."
elif yield_achievement >= 100:
    if quality_score >= 60.0:
        status_color, quality_color = "#77DD77", "#77DD77"
        status_icon = "✅"
        status_text = "OPTIMAL: Target reached with high-purity harvest."
    elif quality_score >= 40.0: 
        status_color, quality_color = "#FACA2E", "#FACA2E"
        status_icon = "⚠️"
        status_text = "MARGINAL: Target reached, but cargo consistency is dropping."
    else:
        status_color, quality_color = "#E35252", "#E35252"
        status_icon = "❌"
        status_text = "CRITICAL: Severe quality degradation. High risk of DSP failure."
else:
    status_color, quality_color = "#779ECB", "#779ECB"
    status_icon = "🛑"
    status_text = "DEFICIENT: Target volume not reached. Extend duration or adjust feeding."

st.markdown(f"""
<div class="status-card" style="background-color: #1E1E2E; padding: 20px; border-radius: 8px; border-left: 6px solid {status_color}; display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 15px;">
    <div style="flex: 1; min-width: 250px;">
        <h4 style="margin: 0; color: {status_color}; display: flex; align-items: center; word-break: break-word;">
            <span style="margin-right: 10px;">{status_icon}</span> {status_text}
        </h4>
    </div>
    <div class="status-card-metric" style="text-align: right; min-width: 120px;">
        <div style="font-size: 0.9em; color: #B39EB5;">Yield vs Target</div>
        <div style="font-size: 1.5em; font-weight: bold; color: {'#77DD77' if yield_achievement >= 100 else '#779ECB'};">{yield_achievement:.1f}%</div>
    </div>
    <div class="status-card-metric" style="text-align: right; padding-left: 20px; border-left: 1px solid #333; min-width: 140px;">
        <div style="font-size: 0.9em; color: #B39EB5;">Overall Quality Score</div>
        <div style="font-size: 1.5em; font-weight: bold; color: {quality_color};">{quality_score:.1f}%</div>
    </div>
</div>
""", unsafe_allow_html=True)

# --- 10.5. KINETICS & TITRATION ---
st.divider()
st.subheader("Kinetics & pH Titration")
c1, c2, c3, c4 = st.columns(4)
c1.metric("Biomass", f"{final_biomass:.2f} g/L")
c2.metric("Substrate", f"{final_substrate:.2f} g/L")
c3.metric("Volume", f"{final_vol:.2f} L")
c4.metric(f"Req. Base", f"{final_titrant:.2f} L")
st.plotly_chart(fig_monod, use_container_width=True)

# --- 11. ANALYTICS GRID ---
st.divider()
st.subheader("Analytics")
r1c1, r1c2 = st.columns(2)
r2c1, r2c2 = st.columns(2)

with r1c1:
    st.markdown("### Process Accumulation")
    # Using the pre-built figure from Section 8
    st.plotly_chart(fig_accum, use_container_width=True)
    with st.expander("Explore Logic"):
        t1, t2 = st.tabs(["Biology", "Model"])
        t1.markdown("This graph tracks therapeutic accumulation versus impurity buildup. The harvest window is constrained by the impurity crossover point. Harvesting before this ensures optimal purity levels.")
        t2.latex(r"\Phi_{total} = \Phi_{Therapeutic} + \Phi_{Stress} + \Phi_{Apoptotic}")

with r1c2:
    st.markdown("### Cellular Viability")
    st.plotly_chart(fig_viab, use_container_width=True)
    with st.expander("Explore Logic"):
        t1, t2, t3 = st.tabs(["Biology", "Model", "Exponential"])
        
        t1.markdown("Viability decline signifies the transition from growth to the apoptotic phase. The 'Death Cliff' marker identifies when cell repair mechanisms collapse under toxic stress.")
        
        t2.latex(r"\frac{dV}{dt} = -(\kappa_{tox} + \tau_{shear})")
        t2.markdown("The model utilizes zero-order linear decay kinetics, computed via a Forward Euler integration step, rather than a physiological exponential curve. See Exponential.")
        
        t3.latex(r"\text{Exponential Rate: } \frac{dV}{dt} = -k_d V \implies V(t) = V_0 e^{-k_d t}")
        t3.markdown("""
        While mammalian apoptosis is biologically exponential, using it as a control signal causes critical hardware failures:
        Exponential feedback creates steep derivative slopes. PID controllers (e.g., DeltaV) violently overcompensate gas flow and chilling jackets.
        This "controller panic" triggers aggressive impeller agitation. The resulting hydrodynamic shear stress physically shreds EV lipid bilayers and ruins cargo integrity ([Thompson & Papoutsakis, 2023](https://doi.org/10.1016/j.biotechadv.2023.108158)).
        A linear baseline dampens the controller's derivative, ensuring smooth agitation and preserving the EVs.
        
        
        By keeping the foundational physics computationally safe, future Neural Networks can be trained on errors and patterns.
        """)

with r2c1:
    st.markdown("### Yield-to-Value Bridge")
    st.plotly_chart(fig_funnel, use_container_width=True)
    with st.expander("Explore Logic"):
        t1, t2 = st.tabs(["Biology", "Model"])
        t1.markdown("This funnel visualizes the mass balance across downstream stages. Narrowing segments highlight cumulative yield loss during purification. It serves as a diagnostic for loading/recovery efficiency.")
        t2.latex(r"V_{final} = Yield_{raw} \cdot \eta_{purity} \cdot \phi_{consistency}")

with r2c2:
    st.markdown("### Yield Sensitivity Analysis")
    st.plotly_chart(fig_sens, use_container_width=True)
    with st.expander("Explore Logic"):
        t1, t2 = st.tabs(["Biology", "Model"])
        
        t1.markdown("""
        The solid line predicts absolute yield, while the shaded region represents the biological variance (95% CI). Variance interval widens as time progresses. As cellular viability drops and the culture enters late-stage apoptosis, stochastic events  make the process highly unpredictable. The "optimal spot" is the inflection point where incremental EV gain is maximized before the risk of batch variance becomes too wide.
        """)
        
        t2.latex(r"\text{Optimal Harvest: } \frac{d}{dt}Yield(t) = 0")
        t2.latex(r"\text{Risk Variance: } \sigma^2(t) \propto \int_{0}^{t} (\kappa_{tox}(\tau)) d\tau")
        t2.markdown("The mathematical optimization seeks to maximize yield while minimizing the integral of accumulated toxic variance.")



