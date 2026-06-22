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
import plotly.graph_objects as go
import plotly.express as px
from fpdf import FPDF
import io

# --- 1. THEME & STYLING ---
st.set_page_config(page_title="EVelution Bioreactor Twin", layout="wide")

# --- 2. CSS INJECTION ---
# Keep all CSS here, strictly separate from your Python logic
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
</style>
""", unsafe_allow_html=True)
    
# --- 2. GLOBAL CONSTANTS ---
fixed_height = 400
fixed_margin = dict(t=40, b=10, l=20, r=20)
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
    
    def run_simulation(self, init_o2, target_temp, target_ph, mixing_homogeneity, duration_hours, s_o2, s_temp, s_ph):
        viability = 100.0
        history = {"Hour": [], "Therapeutic EVs": [], "Stress-Altered EVs": [], "Apoptotic Impurities": [], "Cell Viability (%)": []}
        for hour in range(1, duration_hours + 1):
            rate = BiogenesisEngine.calc_flux(init_o2, target_temp, target_ph, s_o2, s_temp, s_ph, self.base_rate)
            shear = ((mixing_homogeneity - 85.0)**1.5 * 0.1) if mixing_homogeneity > 85 else 0
            viability = max(0, viability - (0.5 + shear))
            history["Hour"].append(hour)
            history["Therapeutic EVs"].append(rate * (viability / 100))
            history["Stress-Altered EVs"].append(rate * ((100 - viability) / 100) * 0.3)
            history["Apoptotic Impurities"].append(self.base_rate * ((100 - viability) / 100) * 5)
            history["Cell Viability (%)"].append(viability)
        return pd.DataFrame(history)

# --- 4. ADVANCED PDF GENERATOR (Text & Data Only) ---
def generate_qms_pdf(cell_line, vol, dur, target, yield_val, purity, consistency, q_score):
    pdf = FPDF()
    pdf.add_page()
    
    # --- PAGE 1: EXECUTIVE SUMMARY ---
    pdf.set_font("Helvetica", "B", 16)
    pdf.set_text_color(119, 158, 203) 
    pdf.cell(0, 10, "EVelution-bio: Digital Twin Projection Report", ln=True)
    pdf.set_font("Helvetica", "I", 10)
    pdf.set_text_color(128, 128, 128)
    pdf.cell(0, 5, "Multi-Machinery Model (MMModel) Bioprocess Simulation", ln=True)
    pdf.line(10, 28, 200, 28)
    pdf.ln(10)
    
    pdf.set_font("Helvetica", "B", 12)
    pdf.set_text_color(0, 0, 0)
    pdf.cell(0, 10, "1. Bioreactor Configuration", ln=True)
    pdf.set_font("Helvetica", "", 11)
    pdf.cell(0, 8, f"Host Cell Line: {cell_line}", ln=True)
    pdf.cell(0, 8, f"Working Volume: {vol} L", ln=True)
    pdf.cell(0, 8, f"Culture Duration: {dur} Hours", ln=True)
    pdf.ln(5)
    
    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(0, 10, "2. Projected Downstream Outcomes", ln=True)
    pdf.set_font("Helvetica", "", 11)
    pdf.cell(0, 8, f"Target Yield: {target:.2e} EVs", ln=True)
    pdf.cell(0, 8, f"Projected Yield: {yield_val:.2e} EVs ({(yield_val/target)*100:.1f}% of Target)", ln=True)
    pdf.cell(0, 8, f"Downstream Purity Recovery: {purity*100:.1f}%", ln=True)
    pdf.cell(0, 8, f"Cargo Consistency: {consistency*100:.1f}%", ln=True)
    pdf.ln(5)
    
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
        * **Yield Performance**: Total therapeutic EV quantity projected post-DSP. Calculated by integrating hourly production flux ($\Phi$) over the culture duration, scaled by volume.
        * **Harvest Concentration**: Instantaneous density of EVs in the bioreactor at the moment of harvest ($\Phi_{final}$). Dictates the upstream input load for downstream purification.
        * **Downstream Purity**: Represents the success rate of the purification workflow (TFF/Chromatography). Modeled as a recovery coefficient ($\eta_{purity}$).
        * **Cargo Consistency**: Quality index denoting the percentage of EVs that contain the active therapeutic payload ($\phi_{consistency}$).
        """)
    with tab_math:
        st.latex(r" \text{ Yield Performance: } Y_{perf} = \left( \sum_{t=1}^{T} \Phi(t) \cdot V_{react} \right) \cdot \eta_{purity} \cdot \phi_{consistency}")
        st.latex(r" \text{ Harvest Conc: } \Phi_{final} = \Phi(t_{harvest})")
        st.latex(r" \text{ Purity: } \eta_{purity} = \frac{EV_{purified}}{EV_{crude}}")
        st.latex(r" \text{ Consistency: } \phi_{consistency} = \frac{EV_{loaded}}{EV_{total}}")


# --- 6. CELL LINE DATABASE ---
cell_line_db = {
    "Human MSCs (Bone Marrow)": {"hypoxia": 1.15, "thermal": 0.85, "ph": 0.95},
    "HEK293T (Suspension)": {"hypoxia": 0.90, "thermal": 1.20, "ph": 0.80},
    "CHO-K1": {"hypoxia": 0.75, "thermal": 0.95, "ph": 1.30},
    "Custom / Empirical DoE": {"hypoxia": 1.00, "thermal": 1.00, "ph": 1.00}
}

# --- 7. SIDEBAR INPUTS ---
with st.sidebar:
    st.header("Cell Line")
    selected_cell = st.selectbox("Select Host Cell Line", list(cell_line_db.keys()), help="Loads validated empirical kinetic parameters.")
    defaults = cell_line_db[selected_cell]
    
    manual_override = st.toggle("Enable Manual Input (DoE Override)")
    s_o2 = st.number_input("Hypoxia Modifier (Ks scalar)", min_value=0.0, max_value=3.0, value=defaults["hypoxia"], step=0.05, disabled=not manual_override)
    s_temp = st.number_input("Thermal Modifier (Ea scalar)", min_value=0.0, max_value=3.0, value=defaults["thermal"], step=0.05, disabled=not manual_override)
    s_ph = st.number_input("pH Modifier (Inhibition scalar)", min_value=0.0, max_value=3.0, value=defaults["ph"], step=0.05, disabled=not manual_override)
    
    st.divider()
    st.header("Experimental")
    vol = st.number_input("Volume (L)", value=50.0)
    o2 = st.slider("Oxygen (%)", 0, 21, 21)
    temp = st.slider("Temp (°C)", 30, 45, 37)
    ph = st.slider("pH", 6.0, 8.0, 7.4)
    mix = st.slider("Mixing (%)", 50, 100, 85)
    dur = st.slider("Duration (h)", 12, 72, 48)
    
    st.divider()
    st.header("Yield Target")
    st.caption("Set target load for downstream processing:")
    t_col1, t_col2 = st.columns(2)
    
    with t_col1:
        target_base = st.number_input("Base", min_value=1.0, max_value=9.9, value=1.0, step=0.1)
    with t_col2:
        target_exp = st.number_input("Magnitude (10^x)", min_value=10, max_value=20, value=15, step=1)
        
    target = target_base * (10 ** target_exp)
    st.markdown(f"<div style='text-align: right; color: #77DD77; font-size: 0.9em;'><b>Active Target: {target:.1e} EVs</b></div>", unsafe_allow_html=True)


# --- 8. RUN SIMULATION & GENERATE GRAPHS ---
model = FedBatchBioreactorModel()
df = model.run_simulation(o2, temp, ph, mix, dur, s_o2, s_temp, s_ph)

avg_viability = df["Cell Viability (%)"].mean() / 100.0
impurity_ratio = df["Apoptotic Impurities"].iloc[-1] / (df["Therapeutic EVs"].iloc[-1] + 1e-9)

dynamic_purity = max(0.2, 0.78 * (1.0 / (1.0 + (impurity_ratio * 0.05)))) 
dynamic_consistency = max(0.3, 0.62 * (avg_viability ** 0.5)) 

total_prod = df["Therapeutic EVs"].sum()
true_val = total_prod * vol * 1000 * dynamic_purity * dynamic_consistency
yield_achievement = (true_val / target) * 100
quality_score = (dynamic_purity * dynamic_consistency) * 100 

# Build Figure 1: Accumulation
fig_accum = px.line(df, x="Hour", y=["Therapeutic EVs", "Stress-Altered EVs", "Apoptotic Impurities"], log_y=True,
              color_discrete_map={"Therapeutic EVs": C_GREEN, "Stress-Altered EVs": C_PURPLE, "Apoptotic Impurities": C_BLUE})
if 'historical_df' in st.session_state and st.session_state['historical_df'] is not None:
    hist_df = st.session_state['historical_df']
    if "Hour" in hist_df.columns and "Therapeutic EVs" in hist_df.columns:
        fig_accum.add_trace(go.Scatter(
            x=hist_df["Hour"], y=hist_df["Therapeutic EVs"],
            mode='lines', name='Golden Batch Benchmark',
            line=dict(color='rgba(227, 82, 82, 0.7)', width=3, dash='dash')
        ))
fig_accum.update_layout(height=fixed_height, margin=fixed_margin, hovermode="x unified")

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

# Build Figure 4: Sensitivity
sens_range = range(12, 96, 6)
sens_data = [model.run_simulation(o2, temp, ph, mix, d, s_o2, s_temp, s_ph)["Therapeutic EVs"].sum() * vol * 1000 * 0.78 * 0.62 for d in sens_range]
fig_sens = px.line(x=sens_range, y=sens_data, color_discrete_sequence=[C_GREEN])
idx = np.argmax(sens_data)
fig_sens.add_trace(go.Scatter(
    x=[list(sens_range)[idx]], y=[sens_data[idx]], 
    mode='markers+text', text=['Optimal Harvest'], textposition='top right',
    textfont=dict(color=C_STAR, size=12, family="sans serif"), marker=dict(size=14, color=C_STAR, symbol='star')
))
fig_sens.update_layout(height=fixed_height, margin=fixed_margin, showlegend=False, hovermode="x unified")

# Dictionary to pass all graphs easily
report_figs = {'accum': fig_accum, 'viab': fig_viab, 'funnel': fig_funnel, 'sens': fig_sens}


# --- 9. SIDEBAR DATA EXPORT ---
with st.sidebar:
    st.divider()
    st.header("Data & Benchmarking")
    
    st.markdown("**Export Current Simulation**")
    pdf_bytes = generate_qms_pdf(
        cell_line=selected_cell, vol=vol, dur=dur, target=target,
        yield_val=true_val, purity=dynamic_purity, consistency=dynamic_consistency, q_score=quality_score
    )
    
    st.download_button(
        label="Download PDF",
        data=pdf_bytes,
        file_name="EVelution_Batch_Projection.pdf",
        mime="application/pdf",
        use_container_width=True
    )
    
    st.markdown("<br>", unsafe_allow_html=True)
    
    st.markdown("**Historical Benchmarking**")
    st.info("Drop a CSV of your best historical run here. The engine will project it as a dashed 'ghost line' on your main graph so you can optimize against it in real-time.")
    
    uploaded_file = st.file_uploader(
        "Upload History File (.csv)", 
        type="csv",
        help="Ensure your CSV contains 'Hour' and 'Therapeutic EVs' columns for the graph overlay to function.",
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
    <div style="font-size: 0.75em; color: #808080; margin-top: 20px; padding-top: 15px; border-top: 1px solid #333; line-height: 1.4; text-align: justify;">
        <b>DATA GOVERNANCE AND LIABILITY DISCLAIMER</b><br>
        Data uploaded to this application is processed strictly in volatile memory (RAM) for ephemeral visualization. The architecture contains no persistent storage, logging, or external transmission protocols. By utilizing this software, end-users assume absolute liability for the safeguarding of proprietary technical trade secrets, in compliance with the Swedish Trade Secrets Act (SFS 2018:558, as amended 2026), the EU Trade Secrets Directive (2016/943), and the US Defend Trade Secrets Act (DTSA). Furthermore, operations align with global electronic record frameworks (FDA 21 CFR Part 11, EMA Annex 11) and international data sovereignty laws (GDPR) via zero-retention processing. Under international safe harbor and ephemeral data processing doctrines, this application acts purely as a localized execution environment, legally exempting the provider from data hosting, processor, and controller liabilities. Users are legally obligated to execute the internal data purge protocol upon session termination.
    </div>
    """, unsafe_allow_html=True)


# --- 10. METRICS & BATCH EVALUATION ---
m1, m2, m3, m4 = st.columns(4)
m1.metric("Yield Performance", f"{true_val:.2e}")
m2.metric("Harvest Conc", f"{df['Therapeutic EVs'].iloc[-1]:.2e} ev/mL")
m3.metric("Downstream Purity", f"{dynamic_purity*100:.1f}%")
m4.metric("Cargo Consistency", f"{dynamic_consistency*100:.1f}%")

st.markdown("### Batch Success Evaluation")

is_startup = (
    selected_cell == "Human MSCs (Bone Marrow)" and
    vol == 50.0 and o2 == 21 and temp == 37 and ph == 7.4 and mix == 85 and dur == 48 and
    target_base == 1.0 and target_exp == 15 and not manual_override
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
        status_text = "MARGINAL: Target reached, but cargo consistency is dropping. Monitor parameters."
    else:
        status_color, quality_color = "#E35252", "#E35252"
        status_icon = """<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="#E35252" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" style="vertical-align: middle; margin-right: 8px;"><path d="M10.29 3.86L1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"></path><line x1="12" y1="9" x2="12" y2="13"></line><line x1="12" y1="17" x2="12.01" y2="17"></line></svg>"""
        status_text = "CRITICAL: Severe quality degradation. High risk of DSP failure."
else:
    status_color, quality_color = "#779ECB", "#779ECB"
    status_icon = "🛑"
    status_text = "DEFICIENT: Target volume not reached. Extend duration or adjust feeding."

st.markdown(f"""
<div style="background-color: #1E1E2E; padding: 20px; border-radius: 8px; border-left: 6px solid {status_color}; display: flex; align-items: center; justify-content: space-between;">
    <div style="flex: 1;">
        <h4 style="margin: 0; color: {status_color}; display: flex; align-items: center;">
            {status_icon} {status_text}
        </h4>
    </div>
    <div style="text-align: right; min-width: 120px;">
        <div style="font-size: 0.9em; color: #B39EB5;">Yield vs Target</div>
        <div style="font-size: 1.5em; font-weight: bold; color: {'#77DD77' if yield_achievement >= 100 else '#779ECB'};">{yield_achievement:.1f}%</div>
    </div>
    <div style="text-align: right; padding-left: 20px; border-left: 1px solid #333; min-width: 140px;">
        <div style="font-size: 0.9em; color: #B39EB5;">Overall Quality Score</div>
        <div style="font-size: 1.5em; font-weight: bold; color: {quality_color};">{quality_score:.1f}%</div>
    </div>
</div>
""", unsafe_allow_html=True)


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
        This "controller panic" triggers aggressive impeller agitation. The resulting **hydrodynamic shear stress** physically shreds EV lipid bilayers and ruins cargo integrity ([Thompson & Papoutsakis, 2023](https://doi.org/10.1016/j.biotechadv.2023.108158)).
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
        t1.markdown("Determines the 'sweet spot' for harvest duration. The inflection point occurs where incremental EV gain is offset by culture necrosis and byproduct toxicity.")
        t2.latex(r"\frac{d}{dt}Yield(t) = 0 \quad at \quad t_{optimal}")

