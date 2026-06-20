# EVelution-Bioreactor-Twin
EVelution is a series of bioinformatics simulations. Bioreactor Twin is a digital twin designed for optimizing extracellular vesicle yield in industrial fed-batch bioreactors.


# EVelution: Fed-Batch Bioreactor Digital Twin
**Powered by the Multi-Machinery Model of EV Biogenesis**

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Streamlit App](https://img.shields.io/badge/Streamlit-Live_Demo-FF4B4B.svg)](#live-demo)
[![License: Copyright](https://img.shields.io/badge/License-Copyright_Evgenia_Trudova-red.svg)](#copyright--terms-of-use)

## Overview
**EVelution** is a time-series bioinformatics simulation and digital twin designed for optimizing Extracellular Vesicle (EV) yield in industrial fed-batch bioreactors. 

Bridging the gap between molecular biology and bioprocess engineering, this tool replaces upstream trial-and-error with mathematical precision. By simulating the environmental stressors that trigger EV shedding (hypoxia, thermal shifts, and acidic environments), EVelution predicts the exact "Sweet Spot" where therapeutic EV production is maximized before cellular viability crashes.

## The Multi-Machinery Model
Current biomanufacturing often assumes EV secretion is a static baseline. The **Multi-Machinery Model** challenges this, proposing that EV biogenesis is a dynamically regulated survival network. EVelution separates EV production into:
1. **The Gradient Baseline:** Constitutive secretion of *Therapeutic EVs* based on metabolic health.
2. **Acute Stress Machinery:** Explosive shedding of *Stress-Altered EVs* triggered only when specific environmental thresholds (O2, Temp, pH) are crossed.

## Key Features
* **48-Hour Time-Series Kinetics:** Simulates dynamic environmental shifts, oxygen depletion, and EV accumulation over a standard fed-batch run.
* **Spatial Zonal Mixing:** Accounts for hardware limitations (impeller speed) by simulating hypoxic/acidic "Dead Zones" alongside well-mixed zones.
* **Cargo Consistency Tracking:** Mathematically differentiates between desired therapeutic cargo and stress-altered survival EVs.
* **Downstream Purity Penalties:** Models the exponential spike in apoptotic bodies as cell viability drops, calculating the direct impact on downstream filtration and Cost of Goods (COGs).
* **Machine-Learning Ready Architecture:** Built with decoupled, parameterized biological constants ready for `scipy.optimize` fitting against physical Nanoparticle Tracking Analysis (NTA) lab data.

