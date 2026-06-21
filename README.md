# EVelution-Bioreactor-Twin
**Powered by the Multi-Machinery Model of EV Biogenesis**



[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Streamlit App](https://img.shields.io/badge/Streamlit-Live_Demo-FF4B4B.svg)](#live-demo)
[![License: Copyright](https://img.shields.io/badge/License-Copyright_Evgenia_Trudova-red.svg)](#copyright--terms-of-use)


<a href="https://evelution-bio.streamlit.app/">
  <img width="1852" height="777" alt="Evelution Bio App" src="https://github.com/user-attachments/assets/ff37be63-63e1-4fee-8407-3510ed58de5c" />
</a>


EVelution-bio is a digital twin and predictive bioinformatics tool designed for the optimization of extracellular Vesicle biomanufacturing. EVelution-bio bridges upstream cell biology and downstream process engineering, helping reduce the costs of physical trial-and-error. By simulating the precise environmental stressors  (hypoxia, thermal shifts, and hydrodynamic shear), the application calculates where therapeutic yield is maximized before cellular viability collapses.

The bioprocessing industry currently suffers from a critical blind spot: optimizing scale-up runs based on raw NTA (Nanoparticle Tracking Analysis) scatter counts. EVelution introduces the Yield-to-Value Bridge, a mathematical framework that deconstructs raw particle counts by severely penalizing the batch for empty cargo and apoptotic debris. Instead of showing a deceptive "total count," EVelution tracks the True Functional Value—isolating only the structurally intact, therapeutically aligned EVs that will survive downstream filtration and chromatography.

Current scale-up strategies often falsely assume EV secretion is a static baseline. The underlying Multi-Machinery Model challenges this by mathematically modeling biogenesis as a dynamically regulated, multi-pathway survival network.

Uses: 
The "Reality Check"
Bioprocess engineers often develop high-yield protocols in wet lab flasks that fail when scaled up. 
Engineers can input their target yield and the specific bioreactor volume to see if their process is physically achievable. The app functions as a "sanity check" to determine if they need to increase tank volume or refine their feeding strategy, preventing the waste.

QC Troubleshooting
When a batch fails downstream quality control, the root cause is often "invisible" because the NTA machine only counts particles—it doesn't distinguish between therapeutic EVs and apoptotic debris.Engineers can upload historical CSV run data into the app. By overlaying physical data with the digital twin's kinetic predictions, they can pinpoint exactly which hour the "Death Cliff" occurred. This identifies whether the batch was ruined by an oxygen-starved "dead zone" or if it was ruined by over-agitation creating hydrodynamic shear.




