**Agent-Based Contamination Modelling in Drinking Water Distribution Networks (DWDN)**

**Description**

This repository contains the integrated modelling framework developed for the paper on agent-based behavioural response during microbial contamination events in drinking water distribution networks (DWDNs).
The framework couples:
  - Agent-Based Modelling (ABM) of consumer behavioural response
  - Stochastic household water demand generation
  - Hydraulic and multi-species water quality simulation (EPANET-MSX)
  - Quantitative Microbial Risk Assessment (QMRA)

The model explicitly represents the socio-technical feedback between advisory-induced demand reduction and contaminant transport, and evaluates resulting infection risk under different detection and advisory timing scenarios.

**Overview**

The framework simulates contamination events in a DWDN under multiple scenarios combining:
  - Detection timing (e.g., early vs moderate)
  - Advisory timing (immediate vs post source identification)
  - Heterogeneous household behavioural response (partial/full adherence, delayed adoption, awareness diffusion)

Behavioural changes modify nodal water demands, which in turn alter:
  - Flow velocities
  - Residence times
  - Chlorine decay
  - Pathogen transport dynamics

Exposure is quantified at the individual level using a QMRA approach, computing dose and infection probability over time.
The framework enables systematic comparison between:
  - Baseline (RAW demand, no behavioural adaptation)
  - ABM-driven demand response scenarios

Main execution scripts:
  - ABM_CCWI_2026_Create_STREAM_demands_LTOWN.m – Defines ABM-based behavioral parameters and STREaM demands
  - ABM_CCWI_2026_assign_STREAM_demands_LTOWN.m – Assigns the newly generated STREaM demands and creates a new .inp
  - ABM_CCWI_2026_setup_MSX_LTOWN.m – Setup and run the MSX
  - ABM_CCWI_2026_risk_evaluation_LTOWN.m – Calculates the infection risk using QMRA
  - ABM_CCWI_2026_visualization.m – Visualization of results

**Requirements**

The framework was developed and tested using:
  - [MATLAB](https://www.mathworks.com/)
  - [EPANET-MATLAB Toolkit](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit)
  - [Symbolic Math Toolbox](https://www.mathworks.com/products/symbolic.html)

**Instructions**
Run the main execution scripts with the exact same order to perform the simulations of the paper and see the results

**Contributors**
- [Sotiris Paraskevopoulos](https://github.com/Sotireas), [ITI, Centre for Research and Technology Hellas (CERTH)](https://www.certh.gr/46A8F1C8.el.aspx)
- [Stelios Vrachimis](https://github.com/SteliosVr), [KIOS Research and Innovation Center of Excellence, University of Cyprus](https://www.kios.ucy.ac.cy/)
- [Marios Kyriakou](https://github.com/Mariosmsk), [KIOS Research and Innovation Center of Excellence, University of Cyprus](https://www.kios.ucy.ac.cy/)

