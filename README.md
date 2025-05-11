# Dispersion engineering of 3D integrated ultra-low-loss silicon nitride waveguides

This repo contains codes for an undergraduate thesis project which aims at transforming the dispersion of microresonators made of 100nm-thick ultra-low-loss silicon nitride waveguides from normal to anomalous by leveraging the inter-resonator mode coupling.

------
## Abstract of the thesis
Microresonators featuring ultra-high quality factors can be realized with ultra-thin silicon nitride waveguides exhibiting exceptionally low loss. However, their inherent normal group velocity dispersion brings fundamental limitations to on-chip nonlinear photonic applications, including the Kerr frequency comb in terms of phase mismatch issues, narrow bandwidth, and incompatibility with bright soliton generation. Anomalous dispersion can be achieved in microresonators without degradation in quality factors by leveraging the interresonator mode coupling between them. Following a review of the physical principles and design methods of 2D coupled microresonators for dispersion engineering, we propose and analyze three novel 3D integrated configurations based on the coupled-mode theory and the finite-difference method in an attempt to achieve broadband near-zero anomalous dispersion in ultra-low-loss silicon nitride microresonators. The designs include: vertically stacked concentric ring resonators, double-layer racetrack resonators with inter-layer mode redistribution, and horizontally offset racetrack resonators with vertical coupling. The most feasible structure among them is the horizontally offset structure, which provides more flexibility and controllability to dispersion engineering, apart from a major increase in integration density. Given the geometry of the microresonators, an at least 7% increase in the anomalous dispersion spectral range can be achieved compared to 2D integrated structures, which can reach 200% when the lengths of the two cavities draw near.

Key words: integrated photonics, silicon nitride, dispersion, 3D integration

------
## File/Folder Structure

The files and folders of this project are organized as follows.

```yaml
SiNDispersionEngineering/
# The following python files are used to simulate mode-conserved coupling based on Coupled Mode Theory (CMT)
├── main.py
├── Waveguides.py
├── Coupled_Waveguides.py
├── Parameter_sweeper.py
├── Data_analyzer.py
├── Functions.py
├── Finding waveguide width satisfying the phase-match condition.ipynb
# These json files are used to store default parameters for plotting
├── Param_plot_curve.json
├── Param_plot_image.json
├── Param_plot_field_profile.json
# This folder contains Lumerical scripts used for Finite Difference Eigenmode solver (FDE) calculation
├── Lumerical Scripts
│   ├── EME_simulation_of_vertical_coupled_rings.lsf
│   ├── Export_index.lsf
│   ├── Find_Dispersion_of_2D_concentric_Rings.lsf
│   ├── Find_modes_of_horizontally_coupled_double_rings_func.lsf
│   ├── Find_modes_of_single_ring_at_multi_wavls_func.lsf
│   ├── Find_modes_of_single_ring_at_multi_wavls.lsf
│   ├── Find_modes_of_vertically_coupled_double_rings_func.lsf
│   ├── Find_modes_of_vertically_coupled_double_rings.lsf
│   ├── main_2D.lsf
│   ├── main_3D.lsf
│   ├── Scanning_Lx_and_gapx_given_gapy_3D_func.lsf
│   ├── Scanning_Lx_for_given_gapx_2D_func.lsf
│   ├── Scanning_Lx_for_given_gapx_2D.lsf
│   ├── setup_RDL.lsf
│   └── Vertical_coupled_rings_EME.lsf
# The following Jupyter Notebooks and python files are used to analyze the data and plot figures
├── Analysis of 2D concentric.ipynb
├── Analysis of 2D parallel.ipynb
├── Analysis of 3D concentric.ipynb
├── Analysis of 3D RDL.ipynb
├── Analysis of 3D offset.py
├── Extract_EME_res.ipynb
# This folder contains parameters used in FDE simulation
├── config
│   ├── Param_800x400.csv
│   ├── Param_L_inner_2_8.csv
│   ├── Param_L_inner_8.csv
│   ├── Param_straight_2_8.csv
│   └── Param_vertical.csv
# This folder is used to store results and figures
├── results
│   ├── 2D concentric rings
│   │    ├── ...
│   ├── 2D parallel rings
│   │    ├── ...
│   ├── 3D concentric rings
│   │    ├── ...
│   ├── 3D RDL
│   │    ├── ...
│   ├── 3D offset racetracks
│   │    └── ...
# This folder contains Jupyter Notebooks used to analyze the approximate analytical solution of bent strip waveguides
├── Analytical Solutions
│   ├── Find_Mode_Single_Ring_Rect_Approx.ipynb
│   ├── Find_Mode_Single_Ring.ipynb

```
------

**Author:** Weihao Xu

**Date:** May. 11th, 2025

**Email:** harryxwh@gmail.com

**End of Documentation**
