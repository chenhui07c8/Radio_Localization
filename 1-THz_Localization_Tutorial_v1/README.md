Sample code for THz localization


Notes: 


1. This code considers SWM (near-field), which requires a high computational cost. We are preparing a simplified far-field model (which is widely used in mmWave systems) with performance loss and reduced complexity. *(coming soon)*
2. Constrained CRB (CCRB) not implemented in the current version, meaning orientation estimation applies to max two orientation unknowns (e.g., alpha & beta). We will provide a CCRB sample code later.
3. Spatial non-stationarity (SNS) is not considered; the analysis of model mismatch can be found in this work:\
    "Channel Model Mismatch Analysis for XL-MIMO Systems from a Localization Perspective." *(coming soon, 30-May)*
4. Hardware impairment (HWI) is not considered. Localization under HWI can be found in this work:\
    "MCRB-based Performance Analysis of 6G Localization under Hardware Impairments."
5. Doppler effect not considered. Localization under UE mobility can be found in this work:\
    "Doppler-Enabled Single-Antenna Localization and Mapping without Synchronization."  *(coming soon, 30-May)*



Please let me know if you have any questions or suggestions. I am happy for all types of discussions and collaborations on the topics related to **5G/6G Radio Localization**. :)
\
Email: hui.chen@chalmers.se; hui.chen@kaust.edu.sa
