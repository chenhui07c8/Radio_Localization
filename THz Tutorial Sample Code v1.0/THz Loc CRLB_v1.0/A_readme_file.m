% Author: hui.chen@kaust.edu.sa; hui.chen@chalmers.se
% THz Loc CRLB_v0.51
% This sample code provides a AOSA-based SWM channel model with CRLB analysis on the position and orientation.
% All the figures in the paper "A Tutorial on Terahertz-Band Localization
% for 6G Communication Systems" are generated using this code.

% Funciton Naming Rules:
% Sim: simulation. Draft: test for different realizations; Debug: debug codes
% get_xx_xx_xx.m: Block of codes for simplicity
% getAaaaAaaAa.m: Algorithms or functions to obtain some results
% the unit of the angles is [degree], although [rad] is used in the paper.

% CRLB_v0.0 210325
% -3D deterministic localization model with BS-UE and BS-RIS-UE channels.
% -Two RIS optimizer, one far-field, one near-field.
% -4 CRLB considerations: 2D4V/2D8V, syn/asyn. 
%     4V: x, y, alpha (gain), beta (clock offset)
%     8V: x, y, alpha, rho_BM, xi_BM, rho_BRM, xi_BRM, beta

% CRLB_v0.1 210414
% -Add subarray (analog) beamforming gain 'A' and antenna gain 'G'.

% CRLB_v0.2 210512
% -Derivative of the equivalent array response AeqBM, AeqBR, AeqMB, AeqMR
% -Add AOSA function with customized beamforming angles at SA level.

% CRLB_v0.3 210527
% -Update parameter notations (based on the paper)
% -Add AOSA simulations

% CRLB_v0.4 210613
% -add multiple transmissions, symbol block size: N*K*T. (N is the number of receivers)
% -add four CRLB modes

% CRLB_v0.41
% -change exp(1j*c.xi) to exp(1j*2*pi/lambdac*c.xi), to remove singular matrix.

% CRLB_v0.42
% -add SPP + SSP realizations (SA distance/SA angle/AOSA distance) for CRLB
% -rewrite derivative in matrix form
% -add downlink CRLB

% CRLB_v0.43
% -add deterministic NLOS channels with N_C clusters
% measurement vector = [rhoBCM, xiBCM, tauBC, tauBC, thetaBC, tauCM, thetaCM, phiCM, OriM, beta]
% state vector = [rhoBCM, xiBCM, Pc, OriM, beta]

% CRLB_v0.44 210703
% -several bugs fixed

% CRLB_v0.44 210722
% -modify RIS channel model as product-distance path loss model

% CRLB_v0.50 210729
% -Add sub-RIS (similar concept as AOSA) to support large RIS
% -get_FIM(c) returns the FIMs of both uplink & downlink
% -make the whole code more structured.

% CRLB_v0.50 220313
% -Be ready for github.

% CRLB_v1.0 220522
% -Uploaded to Github
% 



% TODO:
% 0. Using a near-field (SWM) model is generally much slower than a far-field
% model, we are working on optimizing the code for a better user experience.
% 1. Far-field model and analog array.
% 2. Multiple BS/RIS
% 3. Different types of NLOS (virtual anchor for reflector or scattering points)
% 4. Other types of CRB, e.g., constrained CRB for 3D orientation,
% misspecified CRB for HWI and model mismatch, posterior CRB for tracking
% 5. Model with hardware impairments (HWI)





% Please consider cite this paper if you think this code is helpful  :)

% Chen, H., Sarieddeen, H., Ballal, T., Wymeersch, H., Alouini, M. S., and
% Al-Naffouri, T. Y. "A tutorial on terahertz-band localization for 
% 6G communication systems". IEEE Communications Surveys & Tutorials. 2022.





















