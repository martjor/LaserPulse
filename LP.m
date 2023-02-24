% File Name:            laser_pulse.m
% Author:               Jorge Augusto Martinez-Ortiz
% Date Created:         01/26/2023
% Description:          Analysis of laser pulse in GEC RF Cell experiments

%% LOAD DATA
clear
path = '/Users/jorgeaugustomartinezortiz/Library/CloudStorage/OneDrive-BaylorUniversity/CASPER/projects/plasma_crystal/crystal_analysis/torsions_datasets/laser_pulse_20221006_crop_0_999.csv';

T = pcryReadTable(path);

%% Try some processing
frame = 100;
F = pcryGetFrame(T,frame);
R = pcryNorm2d([F.x F.y],[F.x F.y]);

[corrfun r rw] = twopointcorr(F.x,F.y,3);
figure
plot(r,corrfun)

b = 27;

den = sum(exp(-(R.^2/(b*b))),2);

figure
pcryVoronoi(F.x,F.y,den)
