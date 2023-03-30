% File Name:            LP.m
% Author:               Jorge Augusto Martinez-Ortiz
% Date Created:         01/26/2023
% Description:          Analysis of laser pulse in GEC RF Cell experiments.
%                       We simply use this script as a test to determine
%                       the method we will use to do our analysis
%
%                       03.13.2023: After meeting with Parker we determined
%                       that we will be coloring the cells using an
%                       equation similar to that of the configurational
%                       energy.

%% LOAD DATA
clear
path = '/Users/jorgeaugustomartinezortiz/Library/CloudStorage/OneDrive-BaylorUniversity/CASPER/projects/plasma_crystal/crystal_analysis/torsions_datasets/laser_pulse_20221006_crop_0_999.csv';

T = pcryReadTable(path);
%% CALCULATE SOME IMPORTANT PARAMETERS OF THIS CRYSTAL
numFrames = max(T.frame) + 1;               % We also account for the 0th frame
numParticles = zeros(numFrames,1);
height = max(T.y) - min(T.y);
width = max(T.x) - min(T.y);
area = height * width;

for frame = 1:numFrames
    currFrame = pcryGetFrame(T,frame);
    numParticles(frame) = size(currFrame,1);               % Number of entries in the current frame
end

numParticlesMean = round(mean(numParticles));
areaPerParticle = area / numParticlesMean;

%% METHOD 1: Measuring Particle Proximity Based On Cell Area
numFramesMovie = 100;
MovArea(numFramesMovie) = struct('cdata',[],'colormap',[]);
MovDist(numFramesMovie) = struct('cdata',[],'colormap',[]);

%% METHOD 1: Area
% I calculate the area of the section of the crystal I am looking at and
% divide it by the mean number of particles. This gives me the area per
% particle. Then, I calculate the area of a particle's voronoi cell on each
% frame and I calculate the difference between this area and the mean area
% per particle. This allow me to color the cells based on whether they're

f1 = figure;
ax1 = axes(f1);
axis equal
xlim([min(T.x) max(T.x)]);
ylim([min(T.y) max(T.y)]);
areaPercent = 0.5;
clim([-areaPercent areaPercent] * areaPerParticle);
for i = 1:numFramesMovie
    F = pcryGetFrame(T,i);

    % Voronoi Cells for the area method
    [v,c] = voronoin([F.x F.y]);
    area = pcryAreaVoronoi(v,c);

    cla(ax1)
    pcryVoronoi(F.x,F.y,area-areaPerParticle);
    title(sprintf("Frame %i",i))
    
    cdata = print('-RGBImage','-r120');
    MovArea(i) = im2frame(cdata);
end

%% Distance method
% I use the pair correlation function to calculate the mean interparticle
% distance of the particles in the crystal. On each frame, I calculate the
% distance between particles and for each particle, I generate a score
% based on the distance between particles. The score decreases as e^-(r^2).
f2 = figure;
ax2 = axes(f2);
axis equal
xlim([min(T.x) max(T.x)]);
ylim([min(T.y) max(T.y)]);
clim([-2 2]);
binDist = 3;
[corrfun r rw] = twopointcorr(F.x,F.y,binDist);
meanDist = 27;

% Calculate what should be the expected value for a single dust particle
numCellParticles = 6;
meanScore = numCellParticles * exp(-1);

for i = 1:numFramesMovie
    F = pcryGetFrame(T,i);

    % Distance between points
    R = pcryNorm2d([F.x F.y],[F.x F.y]);

    % Apply scoring
    score = sum(exp(-(R/meanDist).^2),2);

    % Update plot
    cla(ax2)
    pcryVoronoi(F.x,F.y,score-meanScore);
    title(sprintf("Frame %i",i));

    cdata = print('-RGBImage','-r120');
    MovDist(i) = im2frame(cdata);
end

%% Write video
v = VideoWriter('movie.avi');
v.Quality = 100;
open(v);
writeVideo(v,MovArea);
close(v);