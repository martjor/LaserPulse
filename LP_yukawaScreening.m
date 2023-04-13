% File Name:                LP_yukawaScreening
% Author:                   Jorge A. Martinez-Ortiz
% Date Created:             03.15.2023
% Description:              I created this script following the
%                           conversation that I had with partker. We
%                           decided that we will adopt a method similar to
%                           the distance or scoring method I used in LP.m
%                           Coincidentally, using the decaying exponential
%                           is very commonly used not just to calculate the
%                           shielded electrostatic interactions, but also
%                           the configurational energy. We decided that
%                           sticking to this woul allow us to justify our
%                           method on physics grounds more easily. 
%
%                           We also decided that we will calculate the mean
%                           'score' or 'energy' on each particle across the
%                           duration of the video and we will subtract it
%                           from that of a specific frame 

% Load data
fileName = 'Results1.csv';
%fileName = '/Users/jorgeaugustomartinezortiz/Library/CloudStorage/OneDrive-BaylorUniversity/CASPER/projects/plasma_crystal/crystal_analysis/torsions_datasets/laser_pulse_20221006_crop_0_999.csv';
T = pcryReadTable(fileName,'fiji');

particle = unique(T.particle);
frame = unique(T.frame);

numParticles = numel(particle);
numFrames = numel(frame);
%% GENERATE PAIR CORRELATION FUNCTION
% NOTE: This step is not necessary if you already now your value fo b (mean
% interparticle distance).
dr = 3;
g = {};
r = {};
for i = 1:numFrames
    f = pcryGetFrame(T,frame(i));
    [g{i},r{i},~] = twopointcorr(f.x,f.y,dr);
    fprintf("Frame %i\n",frame(i));
end

[gAvg,rTotal] = corrfun(g,r);
figure
plot(rTotal,gAvg)

%% CALCULATE BACKGROUND CONFIGURATIONAL ENERGY
b = 27;
U = zeros(numParticles,1);
count = zeros(numParticles,1);

for i = 1:numFrames
    % Get Current frame
    f = pcryGetFrame(T,frame(i));

    % Get index to particles
    idx = f.particle;

    % Add Energies to the corresponding particles and count
    U(idx) = U(idx) + energyConfig(f.x,f.y,b);
    count(idx) = count(idx) + 1;

    fprintf("Frame (%i/%i)\n",i,numFrames);
end

% Calculate background energy for each particle
U = U ./ count;

%% GENERATE VIDEO
% Filename where the video will be recorded
fileName = 'PercentDeviation.avi';

% Frames that will be animated/recorded
frameVec = 1:(numFrames-1);

% Create animation figure
animFunc = @(frame) plotEnergyDev(T,U,b,frame);

% Instantiate Animator object
Anim = pcryAnimator(animFunc);

% Generate figure container and specify axes properties
figure
ax = axes;
xlabel('x [Pixels]')
ylabel('y [Pixels]')
axis equal
xlim([min(T.x) max(T.x)]);
ylim([min(T.y) max(T.y)]);
colorbar
colormap spring
clim(ax,[-1 1])

% Animate first to verify that you like the video
Anim.animate(frameVec);

% Record video
%Anim.recordVideo(frameVec,fileName);

%% FUNCTION DEFINITIONS (NO NEED TO EXECUTE)

% Function that will generate a plot
% This function will be pased as the argument to instantiate the animator
% object. 
% The function should take specify the input variables and the way they
% will be accessed whenever the frame number is provided
function plotEnergyDev(T,U,b,frame)
    % Get current frame
    F = pcryGetFrame(T,frame);
    x = F.x;
    y = F.y;

    % Get indices of particles in frame
    % The +1 is very important here, since matlab indices go start at 1
    % particle indices start at 0. This step is necessary to match the
    % right energies to the right particle. It probably is a good idea to
    % pre-process the data to change the indices of the particles to start
    % at 1 altogether
    idx = F.particle;

    % Calculate the configurational energy of the particles
    E = energyConfig(x,y,b);

    % Gather the mean configurational energy of the particles
    EAvg = U(idx);

    % Calculate the percent diff
    diff = (E-EAvg) ./ EAvg;

    % Generate voronoi plot
    pcryVoronoi(x,y,diff);
    title(sprintf("Frame %i",frame));
end

function E = energyConfig(x,y,b)
    % Calculate distances between particles
    r = pcryNorm2d([x y],[x y]);
    ru = triu(r);
    rl = tril(r);

    % Distances between particles not including the main diagonal which is
    % only zeros
    r = rl(:,1:(end-1)) + ru(:,2:end);

    % Calculate proportionality factor of configurational energy
    E = sum(exp(-r/b)./r,2);
end

% Averages the pair correlation function across multiple frames
% This ended up being not very useful
function [gAvg,rBins] = corrfun(g,r)
    numFrames = numel(g);
    rBins = [];

    % Find the unique values for r
    for i = 1:numFrames
        rBins = union(rBins,r{i});
    end
    numBins = numel(rBins);

    % Find the average value of g(r) at every r
    gAvg = zeros(numBins,1);
    for i = 1:numBins
        % Current r value
        rCurr = rBins(i);
        
        % Reset counter and current value of g
        count = 0;
        gTemp = 0;
        
        % Iterate through the pair-correlation-fcn of every frame
        for j = 1:numFrames
            % Create an index to the bin corresponding to the current val
            idx = rCurr == r{j};

            % If the bin value is in this frame, process
            if any(idx)
                gCurr = g{j};
                gTemp = gTemp + gCurr(idx);
                count = count + 1;
            end
        end

        gAvg(i) = gTemp / count;
        fprintf("Processing (%.1f)\n",i/numel(rBins)*100);
    end
end