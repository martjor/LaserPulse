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
fileName = '/Users/jorgeaugustomartinezortiz/Library/CloudStorage/OneDrive-BaylorUniversity/CASPER/projects/plasma_crystal/crystal_analysis/torsions_datasets/laser_pulse_20221006_crop_0_999.csv';
T = pcryReadTable(fileName);

particle = unique(T.particle);
frame = unique(T.frame);

numParticles = numel(particle);
numFrames = numel(frame);

% Fill gaps for particles that disappear
TNan = pcryFillNaN(T);

%% GENERATE PAIR CORRELATION FUNCTION
dr = 3;
g = {};
r = {};
for i = 1:numFrames
    f = pcryGetFrame(T,frame(i));
    [g{i},r{i},~] = twopointcorr(f.x,f.y,dr);
    fprintf("Frame %i\n",frame(i));
end

%% APPROXIMATE INTERPARTICLE DISTANCE
[gAvg,rTotal] = corrfun(g,r);
figure
plot(rTotal,gAvg)

b = 27;

%% CALCULATE BACKGROUND CONFIGURATIONAL ENERGY
U = zeros(numParticles,1);
count = zeros(numParticles,1);

for i = 1:numFrames
    % Get Current frame
    f = pcryGetFrame(T,frame(i));

    % Get index to particles
    idx = f.particle + 1;

    % Add Energies to the corresponding particles and count
    U(idx) = U(idx) + Econfig(f.x,f.y,b);
    count(idx) = count(idx) + 1;

    fprintf("Frame (%i/%i)\n",i,numFrames);
end

% Calculate background energy for each particle
U = U ./ count;

%% DISPLAY

figure
ax = axes;
xlim([min(T.x) max(T.x)]);
ylim([min(T.y) max(T.y)]);
clim(ax,[-1 1])

M(numFrames) = struct('cdata',[],'colormap',[]);
for i = 1:numFrames
    f = pcryGetFrame(T,frame(i));
    idx = f.particle + 1;
    E = Econfig(f.x,f.y,27);
    EAvg = U(idx);
    div = (E-EAvg) ./ EAvg;
    
    cla
    pcryVoronoi(f.x,f.y,log(abs(div)+1))
    title(sprintf("Frame %i",frame(i)));
    drawnow

    cdata = print('-RGBImage','-r120');
    M(i) = im2frame(cdata);
end

%%
writeVid('PercentDeviation.avi',M)



function E = Econfig(x,y,b)
    % Calculate distances between particles
    r = pcryNorm2d([x y],[x y]);
    ru = triu(r);
    rl = tril(r);

    r = rl(:,1:(end-1)) + ru(:,2:end);

    % Subtract 1 to compensate for self-interaction
    E = sum(exp(-r/b)./r,2);
end


function [gAvg,rTotal] = corrfun(g,r)
    rTotal = [];

    % Find the unique values for r
    for i = 1:numel(r)
        rTotal = union(rTotal,r{i});
    end

    % Find the average value of g(r) at every r
    gAvg = zeros(numel(rTotal),1);
    for i = 1:numel(rTotal)
        rCurr = rTotal(i);
        gCurr = g{i};
        
        count = 0;
        gTemp = 0;
        
        for j = 1:numel(r)
            idx = rCurr == r{j};

            if any(idx)
                gTemp = gTemp + gCurr(i);
                count = count + 1;
            end
        end

        gAvg(i) = gTemp / count;
        fprintf("Processing (%.1f)\n",i/numel(rTotal)*100);
    end
end

function writeVid(name,M)
    v = VideoWriter(name);
    v.Quality = 100;
    open(v);
    writeVideo(v,M);
    close(v);
end

