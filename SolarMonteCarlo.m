clc; clear; close all; tic; addpath("Common");
% =================================================
%
%         Monte Carlo Simulation of Solar Visibility
%
%                   By Ian Byres
%
%  ================================================
%% Notes

% Uses Monte Carlo method to vary classical orbital elements of orbit based
% upon mean value and 3sigma dispersion. Computes respective orbits using
% Runge Kutta numerical integrator and computes time in Earths shadow vs
% period for each orbit to find gaussian distribution of time in sun vs
% shadow.

% All units converted first and computations done in [km,kg,s]

% Works for circular and elliptical orbits

% References:
%   1. Curtis, Howard D. Orbital Mechanics for Engineering Students. 
%      4th ed., Butterworth-Heinemann, 2021. 

%% ===================== EDIT HERE ====================================

% Simulation timing (date is used to find position of sun relative to)
startDate = [2024, 12, 17, 0, 0, 0]; % [Year, month, day, hour, min, sec]
tInt = [0.5, 0]; % Integrating time interval in [minutes,seconds]

% Number of orbits to simulate
numOrbits = 2000;

% Nominal values for orbital elements
meanEcc = 0.5; % Nominal Eccentricity
meanPerAlt = 600; % Nominal altitude of perigee [km]
meanInc = 66; % Nominal orbital inclination [deg]
meanRAAN = 12; % Nominal right ascension of ascending node [deg]
meanArgPer = 8; % Nominal argument of periapsis [deg]

% Estimated 3 standard deviation range of orbital elements
sigma3Ecc = 0.1; % 3sigma Eccentricity
sigma3PerAlt = 50; % 3sigma altitude of perigee [km]
sigma3Inc = 5; % 3sigma orbital inclination [deg]
sigma3RAAN = 5; % 3sigma right ascension of ascending node [deg]
sigma3ArgPer = 5; % 3sigma argument of periapsis [deg]

% Postprocessing options
plotDispersedFlag = 1; % 1 to plot generated dispersed orbits along with nominal (Up to first 50 orbits), 0 otherwise

%% ===================== DO NOT EDIT =================================

%% Preparing Inputs

% Planet Specifications (Earth)
EarthRadAvg = 6378.1366; % [km]
muEarth = (3.986004418e14)/(1e9); % [km^3/s^2]

% Preparing integration timing
tInt = (tInt(1)*60)+tInt(2);

% Postprocessing options
if plotDispersedFlag == 1
    additionalOrbPlotted = min(numOrbits,50);
    dispOrbitsCell = cell(1,additionalOrbPlotted);
else
    additionalOrbPlotted = 0;
end
%% Dispersing Orbital Elements

% Dispersion of orbital elements
disEcc = meanEcc + randn(numOrbits,1) .* (sigma3Ecc/3);
disPerAlt = meanPerAlt + randn(numOrbits,1) .* (sigma3PerAlt/3);
disInc = meanInc + randn(numOrbits,1) .* (sigma3Inc/3);
disRAAN = meanRAAN + randn(numOrbits,1) .* (sigma3RAAN/3);
disArgPer = meanArgPer + randn(numOrbits,1) .* (sigma3ArgPer/3);

% Restricting range of dispersion
disEcc(disEcc >= 1) = 0.999999;
disEcc(disEcc < 0) = 0;
disPerAlt(disPerAlt < 100) = 100;

% Computing semi major axis, angular momentum, and period
meanSMA = (meanPerAlt + EarthRadAvg)/(1 - meanEcc); % Nominal Semi major axis [km]
meanT = 2*pi*sqrt((meanSMA^3)/muEarth); % Nominal period [s]
meanh = sqrt(muEarth * meanSMA * (1-(meanEcc^2))); % Nominal specific angular momentum

disSMA = (disPerAlt + EarthRadAvg)./(1 - disEcc); % Semi major axis [km]
dish = sqrt(muEarth .* disSMA .* (1-(disEcc.^2))); % Specific angular momentum [km^2/s]
disT = 2*pi.*sqrt((disSMA.^3)./muEarth); % Orbital period [s]

% Finding initial state vector
r0 = zeros(numOrbits,3);
v0 = zeros(numOrbits,3);
for i = 1:1:numOrbits
    [r0(i,:),v0(i,:)] = sv_from_coe([dish(i),disEcc(i),disRAAN(i),disInc(i),disArgPer(i),0],muEarth);
end

%% Computing Sun Position

maxOdeTime = 1:tInt:(max(disT)+1000);
ephTime = repmat(startDate,length(maxOdeTime),1);
tempvar = max(disT)+1000;
for i = 0:1:(size(ephTime,1)-1)
    ephTime(i+1,6) = ephTime(i+1,6) + (tInt*i);
end
ephDate = juliandate(ephTime);
ephSunPos = planetEphemeris(ephDate,'Earth','Sun');
ephSunPosMag = vecnorm(ephSunPos,2,2);

%% Determining orbits and shadowing

sunlitDuration = zeros(numOrbits,1);
for i = 1:1:numOrbits
    % Computing orbit
    initCond = [r0(i,:),v0(i,:)];
    odeTime = (0:tInt:disT(i))';
    [~,Xorb] = ode45((@(t,X) OrbitPropagator(t,X,muEarth)), odeTime, initCond);
    r = Xorb(:,1:3);
    % Determining shadowing
    sunlit = ones(length(odeTime),1);
    rMag = vecnorm(r,2,2);
    sunPosVec = ephSunPos(1:size(r,1),:);
    sunPosMag = ephSunPosMag(1:size(r,1),:);

    sunAngle = acos(  dot(sunPosVec,r,2) ./ (rMag.*sunPosMag)  );
    angle1 = acos(  EarthRadAvg./rMag );
    angle2 = acos(  EarthRadAvg/sunPosMag );

    sunlit(  (angle1+angle2) <= sunAngle   ) = 0;
    sunlitDuration(i) = sum(sunlit) * tInt;

    if (plotDispersedFlag == 1) && (i <= additionalOrbPlotted)
        dispOrbitsCell{i} = r;
    end

end

% Percent of orbit in sunlight
sunlitPercent = (sunlitDuration./disT) * 100;

% Nominal Orbit Simulation
[r0Nom,v0Nom] = sv_from_coe([meanh,meanEcc,meanRAAN,meanInc,meanArgPer,0],muEarth);
initCondNom = [r0Nom,v0Nom];
odeTimeNom = (0:tInt:meanT)';
[~,XorbNom] = ode45((@(t,X) OrbitPropagator(t,X,muEarth)), odeTimeNom, initCondNom);


%% Displaying Table Results

toc;
fprintf("Note: Variance has units squared. For example the variance of period would be [hr^2] \n");
% Computing COE Stats
statsEcc = [meanEcc;sigma3Ecc;sigma3Ecc/3;mean(disEcc);std(disEcc);var(disEcc)];
statsPerAlt = [meanPerAlt;sigma3PerAlt;sigma3PerAlt/3;mean(disPerAlt);std(disPerAlt);var(disPerAlt)];
statsSMA = [meanSMA;"NA";"NA";mean(disSMA);std(disSMA);var(disSMA)];
statsInc = [meanInc;sigma3Inc;sigma3Inc/3;mean(disInc);std(disInc);var(disInc)];
statsRAAN = [meanRAAN;sigma3RAAN;sigma3RAAN/3;mean(disRAAN);std(disRAAN);var(disRAAN)];
statsArgPer = [meanArgPer;sigma3ArgPer;sigma3ArgPer/3;mean(disArgPer);std(disArgPer);var(disArgPer)];
statsT = [meanT/3600;"NA";"NA";mean(disT)/3600;std(disT)/3600;var(disT)/3600];

% Displaying sunlit stats table
sunlitTableResults = [["Mean Value";"Standard Deviation";"Variance"],...
                      [mean(sunlitPercent);std(sunlitPercent);var(sunlitPercent)],...
                      [mean(sunlitDuration)/3600;std(sunlitDuration)/3600;var(sunlitDuration)/3600]];
fprintf("\n                                   <strong>Sunlit Statistics Table</strong> \n \n");
disp(array2table(sunlitTableResults,"VariableNames",[" ","Percent of Orbit Period Sunlit","Total Time in Orbit Sunlit [hr]"]));

% Displaying COE stats table
coestats = [statsEcc,statsPerAlt,statsSMA,statsInc,statsRAAN,statsArgPer,statsT];
coestatsTopRow = [" ","Eccentricity","Perigee Altitude [km]","Semi-Major Axis [km]","Inclination [deg]",...
               "Right Ascension of Ascending Node [deg]","Argument of Perigee [deg]","Period [hr]"];
coestatsLeftColumn = ["Nominal Value";"Inputted 3 Standard Deviations";"Inputted Standard Deviation";...
                   "Sim Results Mean";"Sim Results Standard Deviation";"Sim Results Variance"];
coetableArray = [coestatsLeftColumn,coestats];
fprintf("\n                                                                          "+...
    "<strong>Orbital Elements Statistics Table</strong> \n \n");
disp(array2table(coetableArray,"VariableNames",coestatsTopRow));

% Plot of Orbit

% Setting up figure
f1 = figure(1);
f1.Name = 'Orbit Plot';
hold on;

% Displaying Earth
[XEIMG, map] = rgb2ind(imread('Earth.jpg'),128);
[xEIMG,yEIMG,zEIMG] = sphere(50);
xEIMG = EarthRadAvg*xEIMG;
yEIMG = EarthRadAvg*yEIMG;
zEIMG = EarthRadAvg*zEIMG;
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.Cdata = flipud(XEIMG);
surface(xEIMG,yEIMG,zEIMG,props);
colormap(map);

% Displaying Sun Vector
sunPointingVector = (ephSunPos(1,:)/ephSunPosMag(1)) * max(vecnorm(XorbNom(:,1:3),2,2));
quiver3(0,0,0,sunPointingVector(1),sunPointingVector(2),sunPointingVector(3),'Color','#fcc203')

% Plotting Nominal Trajectory
plot3(XorbNom(:,1),XorbNom(:,2),XorbNom(:,3),'linewidth',2.25,'Color','#14c400');

% Plotting Dispersed Trajectories
if plotDispersedFlag == 1
    for i = 1:1:length(dispOrbitsCell)
        plot3(dispOrbitsCell{i}(:,1),dispOrbitsCell{i}(:,2),dispOrbitsCell{i}(:,3),'linewidth',1.25,'Color','#a82a35');
    end
end

% Formatting Plot
title('Trajectory Propagations');
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
if plotDispersedFlag == 1
    legend('','Sun Pointing Vector at Start of Nominal Orbit','Nominal Trajectory','Dispersed Trajectories');
else
    legend('','Sun Pointing Vector','Nominal Trajectory');
end
set(gca,'Color','#0a0a0d');
axis equal;
grid off;
view([18 29]);
f1.Color = "#0a0a0d";
f1.WindowState = 'maximized';

%% Probability Distribution Plots for Sunlight

% Setting up figure
f2 = figure(2);
f2.Name = 'Probability Distributions Sunlight';
f2.WindowState = 'maximized';
hold on;

% Creating Percentage Plot
subplot(1,2,1);
normalpd1 = fitdist(sunlitPercent,'Normal');
plot(normalpd1);

title('Probability Curve of Percentage of Orbit Being Sunlit');
xlabel('Percentage of Orbit Duration Sunlit');
ylabel('Probability Density');
if (mean(sunlitPercent) + 4*std(sunlitPercent)) > 100
xlim([-inf 100]);
end

dim = [0.40 0.7 0.2 0.2];
str = {'Mean: ' num2str(mean(sunlitPercent)),'Standard Deviation: ' std(sunlitPercent),'Variance: ' var(sunlitPercent)};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

% Creating Time Plot
subplot(1,2,2);
normalpd2 = fitdist(sunlitDuration./3600,'Normal');
plot(normalpd2);

title('Probability Curve of Sunlit Duration in Orbit');
xlabel('Total Time of Orbit Being Sunlit [hr]');
ylabel('Probability Density');

dim = [0.83 0.7 0.2 0.2];
str = {'Mean: ' num2str(mean(sunlitDuration)/3600),'Standard Deviation: ' std(sunlitDuration)/3600,'Variance: ' var(sunlitDuration)/3600};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

%% Probability Distribution Plots for COE
% Setting up figure
f3 = figure(3);
f3.Name = 'Probability Distributions COE';
f3.WindowState = 'maximized';
hold on;

% Eccentricity
subplot(2,3,1);
plot(fitdist(disEcc,'Normal'));
title('Probability Curve of Eccentricity');
xlabel('Eccentricity');
ylabel('Probability Density');
xlim([0 1]);

% Semi-Major Axis
subplot(2,3,2);
plot(fitdist(disSMA,'Normal'));
title('Probability Curve of Semi-Major Axis');
xlabel('Semi-Major Axis [km]');
ylabel('Probability Density');

% Inclination
subplot(2,3,3);
plot(fitdist(disInc,'Normal'));
title('Probability Curve of Inclination');
xlabel('Inclination [deg]');
ylabel('Probability Density');

% Right Ascension of Ascending Node
subplot(2,3,4);
plot(fitdist(disRAAN,'Normal'));
title('Probability Curve of Right Ascension of Ascending Node');
xlabel('Right Ascension of Ascending Node [deg]');
ylabel('Probability Density');

% Argument of Perigee
subplot(2,3,5);
plot(fitdist(disArgPer,'Normal'));
title('Probability Curve of Argument of Perigee');
xlabel('Argument of Perigee [deg]');
ylabel('Probability Density');

% Orbital Period
subplot(2,3,6);
plot(fitdist(disT./3600,'Normal'));
title('Probability Curve of Orbital Period');
xlabel('Orbital Period [hr]');
ylabel('Probability Density');




%% Notes
%{
Todo: 
- Animation

%}
%% Unpurturbed Orbit Function
function Xdot = OrbitPropagator(~,X,muEarth)
    % Pulling position from input state vector
    r = [X(1);X(2);X(3)];

    % Gravitational force from Earth
    a = -(muEarth/(vecnorm(r)^3)).*r;

    % Sending new velocity & acceleration data into loop
    Xdot = [X(4);X(5);X(6);a(1);a(2);a(3)];
end














