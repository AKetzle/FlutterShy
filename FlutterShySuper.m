%% FlutterShy Super (Supersonic Flutter Prediction)
% Alexander Ketzle, Written for the Mississippi State University Space Cowboys and the benefit of the rocketry community
% First written April 2026
% Last updated: Oct 9 2026
% Based on the methods by J. P. Kearns, 1962 and Theodorsen and Garrick, 1940
clc, clear, close all;

%% User Inputs

inputFile = "InputExample.txt"; % path to FlutterShy Input File

inputTable = cell2table(readcell(inputFile,"Delimiter",' = '));
inputTable = rows2vars(inputTable,"VariableNamesSource",1);

% fin parameters
c = cell2mat(inputTable.c); % fin average chord, ft
m = cell2mat(inputTable.m); % mass per unit span, slug / ft (span)
parameters.x_bar = cell2mat(inputTable.x_bar); % chord-normalized distance from c.g. to elastic axis, % chord
parameters.r_bar = cell2mat(inputTable.r_bar); % chord-normalized radius of gyration, % chord
parameters.freq_h = cell2mat(inputTable.freq_h); % bending frequency, rad/s
parameters.freq_alpha = cell2mat(inputTable.freq_alpha); % torsion frequency, rad/s
parameters.a_h = cell2mat(inputTable.a_h);
parameters.g_h = cell2mat(inputTable.g_h);
parameters.g_alpha = cell2mat(inputTable.g_alpha);

% simulation parameters
site_altitude = cell2mat(inputTable.site_altitude); % altitude of launch site above sea level (MUST MATCH RAS SIM), feet
RAS_Filepath = cell2mat(inputTable.RAS_Filepath); % filepath to the RASAERO II flight sim file (.csv)

% advanced simulation controls
if (~ismember('invkstepsize',inputTable.Properties.VariableNames))
    parameters.invkstepsize = 0.0001;
else
    parameters.invkstepsize = cell2mat(inputTable.invkstepsize); % increasing resolution exponentially increases calculation time
end

if (~ismember('invkMax',inputTable.Properties.VariableNames))
    parameters.invkMax = 8;
else
    parameters.invkMax = cell2mat(inputTable.invkMax); % max 1/k value to calc to
end

if (~ismember('machGate',inputTable.Properties.VariableNames))
    parameters.machGate = 1.01;
else
    parameters.machGate = cell2mat(inputTable.machGate); % Don't change this unless you know what you're doing
end

if (~ismember('subsonicCorrection',inputTable.Properties.VariableNames))
    parameters.subsonicCorrection = 'none';
else
    parameters.subsonicCorrection = cell2mat(inputTable.subsonic_model);
end
%% Calculation

parameters.b = c / 2; % average semi-chord, ft

RasData = readRASData(RAS_Filepath);
Altitude = RasData.Altitude;
parameters.velocity = RasData.Velocity;
parameters.mach = RasData.Mach;
Time = RasData.Time;
ApogeeTime = RasData.ApogeeTime;

% Note to self: Just use atmos.m, not the calibrated one, it's closer to what RAS is saying
[parameters.rho,~,~,parameters.a,~] = atmos(Altitude + site_altitude); % get the atmospheric properties at each RAS time step
parameters.mu = m ./ (pi() .* parameters.rho .* parameters.b.^2); % mass ratio parameter

FlutterShyResults = FlutterShy(parameters);

%% Report & Plotting

[fsmin, fsminidx] = min(FlutterShyResults.fs_flutter); % minimum flutter f.s.
minfs_fluttervel = FlutterShyResults.V_f2(fsminidx); % flutter velocity at min f.s.
minfs_rasvel = parameters.velocity(fsminidx); % RAS velocity at min f.s.
minfs_time = Time(fsminidx); % time of min flutter f.s.

[vfmin, vfminidx] = min(FlutterShyResults.V_f2); % minimum flutter velocity
minfluttervel_rasvel = parameters.velocity(vfminidx); % ras velocity at minimum flutter velocity
minfluttervel_fs = FlutterShyResults.fs_flutter(vfminidx); % f.s. of min flutter velocity
minfluttervel_time = Time(vfminidx); % time of minimum flutter velocity

[rasvelmax, maxrasvelidx] = max(parameters.velocity); % maximum rasaero velocity
maxrasvel_fluttervel = FlutterShyResults.V_f2(maxrasvelidx); % flutter vel at max ras velocity
maxrasvel_fs = FlutterShyResults.fs_flutter(maxrasvelidx); % flutter f.s. at amx ras velocity


q = 0.5 .* parameters.rho .* parameters.velocity.^2;
% the report

fprintf("Minimum Flutter F.S.:                     %g\n" + ...
        "Flutter Velocity at Minimum F.S.:         %g ft/s\n" + ...
        "RAS Velocity at Minimum F.S.:             %g ft/s\n" + ...
        "Time of Minimum F.S.:                     %g seconds\n" + ...
        "====================================================\n" + ...
        "Minimum Flutter Velocity:                 %g ft/s\n" + ...
        "RAS Velocity at Minimum Flutter Velocity: %g ft/s\n" + ...
        "F.S. at Minimum Flutter Velocity:         %g\n" + ...
        "Time of Minimum Flutter Velocity:         %g seconds\n" + ...
        "====================================================\n" + ...
        "Maximum RAS Velocity:                     %g ft/s\n" + ...
        "Flutter Velocity at Maximum RAS Velocity: %g ft/s\n" + ...
        "Flutter F.S. at Maximum RAS Velocity:     %g \n",fsmin,minfs_fluttervel,minfs_rasvel,minfs_time,vfmin,minfluttervel_rasvel,minfluttervel_fs,minfluttervel_time,rasvelmax,maxrasvel_fluttervel,maxrasvel_fs);


figure(Name="Flutter Factor of Safety vs. Time");
plot(Time,FlutterShyResults.fs_flutter,LineWidth=1.5);
yline(1.5,LineWidth=1.5);
yline(1,LineWidth=1.5);
xlim([0,ApogeeTime])
ylim([0,inf])
xlabel("Time (s)");
title("Flutter Factor of Safety vs. Time")
ylabel("V_f / Sim Velocity");
legend("Velocity Ratio");
fontsize(16,"points");

figure(Name="Flutter Velocity vs. Sim Velocity");
plot(parameters.velocity,FlutterShyResults.V_f);
xlabel("Simulation Velocity (ft/s)");
ylabel("Flutter Velocity (ft/s)");
title("Flutter Velocity vs. Sim Velocity");
fontsize(16,"points");

figure(Name="Flutter Mach vs. Sim Mach");
plot(parameters.mach,FlutterShyResults.M_f);
xline(1);
xline(0.8);
xline(1.2);
xlabel("Sim Mach");
ylabel("Flutter Mach");
title("Flutter Mach vs. Sim Mach");
fontsize(16,"points");

figure(Name="Mach Numbers vs. Time");
plot(Time,parameters.mach,Time,FlutterShyResults.M_f)
yline(1.2)
legend("Sim Mach Number","Flutter Mach Number");
xlabel("Time (s)");
ylabel("Mach Number");
title("Mach Numbers vs. Time");
xlim([0,ApogeeTime])
fontsize(16,"points");

figure(Name="Flutter Factor of Safety vs. Dynamic Pressure");
plot(q,FlutterShyResults.fs_flutter)
ylim([0,inf])
xlabel("Dynamic Pressure (lbf/ft^2)");
ylabel("Flutter Factor of Safety");
title("Flutter Factor of Safety vs. Dynamic Pressure");
fontsize(16,"points");

figure(Name="Flutter Factor of Safety vs. Altitude");
plot(Altitude,FlutterShyResults.fs_flutter);
ylim([0,inf])
yline(1.5)
yline(1)
xlabel("Altitude (ft)");
ylabel("Flutter Factor of Safety");
title("Flutter Factor of Safety vs. Altitude");
fontsize(16,"points");

figure(Name="Flutter FS vs. Sim Mach");
plot(parameters.mach,FlutterShyResults.fs_flutter);
xline(1);
xline(0.8);
xline(1.2);
xlabel("Sim Mach");
ylabel("Flutter FS");
title("Flutter FS vs. Sim Mach");
fontsize(16,"points");

figure(Name="Flutter Velocity vs. Time");
plot(Time,FlutterShyResults.V_f,LineWidth=1.5);
xlim([0,ApogeeTime])
ylim([0,inf])
xlabel("Time (s)");
title("Flutter Velocity vs. Time")
ylabel("V_f (ft/s)");
fontsize(16,"points");

%% Helper Functions

function FlutterShyResults = FlutterShy(parameters)
    mu = parameters.mu;
    r_bar = parameters.r_bar;
    Mach = parameters.mach;
    x_bar = parameters.x_bar;
    b = parameters.b;
    freq_h = parameters.freq_h;
    freq_alpha = parameters.freq_alpha;
    g_h = parameters.g_h;
    g_alpha = parameters.g_alpha;
    machGate = parameters.machGate;
    invkstepsize = parameters.invkstepsize;
    invkMax = parameters.invkMax;
    Velocity = parameters.velocity;
    a_h = parameters.a_h;
    a = parameters.a;

    % Supersonic
    V_f_sup = kearnsSupersonic(mu, r_bar, Mach, x_bar, b, freq_h, freq_alpha, machGate);

    % Subsonic
    i = sqrt(-1);
    invkrange = [invkstepsize,invkMax];
    n = uint32(((invkrange(2) - invkrange(1)) ./ invkstepsize) + 1);
    invk = linspace(invkrange(1),invkrange(2),n);
    k = 1 ./ invk;
    Ch_k = besselh(1,2,k) ./ (besselh(1,2,k) + (i .* besselh(0,2,k))); % Theodorsen Function; Less lines to compute than using the other bessel functions
    F = real(Ch_k);
    G = imag(Ch_k);
    G2_k = 2 .* G ./ k;

    V_f_sub = zeros(size(mu));
    iters = size(mu,1);
    for j = 1:iters
        V_f_sub(j) = TR496TR685(freq_alpha, freq_h, a_h, x_bar, r_bar, b, mu(j), F, g_h, g_alpha, k, G2_k,invk);
        V_f_sub(j) = V_f_sub(j) .*(Mach(j)<=machGate);
    end

    M_f_sub1 = V_f_sub ./ a;
    % note: change 1/k correction to inline math eq for accuracy's sake
    % note: the supersonic correction may actually be hurting here. needs verification. possibly remove?
    if strcmp(parameters.subsonicCorrection,'none')
        M_f_sub = M_f_sub1;
    elseif strcmp(parameters.subsonicCorrection,'tr685')
        M_f_sub = sqrt(M_f_sub1.^2 .* (sqrt(1 - (M_f_sub1.^4 ./ 4)) - (M_f_sub1.^2 ./ 2))); % subsonic vel calc from tr685
    elseif strcmp(parameters.subsonicCorrection,'mat1')
        M_f_sub = sqrt(M_f_sub1.^2 .* (sqrt(4 + (M_f_sub1.^4)) - (M_f_sub1.^2))) ./ sqrt(2); % first supersonic vel calc derived via matlab
    elseif strcmp(parameters.subsonicCorrection,'mat2')
        M_f_sub = (sqrt(M_f_sub1.^2 .* (sqrt(4 + (M_f_sub1.^4)) - (M_f_sub1.^2))) ./ sqrt(2) .* (M_f_sub1>=1)) + (sqrt(M_f_sub1.^2 .* (sqrt(1 - (M_f_sub1.^4 ./ 4)) - (M_f_sub1.^2 ./ 2))) .* (M_f_sub1<1)); % second supersonic vel calc derived via matlab
    elseif strcmp(parameters.subsonicCorrection,'kearns')
        M_f_sub = kearnsSupersonic(mu, r_bar, M_f_sub1, x_bar, b, freq_h, freq_alpha, machGate) ./ a;
    end
    V_f_sub = M_f_sub .* a;
    M_f_sup = V_f_sup ./ a;
    FlutterShyResults.M_f_sub = M_f_sub;
    FlutterShyResults.V_f = V_f_sub + V_f_sup;
    FlutterShyResults.M_f = M_f_sub + M_f_sup;
    FlutterShyResults.fs_flutter = (FlutterShyResults.V_f ./ abs(Velocity)) .* (Velocity > 50);
    FlutterShyResults.fs_flutter(FlutterShyResults.fs_flutter == 0) = NaN;
    FlutterShyResults.V_f2 = FlutterShyResults.V_f .* (Velocity > 50);
    FlutterShyResults.V_f2(FlutterShyResults.V_f2 == 0) = NaN;
end


%% Todo list
%{
NOT STARTED - Accept multiple unit systems
STARTED - fix subsonic to use vector mu so that for loop can go away - NOTE: this has evolved
NOT STARTED - ensure you can have multiple input vectors so this can be used in optimization stuff
STARTED - tear things apart from the file so you can use this without a ras sim and just a set of conditions - NOTE: almost there
STARTED - allow input of fin physical parameters to calculate needed parameters - only thing keeping this from complete is natural frequencies
NOT STARTED - 3D plots? might be neat
STARTED - Check/Ensure GNU Octave Compatibility - GNU Octave doesn't implement readcell().

COMPLETE - Add toggle between subsonic compressibility corrections for user choice
COMPLETE - allow input of file for analysis parameters
%}