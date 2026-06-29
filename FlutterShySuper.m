%% FlutterShy Super (Supersonic Flutter Prediction)
% Alexander Ketzle, Written for the Mississippi State University Space Cowboys and the benefit of the rocketry community
% First written April 2026
% Last updated: May 22 2026
% Based on the methods by J. P. Kearns, 1962 and Theodorsen and Garrick, 1940
clc, clear, close all;

%% User Inputs

inputFile = "FILE PATH HERE"; % path to FlutterShy Input File

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
parameters.invkstepsize = cell2mat(inputTable.invkstepsize); % increasing resolution exponentially increases calculation time
parameters.invkMax = cell2mat(inputTable.invkMax); % max 1/k value to calc to
parameters.machGate = cell2mat(inputTable.machGate); % Don't change this unless you know what you're doing
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
    V_f_sub = zeros(size(mu));
    iters = size(mu,1);
    for i = 1:iters
        V_f_sub(i) = TR496TR685(freq_alpha, freq_h, a_h, x_bar, r_bar, b, mu(i), invkstepsize, invkMax, g_h, g_alpha);
        V_f_sub(i) = V_f_sub(i) .*(Mach(i)<=machGate);
    end
    
    % V_f_sub = TR496TR685(freq_alpha, freq_h, a_h, x_bar, r_bar, b, mu, invkstepsize, invkMax, g_h, g_alpha).*(RAS_Mach<=machGate);
    
    
    
    M_f_sub1 = V_f_sub ./ a;
    % note: change 1/k correction to inline math eq for accuracy's sake
    % note: the supersonic correction may actually be hurting here. needs verification. possibly remove?
    %M_f_sub = sqrt(M_f_sub1.^2 .* (sqrt(1 - (M_f_sub1.^4 ./ 4)) - (M_f_sub1.^2 ./ 2))); % subsonic vel calc from tr685
    %M_f_sub = sqrt(M_f_sub1.^2 .* (sqrt(4 + (M_f_sub1.^4)) - (M_f_sub1.^2))) ./ sqrt(2); % supersonic vel calc derived via matlab
    M_f_sub = (sqrt(M_f_sub1.^2 .* (sqrt(4 + (M_f_sub1.^4)) - (M_f_sub1.^2))) ./ sqrt(2) .* (M_f_sub1>=1)) + (sqrt(M_f_sub1.^2 .* (sqrt(1 - (M_f_sub1.^4 ./ 4)) - (M_f_sub1.^2 ./ 2))) .* (M_f_sub1<1));
    V_f_sub = M_f_sub .* a;
    M_f_sup = V_f_sup ./ a;
    FlutterShyResults.V_f = V_f_sub + V_f_sup;
    FlutterShyResults.M_f = M_f_sub + M_f_sup;
    FlutterShyResults.fs_flutter = (FlutterShyResults.V_f ./ abs(Velocity)) .* (Velocity > 50);
    FlutterShyResults.fs_flutter(FlutterShyResults.fs_flutter == 0) = NaN;
    FlutterShyResults.V_f2 = FlutterShyResults.V_f .* (Velocity > 50);
    FlutterShyResults.V_f2(FlutterShyResults.V_f2 == 0) = NaN;
end

function [Uf] = TR496TR685(freq_alpha, freq_h, a_h, x_bar, r_bar, b, mu, invkstepsize, invkMax, g_h, g_alpha)
    %{
    Calculates flutter velocity based on the sqrt(X) vs 1/k method.
    Originally found in NACA TR496: https://ntrs.nasa.gov/citations/19930090935
    again in NACA TR685: https://ntrs.nasa.gov/citations/19930091762
    and also in Y.C. Fung's "An Introduction to the Theory of
    Aeroelasticity"
    Flutter condition is when the real and imaginary portions of sqrt(X)
    plotted against 1/k cross

    freq_alpha - Natural pitching (pure rotation about E.A.) frequency of fin, rad/s
    freq_h - Natural plunge (pure bending about E.A.) frequency of fin, rad/s
    a_h - location of elastic axis (E.A.) of fin behind mid-chord divided by semichord length, unitless
    x_alpha - location of c.g. behind E.A. as ratio of semichord, unitless
    r_alpha - reduced radius of gyration around E.A. divided by semichord, unitless
    b - semichord length (half of length of fin), feet (can be other length unit, defines velocity as [unit(b)]/s)
    mu - nondimensional mass ratio, unitless
    g_h - plunge damping coefficient, unitless, typ 0.005 (per FinSim)
    g_alpha - pitch damping coefficient, unitless, typ 0.005 (per FinSim)
    invkstepsize - size of difference between discrete points of 1/k, unitless recommended between 0.0001 - 0.000001
    invkMax - maximum value of 1/k to calculate to, unitless, typ 6-14, adjust higher if the parabola is not closed
    %}
    % set up constants
    i = sqrt(-1);
    invkrange = [invkstepsize,invkMax];
    n = uint32(((invkrange(2) - invkrange(1)) ./ invkstepsize) + 1);
    invk = linspace(invkrange(1),invkrange(2),n);
    k = 1 ./ invk;
    Ch_k = besselh(1,2,k) ./ (besselh(1,2,k) + (i .* besselh(0,2,k))); % Theodorsen Function; Less lines to compute than using the other bessel functions
    F = real(Ch_k);
    G = imag(Ch_k);
    
    % some repeated calculations done here
    freqratiosq = freq_h.^2 ./ freq_alpha.^2;
    musq = mu.^2;
    rbarsq = r_bar.^2;
    mu_rbarsq = mu .* rbarsq;
    G2_k = 2 .* G ./ k;

    % calculate the determinant elements
    A_R = -(mu + 1) - (G2_k);
    A_I = 2 .* F ./ k;
    B_R = -((mu .* x_bar) - a_h) + (A_I ./ k) - ((0.5 - a_h) .* G2_k);
    B_I = (1 ./ k) .* (1 + (G2_k) + ((0.5 - a_h) .* 2 .* F));
    D_R = -((mu .* x_bar) - a_h) + ((0.5 + a_h) .* G2_k);
    D_I = -(0.5 + a_h) .* A_I;
    E_R = -((mu_rbarsq) + a_h.^2 + 0.125) + ((0.25 - a_h.^2) .* G2_k) - ((0.5 + a_h) .* A_I ./ k);
    E_I = (1 ./ k) .* ((0.5 - a_h) - ((0.5 + a_h) .* G2_k) - ((0.25 - a_h.^2) * 2 .* F));
    
    % calculate the real components of the determinant
    delta_R_A = (1 - (g_h * g_alpha)) .* musq .* rbarsq .* freqratiosq;
    delta_R_B = (mu .* freqratiosq .* (E_R - (g_h .* E_I))) + (mu_rbarsq .* (A_R - (g_alpha .* A_I)));
    delta_R_C = (A_R .* E_R) - (B_R .* D_R) - (A_I .* E_I) + (B_I .* D_I);
    
    % calculate the imaginary components of the determinant
    delta_I_A = (g_h + g_alpha) .* musq .* rbarsq .* freqratiosq;
    delta_I_B = (mu .* freqratiosq .* ((g_h .* E_R) + E_I)) + (mu_rbarsq .* (A_I + (g_alpha .* A_R)));
    delta_I_C = (A_I .* E_R) - (B_R .* D_I) + (A_R .* E_I) - (B_I .* D_R);
    
    % calculate the real component roots
    realsqrt = sqrt(delta_R_B.^2 - (4 .* delta_R_A .* delta_R_C)); % note: maybe apply complex check here
    X_R1 = (-delta_R_B - realsqrt) ./ (2 .* delta_R_A);
    X_R2 = (-delta_R_B + realsqrt) ./ (2 .* delta_R_A);
    iscomplex = (imag(X_R1) ~= 0);
    X_R1(iscomplex) = NaN;
    iscomplex = (imag(X_R2) ~= 0);
    X_R2(iscomplex) = NaN;

    % calculate the imaginary component roots
    imagsqrt = sqrt(delta_I_B.^2 - (4 .* delta_I_A .* delta_I_C)); % note: maybe also apply a complex check here? idk
    X_I11 = (((-delta_I_B - imagsqrt) ./ (2 .* delta_I_A)) .* (~delta_I_A == 0));
    xnan = isnan(X_I11);
    X_I11(xnan) = 0;
    X_I12 = ((-delta_I_C ./ delta_I_B) .* (delta_I_A == 0));
    xnan = isnan(X_I12);
    X_I12(xnan) = 0;
    X_I1 = X_I11 + X_I12;
    X_I21 = ((-delta_I_B + imagsqrt) ./ (2 .* delta_I_A)).*(~delta_I_A == 0);
    xnan = isnan(X_I21);
    X_I21(xnan) = 0;
    X_I22 = (-delta_I_C ./ delta_I_B) .* (delta_I_A == 0);
    xnan = isnan(X_I22);
    X_I22(xnan) = 0;
    X_I2 = X_I21 + X_I22;

    % calculate sqrt(X) and sanitize bad output
    rt_X_R1 = sqrt(X_R1);
    rt_X_R2 = sqrt(X_R2);
    rt_X_I1 = sqrt(X_I1);
    iscomplex = (imag(rt_X_I1) ~= 0);
    rt_X_I1(iscomplex) = NaN;
    rt_X_I2 = sqrt(X_I2);
    iscomplex = (imag(rt_X_I2) ~= 0);
    rt_X_I2(iscomplex) = NaN;
    
    % calculate the intercepts of the real and imaginary components
    XRatio1 = (abs(1 - abs(rt_X_R1 ./ rt_X_I1)) .* (~isnan(rt_X_I1))) + (abs(1 - abs(rt_X_R1 ./ rt_X_I2)) .* (isnan(rt_X_I1)));
    XRatio2 = (abs(1 - abs(rt_X_R2 ./ rt_X_I1)) .* (~isnan(rt_X_I1))) + (abs(1 - abs(rt_X_R2 ./ rt_X_I2)) .* (isnan(rt_X_I1)));
    [~,idx1] = find(XRatio1 == min(XRatio1,[],2,"omitnan"));
    [~,idx2] = find(XRatio2 == min(XRatio2,[],2,"omitnan"));
    r1 = min(XRatio1,[],2,"omitnan");
    r2 = min(XRatio2,[],2,"omitnan");

    % calculate flutter velocity
    Uf = freq_alpha .* b .* (((invk(idx1).' ./ rt_X_R1(idx1)).*(r1 < r2)) + ((invk(idx2).' ./ rt_X_R2(idx2)).*(r2 < r1)));
end

function V_f_sup = kearnsSupersonic(mu, r_bar, mach, x_bar, b, freq_h, freq_alpha, machGate)
    sqrt1 = mu .* r_bar.^2 .* sqrt((mach>machGate).*(mach.^2 - 1)) ./ (x_bar .* b);
    sqrt2N = (1 - (freq_h ./ freq_alpha).^2).^2 + (4 .* (x_bar ./ r_bar).^2 .* (freq_h ./ freq_alpha).^2);
    sqrt2D = 2 * (1 + (freq_h ./ freq_alpha).^2);
    V_f_sup = freq_alpha .* b .* sqrt(sqrt1 .* (sqrt2N ./ sqrt2D));
end

function fin = calculateFinProperties(fin)
    t = fin.thickness;
    cr = fin.rootchord;
    ct = fin.tipchord;
    fin.c = (ct + cr) / 2; % fin avg. chord
    h = fin.span;
    fin.b = fin.c./2; % fin avg. semichord
    % REFERENCE r_bar = sqrt(h * J0 / (b^2 * fin Volume))
    if strcmp(fin.airfoil,'rectangular') % rectangular/flat cross-section
        fin.planformArea = fin.c .* h;
        fin.volume = fin.planformArea .* t;
        fin.J0 = fin.c .* t .* (fin.c.^2 + t.^2) ./ 12; % polar moment of inertia - rectangular airfoil simplification
        fin.r_bar = sqrt((1 + (t ./ fin.c).^2) ./ 3); % reduced radius of gyration
    elseif strcmp(fin.airfoil,'hexagonal')
        hc = fin.chamferHeight;
        Achamf = t .* hc ./ 2;
        xbarchamf = hc ./ 3;
        Ixchamf = hc .* t.^3 ./ 48;
        Iychamf = t .* hc.^3 ./ 36;
        rectbase = fin.c - (2 .* hc); % avg rectangle section width
        %Arect = rectbase .* t;
        Ixrect = rectbase .* t.^3 ./ 12;
        Iyrect = t .* rectbase.^3 ./ 12;
        d = xbarchamf + (rectbase ./ 2);
        Ix = Ixrect + (2 .* (Ixchamf + (Achamf .* d)));
        Iy = Iyrect + (2 .* Iychamf);
        fin.volume = (2 .* Achamf .* h) + (((cr - (2 .* hc)) + (ct - (2 .* hc))) .* h .* t ./ 2);
        fin.J0 = Ix + Iy;
        fin.r_bar = sqrt(h .* fin.J0 ./ (fin.b.^2 .* fin.volume));
    elseif strcmp(fin.airfoil,'diamond')
        hc = fin.b;
        sweep = fin.sweep;
        Achamf = t .* hc ./ 2;
        xbarchamf = hc ./ 3;
        Ixchamf = hc .* t.^3 ./ 48;
        Iychamf = t .* hc.^3 ./ 36;
        d = xbarchamf;
        Ix = 2 .* (Ixchamf + (Achamf .* d));
        Iy = 2 .* Iychamf;
        fin.midspan = sqrt(h.^2 + (sweep + (ct ./ 2) - (cr ./ 2)).^2);
        fin.volume = 0.25 .* t .* h .* (cr + ct);
        fin.J0 = Ix + Iy;
        fin.r_bar = sqrt(h .* fin.J0 ./ (fin.b.^2 .* fin.volume));
    end
end

function RasData = readRASData(filepath)
    dataMatrix = readmatrix(filepath,"NumHeaderLines",1);
    rasAlt1 = dataMatrix(:,23);
    [RasData.Apogee,apoindex] = max(rasAlt1);
    fprintf("Maximum RASAero II Altitude: %g ft\n",RasData.Apogee);
    RasData.Altitude = dataMatrix(1:apoindex,23);
    RasData.Time = dataMatrix(1:apoindex,1);
    RasData.ApogeeTime = RasData.Time(apoindex);
    RasData.Mach = dataMatrix(1:apoindex,4);
    RasData.Velocity = dataMatrix(1:apoindex,18);
end
%% Todo list
%{ 
NOT STARTED - Accept multiple unit systems
STARTED - fix subsonic to use vector mu so that for loop can go away - NOTE: this has evolved
NOT STARTED - ensure you can have multiple input vectors so this can be used in optimization stuff
STARTED - tear things apart from the file so you can use this without a ras sim and just a set of conditions - NOTE: almost there
STARTED - allow input of fin physical parameters to calculate needed parameters
NOT STARTED - 3D plots? might be neat
NOT STARTED - Check/Ensure GNU Octave Compatibility

COMPLETE - allow input of file for analysis parameters
%}