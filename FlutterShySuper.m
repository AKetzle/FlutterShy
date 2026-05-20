%% FlutteryShy Super (Supersonic Flutter Prediction)
% Alexander Ketzle, Written for the Mississippi State University Space Cowboys and the benefit of the rocketry community
% Written April 2026
% Based on the methods by J. P. Kearns, 1962 and Theodorsen and Garrick, 1940
clc, clear, close all;

%% User Inputs

% fin parameters
c = 0; % fin average chord, ft
m = 0; % mass per unit span, slug / ft (span)
x_bar = 0.5; % chord-normalized distance from c.g. to leading edge, per chord
r_bar = 0; % semichord-normalized radius of gyration, % semichord
freq_h = 0; % bending frequency, rad/s
freq_alpha = 0; % torsion frequency, rad/s
a_h = 0.5; % chord-normalized distance from elastic axis to leading edge, per chord
g_h = 0.000; % bending damping ratio
g_alpha = 0.000; % torsion damping ratio

% simulation parameters
site_altitude = 0; % altitude of launch site above sea level (MUST MATCH RAS SIM), feet
RAS_Filepath = "FILE PATH HERE"; % filepath to the RASAERO II flight sim file (.csv)

% advanced simulation controls
invkstepsize = 0.0001; % increasing resolution exponentially increases calculation time
invkMax = 8; % max 1/k value to calc to
machGate = 1.01; % Don't change this unless you know what you're doing
%% Calculation

b = c / 2; % average semi-chord, ft
RAS_DATA = readmatrix(RAS_Filepath,"NumHeaderLines",1);
RAS_Alt1 = RAS_DATA(:,23);
[apogee, apoidx] = max(RAS_Alt1);
fprintf("Maximum RASAero II Altitude: %g ft\n",apogee);
RAS_Alt = RAS_DATA(1:apoidx,23);
RAS_Time = RAS_DATA(1:apoidx,1);
apoTime = RAS_Time(apoidx);
RAS_Mach = RAS_DATA(1:apoidx,4);
RAS_Vel = RAS_DATA(1:apoidx,18);
% Note to self: Just use atmos.m, not the calibrated one, it's closer to what RAS is saying
[rho,~,T,a,~] = atmos(RAS_Alt + site_altitude); % get the atmospheric properties at each RAS time step
mu = m ./ (pi() .* rho .* b.^2); % mass ratio parameter
q = 0.5 .* rho .* RAS_Vel.^2;

% Supersonic
V_f_sup = kearnsSupersonic(mu, r_bar, RAS_Mach, x_bar, b, freq_h, freq_alpha, machGate);

% Subsonic
V_f_sub = zeros(size(mu));
iters = size(mu,1);
for i = 1:iters
    V_f_sub(i) = TR496TR685(freq_alpha, freq_h, a_h, x_bar, r_bar, b, mu(i), invkstepsize, invkMax, g_h, g_alpha);
    V_f_sub(i) = V_f_sub(i) .*(RAS_Mach(i)<=machGate);
end
M_f_sub1 = V_f_sub ./ a;
% note: change 1/k correction to inline math eq for accuracy's sake
%M_f_sub = sqrt(M_f_sub1.^2 .* (sqrt(1 - (M_f_sub1.^4 ./ 4)) - (M_f_sub1.^2 ./ 2))); % subsonic vel calc from tr685
%M_f_sub = sqrt(M_f_sub1.^2 .* (sqrt(4 + (M_f_sub1.^4)) - (M_f_sub1.^2))) ./ sqrt(2); % supersonic vel calc derived via matlab
M_f_sub = (sqrt(M_f_sub1.^2 .* (sqrt(4 + (M_f_sub1.^4)) - (M_f_sub1.^2))) ./ sqrt(2) .* (M_f_sub1>=1)) + (sqrt(M_f_sub1.^2 .* (sqrt(1 - (M_f_sub1.^4 ./ 4)) - (M_f_sub1.^2 ./ 2))) .* (M_f_sub1<1));
M_f_sup = V_f_sup ./ a;
V_f = V_f_sub + V_f_sup;
M_f = M_f_sub + M_f_sup;
%% Plotting
close all;
fs_flutter = (V_f ./ abs(RAS_Vel)) .* (RAS_Vel > 50);
fs_flutter(fs_flutter == 0) = NaN;
V_f2 = V_f .* (RAS_Vel > 50);
V_f2(V_f2 == 0) = NaN;

[fsmin, fsminidx] = min(fs_flutter); % minimum flutter f.s.
minfs_fluttervel = V_f2(fsminidx); % flutter velocity at min f.s.
minfs_rasvel = RAS_Vel(fsminidx); % RAS velocity at min f.s.
minfs_time = RAS_Time(fsminidx); % time of min flutter f.s.

[vfmin, vfminidx] = min(V_f2); % minimum flutter velocity
minfluttervel_rasvel = RAS_Vel(vfminidx); % ras velocity at minimum flutter velocity
minfluttervel_fs = fs_flutter(vfminidx); % f.s. of min flutter velocity
minfluttervel_time = RAS_Time(vfminidx); % time of minimum flutter velocity

[rasvelmax, maxrasvelidx] = max(RAS_Vel); % maximum rasaero velocity
maxrasvel_fluttervel = V_f2(maxrasvelidx); % flutter vel at max ras velocity
maxrasvel_fs = fs_flutter(maxrasvelidx); % flutter f.s. at amx ras velocity

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
plot(RAS_Time,fs_flutter,LineWidth=1.5);
yline(1.5,LineWidth=1.5);
yline(1,LineWidth=1.5);
xlim([0,apoTime])
ylim([0,inf])
xlabel("Time (s)");
title("Flutter Factor of Safety vs. Time")
ylabel("V_f / RAS Velocity");
legend("Velocity Ratio");
fontsize(16,"points");

figure(Name="Flutter Velocity vs. Sim Velocity");
plot(RAS_Vel,V_f);
xlabel("Simulation Velocity (ft/s)");
ylabel("Flutter Velocity (ft/s)");
title("Flutter Velocity vs. Sim Velocity");
fontsize(16,"points");

figure(Name="Flutter Mach vs. Sim Mach");
plot(RAS_Mach,M_f);
xline(1);
xline(0.8);
xline(1.2);
xlabel("Sim Mach");
ylabel("Flutter Mach");
title("Flutter Mach vs. Sim Mach");
fontsize(16,"points");

figure(Name="Mach Numbers vs. Time");
plot(RAS_Time,RAS_Mach,RAS_Time,M_f)
yline(1.2)
legend("Sim Mach Number","Flutter Mach Number");
xlabel("Time (s)");
ylabel("Mach Number");
title("Mach Numbers vs. Time");
xlim([0,apoTime])
fontsize(16,"points");

figure(Name="Flutter Factor of Safety vs. Dynamic Pressure");
plot(q,fs_flutter)
ylim([0,inf])
xlabel("Dynamic Pressure (lbf/ft^2)");
ylabel("Flutter Factor of Safety");
title("Flutter Factor of Safety vs. Dynamic Pressure");
fontsize(16,"points");

figure(Name="Flutter Factor of Safety vs. Altitude");
plot(RAS_Alt,fs_flutter);
ylim([0,inf])
yline(1.5)
yline(1)
xlabel("Altitude (ft)");
ylabel("Flutter Factor of Safety");
title("Flutter Factor of Safety vs. Altitude");
fontsize(16,"points");

figure(Name="Flutter FS vs. Sim Mach");
plot(RAS_Mach,fs_flutter);
xline(1);
xline(0.8);
xline(1.2);
xlabel("Sim Mach");
ylabel("Flutter FS");
title("Flutter FS vs. Sim Mach");
fontsize(16,"points");

figure(Name="Flutter Velocity vs. Time");
plot(RAS_Time,V_f,LineWidth=1.5);
xlim([0,apoTime])
ylim([0,inf])
xlabel("Time (s)");
title("Flutter Velocity vs. Time")
ylabel("V_f (ft/s)");
fontsize(16,"points");

%% Helper Functions

function [Uf, flutterPoint] = TR496TR685(freq_alpha, freq_h, a_h, x_alpha, r_alpha, b, mu, invkstepsize, invkMax, g_h, g_alpha)
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
    i = sqrt(-1);
    invkrange = [invkstepsize,invkMax];
    n = uint32(((invkrange(2) - invkrange(1)) ./ invkstepsize) + 1);
    invk_set = linspace(invkrange(1),invkrange(2),n);
    k_set = 1 ./ invk_set;
    
    k = k_set;
    k_inv = invk_set;
    
    Ch_k = besselh(1,2,k) ./ (besselh(1,2,k) + (i .* besselh(0,2,k))); % Theodorsen Function; Less lines to compute than using the other bessel functions
    F = real(Ch_k);
    G = imag(Ch_k);
    
    A_R = -(mu + 1) - (2 .* G ./ k);
    A_I = 2 .* F ./ k;
    B_R = -((mu .* x_alpha) - a_h) + (2 .* F ./ k.^2) - ((0.5 - a_h) .* 2 .* G ./ k);
    B_I = (1 ./ k) .* (1 + (2 .* G ./ k) + ((0.5 - a_h) .* 2 .* F));
    D_R = -((mu .* x_alpha) - a_h) + ((0.5 + a_h) .* 2 .* G ./ k);
    D_I = -(0.5 + a_h) .* 2 .* F ./ k;
    E_R = -((mu .* r_alpha.^2) + a_h.^2 + 0.125) + ((0.25 - a_h.^2) .* 2 .* G ./ k) - ((0.5 + a_h) .* 2 .* F ./ k.^2);
    E_I = (1 ./ k) .* ((0.5 - a_h) - ((0.5 + a_h) * 2 .* G ./ k) - ((0.25 - a_h.^2) * 2 .* F));
    
    delta_R_A = (1 - (g_h * g_alpha)) * mu.^2 * r_alpha.^2 * freq_h.^2 ./ freq_alpha.^2;
    delta_R_B = (mu .* freq_h.^2 ./ freq_alpha.^2 .* (E_R - (g_h .* E_I))) + (mu .* r_alpha.^2 .* (A_R - (g_alpha .* A_I)));
    delta_R_C = (A_R .* E_R) - (B_R .* D_R) - (A_I .* E_I) + (B_I .* D_I);
    
    delta_I_A = (g_h + g_alpha) .* mu.^2 .* r_alpha.^2 .* freq_h.^2 ./ freq_alpha.^2;
    delta_I_B = (mu .* freq_h.^2 ./ freq_alpha.^2 .* ((g_h .* E_R) + E_I)) + (mu .* r_alpha.^2 .* (A_I + (g_alpha .* A_R)));
    delta_I_C = (A_I .* E_R) - (B_R .* D_I) + (A_R .* E_I) - (B_I .* D_R);
    
    X_R1 = (-delta_R_B - sqrt(delta_R_B.^2 - (4 .* delta_R_A .* delta_R_C))) ./ (2 .* delta_R_A);
    X_R2 = (-delta_R_B + sqrt(delta_R_B.^2 - (4 .* delta_R_A .* delta_R_C))) ./ (2 .* delta_R_A);

    if(~delta_I_A == 0) % change this to an inline math eq when you figure out mu's deal
        X_I1 = (-delta_I_B - sqrt(delta_I_B.^2 - (4 .* delta_I_A .* delta_I_C))) ./ (2 .* delta_I_A);
        X_I2 = (-delta_I_B + sqrt(delta_I_B.^2 - (4 .* delta_I_A .* delta_I_C))) ./ (2 .* delta_I_A);
    else
        X_I1 = -delta_I_C ./ delta_I_B;
        X_I2 = -delta_I_C ./ delta_I_B;
    end

    iscomplex = (imag(X_R1) ~= 0);
    X_R1(iscomplex) = NaN;
    iscomplex = (imag(X_R2) ~= 0);
    X_R2(iscomplex) = NaN;

    rt_X_R1 = sqrt(X_R1);
    rt_X_R2 = sqrt(X_R2);
    rt_X_I1 = sqrt(X_I1);
    iscomplex = (imag(rt_X_I1) ~= 0);
    rt_X_I1(iscomplex) = NaN;
    rt_X_I2 = sqrt(X_I2);
    iscomplex = (imag(rt_X_I2) ~= 0);
    rt_X_I2(iscomplex) = NaN;
    
    solutionmatrix = [k; k_inv; X_R1; X_R2; X_I1; rt_X_R1; rt_X_R2; rt_X_I1; rt_X_I2];
    
    XRatio1 = (abs(1 - abs(solutionmatrix(6,:) ./ solutionmatrix(8,:))) .* (~isnan(rt_X_I1))) + (abs(1 - abs(solutionmatrix(6,:) ./ solutionmatrix(9,:))) .* (isnan(rt_X_I1)));
    XRatio2 = (abs(1 - abs(solutionmatrix(7,:) ./ solutionmatrix(8,:))) .* (~isnan(rt_X_I1))) + (abs(1 - abs(solutionmatrix(7,:) ./ solutionmatrix(9,:))) .* (isnan(rt_X_I1)));
    
    % if (~isnan(rt_X_I1)) % change this to an inline math eq
    %     XRatio1 = abs(1 - abs(solutionmatrix(6,:) ./ solutionmatrix(8,:)));
    %     XRatio2 = abs(1 - abs(solutionmatrix(7,:) ./ solutionmatrix(8,:)));
    % else
    %     XRatio1 = abs(1 - abs(solutionmatrix(6,:) ./ solutionmatrix(9,:)));
    %     XRatio2 = abs(1 - abs(solutionmatrix(7,:) ./ solutionmatrix(9,:)));
    % end
    
    [~,idx1] = find(XRatio1 == min(XRatio1,[],"omitnan"));
    [~,idx2] = find(XRatio2 == min(XRatio2,[],"omitnan"));
    r1 = XRatio1(:,idx1);
    r2 = XRatio2(:,idx2);
    if(r1 < r2)
        flutterPoint = solutionmatrix(:,idx1);
        Uf = freq_alpha .* b .* flutterPoint(2) ./ flutterPoint(6);
    else
        flutterPoint = solutionmatrix(:,idx2);
        Uf = freq_alpha .* b .* flutterPoint(2) ./ flutterPoint(7);
    end
end

function V_f_sup = kearnsSupersonic(mu, r_bar, RAS_Mach, x_bar, b, freq_h, freq_alpha, machGate)
    sqrt1 = mu .* r_bar.^2 .* sqrt((RAS_Mach>machGate).*(RAS_Mach.^2 - 1)) ./ (x_bar .* b);
    sqrt2N = (1 - (freq_h ./ freq_alpha).^2).^2 + (4 .* (x_bar ./ r_bar).^2 .* (freq_h ./ freq_alpha).^2);
    sqrt2D = 2 * (1 + (freq_h ./ freq_alpha).^2);
    V_f_sup = freq_alpha .* b .* sqrt(sqrt1 .* (sqrt2N ./ sqrt2D));
end

%% Todo list
%{ 
Accept multiple unit systems
fix subsonic to use vector mu so that for loop can go away
ensure you can have multiple input vectors so this can be used in optimization stuff
tear things apart from the file so you can use this without a ras sim and just a set of conditions
%}