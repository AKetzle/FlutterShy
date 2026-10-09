function [Uf] = TR496TR685(freq_alpha, freq_h, a_h, x_bar, r_bar, b, mu, F, g_h, g_alpha,k, G2_k,invk)
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

    % some repeated calculations done here
    freqratiosq = freq_h.^2 ./ freq_alpha.^2;
    musq = mu.^2;
    rbarsq = r_bar.^2;
    mu_rbarsq = mu .* rbarsq;

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
    %XRatio1 = (abs(1 - abs(rt_X_R1 ./ rt_X_I1)) .* (~isnan(rt_X_I1))) + (abs(1 - abs(rt_X_R1 ./ rt_X_I2)) .* (isnan(rt_X_I2)));
    %XRatio2 = (abs(1 - abs(rt_X_R2 ./ rt_X_I1)) .* (~isnan(rt_X_I1))) + (abs(1 - abs(rt_X_R2 ./ rt_X_I2)) .* (isnan(rt_X_I2)));
    if (~isnan(rt_X_I1))
        XRatio1 = abs(1 - abs(rt_X_R1 ./ rt_X_I1));
        XRatio2 = abs(1 - abs(rt_X_R2 ./ rt_X_I1));
    else
        XRatio1 = abs(1 - abs(rt_X_R1 ./ rt_X_I2));
        XRatio2 = abs(1 - abs(rt_X_R2 ./ rt_X_I2));
    end
    [~,idx1] = find(XRatio1 == min(XRatio1,[],2,"omitnan"));
    [~,idx2] = find(XRatio2 == min(XRatio2,[],2,"omitnan"));
    r1 = min(XRatio1,[],2,"omitnan");
    r2 = min(XRatio2,[],2,"omitnan");

    % calculate flutter velocity
    Uf = freq_alpha .* b .* (((invk(idx1).' ./ rt_X_R1(idx1)).*(r1 < r2)) + ((invk(idx2).' ./ rt_X_R2(idx2)).*(r2 < r1)));
end
