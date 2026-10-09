function V_f_sup = kearnsSupersonic(mu, r_bar, mach, x_bar, b, freq_h, freq_alpha, machGate)
    sqrt1 = mu .* r_bar.^2 .* sqrt((mach>machGate).*(mach.^2 - 1)) ./ (x_bar .* b);
    sqrt2N = (1 - (freq_h ./ freq_alpha).^2).^2 + (4 .* (x_bar ./ r_bar).^2 .* (freq_h ./ freq_alpha).^2);
    sqrt2D = 2 * (1 + (freq_h ./ freq_alpha).^2);
    V_f_sup = freq_alpha .* b .* sqrt(sqrt1 .* (sqrt2N ./ sqrt2D));
end