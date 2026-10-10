function fin = calculateFinProperties(fin)
    t = fin.thickness;
    cr = fin.rootchord;
    ct = fin.tipchord;
    fin.c = (ct + cr) / 2; % fin avg. chord
    h = fin.span;
    fin.b = fin.c./2; % fin avg. semichord
    beta1 = 1.87510407;
    %sweep = fin.sweep;
    if strcmp(fin.airfoil,'rectangular') % rectangular/flat cross-section
        fin.planformArea = fin.c .* h;
        fin.volume = fin.planformArea .* t;
        fin.J0 = fin.c .* t .* (fin.c.^2 + t.^2) ./ 12; % polar moment of inertia - rectangular airfoil simplification
        fin.r_bar = sqrt((1 + (t ./ fin.c).^2) ./ 3); % reduced radius of gyration
        Ix = fin.c * t^3 / 12; % fin x area moment of inertia
        fin.foilArea = fin.c .* t;
    elseif strcmp(fin.airfoil,'hexagonal')
        hc = fin.chamferHeight;
        Achamf = t .* hc ./ 2;
        xbarchamf = hc ./ 3;
        Ixchamf = hc .* t.^3 ./ 48;
        Iychamf = t .* hc.^3 ./ 36;
        rectbase = fin.c - (2 .* hc); % avg rectangle section width
        Ixrect = rectbase .* t.^3 ./ 12;
        Iyrect = t .* rectbase.^3 ./ 12;
        d = xbarchamf + (rectbase ./ 2);
        Ix = Ixrect + (2 .* (Ixchamf + (Achamf .* d)));
        Iy = Iyrect + (2 .* Iychamf);
        fin.volume = (2 .* Achamf .* h) + (((cr - (2 .* hc)) + (ct - (2 .* hc))) .* h .* t ./ 2);
        fin.J0 = Ix + Iy;
        fin.r_bar = sqrt(h .* fin.J0 ./ (fin.b.^2 .* fin.volume));
        fin.foilArea = (2 .* Achamf) + (fin.c .* t);
    elseif strcmp(fin.airfoil,'diamond')
        hc = fin.b;
        Achamf = t .* hc ./ 2;
        xbarchamf = hc ./ 3;
        Ixchamf = hc .* t.^3 ./ 48;
        Iychamf = t .* hc.^3 ./ 36;
        d = xbarchamf;
        Ix = 2 .* (Ixchamf + (Achamf .* d));
        Iy = 2 .* Iychamf;
        %fin.midspan = sqrt(h.^2 + (sweep + (ct ./ 2) - (cr ./ 2)).^2);
        fin.volume = 0.25 .* t .* h .* (cr + ct);
        fin.J0 = Ix + Iy;
        fin.r_bar = sqrt(h .* fin.J0 ./ (fin.b.^2 .* fin.volume));
        fin.foilArea = 2 .* Achamf;
    elseif strcmp(fin.airfoil,'biconvex')
        %fin.midspan = sqrt(h.^2 + (sweep + (ct ./ 2) - (cr ./ 2)).^2);
        fin.volume = t .* h .* (cr + ct) ./ 3;
        fin.J0 = fin.c .* t .* ((t.^2 ./ 6) + (fin.c.^2 ./ 12));
        fin.r_bar = sqrt(h .* fin.J0 ./ (fin.b.^2 .* fin.volume));
        %TODO: MAKE THIS COMPATIBLE WITH FREQUENCY EQUATIONS
    end
    fin.m = fin.volume .* fin.density ./ h;
    if (isfield(fin,'E'))
        fin.freq_h = beta1^2 .* sqrt(fin.E .* Ix ./ (fin.m .* h.^4));
        if(isfield(fin,'poisson'))
            fin.G = fin.E ./ (2 .* (1 + fin.poisson));
        end
    end
    if(isfield(fin,'G'))
        if strcmp(fin.airfoil,'rectangular')
            K_tors = (fin.c ./ 2) .* (t ./ 2).^3 .* ((16 / 3) - (3.36 .* (t ./ fin.c) .*(1 - (t.^4 ./ (12 .* fin.c.^4)))));
        else
            K_tors = 4 .* Ix ./ (1 + (16 .* (Ix ./ (fin.foilArea .* fin.c^2))));
        end
        fin.freq_alpha = pi() .* sqrt(K_tors .* fin.G ./ (h^2 .* fin.density .* fin.J0)) / 2;
    end
end