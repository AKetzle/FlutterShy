function fin = calculateFinProperties(fin)
    t = fin.thickness;
    cr = fin.rootchord;
    ct = fin.tipchord;
    fin.c = (ct + cr) / 2; % fin avg. chord
    h = fin.span;
    fin.b = fin.c./2; % fin avg. semichord
    %sweep = fin.sweep;
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
    elseif strcmp(fin.airfoil,'biconvex')
        %fin.midspan = sqrt(h.^2 + (sweep + (ct ./ 2) - (cr ./ 2)).^2);
        fin.volume = t .* h .* (cr + ct) ./ 3;
        fin.J0 = fin.c .* t .* ((t.^2 ./ 6) + (fin.c.^2 ./ 12));
        fin.r_bar = sqrt(h .* fin.J0 ./ (fin.b.^2 .* fin.volume));
    end
    fin.m = fin.volume .* fin.density / h;
end