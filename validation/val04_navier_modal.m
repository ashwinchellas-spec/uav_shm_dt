function R = val04_navier_modal()
%VAL04_NAVIER_MODAL  Free vibration of SS-1 cross-ply plates vs Navier FSDT.
%   (a) symmetric [0/90/90/0], E1/E2 = 40, a/h = 10  - first 6 flexural modes
%   (b) antisymmetric [0/90], E1/E2 = 40, a/h = 10   - tests B-coupling and I1
%   Mesh 30x30 HR, consistent mass incl. rotary inertia I2.  Only modes
%   whose kinetic energy is >50% transverse are compared (flexural modes).
%   Non-dimensional frequency  w_bar = omega a^2 / h sqrt(rho / E2).
m = material_library('reddy_e40'); a = 1; h = 0.1;
R = [];
cases = {[0 90 90 0], 'symmetric [0/90/90/0]'; [0 90], 'antisymmetric [0/90]'};
for c = 1:2
    th = cases{c,1}; n = numel(th);
    L.mat = repmat(m, 1, n); L.theta = th; L.t = h/n*ones(1,n);
    sec = laminate_section(L, 5/6);
    nav = navier_fsdt(sec, a, a, 'modal', 6);
    mesh = rect_mesh(a, a, 30, 30);
    model = apply_bc(assemble_system(mesh, sec, struct('shear','hr')), 'SSSS_crossply');
    [~, Phi, om] = solve_modal(model, 30);
    % keep flexural modes only: transverse (w) share of the kinetic energy.
    % The FE model also has in-plane modes (e.g. m = 0 shear modes,
    % w_bar = pi sqrt(G12/E2) a/h = 24.3) that the flexural Navier series excludes.
    wd = false(model.ndof,1); wd(3:5:end) = true;
    Mw = model.M; Mw(~wd,:) = 0; Mw(:,~wd) = 0;
    fw = sum(Phi .* (Mw*Phi), 1);           % phi' M phi = 1 -> fraction
    om = om(fw > 0.5);
    sc = a^2/h*sqrt(m.rho/m.E2);
    nchk = 6; if c == 2, nchk = 3; end
    for k = 1:nchk
        R = [R, vres('V04 Navier modal', sprintf('%s mode %d (m,n)=(%d,%d)', cases{c,2}, k, ...
             nav.mn(k,1), nav.mn(k,2)), om(k)*sc, nav.omega(k)*sc, 1.5e-2)]; %#ok<AGROW>
    end
end
end
