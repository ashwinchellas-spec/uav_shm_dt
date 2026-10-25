function env = flight_envelope(ac, g)
%FLIGHT_ENVELOPE  Manoeuvre + gust V-n diagram and the design load factors.
%
%   env = flight_envelope(cfg.ac, cfg.g)
%   Stall boundary   n_s(V) = rho V^2 S CLmax / (2 W)
%   Manoeuvre corner V_A = V_S1 sqrt(n_pos)
%   Gust lines       n = 1 +/- dn(V),  Pratt-Walker (gust_load_factor.m)
%   Gust loads are capped by the stall boundary (the wing cannot produce
%   more lift than CLmax at that speed).
%
%   env.n_limit_pos / env.n_limit_neg : design limit load factors
%   env.critical : text description of the governing condition
if nargin < 2, g = 9.81; end
W = ac.mass * g;
q = @(V) 0.5 * ac.rho * V.^2;
env.Vs1 = sqrt(2*W / (ac.rho * ac.S * ac.CLmax));
env.Va  = env.Vs1 * sqrt(ac.n_pos);
nstall_p = @(V) q(V) * ac.S * ac.CLmax / W;
nstall_n = @(V) q(V) * ac.S * ac.CLmin / W;
[dn_c, gi] = gust_load_factor(setfield(ac, 'V', ac.Vc), ac.Ude_c); %#ok<SFLD>
[dn_d, ~ ] = gust_load_factor(setfield(ac, 'V', ac.Vd), ac.Ude_d); %#ok<SFLD>
env.gust = gi;
env.dn_c = dn_c; env.dn_d = dn_d;
% candidate points (capped by stall)
pts = struct('name', {}, 'V', {}, 'n', {});
pts(end+1) = struct('name', 'Manoeuvre +n (V_A..V_D)', 'V', env.Va, 'n', ac.n_pos);
pts(end+1) = struct('name', 'Manoeuvre -n',            'V', ac.Vc,  'n', ac.n_neg);
pts(end+1) = struct('name', 'Gust +Ude at V_C', 'V', ac.Vc, 'n', min(1 + dn_c, nstall_p(ac.Vc)));
pts(end+1) = struct('name', 'Gust +Ude at V_D', 'V', ac.Vd, 'n', min(1 + dn_d, nstall_p(ac.Vd)));
pts(end+1) = struct('name', 'Gust -Ude at V_C', 'V', ac.Vc, 'n', max(1 - dn_c, nstall_n(ac.Vc)));
pts(end+1) = struct('name', 'Gust -Ude at V_D', 'V', ac.Vd, 'n', max(1 - dn_d, nstall_n(ac.Vd)));
env.points = pts;
nv = [pts.n];
[env.n_limit_pos, ip] = max(nv);
[env.n_limit_neg, in] = min(nv);
env.critical = pts(ip).name; env.critical_V = pts(ip).V;
env.critical_neg = pts(in).name;
env.nstall_p = nstall_p; env.nstall_n = nstall_n;
end
