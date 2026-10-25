function M = run_morris(r, outfile)
%RUN_MORRIS  Morris elementary-effects screening of the healthy wing.
%
%   M = run_morris(r, outfile)
%   Factors : the 12 manufacturing / environmental / boundary variables of
%             the 'shm' tier (ranges from doe_config) - no damage.
%   Outputs : f1, f2 (1T), f3 (2B), f5, f6 [Hz], tip deflection & twist
%             under the 1 g load, minimum reserve factor RF_min.
%   Cost    : r (k + 1) = 13 r simulations.
%   Use     : ranks which variabilities matter, and gives the size of the
%             "environmental noise" that a damage signal must exceed.
spec = doe_config('shm');
k = 12;                                    % E1 ... kroot
cfg0 = wing_baseline_config();
ctx = doe_context(spec, cfg0);
[Uk, delta] = morris_design(k, r, 4, 11);
U = [Uk, 0.5*ones(size(Uk,1), spec.k - k)];
P = doe_map(spec, U);
names = {'f1', 'f2', 'f3', 'f5', 'f6', 'w_tip', 'twist', 'RF_min'};
Y = nan(size(P,1), numel(names));
for i = 1:size(P,1)
    [cfg, dm] = sample_to_case(spec, P(i,:), 'healthy', cfg0);
    o = run_sample(cfg, dm, ctx);
    Y(i,:) = [o.freq([1 2 3 5 6]), o.tipw(2,1), o.twist(1), o.RF(5)];
end
S = morris_indices(Uk, Y, delta, k);
M.S = S; M.U = Uk; M.P = P(:,1:k); M.Y = Y; M.delta = delta;
M.factors = {spec.vars(1:k).name}; M.outputs = names; M.r = r;
if nargin > 1, save(outfile, '-struct', 'M', '-v7'); end
end
