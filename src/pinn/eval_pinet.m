function M = eval_pinet(out, D, F, db, rows, thr_rows)
%EVAL_PINET  Test metrics of a trained damage network.
%   M = eval_pinet(out, D, F, db, rows, thr_rows)
%   rows     : evaluation rows (test set)
%   thr_rows : rows used to set the detection threshold (healthy
%              VALIDATION rows -> 1 % false alarms)
%   M.acc, M.macroF1, M.conf (7x7, rows = true), M.auc (1x7), M.dr (1x7),
%   M.fa, M.loc_med / M.loc_hit50 (classes 2-5, all and significant),
%   M.ext_r2, M.ext_rho (Spearman), M.rf_rmse, M.rf_r2, M.phys (mean physics residual).
y = F.y_class; nc = numel(F.class_names);
[~, yp] = max(out.pc, [], 2);
r = rows;
M.conf = zeros(nc);
for i = find(r)', M.conf(y(i), yp(i)) = M.conf(y(i), yp(i)) + 1; end
M.acc = sum(diag(M.conf)) / sum(M.conf(:));
prec = diag(M.conf)' ./ max(sum(M.conf, 1), 1); rec = diag(M.conf)' ./ max(sum(M.conf, 2)', 1);
M.f1 = 2*prec.*rec ./ max(prec + rec, eps); M.macroF1 = mean(M.f1);
score = 1 - out.pc(:,1);
thr = prctile_simple(score(thr_rows & y == 1), 99);
M.fa = mean(score(r & y == 1) > thr);
M.auc = nan(1, nc); M.dr = nan(1, nc);
for c = 2:nc
    q = r & (y == c | y == 1);
    M.auc(c) = roc_auc(score(q), y(q) == c);
    M.dr(c) = mean(score(r & y == c) > thr);
end
% localisation (classes 2-5: localisable physics), all and significant
sid = F.sample_id;
sig = db.P(sid,17) > 0.5 & db.P(sid,15) > 0.03;
e = sqrt(sum((out.loc - F.y_loc).^2, 2));
M.loc_med = nan(2, nc); M.loc_hit50 = nan(2, nc);
for c = 2:5
    q = r & y == c;              M.loc_med(1,c) = median(e(q));  M.loc_hit50(1,c) = mean(e(q) < 0.05);
    q = r & y == c & sig;        M.loc_med(2,c) = median(e(q));  M.loc_hit50(2,c) = mean(e(q) < 0.05);
end
% extent (damage amount) and residual strength
q = r & ismember(y, 2:6);
M.ext_r2 = 1 - sum((out.s(q) - D.T.ext(q)).^2) / sum((D.T.ext(q) - mean(D.T.ext(q))).^2);
M.ext_rho = spearman_rho(out.s(q), D.T.ext(q));
rfp = out.rf*D.rf_sd + D.rf_mu; rft = F.y_RF;
M.rf_rmse = sqrt(mean((rfp(r) - rft(r)).^2));
M.rf_r2 = 1 - sum((rfp(r) - rft(r)).^2) / sum((rft(r) - mean(rft(r))).^2);
qd = r & y >= 2 & y <= 6 & rft < prctile_simple(rft(r & y == 1), 1);   % clearly weakened wings
M.rf_rmse_weak = sqrt(mean((rfp(qd) - rft(qd)).^2)); M.n_weak = sum(qd);
M.phys = mean(mean((out.xh(r,:) - D.T.dXn(r,:)).^2, 2));
end
