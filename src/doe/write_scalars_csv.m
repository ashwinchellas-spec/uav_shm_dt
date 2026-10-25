function write_scalars_csv(db, file)
%WRITE_SCALARS_CSV  One row per sample: inputs, class, scalar outputs.
%   Shared by merge_database and repair_database.
N = numel(db.cls); nm = size(db.freq, 2);
fid = fopen(file, 'w');
hdr = [{'id', 'class'}, db.var_names, arrayfun(@(k) sprintf('f%d_Hz', k), 1:nm, 'UniformOutput', false), ...
       {'tipw_L1_mid_m', 'twist_L1_deg', 'tipw_L2_mid_m', 'RF_TsaiWu', 'RF_Hashin', 'RF_core', ...
        'RF_ILSS', 'RF_min', 'DAF_gust', 'DAF_sharp', 'mass_kg', 'ok'}];
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:N
    fprintf(fid, '%d,%s,', i, db.classes{db.cls(i)});
    fprintf(fid, '%.8g,', db.P(i,:), db.freq(i,:), db.tipw(i,2,1), db.twist(i,1), db.tipw(i,2,2), ...
            db.RF(i,1:5), db.gust(i,1), db.gust(i,5), db.mass(i));
    fprintf(fid, '%d\n', db.ok(i));
end
fclose(fid);
end
