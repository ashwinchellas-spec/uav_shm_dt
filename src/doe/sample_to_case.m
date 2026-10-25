function [cfg, damages] = sample_to_case(spec, p, cls, cfg0)
%SAMPLE_TO_CASE  Physical DOE sample -> model configuration + damage list.
%   [cfg, damages] = sample_to_case(spec, p, cls, cfg0)
%   p   : 1 x k physical values (order of spec.vars)
%   cls : class name ('healthy', 'debond', ...)
%   cfg0: baseline configuration (wing_baseline_config)
v = cell2struct(num2cell(p(:)'), {spec.vars.name}, 2);
cfg = cfg0;
cfg.fe.bc = 'elastic_root';
cfg.fe.kroot = v.kroot;
damages = struct('type', {}, 'center', {}, 'radius', {}, 'severity', {}, ...
                 'profile', {}, 'face', {});
L = cfg.geom.span;
switch spec.tier
    case 'shm'
        cfg.unc = struct('E1', v.E1, 'E2', v.E2, 'G12', v.G12, 'Gc', v.Gc, 'Ec', v.Ec, ...
            'rho_f', v.rho_f, 'rho_c', v.rho_c, 'tply', v.tply, 'tcore', v.tcore, ...
            'dtheta', v.dtheta, 'T', v.T);
        c = cfg.geom.c_root;
        if ~strcmp(cls, 'healthy')
            face = 'bottom'; if v.dmg_face > 0.5, face = 'top'; end
            prof = 'uniform';
            if any(strcmp(cls, {'skin_fiber', 'core_crush'})), prof = 'gaussian'; end
            if strcmp(cls, 'root_loosening')
                % loosened band at the root: centre y, half-width 0.1c..0.5c
                hw = c*(0.1 + 0.4*(v.dmg_r - 0.015)/0.045);
                damages(1) = struct('type', cls, 'center', [0, v.dmg_y*c], 'radius', [0, hw], ...
                    'severity', v.dmg_D, 'profile', 'uniform', 'face', face);
            else
                damages(1) = struct('type', cls, 'center', [v.dmg_x*L, v.dmg_y*c], ...
                    'radius', [v.dmg_r, v.dmg_r*v.dmg_ar], 'severity', v.dmg_D, ...
                    'profile', prof, 'face', face);
            end
        end
    case 'design'
        cfg.lam.face_angles = [v.th1 v.th2 v.th3 v.th4];
        im = min(numel(spec.face_mats), floor(v.face_mat*numel(spec.face_mats)) + 1);
        cfg.lam.face_mat = spec.face_mats{im};
        cfg.lam.t_ply = v.tply_mm*1e-3;
        cfg.lam.t_core = v.tcore_mm*1e-3;
        % foam grade: G scales with density^2 (Gibson & Ashby, open/closed-cell foams)
        cfg.unc = struct('Gc', v.core_grade, 'Ec', v.core_grade, 'rho_c', sqrt(v.core_grade));
        cfg.geom.taper = v.taper;
        cfg.geom.sweep_le = v.sweep;
end
end
