function model = assemble_system(mesh, sections, opts)
%ASSEMBLE_SYSTEM  Global sparse stiffness and mass matrices.
%
%   model = assemble_system(mesh, sections, opts)
%   mesh      : from wing_mesh / rect_mesh (mesh.sec_id gives each
%               element's section index into 'sections')
%   sections  : struct array (or cell) of laminate_section outputs
%   opts.shear: 'hr' (default) | 'sri' | 'full'
%
%   model.K, model.M  : ndof x ndof sparse, ndof = 5*nn
%   model.rec         : per-element recovery operators (Bm, Bb, Rq: 3/3/2 x 20 x ne)
%   model.edofs       : ne x 20 global DOF numbers of each element
%
%   Elements with identical shape and section reuse a cached element
%   matrix, so uniform meshes assemble very fast (important for the
%   thousands of DOE runs in Milestone 3).

if nargin < 3, opts = struct(); end
if ~isfield(opts, 'shear'), opts.shear = 'hr'; end
if iscell(sections), sections = [sections{:}]; end

ne = mesh.ne; nn = mesh.nn; ndof = 5*nn;
edofs = zeros(ne, 20);
for a = 1:4
    edofs(:, 5*a-4:5*a) = 5*(mesh.elems(:,a) - 1) + (1:5);
end

II = zeros(400*ne, 1); JJ = II; VK = II; VM = II;
Bm = zeros(3,20,ne); Bb = Bm; Rq = zeros(2,20,ne); area = zeros(ne,1);
cache = containers.Map('KeyType','char','ValueType','any');

for e = 1:ne
    xy  = mesh.nodes(mesh.elems(e,:), :);
    sid = mesh.sec_id(e);
    rel = xy - xy(1,:);
    key = sprintf('%d|%.9g,', sid, round(rel(:)'*1e9)/1e9);
    if isKey(cache, key)
        c = cache(key); Ke = c.Ke; Me = c.Me; r = c.r;
    else
        [Ke, Me, r] = elem_fsdt_q4(xy, sections(sid), opts.shear);
        c.Ke = Ke; c.Me = Me; c.r = r; cache(key) = c;
    end
    idx = (e-1)*400 + (1:400);
    d = edofs(e,:);
    [jj, ii] = meshgrid(d, d);
    II(idx) = ii(:); JJ(idx) = jj(:);
    VK(idx) = Ke(:); VM(idx) = Me(:);
    Bm(:,:,e) = r.Bm; Bb(:,:,e) = r.Bb; Rq(:,:,e) = r.Rq; area(e) = r.area;
end

model.K = sparse(II, JJ, VK, ndof, ndof);
model.M = sparse(II, JJ, VM, ndof, ndof);
model.K = (model.K + model.K') / 2;
model.M = (model.M + model.M') / 2;
model.edofs = edofs;
model.rec.Bm = Bm; model.rec.Bb = Bb; model.rec.Rq = Rq; model.rec.area = area;
model.mesh = mesh; model.sections = sections; model.opts = opts;
model.ndof = ndof;
model.n_unique_elems = cache.Count;
end
