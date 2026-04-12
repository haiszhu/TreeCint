function test_ccpvdz
% H2O / cc-pVDZ full ERI pipeline using the standalone int2-bdmk-scf executable.
%
% Differences from h2o_ccpvdz_vmunu_vijkl.m:
%   - BDMK runs via int2-bdmk-scf (Fortran, OMP-parallel, no MEX call)
%   - Vijkl assembled with two matrix multiplies instead of a Norb^2 loop
%   - OMP thread count set via omp_num_threads parameter
%   - Paths derived from the script location (no hardcoded user paths)

% ---- paths ---------------------------------------------------------------
tree_root  = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(fullfile(tree_root, 'matlab'));
addpath(fullfile(tree_root, 'external', 'treefun'));

int2bdmk_exe = fullfile(tree_root, 'bin', 'int2-bdmk-scf');
if ~isfile(int2bdmk_exe)
  error('int2-bdmk-scf not found at %s.\nRun: make int2-bdmk-scf', int2bdmk_exe);
end

% ---- parameters ----------------------------------------------------------
omp_num_threads = 4;      % OMP threads for int2-bdmk-scf

treefun_order = 8;
treefun_eps   = 1e-3;
isdf_eps      = 1e-3;
nd            = 290;

% treefun_order = 8;
% treefun_eps   = 1e-6;
% isdf_eps      = 1e-6;
% nd            = 290;

% ---- molecule ------------------------------------------------------------
rad    = 15;
geom   = sprintf(['O    0    0.       0.\n', ...
                  'H    0    -0.757   0.587\n', ...
                  'H    0    0.757    0.587\n']);
molname = 'h2o';
basmod  = 'cc-pvdz.dat';
basis   = fullfile(tree_root, 'basis', basmod);
mol     = gto(geom, basis);
eval_name = 'GTOval_sph';

% ---- Stage 1: adaptive octree -------------------------------------------
disp('=========Start treefun=======');
opts = struct('balance',   true, ...
              'tol',       treefun_eps, ...
              'checkpts',  mol.checkpts, ...
              'ifcoeffs',  false);
func = @(x,y,z) mol.eval_gto(eval_name, cat(4, x, y, z));
tic;
f = treefun3(func, [-rad rad -rad rad -rad rad], treefun_order, opts);
t_treefun = toc;
disp("    treefun time    : " + t_treefun + " s");
disp("    treefun order   : " + treefun_order);
disp("    treefun boxes   : " + size(f.domain,2));
disp('=========End treefun=======');

% ---- Stage 2: upsample to BDMK quadrature grid --------------------------
disp('=========Start upsample adaptive grid=======');
Norb  = mol.nao_nr;
ndim  = 3;
ratio = 0.5 / rad;
ipoly = 0;
f2    = f;
f2.n  = floor(1.5 * f.n);
tic;
[src, nleafbox, srcleaf, wtsleaf, norder, npbox, nboxes, nlevels, ...
 ltree, itree, iptr, centers, boxsize] = treefun2bdmk(f2, ndim, ratio, ipoly);
DS1 = src / ratio;
t_grid = toc;
disp("    upsample time   : " + t_grid + " s");
disp('=========End upsample adaptive grid=======');

% ---- Stage 3: evaluate orbital products at quadrature points -------------
disp('=========Start Norb^2 basis eval=======');
tic;
src0  = src / ratio;
npts  = numel(src0(1,:));
fvals0 = squeeze(mol.eval_gto(eval_name, ...
           cat(4, squeeze(src0(1,:,:)), squeeze(src0(2,:,:)), squeeze(src0(3,:,:)))));
fvals0 = reshape(fvals0, [npts, Norb]);
fvals_ij = zeros(npts, Norb*(Norb+1)/2);
tmpidx = 0;
for i = 1:Norb
  for j = 1:i
    tmpidx = tmpidx + 1;
    fvals_ij(:, tmpidx) = fvals0(:,i) .* fvals0(:,j);
  end
end
t_eval = toc;
disp("    basis eval time : " + t_eval + " s");
disp('=========End Norb^2 basis eval=======');

% ---- Stage 4: interpolative decomposition (ISDF) -------------------------
disp('=========Start interpolative decomposition=======');
tic;
A   = fvals_ij;
[SK, RD, T] = id_libid(A, nd);
Ask = A(:, SK);
idcoefs = zeros(nd, size(A,2));
idcoefs(:, SK) = eye(nd);
idcoefs(:, RD) = T;

% Write ISDF vectors to HDF5 (read by int2-bdmk-scf)
eps_string       = regexprep(sprintf('%.0e', isdf_eps), 'e-0*', 'e-');
isdf_filename    = fullfile(pwd, sprintf('myisdf_%s.h5', eps_string));
dims_rd          = [2, npts, nd];
interp_vecs      = zeros(dims_rd);
interp_vecs(1,:,:) = Ask;
file_id = H5F.create(isdf_filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
ipv_dspace = H5S.create_simple(3, fliplr(size(interp_vecs)), []);
ipv_dcpl   = H5P.create('H5P_DATASET_CREATE');
ipv_dset   = H5D.create(file_id, '/interpolating_vectors', 'H5T_NATIVE_DOUBLE', ipv_dspace, ipv_dcpl);
H5D.write(ipv_dset, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', interp_vecs);
H5D.close(ipv_dset); H5F.close(file_id);
t_id = toc;
disp("    ID time         : " + t_id + " s  (nd=" + nd + ")");
disp('=========End interpolative decomposition=======');

% ---- Stage 5: write tree metadata to HDF5 --------------------------------
disp('=========Start write tree HDF5=======');
ikernel = 1;
beta    = 6.0;
treefun_filename = fullfile(pwd, sprintf('treefun_%s_%s.h5', molname, erase(basmod,'.dat')));
write_treefun_h5(molname, basmod, ndim, treefun_eps, ikernel, beta, ...
                 ipoly, norder, npbox, nboxes, nlevels, ltree, itree, ...
                 iptr, centers, boxsize, DS1, nleafbox, wtsleaf, ratio);
disp('=========End write tree HDF5=======');

% ---- Stage 6: BDMK via int2-bdmk-scf ------------------------------------
disp('=========Start int2-bdmk-scf=======');
output_filename = fullfile(pwd, sprintf('bdmk_%s.h5', eps_string));
if isfile(output_filename); delete(output_filename); end

cmd = sprintf('export OMP_NUM_THREADS=%d; "%s" "%s" "%s" "%s"', ...
              omp_num_threads, int2bdmk_exe, ...
              treefun_filename, isdf_filename, output_filename);
tic;
[status, out] = system(cmd);
t_bdmk = toc;
if status ~= 0
  error('int2-bdmk-scf failed (status=%d):\n%s', status, out);
end
fprintf('%s\n', out);
disp("    BDMK time       : " + t_bdmk + " s  (" + omp_num_threads + " threads)");
disp('=========End int2-bdmk-scf=======');

% ---- Stage 7: read Vmunu from HDF5 ---------------------------------------
Vmunu = h5read(output_filename, '/Vmunu');

% ---- Stage 8: Vijkl via two matrix multiplies (no Norb^2 loop) ----------
% Original computeVijkl loops over (i,k) calling Norb^2 small DGEMMs.
% Here we reshape to a single (nd x Norb^2) collocation matrix C and compute:
%   VijklM = C' * Vmunu * C   [Norb^2 x Norb^2]
% then reshape + permute to [Norb, Norb, Norb, Norb].
%
% Index convention: C(mu, j+(i-1)*Norb) = idcoefs(mu, tri_idx(i,j)).
% After reshape(VijklM, [Norb,Norb,Norb,Norb]) MATLAB gives T(j,i,l,k),
% so permute([2,1,4,3]) recovers Vijkl(i,j,k,l).
disp('=========Start Vijkl computation=======');
tic;

% Build symmetrized collocation matrix [nd x Norb x Norb]
CM3 = zeros(nd, Norb, Norb);
tmpidx = 0;
for i = 1:Norb
  for j = 1:i
    tmpidx = tmpidx + 1;
    CM3(:,j,i) = idcoefs(:, tmpidx);
  end
end
for i = 1:Norb
  for j = i+1:Norb
    CM3(:,j,i) = CM3(:,i,j);
  end
end
CM = reshape(CM3, [nd, Norb^2]);   % [nd x Norb^2]

% Two DGEMM calls: Vmunu*CM then CM'*(result)
VijklM  = CM' * (Vmunu * CM);     % [Norb^2 x Norb^2]
Vijkl   = permute(reshape(VijklM, [Norb, Norb, Norb, Norb]), [2, 1, 4, 3]);

t_vijkl = toc;
disp("    Vijkl time      : " + t_vijkl + " s");
disp('=========End Vijkl computation=======');

% ---- Save results --------------------------------------------------------
eri_mat = ['ERI_' molname '_' erase(erase(basmod,'.dat'),'-') '_' eps_string '.mat'];
save(eri_mat, 'Vijkl');

eri_h5 = ['ERI_' molname '_' erase(erase(basmod,'.dat'),'-') '_' eps_string '.h5'];
file_id        = H5F.create(eri_h5, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
vijkl_dspace   = H5S.create_simple(4, fliplr([Norb,Norb,Norb,Norb]), []);
vijkl_dcpl     = H5P.create('H5P_DATASET_CREATE');
vijkl_dset     = H5D.create(file_id, '/DS1', 'H5T_NATIVE_DOUBLE', vijkl_dspace, vijkl_dcpl);
H5D.write(vijkl_dset, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', Vijkl);
H5D.close(vijkl_dset); H5S.close(vijkl_dspace); H5F.close(file_id);

% ---- Stage 9: Python / PySCF verification --------------------------------
[py_ok, py_metrics] = verify_eri_with_python(eri_h5, geom, tree_root);
assert(py_ok, 'Python ERI verification failed.');
disp('=========Python ERI verification=======');
fprintf('    max_abs = %.16e\n', py_metrics.max_abs);
fprintf('    rel_fro = %.16e\n', py_metrics.rel_fro);
fprintf('    rmse    = %.16e\n', py_metrics.rmse);
disp('=========End Python ERI verification=======');

disp('Pipeline complete.');
end

% --------------------------------------------------------------------------
function [ok, metrics] = verify_eri_with_python(eri_h5_filename, geom, tree_root)
ok      = false;
metrics = struct('max_abs', NaN, 'rel_fro', NaN, 'rmse', NaN);

% Prefer pyenv shim; fall back to system python3
candidates = {'/Users/hzhu/.pyenv/shims/python', 'python3', 'python'};
target_py = '';
for c = candidates
  [s, p] = system(sprintf('command -v %s 2>/dev/null', c{1}));
  if s == 0 && ~isempty(strtrim(p))
    target_py = strtrim(p); break;
  end
end
if isempty(target_py)
  error('No Python executable found. Install Python with numpy, h5py, pyscf, scipy.');
end
fprintf('Python: %s\n', target_py);

verify_script = fullfile(tree_root, 'test', 'isdf', 'H2O', 'verify_eri_h2o_ref.py');
if ~isfile(verify_script)
  error('Verification script missing: %s', verify_script);
end

tmp_mat      = [tempname '.mat'];
cleanup_obj  = onCleanup(@() safe_delete(tmp_mat));
cmd = sprintf('"%s" "%s" "%s" "%s"', target_py, verify_script, eri_h5_filename, tmp_mat);
[status, out] = system(cmd);
if status ~= 0
  error(['Python verification failed:\n  %s\n' ...
         'Install packages: %s -m pip install numpy h5py pyscf scipy\n' ...
         'Output:\n%s'], target_py, target_py, out);
end
result          = load(tmp_mat);
metrics.max_abs = result.max_abs(1);
metrics.rel_fro = result.rel_fro(1);
metrics.rmse    = result.rmse(1);
ok = true;
end

function safe_delete(path)
if exist(path, 'file'); delete(path); end
end
