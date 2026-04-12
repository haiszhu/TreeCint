function mol = gto(atom,basis)
% some functionality of pyscf gto
% 
%
% 02/10/25, Hai

% global ANG_OF NCTR_OF

BLKSIZE = 56;

% constants in Python code
CHARGE_OF   = 1; % indexing +1
PTR_COORD   = 2;
NUC_MOD_OF  = 3;
PTR_ZETA    = 4;
ATM_SLOTS   = 6; 
ATOM_OF     = 1; % indexing +1
ANG_OF      = 2; % ...
NPRIM_OF    = 3; % ...
NCTR_OF     = 4; % ...
% pointer to env
PTR_ENV_START   = 20;
%
NUC_POINT   = 1;

%
BOHR2ANG = 0.529177210903;

%
unit = 'Angstrom';

%
geom = atom;

% split and store
lines = splitlines(geom);
atoms_list = {};
for i = 1:numel(lines)
  line = strtrim(lines{i});
  if isempty(line)
    continue
  end
  parts = strsplit(line);
  if numel(parts) < 4
    continue  % 
  end
  %
  for k = 2:4
    if endsWith(parts{k}, '.')
      parts{k} = extractBefore(parts{k}, strlength(parts{k}));
    end
  end
  symb = parts{1};
  x = str2double(parts{2});
  y = str2double(parts{3});
  z = str2double(parts{4});
  % 
  atoms_list{i,1} = symb;
  atoms_list{i,2} = [x, y, z];
end

% format atom
nAtoms = size(atoms_list, 1);
scale = 1.0;
if strcmpi(unit,'Angstrom')
  scale = 1.0 / BOHR2ANG;
end
checkpts = [];
for i = 1:nAtoms
  symb = atoms_list{i,1};
  coords = atoms_list{i,2}; % in Angstrom
  coords_bohr = coords * scale;
  atoms_list_bohr{i,1} = symb;
  atoms_list_bohr{i,2} = coords_bohr;  % now in Bohr
  %
  checkpts = [checkpts;coords];
end
checkpts = checkpts';

% load(basisfile, symb, optimize...)
all_symbols = atoms_list_bohr(:,1);    
unique_elems = unique(all_symbols, 'stable');   

% create an empty containers.Map for the final basis_dict
basis_dict = containers.Map('KeyType','char', 'ValueType','any');
for ielem = 1:numel(unique_elems)
  symb = unique_elems{ielem};
  % load & store
  basis_dict(symb) = loadNWChemBasis(basis, symb);
end

%
pre_env = zeros(1,PTR_ENV_START);
% ELEMENTS_PROTON = containers.Map({'N','H'}, [7, 1]);
ELEMENTS_PROTON = buildElementsProtonMap();

%
atm = zeros(nAtoms, ATM_SLOTS);
env = [pre_env];
ptr_env = numel(pre_env);
for ia = 1:nAtoms
  % 
  symb   = atoms_list_bohr{ia, 1};
  coords = atoms_list_bohr{ia, 2}; 

  % 
  nuc_charge = ELEMENTS_PROTON(symb);  

  % 
  zeta = 0.0;
  nuclear_model = NUC_POINT;

  % 
  env0 = [coords, zeta];
  len_env0 = numel(env0);  

  % 
  atm0 = zeros(1, ATM_SLOTS);
  atm0(CHARGE_OF)  = nuc_charge;
  atm0(PTR_COORD)  = ptr_env;         
  atm0(NUC_MOD_OF) = nuclear_model;   
  atm0(PTR_ZETA)   = ptr_env + 3;     

  % 
  atm(ia,:) = atm0;
  env = [env, env0];  
  ptr_env = ptr_env + len_env0; 
end

% initialize basis dictionary
if isa(basis_dict, 'struct')
  allSymbols = fieldnames(basis_dict);
  isMapInput = false;
else
  % if basis_dict is a containers.Map
  % allSymbols = keys(basis_dict);
  allSymbols = unique_elems; % respect key order in specified geo...
  isMapInput = true;
end
% 
basdic = containers.Map();
for iSym = 1:numel(allSymbols)
  symb = allSymbols{iSym};

  if isMapInput
    basis_add = basis_dict(symb);
  else
    basis_add = basis_dict.(symb);
  end

  ptr_env0 = ptr_env;
  atom_id = 0;
  bas0 = [];   
  env0 = [];
  for ib = 1:numel(basis_add)
    b = basis_add{ib};  
    angl = b{1};  
    kappa = 0;
    b_coefftmp = b(2:end);  
    b_coeff = cell2mat(b_coefftmp(:));  
    [~, idxSort] = sort(b_coeff(:,1), 'descend');
    b_coeff = b_coeff(idxSort,:);
    es = b_coeff(:,1);
    cs = b_coeff(:,2:end);  
    [nprim, nctr] = size(cs);

    n1 = double((angl*2 + 2) + 1) * 0.5;
    gaussian_intangles = gamma(n1) ./ (2 .* ((2.0 .* es).^n1));
    gto_normangles = 1.0 ./ sqrt(gaussian_intangles);
    for ip = 1:nprim
      cs(ip,:) = cs(ip,:) * gto_normangles(ip);
    end

    ee = es(:) + es(:).';
    n1 = ((angl*2 + 2) + 1) * 0.5;
    ee = gamma(n1) ./ (2 .* (ee.^n1));
    % double check
    scales = zeros(1, nctr);
    for ic = 1:nctr
      Cp = cs(:, ic);
      val = 0;
      for p = 1:nprim
        for q = 1:nprim
          val = val + Cp(p)* ee(p,q)* Cp(q);
        end
      end
      scales(ic) = 1.0 / sqrt(val);
    end
    for ic = 1:nctr
      cs(:,ic) = cs(:,ic) * scales(ic);
    end

    env0 = [env0, [ reshape(es,1,[]), reshape(cs,1,[]) ]]; % env0.append(es) then env0.append(cs.T.reshape(-1))
                                                           % cs is column wise in python... 
    ptr_exp = ptr_env0;                
    ptr_coeff = ptr_exp + nprim;         
    ptr_env0 = ptr_coeff + nprim*nctr;  
    bas0 = [bas0; [atom_id, angl, nprim, nctr, kappa, ptr_exp, ptr_coeff, 0]]; 
  end

  env0 = env0;
  bas0 = bas0;

  ptr_env = ptr_env + numel(env0);
  basdic(symb) = bas0;
  env = [env, env0];
end

%
% [bas] = gatherBasForAtoms(basdic, atoms_list_bohr);
basBlocks = cell(nAtoms, 1);
for ia=1:nAtoms
  atomi = atoms_list_bohr(ia,:);
  symb = atomi{1};  % e.g. 'N' or 'H'
  b = basdic(symb);
  b(:, ATOM_OF) = ia-1;
  basBlocks{ia} = b;
end
bas = vertcat(basBlocks{:});

%
ao_loc = make_loc(bas, 'sph');
natm = size(atm,1);
nbas = size(bas,1);
shls_slice = [0, nbas];

% num of grid points to be evaluated..., for simplicity, start with 1
ngrids = 1;
non0tab = ones( floor((ngrids + BLKSIZE - 1) / BLKSIZE), nbas);

% total number of contracted GTOs for the given mol object
nao_nr = sum((bas(:,ANG_OF) * 2 + 1) .* bas(:,NCTR_OF));

%
mol = [];
mol.atm = atm;
mol.bas = bas;
mol.env = env;
%
mol.ao_loc = ao_loc;
mol.natm = natm;
mol.nbas = nbas;
mol.shls_slice = shls_slice;
%
mol.ngrids = 1;
mol.non0tab = non0tab;
%
mol.nao_nr = nao_nr;

%% 1st tricky part... eval_gto... 02/10/25 Hai here
% eval_name = 'GTOval_sph';
% [x, y, z] = ndgrid(linspace(0,1,5), linspace(0,1,5), linspace(0,1,5));
% coords = [x(:) y(:) z(:)];
% coords = cat(4,x,y,z);
% myeval_gto(eval_name, coords, shls_slice, ao_loc, non0tab, atm, natm, bas, nbas, env);
mol.eval_gto = @(eval_name, coords) myeval_gto( eval_name, coords, shls_slice, ao_loc, non0tab, atm, natm, bas, nbas, env);

%% add checkpts for building adaptive grid
mol.checkpts = checkpts;

%% 2nd, with Norb^2 basis
mol.eval_gto2 = @(eval_name, coords) myeval_gto2( eval_name, coords, shls_slice, ao_loc, non0tab, atm, natm, bas, nbas, env);

%% 3rd, with just i-th basis
mol.eval_gtoi = @(eval_name, coords, i) myeval_gtoi( eval_name, coords, i, shls_slice, ao_loc, non0tab, atm, natm, bas, nbas, env);
end

function fi = myeval_gto(eval_name, coords, shls_slice, ao_loc, non0tab, atm, natm, bas, nbas, env)
% global ANG_OF NCTR_OF 
persistent ANG_OF NCTR_OF nd;
if isempty(nd)
  ANG_OF      = 2; % ...
  NCTR_OF     = 4; % ...
  nd = sum((bas(:,ANG_OF) * 2 + 1) .* bas(:,NCTR_OF));
end
if contains(eval_name, 'GTOval_sph')
  dims = size(coords);
  ndims_var = ndims(coords); % only support two cases
  if ndims_var == 4
    % assume something like n1 x n2 x n3 x 3, and only this
    n1 = dims(1);
    n2 = dims(2);
    n3 = dims(3);
    ngrids = n1*n2*n3;
    x = coords(:,:,:,1);
    y = coords(:,:,:,2);
    z = coords(:,:,:,3);
    xyz = [x(:), y(:), z(:)]';
    xyz = xyz(:);
    fi = zeros([nd ngrids]); 
    atm = reshape(atm',1,[]); % row major in C...
    bas = reshape(bas',1,[]);
    fi = GTOval_sph_generic_mwrap_mex(ngrids, shls_slice, ao_loc, fi, xyz, non0tab, atm, natm, bas, nbas, env);
    fi = reshape(fi,[nd ngrids])';
    fi = reshape(fi,[n1 n2 n3 nd]);
    % keyboard
  elseif ndims_var == 2
    n1 = dims(1);
    n2 = dims(2);
    if n1==3
      xyz = coords;
      ngrids = n2;
      xyz = xyz(:);
      %
      fi = zeros([nd ngrids]); 
      atm = reshape(atm',1,[]); % row major in C...
      bas = reshape(bas',1,[]);
      fi = GTOval_sph_generic_mwrap_mex(ngrids, shls_slice, ao_loc, fi, xyz, non0tab, atm, natm, bas, nbas, env);
      %
      if (n1==3) && (n2==3)
        disp('I am assuming your input coords is xyz by num of grid pts.'); % what can I do...
      end
      fi = reshape(fi,[nd ngrids]);
      %
    end
    if (n1~=3) && (n2==3)
      x = coords(:,1);
      y = coords(:,2);
      z = coords(:,3);
      xyz = [x(:), y(:), z(:)]';
      ngrids = n1;
      xyz = xyz(:);
      %
      fi = zeros([nd ngrids]); 
      atm = reshape(atm',1,[]); % row major in C...
      bas = reshape(bas',1,[]);
      fi = GTOval_sph_generic_mwrap_mex(ngrids, shls_slice, ao_loc, fi, xyz, non0tab, atm, natm, bas, nbas, env);
      %
      fi = reshape(fi,[nd ngrids])';
      %
    end
    
  else
    xyz = [];
    disp('Variable does not match expected tensor dimensions.');
    fi = [];
    % keyboard
  end  
  
end
end

function fi = myeval_gto2(eval_name, coords, shls_slice, ao_loc, non0tab, atm, natm, bas, nbas, env)
% global ANG_OF NCTR_OF 
persistent ANG_OF NCTR_OF Norb nd triuidx triflag;
if isempty(triuidx)
  ANG_OF      = 2; % ...
  NCTR_OF     = 4; % ...
  Norb = sum((bas(:,ANG_OF) * 2 + 1) .* bas(:,NCTR_OF));
  nd = Norb*(Norb+1)/2;
  triuidx = 1:Norb^2;
  triflag = true(Norb,Norb);
  triflag = triu(triflag);
  triuidx = triuidx(triflag(:));
end
if contains(eval_name, 'GTOval_sph')
  dims = size(coords);
  ndims_var = ndims(coords); % only support two cases
  if ndims_var == 4
    % assume something like n1 x n2 x n3 x 3, and only this
    n1 = dims(1);
    n2 = dims(2);
    n3 = dims(3);
    ngrids = n1*n2*n3;
    x = coords(:,:,:,1);
    y = coords(:,:,:,2);
    z = coords(:,:,:,3);
    xyz = [x(:), y(:), z(:)]';
    xyz = xyz(:);
    fi = zeros([Norb ngrids]); 
    atm = reshape(atm',1,[]); % row major in C...
    bas = reshape(bas',1,[]);
    fi = GTOval_sph_generic_mwrap_mex(ngrids, shls_slice, ao_loc, fi, xyz, non0tab, atm, natm, bas, nbas, env);
    fi = reshape(fi,[Norb ngrids]);
    % these will dominate due to ^2 scaling...
    firs = reshape(fi, [Norb, 1, ngrids]);
    fij_full = pagemtimes(firs, permute(firs, [2, 1, 3])); 
    ftmp = reshape(fij_full,[Norb^2 ngrids]);
    ftmp = ftmp(triuidx,:)';
    fi = reshape(ftmp,[n1 n2 n3 nd]);
  
  else % nothing else...  
    xyz = [];
    disp('Variable does not match expected tensor dimensions.');
    fi = [];
    % keyboard
  end  
  
end
end

function fi = myeval_gtoi(eval_name, coords, i, shls_slice, ao_loc, non0tab, atm, natm, bas, nbas, env)
% global ANG_OF NCTR_OF 
persistent ANG_OF NCTR_OF ndi;
if isempty(ndi)
  ANG_OF      = 2; % ...
  NCTR_OF     = 4; % ...
  ndi = 1;
end
if contains(eval_name, 'GTOval_sph')
  dims = size(coords);
  ndims_var = ndims(coords); % only support two cases
  if ndims_var == 4
    % assume something like n1 x n2 x n3 x 3, and only this
    n1 = dims(1);
    n2 = dims(2);
    n3 = dims(3);
    ngrids = n1*n2*n3;
    x = coords(:,:,:,1);
    y = coords(:,:,:,2);
    z = coords(:,:,:,3);
    xyz = [x(:), y(:), z(:)]';
    xyz = xyz(:);
    fi = zeros([ndi ngrids]); 
    atm = reshape(atm',1,[]); % row major in C...
    bas = reshape(bas',1,[]);
    fi = GTOival_sph_generic_mwrap_mex(i, ngrids, shls_slice, ao_loc, fi, xyz, non0tab, atm, natm, bas, nbas, env);
    fi = reshape(fi,[ndi ngrids])';
    fi = reshape(fi,[n1 n2 n3 ndi]);
    % keyboard
  elseif ndims_var == 2
    n1 = dims(1);
    n2 = dims(2);
    if n1==3
      xyz = coords;
      ngrids = n2;
      xyz = xyz(:);
      %
      fi = zeros([ndi ngrids]); 
      atm = reshape(atm',1,[]); % row major in C...
      bas = reshape(bas',1,[]);
      fi = GTOival_sph_generic_mwrap_mex(i, ngrids, shls_slice, ao_loc, fi, xyz, non0tab, atm, natm, bas, nbas, env);
      %
      if (n1==3) && (n2==3)
        disp('I am assuming your input coords is xyz by num of grid pts.'); % what can I do...
      end
      fi = reshape(fi,[ndi ngrids]);
      %
    end
    if (n1~=3) && (n2==3)
      x = coords(:,1);
      y = coords(:,2);
      z = coords(:,3);
      xyz = [x(:), y(:), z(:)]';
      ngrids = n1;
      xyz = xyz(:);
      %
      fi = zeros([ndi ngrids]); 
      atm = reshape(atm',1,[]); % row major in C...
      bas = reshape(bas',1,[]);
      fi = GTOival_sph_generic_mwrap_mex(i, ngrids, shls_slice, ao_loc, fi, xyz, non0tab, atm, natm, bas, nbas, env);
      %
      fi = reshape(fi,[ndi ngrids])';
      %
    end
    
  else
    xyz = [];
    disp('Variable does not match expected tensor dimensions.');
    fi = [];
    % keyboard
  end  
  
end
end