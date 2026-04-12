function new_basis = loadNWChemBasis(basisfile, symb)
% load basis
%
% for example, basis = loadNWChemBasis('cc-pvdz.dat', 'H');

fid = fopen(basisfile, 'r');
if fid < 0
  error('Could not open file: %s', basisfile);
end

fileLines = {};
lineIdx = 0;
tline = fgetl(fid);
while ischar(tline)
  lineIdx = lineIdx + 1;
  fileLines{lineIdx,1} = tline;
  tline = fgetl(fid);
end
fclose(fid);

raw_basis_lines = {};
collecting = false;

for i = 1:numel(fileLines)
  line = strtrim(fileLines{i});
  if isempty(line) || contains(upper(line), 'END')
    if collecting
      break;
    end
    continue
  end

  tokens = strsplit(line);
  if ~isempty(tokens)
    if strcmpi(tokens{1}, symb)
      collecting = true;
      raw_basis_lines{end+1,1} = line;
    elseif collecting
      % 
      if isletter(tokens{1}(1)) && ~strcmpi(tokens{1}, symb)
        break;
      else
        raw_basis_lines{end+1,1} = line;
      end
    end
  elseif collecting
    raw_basis_lines{end+1,1} = line;
  end
end

if ~collecting || isempty(raw_basis_lines)
  warning('No lines found for symbol %s in file %s', symb, basisfile);
  new_basis = {};
  return;
end

new_basis = parseNWChemBlock(raw_basis_lines);
end

function basisData = parseNWChemBlock(lines)
MAXL = 8;
MAPSPDF = struct('S',0, 'P',1, 'D',2, 'F',3, 'G',4, 'H',5, 'I',6);

basis_parsed = cell(1,MAXL);
for i=1:MAXL
  basis_parsed{i} = {};
end

current_l = -1;
sp_mode = false;

for i=1:numel(lines)
  line = strtrim(lines{i});
  if isempty(line)
    continue;
  end

  if isletter(line(1))
    tokens = strsplit(line);
    if numel(tokens)==1
      key = upper(tokens{1});
    else
      key = upper(tokens{end});
    end

    if strcmp(key, 'SP')
      sp_mode = true;
      basis_parsed{1}{end+1} = {0}; 
      basis_parsed{2}{end+1} = {1}; 
      current_l = -1;
    else
      sp_mode = false;
      if isfield(MAPSPDF, key)
        current_l = MAPSPDF.(key);
        basis_parsed{current_l+1}{end+1} = {current_l};
      end
    end
  else
    % numeric
    dat = strrep(line,'D','E');
    parts = strsplit(dat);
    if startsWith(line, '#')
      % It's just a comment line, so ignore
      continue
    end
    floats = cellfun(@str2double, parts, 'UniformOutput',false);
    if any(cellfun(@(x) isnan(x), floats))
      error('Failed to parse line: %s', line);
    end
    row = cell2mat(floats);

    if sp_mode
      if numel(row)>=3
        s_shell = basis_parsed{1}{end};
        p_shell = basis_parsed{2}{end};
        s_shell{end+1} = [row(1), row(2)];
        p_shell{end+1} = [row(1), row(3)];
        basis_parsed{1}{end} = s_shell;
        basis_parsed{2}{end} = p_shell;
      end
    else
      if current_l<0
        continue
      end
      shellcell = basis_parsed{current_l+1}{end};
      shellcell{end+1} = row;
      basis_parsed{current_l+1}{end} = shellcell;
    end
  end
end

% Flatten:
tempData = {};
for L=1:MAXL
  shells = basis_parsed{L};
  for s=1:numel(shells)
    tempData{end+1,1} = shells{s}; 
  end
end

% Remove zero-coeff lines
new_basis = {};
for iSh=1:numel(tempData)
  sh = tempData{iSh};
  Lval = sh{1};
  ec = sh(2:end);
  keepEC = {};
  for iEC=1:numel(ec)
    arr = ec{iEC};
    if any(abs(arr(2:end))>1e-14)
      keepEC{end+1} = arr;
    end
  end
  if ~isempty(keepEC)
    new_basis{end+1,1} = [{Lval}, keepEC]; 
  end
end

basisData = new_basis;
end