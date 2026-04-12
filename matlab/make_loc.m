function ao_loc = make_loc(bas, key)

% Define the columns used in bas
ANG_OF   = 2;  
NCTR_OF  = 4; 
KAPPA_OF = 5; 

% Decide which shape formula to use
if contains(key, 'cart')
  l = bas(:, ANG_OF);
  dims = ((l + 1) .* (l + 2)) ./ 2 .* bas(:, NCTR_OF);
elseif contains(key, 'sph')
  l = bas(:, ANG_OF);
  dims = (2 .* l + 1) .* bas(:, NCTR_OF);
else
  l = bas(:, ANG_OF);
  k = bas(:, KAPPA_OF);
  dims = (4 .* l + 2) .* bas(:, NCTR_OF);
  maskNeg = (k < 0);
  dims(maskNeg) = (2 .* l(maskNeg) + 2) .* bas(maskNeg, NCTR_OF);
  maskPos = (k > 0);
  dims(maskPos) = (2 .* l(maskPos)) .* bas(maskPos, NCTR_OF);
end

% Build the output "ao_loc" array, length(dims)+1
dims = dims;
ao_loc = zeros(1, numel(dims) + 1);
ao_loc(1) = 0;
ao_loc(2:end) = cumsum(dims);
end
