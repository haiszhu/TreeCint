function Vijkl = computeVijkl(nd, Norb, idcoefs, Vmunu, Vijkl)
%

collocation_matrix = zeros(nd,Norb,Norb);
tmpidx = 0;
for i = 1:Norb
  for j = 1:i
    tmpidx = tmpidx + 1;
    collocation_matrix(:,j,i) = idcoefs(:,tmpidx);
  end
end
for i = 1:Norb
  for j = i+1:Norb
    collocation_matrix(:,j,i) = collocation_matrix(:,i,j);
  end
end
collocation_matrix = reshape(collocation_matrix,[nd Norb^2]);

Vijkl = zeros(Norb,Norb,Norb,Norb);
for i = 1:Norb
  collocation_matrix_i = collocation_matrix(:,(i-1)*Norb+(1:Norb));
  for k = 1:Norb
    collocation_matrix_k = transpose(collocation_matrix(:,(k-1)*Norb+(1:Norb)));
    Vijkl(i,:,k,:) = (collocation_matrix_k*(Vmunu*collocation_matrix_i))';
  end
end
end
