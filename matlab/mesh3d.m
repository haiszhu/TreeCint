function xyz = mesh3d(x, nx, y, ny, z, nz, xyz)
% xyz = zeros(3, nx * ny * nz);
ind = 0;
for iz = 1:nz
  for iy = 1:ny
    for ix = 1:nx
      ind = ind + 1;
      xyz(1, ind) = x(ix);
      xyz(2, ind) = y(iy);
      xyz(3, ind) = z(iz);
    end
  end
end
end
