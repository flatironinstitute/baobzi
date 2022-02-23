function [y] = hankelreal(z)
  y = real(besselh(0, z(1) + i * z(2)));
end
