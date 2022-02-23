function [y] = hankelimag(z)
  y = imag(besselh(0, z(1) + i * z(2)));
end
