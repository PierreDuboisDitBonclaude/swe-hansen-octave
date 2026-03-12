function [jn, yn, h1n, h2n] = sph_all_bessel(N, x)
  n  = 0:N;
  nu = n + 0.5;
  pref = sqrt(pi./(2*x));

  jn = pref .* besselj(nu, x);
  yn = pref .* bessely(nu, x);

  h1n = jn + 1i*yn;
  h2n = jn - 1i*yn;

  if x==0
    jn = [1, zeros(1,N)];
    yn(:) = -Inf;
    h1n = jn + 1i*yn;
    h2n = jn - 1i*yn;
  end
end

