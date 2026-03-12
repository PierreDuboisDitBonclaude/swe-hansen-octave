function [mu, w] = clenshaw_curtis_weights(Ntheta)

  if (Ntheta < 2)
    error("Ntheta must be >= 2");
  endif

  N = Ntheta - 1;

  k = (0:N)';
  theta = pi * k / N;
  mu = cos(theta);

  w = zeros(Ntheta,1);

  kk = (2:N)';          % interior indices
  v = ones(N-1,1);

  if (mod(N,2) == 0)
    w0 = 1.0 / (N^2 - 1.0);
    w(1) = w0;
    w(end) = w0;

    for m = 1:(N/2 - 1)
      v = v - 2.0 * cos(2.0*m*theta(kk)) / (4.0*m^2 - 1.0);
    endfor

    v = v - cos(N*theta(kk)) / (N^2 - 1.0);

  else
    w0 = 1.0 / (N^2);
    w(1) = w0;
    w(end) = w0;

    for m = 1:((N+1)/2 - 1)
      v = v - 2.0 * cos(2.0*m*theta(kk)) / (4.0*m^2 - 1.0);
    endfor
  endif

  w(kk) = 2.0 * v / N;

endfunction

