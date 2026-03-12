function T = octave_cmap_bins(clim, N)

  cmin = clim(1);
  cmax = clim(2);

  if (cmax <= cmin)
    error("clim muss cmax > cmin erfüllen");
  endif

  T = zeros(N, 3);

  for k = 1:N
    lower = cmin + (k - 1) * (cmax - cmin) / N;
    upper = cmin + k       * (cmax - cmin) / N;

    T(k, :) = [k, lower, upper];
  endfor
endfunction
