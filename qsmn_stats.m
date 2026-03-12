function R = qsmn_stats(Q)

  Q = squeeze(Q);
  sz = size(Q);
  if numel(sz) ~= 3
    error('Q needs to be size SxMxN after squeeze.');
  end
  S = sz(1); M = sz(2); N = sz(3);

  P = real(0.5 * Q .* conj(Q));         % SxMxN
  P(~isfinite(P)) = -Inf;

  % Flatten + sort
  [vals, idx] = sort(P(:), 'descend');

  good = isfinite(vals) & (vals > -Inf);
  vals = vals(good);
  idx  = idx(good);

  if isempty(vals)
    R = zeros(0,7);
    return;
  end

  [s, m, n] = ind2sub([S M N], idx);


  if mod(M,2)==1
    m = m - (M+1)/2;
  end

  Psum = sum(vals);
  if Psum <= 0 || ~isfinite(Psum)
    Psum = 1;
  end

  Ppct = 100 * (vals / Psum);
  Pcum = cumsum(vals);
  PcumPct = 100 * (Pcum / Psum);

  R = [double(s(:)) double(m(:)) double(n(:)) double(vals(:)) double(Ppct(:)) double(Pcum(:)) double(PcumPct(:))];

end
