function f=gb_near(grid,phys, b,r)
  % CSP Gaussian beam in spherical components (exact Maxwell field via electric Hertz vector)
  % Convention:
  %   theta ... polar angle [0, pi]
  %   phi   ... azimuth [0, 2*pi)
  %   r     ... scalar radius
  %   k     ... wavenumber [rad/m]
  %   b     ... complex-source shift parameter [m]
  %   E0    ... complex amplitude (optional, default 1)
  %   eta0  ... wave impedance (optional, default 376.730313668)
  %
  % Time convention: exp(-1i*omega*t)
  % Scalar CSP kernel:
  %   psi = E0 * exp(1i*k*Rtilde) / Rtilde,
  %   Rtilde = sqrt(x^2 + y^2 + (z - 1i*b)^2)
  % Hertz vector:
  %   Pi_e = xhat * psi
  % Exact fields in Cartesian coordinates:
  %   E = k^2 Pi_e + grad(div(Pi_e))
  %   H = -(1i*k/eta0) curl(Pi_e)
  % Then projected to spherical basis.
  eta0=phys.eta_0_imp;
  E0=1;
  k=phys.k;
  theta=grid.theta_grid;
  phi=grid.phi_grid;

  % Cartesian coordinates
  x = r .* sin(theta) .* cos(phi);
  y = r .* sin(theta) .* sin(phi);
  z = r .* cos(theta);
  zc = z - 1i*b;

  % Complex source distance
  Rtilde = sqrt(x.^2 + y.^2 + zc.^2);

  % Scalar kernel and derivative helpers
  expfac = E0 .* exp(1i*k.*Rtilde);
  psi = expfac ./ Rtilde;
  A = expfac .* (1i*k.*Rtilde - 1) ./ (Rtilde.^3);
  B = expfac .* (-(k.^2).*(Rtilde.^2) - 3i*k.*Rtilde + 3) ./ (Rtilde.^5);

  % Exact Cartesian field components
  Ex = k.^2 .* psi + A + B .* x.^2;
  Ey = B .* x .* y;
  Ez = B .* x .* zc;

  Hx = zeros(size(theta));
  Hy = -(1i*k./eta0) .* A .* zc;
  Hz = +(1i*k./eta0) .* A .* y;

  % Spherical unit-vector projections
  st = sin(theta); ct = cos(theta);
  sp = sin(phi);   cp = cos(phi);

  f.E_r     = Ex .* st .* cp + Ey .* st .* sp + Ez .* ct;
  f.E_theta = Ex .* ct .* cp + Ey .* ct .* sp - Ez .* st;
  f.E_phi   = -Ex .* sp + Ey .* cp;

  f.H_r     = Hx .* st .* cp + Hy .* st .* sp + Hz .* ct;
  f.H_theta = Hx .* ct .* cp + Hy .* ct .* sp - Hz .* st;
  f.H_phi   = -Hx .* sp + Hy .* cp;

  f.E=cat(3,f.E_r,f.E_theta,f.E_phi);
  f.H=cat(3,f.H_r,f.H_theta,f.H_phi);
  f.r=r;
  f.phys=phys;
  f.grid=grid;

end
