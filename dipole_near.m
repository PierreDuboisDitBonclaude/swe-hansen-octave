function dp = dipole_near( r,phys,grid)

    k=phys.k;
    theta=grid.theta_grid;
    phi=grid.phi_grid;

    E_r     = phys.eta_0_imp* cos(theta)/2/pi .*exp(-1i*k*r) .* ( -1i./(k*r.^3) + 1/r.^2 );
    E_theta = phys.eta_0_imp* sin(theta)/4/pi .*exp(-1i*k*r) .* ( -1i./(k*r.^3) + 1/r.^2 + 1i*k/r );
    E_phi   = zeros( size(E_r) );

    H_r     = zeros( size(E_r) );
    H_theta = zeros( size(E_theta) );
    H_phi   = 1/4/pi .* sin(theta) .* ( 1i*k/r + 1/r^2 ) .* exp(-1i*k*r);

    # All field components need to have the dimensions of grid.theta_grid
    # E and H stack them together in the third dimension so E and H dimanesion finally are size(grid.theta_grid) x 3
    # Create a struct like below an the field is ready to use

    E=cat(3,E_r,E_theta,E_phi);
    H=cat(3,H_r,H_theta,H_phi);

    dp.E_r = E_r;
    dp.E_theta = E_theta;
    dp.E_phi=E_phi;

    dp.H_r=H_r;
    dp.H_theta=H_theta;
    dp.H_phi=H_phi;

    dp.E=E;
    dp.H=H;
    dp.r=r;
    dp.grid=grid;
    dp.phys=phys;

endfunction
