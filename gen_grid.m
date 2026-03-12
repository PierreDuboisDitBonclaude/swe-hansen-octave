#Generate Grid
function G=gen_grid(N,name)
    if ( strcmp(name,'gl') )
      ps=2*N+1;
      ts=N+1;
      [mu wmu]=gauss_legendre_weights(ts,-1,1);
    elseif ( strcmp(name,'cc') )
      ps=2*N+1;
      ts=2*N+2;
      [mu wmu]=clenshaw_curtis_weights(ts);
    else
      error('no valid discretization input')
    endif
  dphi=2*pi/ps;
  dmu=2/(ts-1);
  theta=acos(mu);
  phi=(0:ps-1).'*dphi;
  [mu_grid,phi_grid]=meshgrid(mu,phi);
  theta_grid=acos(mu_grid);
  theta_grid_lin=acos(mu_grid(:)).';
  mu_grid_lin=mu_grid(:).';
  phi_grid_lin=phi_grid(:).';

  G.MN0=N+1;
  G.M=2*N+1;
  G.N=N;
  G.ts=ts;
  G.theta=theta;
  G.theta_grid=theta_grid;
  G.theta_grid_lin=theta_grid_lin;
  G.mu=mu;
  G.mu_grid=mu_grid;
  G.mu_grid_lin=mu_grid_lin;
  G.dmu=dmu;
  G.wmu=wmu;
  G.ps=ps;
  G.phi=phi;
  G.phi_grid=phi_grid;
  G.phi_grid_lin=phi_grid_lin;
  G.dphi=dphi;
  G.name=name;
endfunction
