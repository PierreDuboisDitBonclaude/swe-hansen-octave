
function gb=gaussian_far(A,theta_beam,r, phys,grid)
  k=phys.k;
  theta=grid.theta_grid;
  phi=grid.phi_grid;
  eta_0_admit=phys.eta_0_admit;

  theta_beam=theta_beam/180*pi;
  b=(20*log10((1+cos(theta_beam))/2)-A)./(20*k*(1-cos(theta_beam))*log10(e));

  E_theta=exp(-1j*k.*r)./(k.*r).*exp(k*b*cos(theta)).*(1+cos(theta)).*(cos(phi));
  E_phi=exp(-1j*k.*r)./(k.*r).*exp(k*b*cos(theta)).*(1+cos(theta)).*(-sin(phi));
  E_r=zeros(size(E_phi));

  H_theta=-E_phi*eta_0_admit;
  H_phi=E_theta*eta_0_admit;
  H_r=zeros(size(E_phi));

  E=cat(3,E_r,E_theta,E_phi);
  H=cat(3,H_r,H_theta,H_phi);

  gb.taper=A;
  gb.taper_angle=theta_beam;
  gb.E_r=E_r;
  gb.E_theta=E_theta;
  gb.E_phi=E_phi;
  gb.H_r= H_r;
  gb.H_theta= H_theta;
  gb.H_phi= H_phi;
  gb.E=E;
  gb.H=H;
  gb.r=r;
  gb.grid=grid;
  gb.phys=phys;
  gb.b=b;

endfunction

