function phys=gen_phys(f)
  eps_0=8.854187818814*10^-12;
  mu_0=1.2566370612720*10^-6;
  c=(eps_0*mu_0)^-(0.5);
  lambda=c/f;
  k=2*pi/lambda;
  eta_0_admit=(eps_0/mu_0)^0.5;
  eta_0_imp=1/eta_0_admit;

  phys.eps_0=eps_0;
  phys.mu_0=mu_0;
  phys.c=c;
  phys.lambda=lambda;
  phys.k=k;
  phys.eta_0_admit=eta_0_admit;
  phys.eta_0_imp=eta_0_imp;
  phys.f=f;
endfunction

